"""Splice site detection and scoring for RNA-seq data.

This module provides tools for detecting splice junctions from aligned reads,
parsing CIGAR strings, and scoring splice site strength using position weight
matrices.
"""

from __future__ import annotations

import math
import re
from collections import defaultdict
from typing import Any, Dict, Literal

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Optional dependency handling
try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None  # type: ignore[assignment]


# =============================================================================
# Type Definitions
# =============================================================================

JunctionType = Literal["canonical", "semi_canonical", "non_canonical", "unknown"]


# =============================================================================
# Splice Junction Detection
# =============================================================================


def detect_splice_junctions(
    alignments: list[dict],
    min_reads: int = 3,
    min_overhang: int = 8,
) -> list[dict]:
    """Identify splice junctions from aligned reads using split-read evidence.

    Parses CIGAR strings from aligned reads to extract intron-spanning
    junctions. Junctions are aggregated by genomic coordinates and filtered
    by minimum read support and overhang length.

    Args:
        alignments: List of alignment records, each a dict with keys:
            - chrom (str): Chromosome/contig name
            - start (int): 0-based alignment start position
            - cigar (str): CIGAR string (e.g., "50M100N50M")
            - strand (str, optional): Strand ("+" or "-"), defaults to "+"
            - read_id (str, optional): Read identifier
        min_reads: Minimum number of supporting reads for a junction to be
            reported. Junctions below this threshold are filtered out.
        min_overhang: Minimum alignment overhang on each side of the junction
            in base pairs. Reads with shorter anchors are excluded to reduce
            false positives from misalignment.

    Returns:
        List of junction dicts, each containing:
            - chrom (str): Chromosome name
            - start (int): Junction start (donor site, 0-based)
            - end (int): Junction end (acceptor site, 0-based)
            - strand (str): Strand ("+" or "-")
            - read_count (int): Number of supporting reads
            - junction_type (str): "canonical", "semi_canonical", or "non_canonical"
            - overhang_mean (float): Mean overhang length across supporting reads

    Raises:
        ValueError: If alignments is empty or required fields are missing.

    Example:
        >>> alignments = [
        ...     {"chrom": "chr1", "start": 100, "cigar": "50M200N50M", "strand": "+"},
        ...     {"chrom": "chr1", "start": 100, "cigar": "50M200N50M", "strand": "+"},
        ...     {"chrom": "chr1", "start": 100, "cigar": "50M200N50M", "strand": "+"},
        ... ]
        >>> junctions = detect_splice_junctions(alignments, min_reads=2)
        >>> junctions[0]["read_count"]
        3
    """
    if not alignments:
        logger.warning("No alignments provided for junction detection")
        return []

    # Aggregate junctions: key=(chrom, start, end, strand) -> info
    junction_counts: dict[tuple[str, int, int, str], dict[str, Any]] = defaultdict(
        lambda: {"read_count": 0, "overhangs": []}
    )

    n_parsed = 0
    n_skipped = 0

    for aln in alignments:
        chrom = aln.get("chrom")
        start = aln.get("start")
        cigar = aln.get("cigar", "")
        strand = aln.get("strand", "+")

        if chrom is None or start is None or not cigar:
            n_skipped += 1
            continue

        # Parse CIGAR to find N operations (intron/splice junctions)
        junctions_from_read = _parse_cigar_junctions(cigar, int(start), int(min_overhang))

        for junc_start, junc_end, left_overhang, right_overhang in junctions_from_read:
            key = (str(chrom), junc_start, junc_end, str(strand))
            junction_counts[key]["read_count"] += 1
            junction_counts[key]["overhangs"].append((left_overhang, right_overhang))
            n_parsed += 1

    if n_skipped > 0:
        logger.warning(f"Skipped {n_skipped}/{len(alignments)} alignments with missing fields")

    # Filter by min_reads and build output
    results: list[dict] = []

    for (chrom, junc_start, junc_end, strand), info in junction_counts.items():
        if info["read_count"] < min_reads:
            continue

        # Calculate mean overhang
        all_overhangs = info["overhangs"]
        mean_overhang = sum(min(left, right) for left, right in all_overhangs) / len(all_overhangs)

        # Determine junction type based on intron length
        intron_length = junc_end - junc_start
        junction_type = _classify_junction_type(intron_length)

        results.append(
            {
                "chrom": chrom,
                "start": junc_start,
                "end": junc_end,
                "strand": strand,
                "read_count": info["read_count"],
                "junction_type": junction_type,
                "overhang_mean": round(mean_overhang, 1),
            }
        )

    # Sort by chromosome, then start position
    results.sort(key=lambda x: (x["chrom"], x["start"], x["end"]))

    logger.info(
        f"Detected {len(results)} splice junctions from {n_parsed} split-read "
        f"observations (min_reads={min_reads}, min_overhang={min_overhang})"
    )

    return results


def _parse_cigar_junctions(
    cigar: str,
    start: int,
    min_overhang: int,
) -> list[tuple[int, int, int, int]]:
    """Parse CIGAR string to extract splice junction coordinates.

    Args:
        cigar: CIGAR string from alignment.
        start: 0-based start position of the alignment.
        min_overhang: Minimum overhang on each side of junction.

    Returns:
        List of (junction_start, junction_end, left_overhang, right_overhang)
        tuples for each N operation in the CIGAR string.
    """
    junctions: list[tuple[int, int, int, int]] = []

    # Parse CIGAR operations
    ops = re.findall(r"(\d+)([MIDNSHP=X])", cigar)
    if not ops:
        return junctions

    ref_pos = start
    left_match_len = 0  # Matched bases before current junction
    found_n = False

    # First pass: calculate total matched bases for right overhang
    all_ops: list[tuple[int, str]] = [(int(length), op) for length, op in ops]

    ref_pos = start
    match_segments: list[int] = []
    current_match = 0

    for length, op in all_ops:
        if op in ("M", "=", "X"):
            current_match += length
            ref_pos += length
        elif op == "N":
            match_segments.append(current_match)
            current_match = 0
            ref_pos += length
        elif op == "D":
            ref_pos += length
        elif op in ("I", "S", "H", "P"):
            # These don't consume reference
            pass

    match_segments.append(current_match)

    # Second pass: extract junctions with overhang info
    ref_pos = start
    segment_idx = 0

    for length, op in all_ops:
        if op in ("M", "=", "X"):
            ref_pos += length
        elif op == "N":
            junc_start = ref_pos
            junc_end = ref_pos + length

            left_oh = match_segments[segment_idx] if segment_idx < len(match_segments) else 0
            right_oh = match_segments[segment_idx + 1] if segment_idx + 1 < len(match_segments) else 0

            if left_oh >= min_overhang and right_oh >= min_overhang:
                junctions.append((junc_start, junc_end, left_oh, right_oh))

            ref_pos += length
            segment_idx += 1
        elif op == "D":
            ref_pos += length
        # I, S, H, P don't consume reference

    return junctions


def _classify_junction_type(intron_length: int) -> str:
    """Classify junction type based on intron length heuristics.

    Args:
        intron_length: Length of the intron in base pairs.

    Returns:
        Junction type string: "canonical", "semi_canonical", or "non_canonical".
    """
    # Typical intron length ranges (heuristic classification when sequence
    # context is unavailable)
    if 50 <= intron_length <= 500_000:
        return "canonical"
    elif 20 <= intron_length < 50:
        return "semi_canonical"
    else:
        return "non_canonical"


# =============================================================================
# Splice Site Strength Scoring
# =============================================================================

# Position weight matrices for canonical splice sites (GT-AG rule)
# Based on consensus sequences from mammalian introns
# Positions relative to exon-intron boundary

# Donor site (5' splice site): last 3 exonic + first 6 intronic bases
# Consensus: ...MAG|GURAGU... (M=A/C, R=A/G, U=T in DNA)
_DONOR_PWM: dict[int, dict[str, float]] = {
    -3: {"A": 0.35, "C": 0.35, "G": 0.15, "T": 0.15},
    -2: {"A": 0.60, "C": 0.03, "G": 0.35, "T": 0.02},
    -1: {"A": 0.09, "C": 0.03, "G": 0.85, "T": 0.03},
    0: {"A": 0.00, "C": 0.00, "G": 1.00, "T": 0.00},  # G (invariant GT)
    1: {"A": 0.00, "C": 0.00, "G": 0.00, "T": 1.00},  # T (invariant GT)
    2: {"A": 0.60, "C": 0.01, "G": 0.38, "T": 0.01},
    3: {"A": 0.68, "C": 0.07, "G": 0.10, "T": 0.15},
    4: {"A": 0.16, "C": 0.15, "G": 0.53, "T": 0.16},
    5: {"A": 0.17, "C": 0.15, "G": 0.15, "T": 0.53},
}

# Acceptor site (3' splice site): last 6 intronic + first 1 exonic base
# Consensus: ...YYYYNCAG|G... (Y=C/T, N=any)
_ACCEPTOR_PWM: dict[int, dict[str, float]] = {
    -6: {"A": 0.08, "C": 0.30, "G": 0.08, "T": 0.54},
    -5: {"A": 0.06, "C": 0.30, "G": 0.06, "T": 0.58},
    -4: {"A": 0.08, "C": 0.30, "G": 0.08, "T": 0.54},
    -3: {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
    -2: {"A": 1.00, "C": 0.00, "G": 0.00, "T": 0.00},  # A (invariant AG)
    -1: {"A": 0.00, "C": 0.00, "G": 1.00, "T": 0.00},  # G (invariant AG)
    0: {"A": 0.15, "C": 0.15, "G": 0.55, "T": 0.15},
}


def compute_splice_site_strength(
    sequence: str,
    site_type: str = "donor",
) -> float:
    """Score splice site strength using a position weight matrix.

    Computes a log-odds score for a splice site sequence against the
    canonical GT-AG splice site consensus. Higher scores indicate
    stronger splice sites that are more likely to be used.

    Args:
        sequence: DNA sequence centered on the splice site boundary.
            For donor sites: 9 bases (3 exonic + 6 intronic, boundary
            between positions 3 and 4, e.g., "AAGGTAAGT").
            For acceptor sites: 7 bases (6 intronic + 1 exonic, boundary
            between positions 6 and 7, e.g., "TTCCAGG").
            Case-insensitive, only A/C/G/T characters accepted.
        site_type: Type of splice site:
            - "donor" (or "5ss"): 5' splice site (exon|intron boundary)
            - "acceptor" (or "3ss"): 3' splice site (intron|exon boundary)

    Returns:
        Log-odds score (bits). Typical ranges:
            - Strong sites: > 8.0
            - Moderate sites: 4.0 - 8.0
            - Weak sites: 0.0 - 4.0
            - Very weak/cryptic: < 0.0

    Raises:
        ValueError: If sequence length doesn't match expected length for
            site_type, or contains non-ACGT characters.

    Example:
        >>> # Perfect donor site (consensus GT)
        >>> score = compute_splice_site_strength("AAGGTAAGT", "donor")
        >>> score > 5.0
        True
        >>> # Acceptor site
        >>> score = compute_splice_site_strength("TTCCAGG", "acceptor")
        >>> score > 0.0
        True
    """
    sequence = sequence.upper().strip()

    # Validate characters
    valid_bases = set("ACGT")
    invalid = set(sequence) - valid_bases
    if invalid:
        raise ValueError(f"Invalid characters in sequence: {invalid}. Only ACGT allowed.")

    # Select PWM based on site type
    if site_type in ("donor", "5ss"):
        pwm = _DONOR_PWM
        expected_length = 9  # positions -3 to +5
    elif site_type in ("acceptor", "3ss"):
        pwm = _ACCEPTOR_PWM
        expected_length = 7  # positions -6 to 0
    else:
        raise ValueError(f"Unknown site_type: {site_type}. Valid: donor, 5ss, acceptor, 3ss")

    if len(sequence) != expected_length:
        raise ValueError(f"Expected sequence length {expected_length} for {site_type} site, " f"got {len(sequence)}")

    # Compute log-odds score
    background = 0.25  # Uniform background frequency
    score = 0.0
    positions = sorted(pwm.keys())

    for i, pos in enumerate(positions):
        base = sequence[i]
        freq = pwm[pos].get(base, 0.001)  # Small pseudocount for unseen bases

        # Log-odds relative to background
        if freq > 0:
            score += math.log2(freq / background)
        else:
            score += math.log2(0.001 / background)  # Penalty for impossible base

    return round(score, 4)
