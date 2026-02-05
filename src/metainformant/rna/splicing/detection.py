"""Alternative splicing detection and quantification for RNA-seq data.

This module provides tools for identifying splice junctions from aligned
reads, classifying alternative splicing events, computing Percent Spliced
In (PSI) values, testing for differential splicing between conditions,
discovering novel junctions, and scoring splice site strength.

All implementations are pure Python with optional numpy/scipy acceleration.
No external R or specialized splicing tool dependencies required.

Main Functions:
    Junction Detection:
        - detect_splice_junctions: Identify splice junctions from aligned reads
        - find_novel_junctions: Discover unannotated splice junctions

    Event Classification:
        - classify_splicing_events: Classify junctions by splicing event type

    Quantification:
        - compute_psi: Percent Spliced In with confidence interval
        - differential_splicing: Test for differential splicing between groups

    Sequence Analysis:
        - compute_splice_site_strength: Score splice site using position weight matrix

Example:
    >>> from metainformant.rna.splicing import detection
    >>> alignments = [
    ...     {"chrom": "chr1", "start": 100, "end": 200, "cigar": "50M100N50M"},
    ...     {"chrom": "chr1", "start": 100, "end": 200, "cigar": "50M100N50M"},
    ... ]
    >>> junctions = detection.detect_splice_junctions(alignments, min_reads=1)
    >>> psi = detection.compute_psi(inclusion_reads=30, exclusion_reads=10)
"""

from __future__ import annotations

import math
import re
from collections import defaultdict
from typing import Any, Dict, List, Literal, Optional, Sequence, Tuple

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Optional dependency handling
try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None  # type: ignore[assignment]

try:
    from scipy import stats as scipy_stats
    from scipy.special import beta as beta_func
    from scipy.special import betainc

    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    scipy_stats = None  # type: ignore[assignment]


# =============================================================================
# Type Definitions
# =============================================================================

SplicingEventType = Literal[
    "exon_skipping",
    "intron_retention",
    "alt_5ss",
    "alt_3ss",
    "mutually_exclusive",
    "unknown",
]

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
# Splicing Event Classification
# =============================================================================


def classify_splicing_events(
    junctions: list[dict],
    gene_model: dict | None = None,
) -> list[dict]:
    """Classify splice junctions into alternative splicing event types.

    Uses junction coordinates and optionally a gene model to classify each
    junction as one of the standard alternative splicing event types.

    Without a gene model, classification uses coordinate-based heuristics
    comparing junctions that share donor or acceptor sites. With a gene
    model, exon annotations enable precise event typing.

    Args:
        junctions: List of junction dicts from detect_splice_junctions(),
            each with keys: chrom, start, end, strand, read_count.
        gene_model: Optional gene annotation dict with structure:
            {
                "genes": {
                    "gene_id": {
                        "chrom": str,
                        "strand": str,
                        "exons": [{"start": int, "end": int, "exon_id": str}, ...],
                        "introns": [{"start": int, "end": int}, ...]
                    }
                }
            }
            If None, classification uses coordinate-based inference only.

    Returns:
        List of event dicts, each containing:
            - event_type (str): One of "exon_skipping", "intron_retention",
              "alt_5ss", "alt_3ss", "mutually_exclusive", "unknown"
            - junctions (list[dict]): The junction(s) involved in this event
            - chrom (str): Chromosome
            - region_start (int): Start of the affected region
            - region_end (int): End of the affected region
            - confidence (float): Classification confidence score (0-1)

    Example:
        >>> junctions = [
        ...     {"chrom": "chr1", "start": 100, "end": 300, "strand": "+", "read_count": 10},
        ...     {"chrom": "chr1", "start": 100, "end": 200, "strand": "+", "read_count": 5},
        ...     {"chrom": "chr1", "start": 200, "end": 300, "strand": "+", "read_count": 5},
        ... ]
        >>> events = classify_splicing_events(junctions)
    """
    if not junctions:
        return []

    events: list[dict] = []

    # Group junctions by chromosome and strand
    grouped: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for junc in junctions:
        key = (junc["chrom"], junc.get("strand", "+"))
        grouped[key].append(junc)

    for (chrom, strand), chrom_junctions in grouped.items():
        # Sort by start position
        sorted_juncs = sorted(chrom_junctions, key=lambda x: (x["start"], x["end"]))

        # If gene model available, use annotation-based classification
        if gene_model is not None:
            annotation_events = _classify_with_gene_model(sorted_juncs, chrom, strand, gene_model)
            events.extend(annotation_events)
            continue

        # Coordinate-based classification
        coord_events = _classify_by_coordinates(sorted_juncs, chrom, strand)
        events.extend(coord_events)

    logger.info(f"Classified {len(events)} splicing events from {len(junctions)} junctions")

    return events


def _classify_with_gene_model(
    junctions: list[dict],
    chrom: str,
    strand: str,
    gene_model: dict,
) -> list[dict]:
    """Classify junctions using gene model annotations.

    Args:
        junctions: Sorted junctions for one chromosome/strand.
        chrom: Chromosome name.
        strand: Strand.
        gene_model: Gene annotation dictionary.

    Returns:
        List of classified event dicts.
    """
    events: list[dict] = []
    genes = gene_model.get("genes", {})

    # Find relevant genes on this chromosome/strand
    relevant_genes = {
        gid: ginfo for gid, ginfo in genes.items() if ginfo.get("chrom") == chrom and ginfo.get("strand", "+") == strand
    }

    for junc in junctions:
        junc_start = junc["start"]
        junc_end = junc["end"]
        best_event_type: SplicingEventType = "unknown"
        best_confidence = 0.0

        for gid, ginfo in relevant_genes.items():
            exons = sorted(ginfo.get("exons", []), key=lambda e: e["start"])
            introns = ginfo.get("introns", [])

            # Check for intron retention: junction spans exactly one annotated intron
            for intron in introns:
                if abs(junc_start - intron["start"]) <= 5 and abs(junc_end - intron["end"]) <= 5:
                    best_event_type = "intron_retention"
                    best_confidence = 0.9
                    break

            # Check for exon skipping: junction skips one or more internal exons
            if best_event_type == "unknown":
                skipped_exons = [e for e in exons if e["start"] > junc_start and e["end"] < junc_end]
                if skipped_exons:
                    best_event_type = "exon_skipping"
                    best_confidence = 0.85

            # Check for alt 5' or 3' splice sites
            if best_event_type == "unknown":
                for other_junc in junctions:
                    if other_junc is junc:
                        continue
                    # Same acceptor, different donor -> alt 5'ss
                    if abs(other_junc["end"] - junc_end) <= 5 and abs(other_junc["start"] - junc_start) > 5:
                        best_event_type = "alt_5ss"
                        best_confidence = 0.7
                        break
                    # Same donor, different acceptor -> alt 3'ss
                    if abs(other_junc["start"] - junc_start) <= 5 and abs(other_junc["end"] - junc_end) > 5:
                        best_event_type = "alt_3ss"
                        best_confidence = 0.7
                        break

        events.append(
            {
                "event_type": best_event_type,
                "junctions": [junc],
                "chrom": chrom,
                "region_start": junc_start,
                "region_end": junc_end,
                "confidence": best_confidence,
            }
        )

    return events


def _classify_by_coordinates(
    junctions: list[dict],
    chrom: str,
    strand: str,
) -> list[dict]:
    """Classify junctions using coordinate-based heuristics.

    When no gene model is available, infers event types from the spatial
    relationships between junctions sharing donor or acceptor sites.

    Args:
        junctions: Sorted junctions for one chromosome/strand.
        chrom: Chromosome name.
        strand: Strand.

    Returns:
        List of classified event dicts.
    """
    events: list[dict] = []

    # Index junctions by donor (start) and acceptor (end) sites
    by_donor: dict[int, list[dict]] = defaultdict(list)
    by_acceptor: dict[int, list[dict]] = defaultdict(list)

    for junc in junctions:
        by_donor[junc["start"]].append(junc)
        by_acceptor[junc["end"]].append(junc)

    classified_juncs: set[int] = set()

    # Detect exon skipping: a long junction that spans the region of two
    # shorter junctions sharing an internal acceptor/donor
    for i, junc in enumerate(junctions):
        if id(junc) in classified_juncs:
            continue

        junc_start = junc["start"]
        junc_end = junc["end"]

        # Look for shorter junctions that together span the same region
        for j, other in enumerate(junctions):
            if i == j or id(other) in classified_juncs:
                continue

            # Check if other junction is contained within this one
            if other["start"] >= junc_start and other["end"] <= junc_end:
                # Find a complementary junction
                for k, third in enumerate(junctions):
                    if k in (i, j) or id(third) in classified_juncs:
                        continue
                    # Complementary: third starts near other's end, ends near junc's end
                    if abs(third["start"] - other["end"]) <= 10 and abs(third["end"] - junc_end) <= 10:
                        events.append(
                            {
                                "event_type": "exon_skipping",
                                "junctions": [junc, other, third],
                                "chrom": chrom,
                                "region_start": junc_start,
                                "region_end": junc_end,
                                "confidence": 0.7,
                            }
                        )
                        classified_juncs.update([id(junc), id(other), id(third)])
                        break

    # Detect alternative 5' splice sites: same acceptor, different donors
    for acceptor, juncs in by_acceptor.items():
        unclassified = [j for j in juncs if id(j) not in classified_juncs]
        if len(unclassified) >= 2:
            events.append(
                {
                    "event_type": "alt_5ss",
                    "junctions": unclassified,
                    "chrom": chrom,
                    "region_start": min(j["start"] for j in unclassified),
                    "region_end": acceptor,
                    "confidence": 0.6,
                }
            )
            for j in unclassified:
                classified_juncs.add(id(j))

    # Detect alternative 3' splice sites: same donor, different acceptors
    for donor, juncs in by_donor.items():
        unclassified = [j for j in juncs if id(j) not in classified_juncs]
        if len(unclassified) >= 2:
            events.append(
                {
                    "event_type": "alt_3ss",
                    "junctions": unclassified,
                    "chrom": chrom,
                    "region_start": donor,
                    "region_end": max(j["end"] for j in unclassified),
                    "confidence": 0.6,
                }
            )
            for j in unclassified:
                classified_juncs.add(id(j))

    # Remaining junctions are unclassified
    for junc in junctions:
        if id(junc) not in classified_juncs:
            events.append(
                {
                    "event_type": "unknown",
                    "junctions": [junc],
                    "chrom": chrom,
                    "region_start": junc["start"],
                    "region_end": junc["end"],
                    "confidence": 0.0,
                }
            )

    return events


# =============================================================================
# PSI Computation
# =============================================================================


def compute_psi(
    inclusion_reads: int,
    exclusion_reads: int,
    effective_length_ratio: float = 1.0,
) -> dict:
    """Compute Percent Spliced In (PSI) with confidence interval.

    PSI measures the fraction of transcripts that include a particular
    exon or splice junction. Uses a Beta distribution model to estimate
    the PSI value and its confidence interval.

    The formula is:
        PSI = (inclusion / effective_length_ratio) /
              (inclusion / effective_length_ratio + exclusion)

    Args:
        inclusion_reads: Number of reads supporting the inclusion isoform
            (e.g., reads spanning the included exon).
        exclusion_reads: Number of reads supporting the exclusion isoform
            (e.g., reads spanning the skipping junction).
        effective_length_ratio: Ratio of effective lengths between inclusion
            and exclusion isoforms. Accounts for length bias in read
            assignment. Default 1.0 assumes equal effective lengths.

    Returns:
        Dictionary with keys:
            - psi (float): Point estimate of PSI (0.0 to 1.0)
            - ci_low (float): Lower bound of 95% confidence interval
            - ci_high (float): Upper bound of 95% confidence interval
            - total_reads (int): Total supporting reads
            - inclusion_reads (int): Input inclusion count
            - exclusion_reads (int): Input exclusion count
            - method (str): "beta_distribution"

    Raises:
        ValueError: If read counts are negative or effective_length_ratio <= 0.

    Example:
        >>> result = compute_psi(inclusion_reads=30, exclusion_reads=10)
        >>> 0.7 < result["psi"] < 0.8
        True
        >>> result["ci_low"] < result["psi"] < result["ci_high"]
        True
    """
    if inclusion_reads < 0 or exclusion_reads < 0:
        raise ValueError(
            f"Read counts must be non-negative: inclusion={inclusion_reads}, " f"exclusion={exclusion_reads}"
        )
    if effective_length_ratio <= 0:
        raise ValueError(f"effective_length_ratio must be positive, got {effective_length_ratio}")

    total_reads = inclusion_reads + exclusion_reads

    if total_reads == 0:
        logger.warning("Zero total reads, PSI undefined")
        return {
            "psi": 0.0,
            "ci_low": 0.0,
            "ci_high": 1.0,
            "total_reads": 0,
            "inclusion_reads": 0,
            "exclusion_reads": 0,
            "method": "beta_distribution",
        }

    # Adjust inclusion reads by effective length ratio
    adjusted_inclusion = inclusion_reads / effective_length_ratio
    psi = adjusted_inclusion / (adjusted_inclusion + exclusion_reads)

    # Beta distribution parameters (Bayesian with uniform prior)
    # Using adjusted counts as pseudo-observations
    alpha = adjusted_inclusion + 1  # +1 for uniform prior
    beta_param = exclusion_reads + 1

    # Compute 95% confidence interval using Beta distribution
    ci_low, ci_high = _beta_confidence_interval(alpha, beta_param, confidence=0.95)

    return {
        "psi": round(psi, 6),
        "ci_low": round(ci_low, 6),
        "ci_high": round(ci_high, 6),
        "total_reads": total_reads,
        "inclusion_reads": inclusion_reads,
        "exclusion_reads": exclusion_reads,
        "method": "beta_distribution",
    }


def _beta_confidence_interval(
    alpha: float,
    beta_param: float,
    confidence: float = 0.95,
) -> tuple[float, float]:
    """Compute confidence interval from Beta distribution.

    Args:
        alpha: Alpha parameter of Beta distribution.
        beta_param: Beta parameter of Beta distribution.
        confidence: Confidence level (e.g., 0.95 for 95% CI).

    Returns:
        Tuple of (lower_bound, upper_bound).
    """
    tail = (1 - confidence) / 2

    if HAS_SCIPY:
        ci_low = float(scipy_stats.beta.ppf(tail, alpha, beta_param))
        ci_high = float(scipy_stats.beta.ppf(1 - tail, alpha, beta_param))
    else:
        # Pure Python approximation using normal approximation to Beta
        mean = alpha / (alpha + beta_param)
        variance = (alpha * beta_param) / ((alpha + beta_param) ** 2 * (alpha + beta_param + 1))
        std = math.sqrt(variance)

        # z-score for 95% CI
        z = 1.96 if confidence == 0.95 else _norm_ppf(1 - tail)
        ci_low = max(0.0, mean - z * std)
        ci_high = min(1.0, mean + z * std)

    return ci_low, ci_high


def _norm_ppf(p: float) -> float:
    """Approximate normal distribution percent point function (inverse CDF).

    Uses the Abramowitz and Stegun rational approximation.

    Args:
        p: Probability (0 < p < 1).

    Returns:
        z-score corresponding to the given probability.
    """
    if p <= 0 or p >= 1:
        return 0.0

    # Rational approximation constants
    a = [
        -3.969683028665376e01,
        2.209460984245205e02,
        -2.759285104469687e02,
        1.383577518672690e02,
        -3.066479806614716e01,
        2.506628277459239e00,
    ]
    b = [
        -5.447609879822406e01,
        1.615858368580409e02,
        -1.556989798598866e02,
        6.680131188771972e01,
        -1.328068155288572e01,
    ]

    if p < 0.5:
        q = math.sqrt(-2 * math.log(p))
        num = ((((a[0] * q + a[1]) * q + a[2]) * q + a[3]) * q + a[4]) * q + a[5]
        den = ((((b[0] * q + b[1]) * q + b[2]) * q + b[3]) * q + b[4]) * q + 1
        return num / den
    else:
        q = math.sqrt(-2 * math.log(1 - p))
        num = ((((a[0] * q + a[1]) * q + a[2]) * q + a[3]) * q + a[4]) * q + a[5]
        den = ((((b[0] * q + b[1]) * q + b[2]) * q + b[3]) * q + b[4]) * q + 1
        return -(num / den)


# =============================================================================
# Differential Splicing
# =============================================================================


def differential_splicing(
    psi_group1: list[float],
    psi_group2: list[float],
    method: str = "empirical_bayes",
) -> dict:
    """Test for differential splicing between two conditions.

    Compares PSI distributions between two groups using delta-PSI
    (difference in mean PSI) and a statistical test to assess significance.

    Supports empirical Bayes shrinkage to improve estimates with small
    sample sizes, as well as standard non-parametric tests.

    Args:
        psi_group1: PSI values for condition 1 (e.g., control).
            Each value should be between 0 and 1.
        psi_group2: PSI values for condition 2 (e.g., treatment).
            Each value should be between 0 and 1.
        method: Statistical method for testing:
            - "empirical_bayes": Empirical Bayes shrinkage with moderated
              t-test (recommended for small samples)
            - "wilcoxon": Wilcoxon rank-sum (Mann-Whitney U) test
            - "permutation": Permutation test on delta-PSI

    Returns:
        Dictionary with keys:
            - delta_psi (float): Mean PSI difference (group2 - group1)
            - p_value (float): Statistical significance
            - mean_psi_group1 (float): Mean PSI in group 1
            - mean_psi_group2 (float): Mean PSI in group 2
            - sd_group1 (float): Standard deviation in group 1
            - sd_group2 (float): Standard deviation in group 2
            - method (str): Method used for testing
            - significant (bool): Whether delta_psi is significant
              (|delta_psi| >= 0.1 and p_value < 0.05)
            - n_group1 (int): Sample size group 1
            - n_group2 (int): Sample size group 2

    Raises:
        ValueError: If either group is empty or contains values outside [0, 1].

    Example:
        >>> control_psi = [0.8, 0.75, 0.82, 0.78]
        >>> treatment_psi = [0.4, 0.45, 0.38, 0.42]
        >>> result = differential_splicing(control_psi, treatment_psi)
        >>> result["delta_psi"] < -0.3
        True
    """
    if not psi_group1 or not psi_group2:
        raise ValueError("Both groups must have at least one PSI value")

    # Validate PSI ranges
    for i, val in enumerate(psi_group1):
        if not (0.0 <= val <= 1.0):
            raise ValueError(f"PSI values must be in [0,1], got {val} in group1[{i}]")
    for i, val in enumerate(psi_group2):
        if not (0.0 <= val <= 1.0):
            raise ValueError(f"PSI values must be in [0,1], got {val} in group2[{i}]")

    n1 = len(psi_group1)
    n2 = len(psi_group2)

    mean1 = sum(psi_group1) / n1
    mean2 = sum(psi_group2) / n2
    delta_psi = mean2 - mean1

    # Standard deviations
    if n1 > 1:
        var1 = sum((x - mean1) ** 2 for x in psi_group1) / (n1 - 1)
        sd1 = math.sqrt(var1)
    else:
        var1 = 0.0
        sd1 = 0.0

    if n2 > 1:
        var2 = sum((x - mean2) ** 2 for x in psi_group2) / (n2 - 1)
        sd2 = math.sqrt(var2)
    else:
        var2 = 0.0
        sd2 = 0.0

    # Statistical test
    if method == "empirical_bayes":
        p_value = _empirical_bayes_test(psi_group1, psi_group2, mean1, mean2, var1, var2)
    elif method == "wilcoxon":
        p_value = _wilcoxon_test(psi_group1, psi_group2)
    elif method == "permutation":
        p_value = _permutation_test(psi_group1, psi_group2, n_permutations=10000)
    else:
        raise ValueError(f"Unknown method: {method}. Valid: empirical_bayes, wilcoxon, permutation")

    # Significance threshold: |delta_PSI| >= 0.1 and p < 0.05
    significant = abs(delta_psi) >= 0.1 and p_value < 0.05

    return {
        "delta_psi": round(delta_psi, 6),
        "p_value": round(p_value, 8),
        "mean_psi_group1": round(mean1, 6),
        "mean_psi_group2": round(mean2, 6),
        "sd_group1": round(sd1, 6),
        "sd_group2": round(sd2, 6),
        "method": method,
        "significant": significant,
        "n_group1": n1,
        "n_group2": n2,
    }


def _empirical_bayes_test(
    group1: list[float],
    group2: list[float],
    mean1: float,
    mean2: float,
    var1: float,
    var2: float,
) -> float:
    """Moderated t-test with empirical Bayes variance shrinkage.

    Shrinks per-group variances toward a pooled prior, improving
    statistical power for small sample sizes (limma-style approach).

    Args:
        group1: PSI values group 1.
        group2: PSI values group 2.
        mean1: Mean of group 1.
        mean2: Mean of group 2.
        var1: Variance of group 1.
        var2: Variance of group 2.

    Returns:
        p-value from the moderated t-test.
    """
    n1 = len(group1)
    n2 = len(group2)

    if n1 < 2 or n2 < 2:
        # Cannot perform t-test with n=1; use simple threshold
        return 1.0 if abs(mean2 - mean1) < 0.1 else 0.05

    # Pool variance with prior (shrinkage toward common variance)
    pooled_var = ((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2)

    # Empirical Bayes shrinkage: blend sample variance with prior
    prior_df = 3.0  # Prior degrees of freedom (tunable)
    prior_var = pooled_var  # Use pooled as prior

    shrunk_var1 = (prior_df * prior_var + (n1 - 1) * var1) / (prior_df + n1 - 1)
    shrunk_var2 = (prior_df * prior_var + (n2 - 1) * var2) / (prior_df + n2 - 1)

    # Moderated t-statistic
    se = math.sqrt(shrunk_var1 / n1 + shrunk_var2 / n2)

    if se < 1e-10:
        return 1.0 if abs(mean2 - mean1) < 1e-10 else 0.0

    t_stat = (mean2 - mean1) / se

    # Degrees of freedom (Welch-Satterthwaite with inflated df from prior)
    df = n1 + n2 - 2 + 2 * prior_df

    # p-value from t-distribution
    if HAS_SCIPY:
        p_value = float(2 * scipy_stats.t.sf(abs(t_stat), df))
    else:
        # Approximate using normal for large df
        p_value = 2 * (1 - _normal_cdf(abs(t_stat)))

    return min(p_value, 1.0)


def _wilcoxon_test(group1: list[float], group2: list[float]) -> float:
    """Wilcoxon rank-sum (Mann-Whitney U) test for PSI comparison.

    Args:
        group1: PSI values group 1.
        group2: PSI values group 2.

    Returns:
        p-value from the test.
    """
    if HAS_SCIPY:
        try:
            _, p_value = scipy_stats.mannwhitneyu(group1, group2, alternative="two-sided")
            return float(p_value)
        except ValueError:
            return 1.0

    # Pure Python Mann-Whitney U implementation
    n1 = len(group1)
    n2 = len(group2)
    combined = [(val, 0) for val in group1] + [(val, 1) for val in group2]
    combined.sort(key=lambda x: x[0])

    # Assign ranks (handle ties by average rank)
    ranks: list[float] = []
    i = 0
    while i < len(combined):
        j = i
        while j < len(combined) and combined[j][0] == combined[i][0]:
            j += 1
        avg_rank = (i + j + 1) / 2  # 1-based average rank
        for _ in range(i, j):
            ranks.append(avg_rank)
        i = j

    # U statistic
    r1 = sum(ranks[k] for k in range(len(combined)) if combined[k][1] == 0)
    u1 = r1 - n1 * (n1 + 1) / 2

    # Normal approximation for p-value
    mu = n1 * n2 / 2
    sigma = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)

    if sigma < 1e-10:
        return 1.0

    z = (u1 - mu) / sigma
    p_value = 2 * (1 - _normal_cdf(abs(z)))
    return min(p_value, 1.0)


def _permutation_test(
    group1: list[float],
    group2: list[float],
    n_permutations: int = 10000,
) -> float:
    """Permutation test for differential splicing.

    Args:
        group1: PSI values group 1.
        group2: PSI values group 2.
        n_permutations: Number of permutations.

    Returns:
        Permutation p-value.
    """
    import random

    observed_diff = abs(sum(group2) / len(group2) - sum(group1) / len(group1))
    combined = group1 + group2
    n1 = len(group1)
    count_extreme = 0

    for _ in range(n_permutations):
        random.shuffle(combined)
        perm_g1 = combined[:n1]
        perm_g2 = combined[n1:]
        perm_diff = abs(sum(perm_g2) / len(perm_g2) - sum(perm_g1) / len(perm_g1))
        if perm_diff >= observed_diff:
            count_extreme += 1

    p_value = (count_extreme + 1) / (n_permutations + 1)
    return p_value


def _normal_cdf(x: float) -> float:
    """Approximate standard normal CDF using error function approximation.

    Args:
        x: z-score value.

    Returns:
        Cumulative probability P(Z <= x).
    """
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))


# =============================================================================
# Novel Junction Discovery
# =============================================================================


def find_novel_junctions(
    junctions: list[dict],
    known_junctions: list[dict],
    max_distance: int = 10,
) -> list[dict]:
    """Identify unannotated splice junctions not present in a reference set.

    Compares detected junctions against a set of known/annotated junctions
    and reports those that are not within max_distance of any known junction.

    Args:
        junctions: Detected junctions from detect_splice_junctions().
        known_junctions: Reference set of annotated junctions, each with
            keys: chrom, start, end. Strand is optional.
        max_distance: Maximum coordinate distance (in bp) for a detected
            junction to be considered matching a known junction. Both start
            and end must be within this distance.

    Returns:
        List of novel junction dicts. Each dict is a copy of the input
        junction with an additional key:
            - nearest_known_distance (int): Distance to the nearest known
              junction (sum of start and end offsets), or -1 if no known
              junctions exist on the same chromosome.

    Example:
        >>> detected = [
        ...     {"chrom": "chr1", "start": 1000, "end": 2000, "strand": "+", "read_count": 5},
        ...     {"chrom": "chr1", "start": 3000, "end": 4000, "strand": "+", "read_count": 8},
        ... ]
        >>> known = [
        ...     {"chrom": "chr1", "start": 1002, "end": 1998},
        ... ]
        >>> novel = find_novel_junctions(detected, known, max_distance=10)
        >>> len(novel)  # Only the chr1:3000-4000 junction is novel
        1
    """
    if not junctions:
        return []

    if not known_junctions:
        logger.info("No known junctions provided, all junctions are novel")
        results = []
        for junc in junctions:
            novel = dict(junc)
            novel["nearest_known_distance"] = -1
            results.append(novel)
        return results

    # Index known junctions by chromosome for fast lookup
    known_by_chrom: dict[str, list[dict]] = defaultdict(list)
    for kj in known_junctions:
        known_by_chrom[kj["chrom"]].append(kj)

    novel: list[dict] = []

    for junc in junctions:
        chrom = junc["chrom"]
        junc_start = junc["start"]
        junc_end = junc["end"]

        chrom_known = known_by_chrom.get(chrom, [])

        if not chrom_known:
            result = dict(junc)
            result["nearest_known_distance"] = -1
            novel.append(result)
            continue

        # Find nearest known junction
        min_dist = float("inf")
        is_known = False

        for kj in chrom_known:
            start_dist = abs(junc_start - kj["start"])
            end_dist = abs(junc_end - kj["end"])

            if start_dist <= max_distance and end_dist <= max_distance:
                is_known = True
                break

            total_dist = start_dist + end_dist
            if total_dist < min_dist:
                min_dist = total_dist

        if not is_known:
            result = dict(junc)
            result["nearest_known_distance"] = int(min_dist) if min_dist != float("inf") else -1
            novel.append(result)

    logger.info(f"Found {len(novel)} novel junctions out of {len(junctions)} total " f"(max_distance={max_distance})")

    return novel


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
