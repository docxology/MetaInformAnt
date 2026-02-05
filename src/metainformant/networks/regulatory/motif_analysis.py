"""Transcription factor binding motif analysis.

This module provides tools for discovering and scoring transcription factor
binding site motifs in regulatory sequences, including k-mer overrepresentation
analysis, position weight matrix (PWM) construction and scoring, and motif
scanning against known motif libraries.
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Optional dependencies
try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None  # type: ignore[assignment]

# Nucleotide constants
NUCLEOTIDES = ["A", "C", "G", "T"]
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


def _reverse_complement(seq: str) -> str:
    """Compute reverse complement of a DNA sequence.

    Args:
        seq: DNA sequence string.

    Returns:
        Reverse complement string.
    """
    return "".join(COMPLEMENT.get(c, "N") for c in reversed(seq.upper()))


def _information_content(pwm: list[dict[str, float]]) -> float:
    """Compute total information content of a PWM.

    Args:
        pwm: Position weight matrix as list of dicts mapping nucleotides
            to probabilities.

    Returns:
        Total information content in bits.
    """
    ic = 0.0
    for pos in pwm:
        for nt in NUCLEOTIDES:
            freq = pos.get(nt, 0.0)
            if freq > 0:
                ic += freq * math.log2(freq / 0.25)
    return ic


def build_pwm(
    sequences: list[str],
    pseudocount: float = 0.01,
) -> list[dict[str, float]]:
    """Build a position weight matrix from aligned sequences.

    Constructs a PWM from a set of equal-length aligned DNA sequences by
    counting nucleotide frequencies at each position and applying a
    pseudocount for numerical stability.

    Args:
        sequences: List of equal-length DNA sequences (aligned).
        pseudocount: Pseudocount added to each nucleotide frequency to
            avoid log(0). Default is 0.01.

    Returns:
        PWM as list of dicts, one per position, mapping nucleotides
        (A, C, G, T) to probability values.

    Raises:
        ValueError: If sequences are empty or have unequal lengths.
    """
    if not sequences:
        raise ValueError("sequences must not be empty")

    seq_len = len(sequences[0])
    for i, seq in enumerate(sequences):
        if len(seq) != seq_len:
            raise ValueError(
                f"All sequences must have equal length. " f"Sequence {i} has length {len(seq)}, expected {seq_len}"
            )

    logger.info("Building PWM from %d sequences of length %d", len(sequences), seq_len)

    n_seqs = len(sequences)
    pwm: list[dict[str, float]] = []

    for pos in range(seq_len):
        counts: dict[str, float] = {nt: pseudocount for nt in NUCLEOTIDES}

        for seq in sequences:
            nt = seq[pos].upper()
            if nt in counts:
                counts[nt] += 1.0

        total = sum(counts.values())
        freqs = {nt: counts[nt] / total for nt in NUCLEOTIDES}
        pwm.append(freqs)

    ic = _information_content(pwm)
    logger.info("PWM built: %d positions, information content = %.2f bits", seq_len, ic)

    return pwm


def score_motif_match(
    sequence: str,
    pwm: list[dict[str, float]],
) -> float:
    """Score a sequence against a position weight matrix.

    Computes a log-likelihood ratio score for the given sequence relative
    to the PWM, comparing the probability under the motif model to a
    uniform background (0.25 per nucleotide).

    Args:
        sequence: DNA sequence to score. Must have length >= PWM length.
        pwm: Position weight matrix as list of dicts mapping nucleotides
            to probabilities.

    Returns:
        Log-likelihood ratio score (higher is better match). Returns
        negative infinity for sequences shorter than the PWM.

    Raises:
        ValueError: If sequence is shorter than the PWM.
    """
    if len(sequence) < len(pwm):
        raise ValueError(f"Sequence length ({len(sequence)}) must be >= PWM length ({len(pwm)})")

    background = 0.25
    score = 0.0

    for pos, pos_freqs in enumerate(pwm):
        nt = sequence[pos].upper()
        freq = pos_freqs.get(nt, 0.001)  # Small value for unknown nucleotides
        score += math.log2(freq / background)

    return round(score, 6)


def find_tf_binding_motifs(
    sequences: list[str],
    method: str = "overrepresentation",
    k_min: int = 6,
    k_max: int = 10,
    top_n: int = 20,
    background_sequences: list[str] | None = None,
) -> list[dict]:
    """Find enriched sequence motifs in regulatory regions.

    Identifies overrepresented k-mers in the input sequences compared to
    expected frequencies (uniform or from background sequences), then
    clusters similar k-mers and builds position weight matrices for the
    top motifs.

    Args:
        sequences: List of DNA sequences to search for motifs (e.g.,
            promoter regions, enhancer sequences).
        method: Discovery method. Currently supports 'overrepresentation'.
        k_min: Minimum k-mer size to search.
        k_max: Maximum k-mer size to search.
        top_n: Number of top motifs to return.
        background_sequences: Optional background sequences for computing
            expected frequencies. If None, uniform distribution is assumed.

    Returns:
        List of dicts, each containing:
            - motif: The consensus k-mer sequence.
            - consensus: Consensus sequence derived from the PWM.
            - p_value: Statistical significance of overrepresentation.
            - n_occurrences: Number of times the motif appears.
            - pwm: Position weight matrix for the motif.

    Raises:
        ValueError: If sequences are empty or method is unsupported.
    """
    if not sequences:
        raise ValueError("sequences must not be empty")

    if method != "overrepresentation":
        raise ValueError(f"Unsupported method: {method}. Use 'overrepresentation'.")

    logger.info(
        "Finding TF binding motifs: %d sequences, k=%d-%d",
        len(sequences),
        k_min,
        k_max,
    )

    # Count k-mers in foreground
    fg_counts: Counter[str] = Counter()
    total_fg_kmers = 0

    for seq in sequences:
        seq_upper = seq.upper()
        for k in range(k_min, k_max + 1):
            for i in range(len(seq_upper) - k + 1):
                kmer = seq_upper[i : i + k]
                if all(c in "ACGT" for c in kmer):
                    fg_counts[kmer] += 1
                    total_fg_kmers += 1

    # Count k-mers in background (or use uniform expectation)
    if background_sequences:
        bg_counts: Counter[str] = Counter()
        total_bg_kmers = 0
        for seq in background_sequences:
            seq_upper = seq.upper()
            for k in range(k_min, k_max + 1):
                for i in range(len(seq_upper) - k + 1):
                    kmer = seq_upper[i : i + k]
                    if all(c in "ACGT" for c in kmer):
                        bg_counts[kmer] += 1
                        total_bg_kmers += 1
    else:
        bg_counts = None  # type: ignore[assignment]
        total_bg_kmers = 0

    # Score each k-mer by overrepresentation
    scored_kmers: list[tuple[str, float, int, float]] = []

    for kmer, count in fg_counts.items():
        k = len(kmer)
        observed_freq = count / total_fg_kmers

        if bg_counts is not None and total_bg_kmers > 0:
            bg_count = bg_counts.get(kmer, 0)
            expected_freq = (bg_count + 1) / (total_bg_kmers + 4**k)
        else:
            expected_freq = 1.0 / (4**k)

        # Log-odds score
        if expected_freq > 0:
            log_odds = math.log2(observed_freq / expected_freq)
        else:
            log_odds = 0.0

        # Binomial test p-value approximation
        p_value = _binomial_pvalue(count, total_fg_kmers, expected_freq)

        scored_kmers.append((kmer, log_odds, count, p_value))

    # Sort by p-value ascending, then by count descending
    scored_kmers.sort(key=lambda x: (x[3], -x[2]))

    # Cluster similar k-mers and take top representatives
    selected_kmers = _cluster_kmers(scored_kmers, top_n=top_n)

    # Build results with PWMs
    results: list[dict] = []
    for kmer, log_odds, count, p_value in selected_kmers:
        # Collect all occurrences for PWM building
        instances = _collect_motif_instances(sequences, kmer)
        if instances:
            pwm = build_pwm(instances, pseudocount=0.01)
            consensus = _consensus_from_pwm(pwm)
        else:
            pwm = build_pwm([kmer], pseudocount=0.01)
            consensus = kmer

        results.append(
            {
                "motif": kmer,
                "consensus": consensus,
                "p_value": round(p_value, 10),
                "n_occurrences": count,
                "pwm": pwm,
            }
        )

    logger.info("Found %d enriched motifs", len(results))

    return results


def _binomial_pvalue(k: int, n: int, p: float) -> float:
    """Compute approximate binomial p-value using normal approximation.

    Args:
        k: Number of successes.
        n: Number of trials.
        p: Expected probability.

    Returns:
        Approximate one-sided p-value.
    """
    if n == 0 or p <= 0 or p >= 1:
        return 1.0

    mean = n * p
    variance = n * p * (1 - p)

    if variance < 1e-15:
        return 0.0 if k > mean else 1.0

    z = (k - mean) / math.sqrt(variance)
    # One-sided p-value (upper tail)
    p_val = 0.5 * math.erfc(z / math.sqrt(2.0))
    return max(0.0, min(1.0, p_val))


def _cluster_kmers(
    scored_kmers: list[tuple[str, float, int, float]],
    top_n: int,
    similarity_threshold: int = 2,
) -> list[tuple[str, float, int, float]]:
    """Cluster similar k-mers and select representatives.

    Groups k-mers that differ by at most similarity_threshold positions,
    keeping the most significant representative from each cluster.

    Args:
        scored_kmers: List of (kmer, score, count, p_value) tuples.
        top_n: Maximum number of representatives to return.
        similarity_threshold: Maximum Hamming distance for clustering.

    Returns:
        List of representative (kmer, score, count, p_value) tuples.
    """
    selected: list[tuple[str, float, int, float]] = []
    used: set[str] = set()

    for entry in scored_kmers:
        kmer = entry[0]
        if kmer in used:
            continue

        # Check if this kmer is too similar to an already selected one
        too_similar = False
        for prev_kmer, _, _, _ in selected:
            if len(kmer) == len(prev_kmer):
                hamming = sum(1 for a, b in zip(kmer, prev_kmer) if a != b)
                if hamming <= similarity_threshold:
                    too_similar = True
                    break
            # Also check reverse complement
            rc = _reverse_complement(kmer)
            if len(rc) == len(prev_kmer):
                hamming_rc = sum(1 for a, b in zip(rc, prev_kmer) if a != b)
                if hamming_rc <= similarity_threshold:
                    too_similar = True
                    break

        if not too_similar:
            selected.append(entry)
            used.add(kmer)
            used.add(_reverse_complement(kmer))

        if len(selected) >= top_n:
            break

    return selected


def _collect_motif_instances(
    sequences: list[str],
    kmer: str,
    max_mismatch: int = 1,
) -> list[str]:
    """Collect all instances of a motif from sequences allowing mismatches.

    Args:
        sequences: List of DNA sequences to search.
        kmer: The motif to search for.
        max_mismatch: Maximum allowed mismatches.

    Returns:
        List of matched subsequences.
    """
    k = len(kmer)
    instances: list[str] = []
    kmer_upper = kmer.upper()
    rc_kmer = _reverse_complement(kmer_upper)

    for seq in sequences:
        seq_upper = seq.upper()
        for i in range(len(seq_upper) - k + 1):
            subseq = seq_upper[i : i + k]
            if not all(c in "ACGT" for c in subseq):
                continue

            # Check forward match
            mismatches = sum(1 for a, b in zip(subseq, kmer_upper) if a != b)
            if mismatches <= max_mismatch:
                instances.append(subseq)
                continue

            # Check reverse complement match
            mismatches_rc = sum(1 for a, b in zip(subseq, rc_kmer) if a != b)
            if mismatches_rc <= max_mismatch:
                instances.append(subseq)

    return instances


def _consensus_from_pwm(pwm: list[dict[str, float]]) -> str:
    """Derive consensus sequence from a PWM.

    Args:
        pwm: Position weight matrix.

    Returns:
        Consensus sequence string.
    """
    consensus = []
    for pos in pwm:
        best_nt = max(NUCLEOTIDES, key=lambda nt: pos.get(nt, 0.0))
        consensus.append(best_nt)
    return "".join(consensus)


def scan_sequence_for_motifs(
    sequence: str,
    motif_library: dict[str, list[dict[str, float]]],
    threshold: float = 0.8,
    scan_reverse: bool = True,
) -> list[dict]:
    """Scan a sequence for matches to known motifs.

    Slides each motif PWM across the sequence and reports positions where
    the normalized score exceeds the threshold.

    Args:
        sequence: DNA sequence to scan.
        motif_library: Dictionary mapping motif IDs to their PWMs (each PWM
            is a list of dicts mapping nucleotides to probabilities).
        threshold: Minimum normalized score (0-1) to report a match.
            The score is normalized by the maximum possible score for
            each PWM.
        scan_reverse: Whether to also scan the reverse complement strand.

    Returns:
        List of dicts, each containing:
            - motif_id: Identifier of the matched motif.
            - position: 0-based start position of the match.
            - strand: '+' or '-' for forward/reverse complement.
            - score: Normalized match score (0-1).

    Raises:
        ValueError: If sequence is empty or motif_library is empty.
    """
    if not sequence:
        raise ValueError("sequence must not be empty")
    if not motif_library:
        raise ValueError("motif_library must not be empty")

    logger.info(
        "Scanning sequence (length %d) against %d motifs",
        len(sequence),
        len(motif_library),
    )

    seq_upper = sequence.upper()
    rc_seq = _reverse_complement(seq_upper) if scan_reverse else ""

    matches: list[dict] = []

    for motif_id, pwm in motif_library.items():
        motif_len = len(pwm)

        # Compute max possible score for normalization
        max_score = 0.0
        min_score = 0.0
        for pos_freqs in pwm:
            max_freq = max(pos_freqs.get(nt, 0.001) for nt in NUCLEOTIDES)
            min_freq = min(pos_freqs.get(nt, 0.001) for nt in NUCLEOTIDES)
            max_score += math.log2(max_freq / 0.25)
            min_score += math.log2(min_freq / 0.25)

        score_range = max_score - min_score
        if score_range < 1e-10:
            continue

        # Scan forward strand
        for i in range(len(seq_upper) - motif_len + 1):
            subseq = seq_upper[i : i + motif_len]
            if not all(c in "ACGT" for c in subseq):
                continue

            raw_score = 0.0
            for pos, pos_freqs in enumerate(pwm):
                nt = subseq[pos]
                freq = pos_freqs.get(nt, 0.001)
                raw_score += math.log2(freq / 0.25)

            normalized = (raw_score - min_score) / score_range
            if normalized >= threshold:
                matches.append(
                    {
                        "motif_id": motif_id,
                        "position": i,
                        "strand": "+",
                        "score": round(normalized, 6),
                    }
                )

        # Scan reverse strand
        if scan_reverse and rc_seq:
            for i in range(len(rc_seq) - motif_len + 1):
                subseq = rc_seq[i : i + motif_len]
                if not all(c in "ACGT" for c in subseq):
                    continue

                raw_score = 0.0
                for pos, pos_freqs in enumerate(pwm):
                    nt = subseq[pos]
                    freq = pos_freqs.get(nt, 0.001)
                    raw_score += math.log2(freq / 0.25)

                normalized = (raw_score - min_score) / score_range
                if normalized >= threshold:
                    # Convert position to forward strand coordinates
                    fwd_pos = len(seq_upper) - i - motif_len
                    matches.append(
                        {
                            "motif_id": motif_id,
                            "position": fwd_pos,
                            "strand": "-",
                            "score": round(normalized, 6),
                        }
                    )

    # Sort by score descending
    matches.sort(key=lambda m: m["score"], reverse=True)

    logger.info("Found %d motif matches above threshold %.2f", len(matches), threshold)

    return matches
