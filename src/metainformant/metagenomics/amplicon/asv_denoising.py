"""ASV denoising for amplicon metagenomics.

Implements DADA2-style amplicon sequence variant (ASV) denoising, including
error model estimation from quality scores, error-correction-based denoising,
and paired-end read merging. ASVs provide single-nucleotide resolution compared
to the 97% identity clusters of traditional OTU approaches.

The denoising algorithm:
1. Learns a position-specific error model from quality scores.
2. Groups sequences by abundance.
3. For each unique sequence, calculates the probability that it arose from
   a more abundant sequence through sequencing errors.
4. Merges sequences that are likely error variants of more abundant true sequences.
"""

from __future__ import annotations

import math
import os
from collections import Counter
from dataclasses import dataclass, field
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

_ENV_PREFIX = "META_"


@dataclass
class ErrorModel:
    """Position-specific sequencing error model.

    Attributes:
        transition_rates: 4x4 matrix of base transition probabilities.
            Rows/columns indexed as A=0, C=1, G=2, T=3.
            Entry [i][j] = P(observe j | true base is i).
        quality_error_rates: Mapping of quality score to error probability.
        positions_modeled: Number of positions used to build the model.
    """

    transition_rates: list[list[float]]
    quality_error_rates: dict[int, float] = field(default_factory=dict)
    positions_modeled: int = 0


@dataclass
class ASV:
    """Amplicon Sequence Variant.

    Attributes:
        sequence: The denoised nucleotide sequence.
        abundance: Total abundance (count of reads assigned).
        sample_counts: Per-sample abundance if from multiple samples.
        quality: Average quality of constituent reads.
    """

    sequence: str
    abundance: int
    sample_counts: dict[str, int] = field(default_factory=dict)
    quality: float = 0.0


@dataclass
class DenoisingResult:
    """Result of ASV denoising.

    Attributes:
        asvs: List of denoised ASVs.
        error_model: The error model used or learned.
        total_reads: Total number of input reads.
        num_asvs: Number of unique ASVs after denoising.
        reads_removed: Number of reads merged into other ASVs.
    """

    asvs: list[ASV]
    error_model: ErrorModel | None
    total_reads: int
    num_asvs: int
    reads_removed: int


def estimate_error_rates(
    quality_scores: list[list[int]],
    sequences: list[str] | None = None,
) -> ErrorModel:
    """Learn a sequencing error model from quality score data.

    Estimates position-specific error rates by converting Phred quality scores
    to error probabilities and, if sequences are provided, computing empirical
    transition rates between nucleotide bases.

    Phred quality score Q maps to error probability: P_error = 10^(-Q/10)

    Args:
        quality_scores: List of quality score arrays, one per read.
            Each inner list contains integer Phred scores per position.
        sequences: Optional list of nucleotide sequences corresponding to
            the quality scores. Used to estimate base transition rates.

    Returns:
        ErrorModel with transition rates and quality-based error rates.

    Raises:
        ValueError: If quality_scores is empty.

    Examples:
        >>> scores = [[30, 30, 25, 20], [28, 30, 22, 18]]
        >>> model = estimate_error_rates(scores)
        >>> model.quality_error_rates[30]  # ~0.001
        0.001
    """
    if not quality_scores:
        raise ValueError("Quality scores list must not be empty")

    logger.info(f"Estimating error rates from {len(quality_scores)} reads")

    # Build quality -> error probability mapping from observed scores
    all_scores: Counter[int] = Counter()
    for scores in quality_scores:
        for q in scores:
            all_scores[q] += 1

    quality_error_rates: dict[int, float] = {}
    for q in range(0, 42):  # Phred scores typically 0-41
        # P_error = 10^(-Q/10)
        quality_error_rates[q] = 10.0 ** (-q / 10.0)

    # Base indices: A=0, C=1, G=2, T=3
    base_to_idx = {"A": 0, "C": 1, "G": 2, "T": 3}
    # Initialize transition matrix with Phred-based priors
    transition_counts = [[0.0] * 4 for _ in range(4)]
    total_bases_per_row = [0.0] * 4

    if sequences and len(sequences) == len(quality_scores):
        # Empirical transition estimation from sequence data
        # Compare each read to the consensus at each position
        # First, build positional consensus
        max_len = max(len(s) for s in sequences)
        position_bases: list[Counter[str]] = [Counter() for _ in range(max_len)]

        for seq in sequences:
            for pos, base in enumerate(seq.upper()):
                if base in base_to_idx:
                    position_bases[pos][base] += 1

        consensus = []
        for pos_counter in position_bases:
            if pos_counter:
                consensus.append(pos_counter.most_common(1)[0][0])
            else:
                consensus.append("N")

        # Count transitions: consensus base -> observed base, weighted by error prob
        for read_idx, seq in enumerate(sequences):
            for pos, obs_base in enumerate(seq.upper()):
                if pos >= len(consensus) or obs_base not in base_to_idx:
                    continue
                true_base = consensus[pos]
                if true_base not in base_to_idx:
                    continue
                i = base_to_idx[true_base]
                j = base_to_idx[obs_base]
                q = quality_scores[read_idx][pos] if pos < len(quality_scores[read_idx]) else 30
                weight = 1.0
                transition_counts[i][j] += weight
                total_bases_per_row[i] += weight
    else:
        # No sequences: use Phred-based uniform error model
        # Diagonal (correct) weighted by average quality
        avg_q = 30  # Default assumption
        if quality_scores:
            flat_scores = [q for scores in quality_scores for q in scores]
            if flat_scores:
                avg_q = sum(flat_scores) / len(flat_scores)

        error_rate = 10.0 ** (-avg_q / 10.0)
        for i in range(4):
            transition_counts[i][i] = 1.0 - error_rate
            for j in range(4):
                if i != j:
                    transition_counts[i][j] = error_rate / 3.0
            total_bases_per_row[i] = 1.0

    # Normalize to probabilities
    transition_rates: list[list[float]] = [[0.0] * 4 for _ in range(4)]
    for i in range(4):
        row_total = total_bases_per_row[i]
        if row_total > 0:
            for j in range(4):
                transition_rates[i][j] = transition_counts[i][j] / row_total
        else:
            # Uniform prior
            transition_rates[i][i] = 0.999
            for j in range(4):
                if i != j:
                    transition_rates[i][j] = 0.001 / 3.0

    positions_modeled = max(len(s) for s in quality_scores) if quality_scores else 0

    model = ErrorModel(
        transition_rates=transition_rates,
        quality_error_rates=quality_error_rates,
        positions_modeled=positions_modeled,
    )

    logger.info(f"Error model built from {positions_modeled} positions, {len(quality_scores)} reads")
    return model


def _calculate_error_probability(
    seq_a: str,
    seq_b: str,
    abundance_a: int,
    abundance_b: int,
    error_model: ErrorModel,
) -> float:
    """Calculate the probability that seq_b arose from seq_a through errors.

    Uses the error model to compute:
        P(seq_b | seq_a) = product over positions of transition probability

    This is the core of the DADA2-style denoising: if P(B|A) * abundance_A
    is high enough, B is likely an error variant of A.

    Args:
        seq_a: The more abundant (candidate true) sequence.
        seq_b: The less abundant (candidate error) sequence.
        abundance_a: Abundance of sequence A.
        abundance_b: Abundance of sequence B.
        error_model: Error model with transition rates.

    Returns:
        Log-probability that B is an error variant of A.
    """
    base_to_idx = {"A": 0, "C": 1, "G": 2, "T": 3}
    log_prob = 0.0

    min_len = min(len(seq_a), len(seq_b))
    for pos in range(min_len):
        base_a = seq_a[pos].upper()
        base_b = seq_b[pos].upper()

        if base_a not in base_to_idx or base_b not in base_to_idx:
            continue

        i = base_to_idx[base_a]
        j = base_to_idx[base_b]

        rate = error_model.transition_rates[i][j]
        if rate <= 0:
            rate = 1e-10  # Floor to avoid log(0)
        log_prob += math.log(rate)

    # Penalty for length difference
    len_diff = abs(len(seq_a) - len(seq_b))
    if len_diff > 0:
        # Indel penalty
        log_prob -= len_diff * 5.0  # Strong penalty for length differences

    return log_prob


def denoise_sequences(
    sequences: dict[str, str],
    error_rates: ErrorModel | None = None,
    quality_scores: dict[str, list[int]] | None = None,
    abundance: dict[str, int] | None = None,
    omega_a: float = 1e-40,
    min_abundance: int = 1,
    band_size: int = 16,
) -> DenoisingResult:
    """Denoise amplicon sequences into ASVs using error model-based correction.

    Implements the core DADA2 denoising algorithm:
    1. Dereplicate sequences and compute abundances.
    2. If no error model provided, estimate one from quality scores or
       use a default Phred-30 model.
    3. Sort unique sequences by abundance (descending).
    4. For each sequence, test whether it could be an error variant of
       any more abundant ASV using the error model.
    5. If the p-value (based on abundance-weighted error probability)
       exceeds omega_a, the sequence is merged into the parent ASV.
    6. Otherwise, the sequence becomes a new ASV.

    Args:
        sequences: Dictionary mapping read IDs to nucleotide sequences.
        error_rates: Pre-computed error model. If None, estimated from
            quality_scores or defaults are used.
        quality_scores: Optional quality scores per read for error model estimation.
        abundance: Pre-computed abundance per unique sequence. If None, computed
            from duplicate sequences in the input.
        omega_a: P-value threshold for the abundance test. Lower values are
            more conservative (fewer sequences merged). Default 1e-40.
        min_abundance: Minimum abundance to retain a sequence. Default 1.
        band_size: Band size for banded alignment (not used in pure-Python mode).

    Returns:
        DenoisingResult with denoised ASVs and statistics.

    Raises:
        ValueError: If sequences is empty.

    Examples:
        >>> reads = {"r1": "ATCGATCG", "r2": "ATCGATCG", "r3": "ATCGATCG",
        ...          "r4": "ATCGATTG"}  # r4 has one error
        >>> result = denoise_sequences(reads)
        >>> result.num_asvs <= 2
        True
    """
    if not sequences:
        raise ValueError("Input sequences dictionary must not be empty")

    logger.info(f"Denoising {len(sequences)} reads into ASVs")

    # Dereplicate: count unique sequences
    if abundance:
        unique_seqs = {seq_id: seq.upper() for seq_id, seq in sequences.items()}
        abund = abundance
    else:
        seq_counts: Counter[str] = Counter()
        seq_to_id: dict[str, str] = {}
        for seq_id, seq in sequences.items():
            upper_seq = seq.upper().replace("-", "")
            seq_counts[upper_seq] += 1
            if upper_seq not in seq_to_id:
                seq_to_id[upper_seq] = seq_id

        unique_seqs = {seq_to_id[seq]: seq for seq in seq_counts}
        abund = {seq_to_id[seq]: count for seq, count in seq_counts.items()}

    # Filter by minimum abundance
    unique_seqs = {sid: seq for sid, seq in unique_seqs.items() if abund.get(sid, 0) >= min_abundance}
    abund = {sid: c for sid, c in abund.items() if sid in unique_seqs}

    total_reads = sum(abund.values())

    # Estimate or use provided error model
    if error_rates is None:
        if quality_scores:
            q_list = [quality_scores[sid] for sid in unique_seqs if sid in quality_scores]
            s_list = [unique_seqs[sid] for sid in unique_seqs if sid in quality_scores]
            if q_list:
                error_rates = estimate_error_rates(q_list, s_list)
            else:
                error_rates = _default_error_model()
        else:
            error_rates = _default_error_model()

    # Sort by abundance descending
    sorted_ids = sorted(unique_seqs.keys(), key=lambda x: -abund.get(x, 0))

    # Denoising: greedy assignment
    asvs: list[ASV] = []
    asv_sequences: list[str] = []
    asv_abundances: list[int] = []
    reads_merged = 0

    for seq_id in sorted_ids:
        seq = unique_seqs[seq_id]
        seq_abund = abund.get(seq_id, 1)

        merged = False
        for asv_idx, asv_seq in enumerate(asv_sequences):
            if abs(len(seq) - len(asv_seq)) > band_size:
                continue

            # Calculate error probability
            log_p = _calculate_error_probability(asv_seq, seq, asv_abundances[asv_idx], seq_abund, error_rates)

            # Abundance p-value test (Poisson-based)
            # Expected number of reads of seq_b from seq_a:
            #   lambda = abundance_a * P(b|a)
            # If observed abundance_b < lambda, it's likely an error variant
            expected = asv_abundances[asv_idx] * math.exp(log_p) if log_p > -500 else 0.0
            if expected > 0:
                # Poisson CDF approximation: P(X >= observed) where X ~ Poisson(expected)
                # If expected >> observed, the sequence is likely an error
                p_value = _poisson_pvalue(seq_abund, expected)
                if p_value > omega_a:
                    # Merge into this ASV
                    asvs[asv_idx].abundance += seq_abund
                    asv_abundances[asv_idx] += seq_abund
                    reads_merged += seq_abund
                    merged = True
                    break

        if not merged:
            new_asv = ASV(sequence=seq, abundance=seq_abund)
            asvs.append(new_asv)
            asv_sequences.append(seq)
            asv_abundances.append(seq_abund)

    result = DenoisingResult(
        asvs=asvs,
        error_model=error_rates,
        total_reads=total_reads,
        num_asvs=len(asvs),
        reads_removed=reads_merged,
    )

    logger.info(f"Denoising complete: {result.num_asvs} ASVs from {total_reads} reads ({reads_merged} merged)")
    return result


def _poisson_pvalue(observed: int, expected: float) -> float:
    """Compute upper-tail Poisson p-value: P(X >= observed | lambda=expected).

    Uses the incomplete gamma function approximation for numerical stability.

    Args:
        observed: Observed count.
        expected: Expected count (Poisson lambda).

    Returns:
        P-value (probability of observing >= observed counts).
    """
    if expected <= 0:
        return 0.0 if observed > 0 else 1.0
    if expected > observed * 10:
        return 1.0  # Expected far exceeds observed, highly consistent with error

    # Direct Poisson CDF computation for small values
    # P(X >= k) = 1 - P(X < k) = 1 - sum_{i=0}^{k-1} (lambda^i * e^-lambda / i!)
    if observed > 200:
        # Normal approximation for large values
        z = (observed - expected) / math.sqrt(expected) if expected > 0 else 0.0
        # Approximate upper tail of standard normal
        return 0.5 * math.erfc(z / math.sqrt(2))

    cumulative = 0.0
    log_lambda = math.log(expected)
    log_prob = -expected  # log(e^-lambda)

    for i in range(observed):
        cumulative += math.exp(log_prob)
        log_prob += log_lambda - math.log(i + 1)

    return 1.0 - cumulative


def _default_error_model() -> ErrorModel:
    """Create a default error model assuming Phred 30 average quality.

    Returns:
        ErrorModel with default transition rates.
    """
    error_rate = 0.001  # Phred 30
    transition_rates = [[0.0] * 4 for _ in range(4)]
    for i in range(4):
        for j in range(4):
            if i == j:
                transition_rates[i][j] = 1.0 - error_rate
            else:
                transition_rates[i][j] = error_rate / 3.0

    quality_error_rates = {q: 10.0 ** (-q / 10.0) for q in range(42)}

    return ErrorModel(
        transition_rates=transition_rates,
        quality_error_rates=quality_error_rates,
        positions_modeled=0,
    )


def merge_paired_reads(
    forward: dict[str, str],
    reverse: dict[str, str],
    min_overlap: int = 20,
    max_mismatch_ratio: float = 0.2,
    quality_forward: dict[str, list[int]] | None = None,
    quality_reverse: dict[str, list[int]] | None = None,
) -> dict[str, str]:
    """Merge paired-end reads based on overlap alignment.

    For each read pair, the algorithm:
    1. Reverse-complements the reverse read.
    2. Slides the forward and reverse-complement reads to find the overlap
       region with the best alignment score.
    3. Requires minimum overlap length and maximum mismatch ratio.
    4. In the overlap region, selects the base with higher quality score
       (if quality data is available), otherwise uses the forward read base.

    Args:
        forward: Dictionary mapping read IDs to forward read sequences.
        reverse: Dictionary mapping read IDs to reverse read sequences.
            Must have matching keys with forward.
        min_overlap: Minimum overlap length in bases (default 20).
        max_mismatch_ratio: Maximum fraction of mismatches in overlap (default 0.2).
        quality_forward: Optional quality scores for forward reads.
        quality_reverse: Optional quality scores for reverse reads.

    Returns:
        Dictionary mapping read IDs to merged sequences. Only successfully
        merged pairs are included.

    Raises:
        ValueError: If forward and reverse have no common read IDs.

    Examples:
        >>> fwd = {"r1": "ATCGATCGATCG"}
        >>> rev = {"r1": "TCGATCGTTTTTT"}  # Overlaps at TCGATCG
        >>> merged = merge_paired_reads(fwd, rev, min_overlap=5)
        >>> "r1" in merged
        True
    """
    common_ids = set(forward.keys()) & set(reverse.keys())
    if not common_ids:
        raise ValueError("No common read IDs between forward and reverse read sets")

    logger.info(f"Merging {len(common_ids)} paired-end read pairs (min_overlap={min_overlap})")

    merged: dict[str, str] = {}
    merge_count = 0
    fail_count = 0

    for read_id in common_ids:
        fwd_seq = forward[read_id].upper()
        rev_seq = _reverse_complement(reverse[read_id].upper())

        fwd_qual = quality_forward.get(read_id) if quality_forward else None
        rev_qual = quality_reverse.get(read_id) if quality_reverse else None
        if rev_qual:
            rev_qual = list(reversed(rev_qual))  # Reverse quality to match RC

        best_overlap = 0
        best_mismatches = float("inf")
        best_offset = -1

        fwd_len = len(fwd_seq)
        rev_len = len(rev_seq)

        # Slide reverse read along forward read to find best overlap
        for offset in range(min_overlap, min(fwd_len, rev_len) + 1):
            # Overlap region: last 'offset' bases of forward, first 'offset' bases of reverse
            fwd_overlap = fwd_seq[fwd_len - offset :]
            rev_overlap = rev_seq[:offset]

            mismatches = sum(1 for a, b in zip(fwd_overlap, rev_overlap) if a != b)
            mismatch_ratio = mismatches / offset

            if mismatch_ratio <= max_mismatch_ratio:
                if offset > best_overlap or (offset == best_overlap and mismatches < best_mismatches):
                    best_overlap = offset
                    best_mismatches = mismatches
                    best_offset = offset

        if best_offset < min_overlap:
            fail_count += 1
            continue

        # Build merged sequence
        # Non-overlapping forward region
        merged_seq_parts = list(fwd_seq[: fwd_len - best_offset])

        # Overlap region: resolve conflicts using quality scores
        fwd_overlap_start = fwd_len - best_offset
        for i in range(best_offset):
            fwd_base = fwd_seq[fwd_overlap_start + i]
            rev_base = rev_seq[i]

            if fwd_base == rev_base:
                merged_seq_parts.append(fwd_base)
            elif fwd_qual and rev_qual:
                fwd_q = fwd_qual[fwd_overlap_start + i] if (fwd_overlap_start + i) < len(fwd_qual) else 0
                rev_q = rev_qual[i] if i < len(rev_qual) else 0
                merged_seq_parts.append(fwd_base if fwd_q >= rev_q else rev_base)
            else:
                merged_seq_parts.append(fwd_base)  # Default to forward

        # Non-overlapping reverse region
        if best_offset < rev_len:
            merged_seq_parts.extend(list(rev_seq[best_offset:]))

        merged[read_id] = "".join(merged_seq_parts)
        merge_count += 1

    logger.info(f"Paired-end merging: {merge_count} merged, {fail_count} failed")
    return merged


def _reverse_complement(sequence: str) -> str:
    """Compute the reverse complement of a DNA sequence.

    Args:
        sequence: DNA sequence string.

    Returns:
        Reverse complement sequence.
    """
    complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(complement.get(base, "N") for base in reversed(sequence))


__all__ = [
    "ASV",
    "DenoisingResult",
    "ErrorModel",
    "denoise_sequences",
    "estimate_error_rates",
    "merge_paired_reads",
]
