"""ASV denoising for amplicon metagenomics.

Implements DADA2-style amplicon sequence variant (ASV) denoising, including
error model estimation from quality scores, error-correction-based denoising,
and paired-end read merging. ASVs provide single-nucleotide resolution compared
to the 97% identity clusters of traditional OTU approaches.

The denoising algorithm:
1. Builds an error model: a 4x4 base-substitution matrix (empirical when
   sequences are supplied, otherwise Phred-derived) plus the Phred definition
   ``P(error) = 10 ** (-Q / 10)``.
2. Groups sequences by abundance.
3. For each unique sequence, calculates the probability that it arose from
   a more abundant sequence through sequencing errors. When per-read Phred
   quality scores are supplied, the error probability at each position comes from
   the (read-averaged) quality of the candidate variant at that position, as in
   DADA2; otherwise the model's position-independent substitution rates are used.
4. Merges sequences whose abundance is consistent with errors from a more
   abundant sequence: the abundance p-value is the Poisson upper tail
   ``P(X >= abundance)`` with ``lambda = parent_abundance * P(variant | parent)``.
"""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Any, cast

from metainformant.core.sequence import reverse_complement
from metainformant.core.utils.logging import get_logger
from metainformant.metagenomics._special import poisson_upper_tail

logger = get_logger(__name__)

# SciPy is optional for this package: without it the Poisson tail uses the exact
# pure-Python implementation in ``metainformant.metagenomics._special``.
try:
    from scipy import stats as scipy_stats

    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    scipy_stats = cast(Any, None)  # optional dependency; callers guard on HAS_SCIPY

_ENV_PREFIX = "META_"

# A uniformly random base call is wrong 75% of the time, so no base call is
# modelled as less reliable than that (Phred 0 would otherwise imply P(error) = 1).
_MAX_BASE_ERROR = 0.75
_MIN_BASE_ERROR = 1e-12
_BASE_INDEX = {"A": 0, "C": 1, "G": 2, "T": 3}

# Phred 30 (P(error) = 0.001): assumed when no quality information is available
_DEFAULT_PHRED = 30
_DEFAULT_ERROR_RATE = 0.001


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
        quality: Mean Phred quality over all positions of the reads assigned to this
            ASV (including reads merged into it). 0.0 when no quality scores were
            supplied to ``denoise_sequences``.
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


def _base_error_probability(phred: float) -> float:
    """Convert a Phred score to a base-call error probability, ``10 ** (-Q / 10)``.

    The result is bounded to ``[1e-12, 0.75]``: a uniformly random base call is
    wrong 75% of the time, so Phred scores below about 1.25 (which nominally imply
    ``P(error) >= 0.75``) are treated as random calls, and the lower bound keeps
    ``log(1 - p)`` and ``log(p)`` finite.

    Args:
        phred: Phred quality score. Negative or NaN scores are treated as 0.

    Returns:
        Error probability in ``[1e-12, 0.75]``.
    """
    q = float(phred)
    if not q >= 0.0:  # also catches NaN
        q = 0.0
    return min(max(math.pow(10.0, -q / 10.0), _MIN_BASE_ERROR), _MAX_BASE_ERROR)


def _quality_profile(read_quals: Sequence[Sequence[float]], length: int) -> tuple[list[float], float, int]:
    """Summarise the Phred scores of the reads that share one unique sequence.

    Args:
        read_quals: Per-position Phred scores, one list per read.
        length: Length of the unique sequence; scores beyond it are ignored.

    Returns:
        Tuple ``(error_profile, phred_sum, phred_count)``. ``error_profile[pos]`` is the
        mean error probability over the reads at that position (NaN where no read
        has a score); ``phred_sum`` and ``phred_count`` give the total Phred score and
        the number of scores, for averaging quality over several sequences.
    """
    err_sum = [0.0] * length
    err_n = [0] * length
    phred_sum = 0.0
    phred_count = 0
    for quals in read_quals:
        for pos in range(min(length, len(quals))):
            q = float(quals[pos])
            if not q >= 0.0:
                q = 0.0
            err_sum[pos] += _base_error_probability(q)
            err_n[pos] += 1
            phred_sum += q
            phred_count += 1
    profile = [err_sum[pos] / err_n[pos] if err_n[pos] else float("nan") for pos in range(length)]
    return profile, phred_sum, phred_count


def estimate_error_rates(
    quality_scores: Sequence[Sequence[float]],
    sequences: Sequence[str] | None = None,
) -> ErrorModel:
    """Learn a sequencing error model from quality score data.

    The per-base error probability is derived from the Phred scores
    (``P_error = 10 ** (-Q / 10)``) and averaged as probabilities, which is larger
    than the probability at the mean Q whenever the scores vary. When sequences are
    provided, the relative frequencies of the substitutions observed against the
    positional consensus decide *which* base replaces the true one; otherwise the
    three alternatives are equally likely. Each row of ``transition_rates`` is
    therefore ``[1 - error, error * share_j, ...]`` and sums to 1.

    Args:
        quality_scores: List of quality score arrays, one per read.
            Each inner list contains integer Phred scores per position.
        sequences: Optional list of nucleotide sequences corresponding to
            the quality scores. Used to estimate base transition rates.

    Returns:
        ErrorModel with transition rates and quality-based error rates.

    Raises:
        ValueError: If quality_scores is empty, or ``sequences`` is given (and not empty)
            with a different length than ``quality_scores``.

    Examples:
        >>> scores = [[30, 30, 25, 20], [28, 30, 22, 18]]
        >>> model = estimate_error_rates(scores)
        >>> model.quality_error_rates[30]  # ~0.001
        0.001
    """
    if len(quality_scores) == 0:
        raise ValueError("Quality scores list must not be empty")

    if sequences is not None and len(sequences) > 0 and len(sequences) != len(quality_scores):
        raise ValueError(
            f"sequences ({len(sequences)}) and quality_scores ({len(quality_scores)}) must have one entry per read"
        )

    logger.info(f"Estimating error rates from {len(quality_scores)} reads")

    quality_error_rates: dict[int, float] = {}
    for q in range(0, 42):  # Phred scores typically 0-41
        # P_error = 10^(-Q/10)
        quality_error_rates[q] = 10.0 ** (-q / 10.0)

    transition_rates: list[list[float]] = [[0.0] * 4 for _ in range(4)]

    if sequences is not None and len(sequences) > 0:
        # Empirical substitution shares from comparison with the positional consensus
        max_len = max(len(s) for s in sequences)
        position_bases: list[Counter[str]] = [Counter() for _ in range(max_len)]

        for seq in sequences:
            for pos, base in enumerate(seq.upper()):
                if base in _BASE_INDEX:
                    position_bases[pos][base] += 1

        consensus = []
        for pos_counter in position_bases:
            if pos_counter:
                consensus.append(pos_counter.most_common(1)[0][0])
            else:
                consensus.append("N")

        # Per consensus base: Phred-derived error probability and observed substitution counts
        mismatch_counts = [[0.0] * 4 for _ in range(4)]
        error_sum = [0.0] * 4
        observations = [0] * 4
        for read_idx, seq in enumerate(sequences):
            read_quals = quality_scores[read_idx]
            for pos, obs_base in enumerate(seq.upper()):
                if pos >= len(consensus) or obs_base not in _BASE_INDEX:
                    continue
                true_base = consensus[pos]
                if true_base not in _BASE_INDEX:
                    continue
                i = _BASE_INDEX[true_base]
                j = _BASE_INDEX[obs_base]
                phred = read_quals[pos] if pos < len(read_quals) else _DEFAULT_PHRED
                error_sum[i] += _base_error_probability(phred)
                observations[i] += 1
                if i != j:
                    mismatch_counts[i][j] += 1.0

        for i in range(4):
            if observations[i] == 0:
                # Consensus base never seen: fall back to the default Phred-30 row
                transition_rates[i][i] = 1.0 - _DEFAULT_ERROR_RATE
                for j in range(4):
                    if i != j:
                        transition_rates[i][j] = _DEFAULT_ERROR_RATE / 3.0
                continue
            error = error_sum[i] / observations[i]
            row_mismatches = sum(mismatch_counts[i])
            transition_rates[i][i] = 1.0 - error
            for j in range(4):
                if i != j:
                    share = mismatch_counts[i][j] / row_mismatches if row_mismatches > 0 else 1.0 / 3.0
                    transition_rates[i][j] = error * share
    else:
        # No sequences: Phred-based uniform error model using the mean error probability
        flat_errors = [_base_error_probability(q) for scores in quality_scores for q in scores]
        error_rate = sum(flat_errors) / len(flat_errors) if flat_errors else _DEFAULT_ERROR_RATE
        for i in range(4):
            transition_rates[i][i] = 1.0 - error_rate
            for j in range(4):
                if i != j:
                    transition_rates[i][j] = error_rate / 3.0

    positions_modeled = max(len(s) for s in quality_scores)

    model = ErrorModel(
        transition_rates=transition_rates,
        quality_error_rates=quality_error_rates,
        positions_modeled=positions_modeled,
    )

    logger.info(f"Error model built from {positions_modeled} positions, {len(quality_scores)} reads")
    return model


def _substitution_share(row: Sequence[float], true_idx: int, observed_idx: int) -> float:
    """Probability that an error at a base replaces ``true_idx`` with ``observed_idx``.

    Args:
        row: Transition-rate row of the true base.
        true_idx: Index of the true base.
        observed_idx: Index of the observed (erroneous) base.

    Returns:
        ``row[observed_idx]`` renormalised over the three possible substitutions;
        uniform (1/3) when the row has no substitution mass.
    """
    off_diagonal = sum(row[k] for k in range(4) if k != true_idx)
    if off_diagonal <= 0:
        return 1.0 / 3.0
    return row[observed_idx] / off_diagonal


def _calculate_error_probability(
    seq_a: str,
    seq_b: str,
    abundance_a: int,
    abundance_b: int,
    error_model: ErrorModel,
    position_error: Sequence[float] | None = None,
) -> float:
    """Calculate the probability that seq_b arose from seq_a through errors.

    Uses the error model to compute:
        P(seq_b | seq_a) = product over positions of the per-position probability

    Without ``position_error`` each position contributes the model's
    ``transition_rates[a][b]``. With it, a matching position contributes
    ``1 - error`` and a mismatch ``error * share(a -> b)``, where ``error`` is the
    Phred-derived error probability at that position and ``share`` is the model's
    substitution distribution given an error. NaN entries (or positions beyond the
    profile) fall back to the model's transition rate.

    This is the core of the DADA2-style denoising: if P(B|A) * abundance_A
    is high enough, B is likely an error variant of A.

    Args:
        seq_a: The more abundant (candidate true) sequence.
        seq_b: The less abundant (candidate error) sequence.
        abundance_a: Abundance of sequence A.
        abundance_b: Abundance of sequence B.
        error_model: Error model with transition rates.
        position_error: Optional per-position error probabilities of sequence B
            derived from its reads' Phred scores.

    Returns:
        Log-probability that B is an error variant of A.
    """
    log_prob = 0.0

    min_len = min(len(seq_a), len(seq_b))
    for pos in range(min_len):
        i = _BASE_INDEX.get(seq_a[pos].upper())
        j = _BASE_INDEX.get(seq_b[pos].upper())

        if i is None or j is None:
            continue

        error = position_error[pos] if position_error is not None and pos < len(position_error) else math.nan
        if math.isnan(error):
            rate = error_model.transition_rates[i][j]
        elif i == j:
            rate = 1.0 - error
        else:
            rate = error * _substitution_share(error_model.transition_rates[i], i, j)
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
    5. If the abundance p-value, the Poisson upper tail ``P(X >= abundance)`` with
       ``lambda = parent_abundance * P(variant | parent)``, exceeds omega_a, the
       sequence is merged into the parent ASV.
    6. Otherwise, the sequence becomes a new ASV.

    Args:
        sequences: Dictionary mapping read IDs to nucleotide sequences.
        error_rates: Pre-computed error model. If None, estimated from
            quality_scores or defaults are used.
        quality_scores: Optional Phred scores per read ID (one list per read, one score
            per position). They set the error probability of every mismatch when a
            variant is tested against a parent (scores are averaged as error
            probabilities over the reads sharing a sequence), feed the error-model
            estimate when ``error_rates`` is None, and give each ASV its mean quality.
            Reads without scores fall back to the error model's substitution rates.
        abundance: Pre-computed abundance per unique sequence. If None, computed
            from duplicate sequences in the input.
        omega_a: P-value threshold for the abundance test. Lower values are
            more conservative (fewer sequences merged). Default 1e-40.
        min_abundance: Minimum abundance to retain a sequence. Default 1.
        band_size: Maximum length difference between a variant and a parent for the
            two to be compared (no alignment is performed).

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

    # Dereplicate: count unique sequences, remembering which reads each one represents
    members: dict[str, list[str]]
    if abundance:
        unique_seqs = {seq_id: seq.upper() for seq_id, seq in sequences.items()}
        abund = abundance
        members = {seq_id: [seq_id] for seq_id in unique_seqs}
    else:
        seq_counts: Counter[str] = Counter()
        seq_to_id: dict[str, str] = {}
        seq_reads: dict[str, list[str]] = {}
        for seq_id, seq in sequences.items():
            upper_seq = seq.upper().replace("-", "")
            seq_counts[upper_seq] += 1
            seq_reads.setdefault(upper_seq, []).append(seq_id)
            if upper_seq not in seq_to_id:
                seq_to_id[upper_seq] = seq_id

        unique_seqs = {seq_to_id[seq]: seq for seq in seq_counts}
        abund = {seq_to_id[seq]: count for seq, count in seq_counts.items()}
        members = {seq_to_id[seq]: reads for seq, reads in seq_reads.items()}

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

    # Per-unique-sequence quality summaries (mean error probability per position, Phred totals)
    position_errors: dict[str, list[float]] = {}
    phred_totals: dict[str, tuple[float, int]] = {}
    if quality_scores:
        for sid, seq in unique_seqs.items():
            read_quals = [quality_scores[rid] for rid in members.get(sid, [sid]) if rid in quality_scores]
            if read_quals:
                profile, phred_sum, phred_count = _quality_profile(read_quals, len(seq))
                position_errors[sid] = profile
                phred_totals[sid] = (phred_sum, phred_count)

    # Sort by abundance descending
    sorted_ids = sorted(unique_seqs.keys(), key=lambda x: -abund.get(x, 0))

    # Denoising: greedy assignment
    asvs: list[ASV] = []
    asv_sequences: list[str] = []
    asv_abundances: list[int] = []
    asv_phred: list[list[float]] = []  # [sum of Phred scores, number of scores] per ASV
    reads_merged = 0

    for seq_id in sorted_ids:
        seq = unique_seqs[seq_id]
        seq_abund = abund.get(seq_id, 1)
        seq_phred = phred_totals.get(seq_id, (0.0, 0))

        merged = False
        for asv_idx, asv_seq in enumerate(asv_sequences):
            if abs(len(seq) - len(asv_seq)) > band_size:
                continue

            # Calculate error probability
            log_p = _calculate_error_probability(
                asv_seq,
                seq,
                asv_abundances[asv_idx],
                seq_abund,
                error_rates,
                position_error=position_errors.get(seq_id),
            )

            # Abundance p-value test (Poisson-based)
            # Expected number of reads of seq_b from seq_a:
            #   lambda = abundance_a * P(b|a)
            # The sequence is consistent with errors when observing at least
            # ``seq_abund`` reads is not improbable under Poisson(lambda).
            expected = asv_abundances[asv_idx] * math.exp(log_p) if log_p > -500 else 0.0
            if expected > 0:
                # Upper-tail probability P(X >= observed) where X ~ Poisson(expected)
                p_value = _poisson_pvalue(seq_abund, expected)
                if p_value > omega_a:
                    # Merge into this ASV
                    asvs[asv_idx].abundance += seq_abund
                    asv_abundances[asv_idx] += seq_abund
                    asv_phred[asv_idx][0] += seq_phred[0]
                    asv_phred[asv_idx][1] += seq_phred[1]
                    reads_merged += seq_abund
                    merged = True
                    break

        if not merged:
            new_asv = ASV(sequence=seq, abundance=seq_abund)
            asvs.append(new_asv)
            asv_sequences.append(seq)
            asv_abundances.append(seq_abund)
            asv_phred.append([seq_phred[0], float(seq_phred[1])])

    # Mean Phred quality of the reads assigned to each ASV (0.0 when no scores were supplied)
    for asv, (asv_phred_sum, asv_phred_count) in zip(asvs, asv_phred):
        if asv_phred_count > 0:
            asv.quality = asv_phred_sum / asv_phred_count

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
    """Compute the upper-tail Poisson p-value ``P(X >= observed)`` for ``X ~ Poisson(expected)``.

    Evaluated as the survival function ``scipy.stats.poisson.sf(observed - 1, expected)``
    (or the equivalent exact pure-Python tail when SciPy is not installed). The tail
    must not be computed as ``1 - cdf``: for the small ``lambda`` and large counts
    that separate real variants from errors the true tail (1e-90 and below) is far
    smaller than the rounding error of the cdf (about 1e-16), so ``1 - cdf`` returns 0,
    a rounding residue, or a negative number, and the residue can exceed ``omega_a``.

    Args:
        observed: Observed count (the abundance of the candidate variant).
        expected: Expected count (Poisson lambda).

    Returns:
        P(X >= observed); 1.0 for ``observed <= 0`` and 0.0 for ``expected <= 0``
        with ``observed > 0``.
    """
    if observed <= 0:
        return 1.0
    if expected <= 0:
        return 0.0
    if HAS_SCIPY:
        return float(scipy_stats.poisson.sf(observed - 1, expected))
    return poisson_upper_tail(observed, expected)


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
    # Input order, so the result does not depend on the per-process string hash seed
    common_ids = [read_id for read_id in forward if read_id in reverse]
    if not common_ids:
        raise ValueError("No common read IDs between forward and reverse read sets")

    logger.info(f"Merging {len(common_ids)} paired-end read pairs (min_overlap={min_overlap})")

    merged: dict[str, str] = {}
    merge_count = 0
    fail_count = 0

    for read_id in common_ids:
        fwd_seq = forward[read_id].upper()
        rev_seq = reverse_complement(reverse[read_id].upper())

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


__all__ = [
    "ASV",
    "DenoisingResult",
    "ErrorModel",
    "denoise_sequences",
    "estimate_error_rates",
    "merge_paired_reads",
]
