"""K-mer analysis module for DNA sequences.

Provides comprehensive k-mer counting, frequency analysis, spectrum computation,
diversity metrics, sequence complexity profiling, homopolymer and microsatellite
detection, and pairwise k-mer profile comparison. All functions operate on raw
DNA sequence strings and return plain Python data structures.
"""

from __future__ import annotations

import math
from collections import Counter
from typing import Dict, List, Tuple

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


# ---------------------------------------------------------------------------
# Core k-mer counting
# ---------------------------------------------------------------------------


def count_kmers(seq: str, k: int) -> Dict[str, int]:
    """Count all k-mers in a DNA sequence.

    Slides a window of size *k* across the sequence and tallies every
    observed k-mer.  The sequence is upper-cased before counting.

    Args:
        seq: DNA sequence string.
        k: Length of each k-mer.  Must be a positive integer no larger
            than the length of *seq*.

    Returns:
        Dictionary mapping each observed k-mer to its count.

    Example:
        >>> count_kmers("ATCGATCG", 2)
        {'AT': 2, 'TC': 2, 'CG': 2, 'GA': 1}
    """
    if not seq or k <= 0:
        return {}

    seq_upper = seq.upper()
    n = len(seq_upper)
    if k > n:
        return {}

    counts: Counter[str] = Counter()
    for i in range(n - k + 1):
        kmer = seq_upper[i : i + k]
        counts[kmer] += 1

    return dict(counts)


# ---------------------------------------------------------------------------
# Frequencies and spectrum
# ---------------------------------------------------------------------------


def kmer_frequencies(seq: str, k: int) -> Dict[str, float]:
    """Compute normalised k-mer frequencies for a DNA sequence.

    Each k-mer count is divided by the total number of k-mers observed
    (i.e. ``len(seq) - k + 1``).

    Args:
        seq: DNA sequence string.
        k: Length of each k-mer.

    Returns:
        Dictionary mapping each observed k-mer to its frequency in [0, 1].

    Example:
        >>> kmer_frequencies("AAAA", 2)
        {'AA': 1.0}
    """
    counts = count_kmers(seq, k)
    if not counts:
        return {}

    total = sum(counts.values())
    return {kmer: count / total for kmer, count in counts.items()}


def kmer_spectrum(seq: str, k: int) -> Dict[int, int]:
    """Compute the k-mer frequency spectrum (frequency of frequencies).

    The spectrum maps each abundance value *f* to the number of distinct
    k-mers that appear exactly *f* times.  This is a standard metric in
    genome assembly quality assessment.

    Args:
        seq: DNA sequence string.
        k: Length of each k-mer.

    Returns:
        Dictionary mapping frequency *f* to number of k-mers with that
        frequency.

    Example:
        >>> kmer_spectrum("ATCGATCG", 2)
        {2: 3, 1: 1}
    """
    counts = count_kmers(seq, k)
    if not counts:
        return {}

    spectrum: Counter[int] = Counter()
    for count in counts.values():
        spectrum[count] += 1

    return dict(spectrum)


# ---------------------------------------------------------------------------
# Unique and over-represented k-mers
# ---------------------------------------------------------------------------


def find_unique_kmers(seq: str, k: int) -> List[str]:
    """Return all k-mers that appear exactly once in the sequence.

    Args:
        seq: DNA sequence string.
        k: Length of each k-mer.

    Returns:
        Sorted list of k-mers with count == 1.

    Example:
        >>> find_unique_kmers("ATCGATCG", 2)
        ['GA']
    """
    counts = count_kmers(seq, k)
    return sorted(kmer for kmer, count in counts.items() if count == 1)


def find_overrepresented_kmers(seq: str, k: int, threshold: float = 2.0) -> Dict[str, float]:
    """Identify k-mers whose frequency exceeds the expected uniform frequency
    by at least *threshold*-fold.

    Under a uniform random model with alphabet size 4, the expected frequency
    of any k-mer is ``1 / 4**k``.  A k-mer is flagged as over-represented when
    ``observed_frequency / expected_frequency >= threshold``.

    Args:
        seq: DNA sequence string.
        k: Length of each k-mer.
        threshold: Minimum fold-enrichment over expected frequency to report.
            Defaults to 2.0.

    Returns:
        Dictionary mapping over-represented k-mers to their fold-enrichment.

    Example:
        >>> find_overrepresented_kmers("AAAAAA", 2, threshold=1.5)
        {'AA': 16.0}
    """
    freqs = kmer_frequencies(seq, k)
    if not freqs:
        return {}

    expected = 1.0 / (4**k)
    if expected == 0.0:
        return {}

    result: Dict[str, float] = {}
    for kmer, freq in freqs.items():
        fold = freq / expected
        if fold >= threshold:
            result[kmer] = fold

    return dict(sorted(result.items(), key=lambda item: item[1], reverse=True))


# ---------------------------------------------------------------------------
# Diversity / entropy
# ---------------------------------------------------------------------------


def kmer_diversity(seq: str, k: int) -> float:
    """Compute Shannon entropy of the k-mer frequency distribution.

    The entropy is calculated in bits (log base 2).  Higher values indicate
    a more diverse k-mer composition; lower values indicate bias towards a
    small set of k-mers.

    Args:
        seq: DNA sequence string.
        k: Length of each k-mer.

    Returns:
        Shannon entropy in bits.  Returns 0.0 for empty or too-short
        sequences.

    Example:
        >>> round(kmer_diversity("ATCGATCG", 2), 4)
        1.9502
    """
    freqs = kmer_frequencies(seq, k)
    if not freqs:
        return 0.0

    entropy = 0.0
    for freq in freqs.values():
        if freq > 0.0:
            entropy -= freq * math.log2(freq)

    return entropy


# ---------------------------------------------------------------------------
# Profile comparison
# ---------------------------------------------------------------------------


def compare_kmer_profiles(seq1: str, seq2: str, k: int) -> Dict[str, float]:
    """Compare k-mer frequency profiles of two sequences.

    Computes three standard similarity/dissimilarity measures between the
    k-mer frequency vectors of *seq1* and *seq2*:

    * **cosine**: cosine similarity (1.0 = identical direction).
    * **jaccard**: Jaccard index over k-mer sets (presence/absence).
    * **bray_curtis**: Bray-Curtis dissimilarity (0.0 = identical, 1.0 =
      completely different).

    Args:
        seq1: First DNA sequence.
        seq2: Second DNA sequence.
        k: Length of each k-mer.

    Returns:
        Dictionary with keys ``cosine``, ``jaccard``, ``bray_curtis``.

    Example:
        >>> result = compare_kmer_profiles("ATCGATCG", "ATCGATCG", 2)
        >>> result["cosine"]
        1.0
    """
    freqs1 = kmer_frequencies(seq1, k)
    freqs2 = kmer_frequencies(seq2, k)

    if not freqs1 and not freqs2:
        return {"cosine": 1.0, "jaccard": 1.0, "bray_curtis": 0.0}
    if not freqs1 or not freqs2:
        return {"cosine": 0.0, "jaccard": 0.0, "bray_curtis": 1.0}

    all_kmers = set(freqs1.keys()) | set(freqs2.keys())

    # --- Cosine similarity ---
    dot = 0.0
    mag1 = 0.0
    mag2 = 0.0
    for kmer in all_kmers:
        v1 = freqs1.get(kmer, 0.0)
        v2 = freqs2.get(kmer, 0.0)
        dot += v1 * v2
        mag1 += v1 * v1
        mag2 += v2 * v2

    denom = math.sqrt(mag1) * math.sqrt(mag2)
    cosine = dot / denom if denom > 0.0 else 0.0

    # --- Jaccard index (set-based) ---
    set1 = set(freqs1.keys())
    set2 = set(freqs2.keys())
    intersection = len(set1 & set2)
    union = len(set1 | set2)
    jaccard = intersection / union if union > 0 else 0.0

    # --- Bray-Curtis dissimilarity ---
    sum_min = 0.0
    sum_total = 0.0
    for kmer in all_kmers:
        v1 = freqs1.get(kmer, 0.0)
        v2 = freqs2.get(kmer, 0.0)
        sum_min += min(v1, v2)
        sum_total += v1 + v2

    bray_curtis = 1.0 - (2.0 * sum_min / sum_total) if sum_total > 0.0 else 0.0

    return {"cosine": cosine, "jaccard": jaccard, "bray_curtis": bray_curtis}


# ---------------------------------------------------------------------------
# Sliding window analysis
# ---------------------------------------------------------------------------


def sliding_window_kmer_complexity(seq: str, k: int, window_size: int, step: int = 1) -> List[Tuple[int, float]]:
    """Compute k-mer diversity (Shannon entropy) in sliding windows.

    For each window starting at position *i*, the Shannon entropy of the
    k-mer distribution within that window is calculated.

    Args:
        seq: DNA sequence string.
        k: Length of each k-mer.
        window_size: Number of bases in each window.  Must be >= *k*.
        step: Step size for sliding the window.  Defaults to 1.

    Returns:
        List of ``(start_position, entropy)`` tuples.

    Example:
        >>> results = sliding_window_kmer_complexity("ATCGATCGATCG", 2, 6, 3)
        >>> len(results) >= 1
        True
    """
    if not seq or k <= 0 or window_size < k or step <= 0:
        return []

    seq_upper = seq.upper()
    n = len(seq_upper)
    if window_size > n:
        return []

    results: List[Tuple[int, float]] = []
    for start in range(0, n - window_size + 1, step):
        window = seq_upper[start : start + window_size]
        entropy = kmer_diversity(window, k)
        results.append((start, entropy))

    return results


# ---------------------------------------------------------------------------
# Homopolymer detection
# ---------------------------------------------------------------------------


def find_homopolymers(seq: str, min_length: int = 3) -> List[Tuple[str, int, int]]:
    """Find homopolymer runs (consecutive identical bases) in a sequence.

    Args:
        seq: DNA sequence string.
        min_length: Minimum run length to report.  Defaults to 3.

    Returns:
        List of ``(base, start_position, run_length)`` tuples, sorted by
        start position.

    Example:
        >>> find_homopolymers("AAATTTCCG", 3)
        [('A', 0, 3), ('T', 3, 3)]
    """
    if not seq or min_length <= 0:
        return []

    seq_upper = seq.upper()
    results: List[Tuple[str, int, int]] = []

    i = 0
    n = len(seq_upper)
    while i < n:
        base = seq_upper[i]
        run_start = i
        while i < n and seq_upper[i] == base:
            i += 1
        run_length = i - run_start
        if run_length >= min_length:
            results.append((base, run_start, run_length))

    return results


# ---------------------------------------------------------------------------
# Microsatellite / STR detection
# ---------------------------------------------------------------------------


def find_microsatellites(seq: str, min_repeats: int = 3, max_unit_length: int = 6) -> List[Tuple[str, int, int, int]]:
    """Find microsatellite (short tandem repeat) loci in a DNA sequence.

    Searches for tandemly repeated motifs of length 1 to *max_unit_length*
    that occur at least *min_repeats* consecutive times.  Homopolymers
    (unit length 1) are included.

    Args:
        seq: DNA sequence string.
        min_repeats: Minimum number of consecutive repeats to report.
            Defaults to 3.
        max_unit_length: Maximum length of the repeat unit.  Defaults to 6.

    Returns:
        List of ``(unit, start, end, repeat_count)`` tuples sorted by start
        position.  *end* is exclusive (Python slice convention).

    Example:
        >>> find_microsatellites("ACACACACTTTT", min_repeats=3, max_unit_length=2)
        [('AC', 0, 8, 4), ('T', 8, 12, 4)]
    """
    if not seq or min_repeats <= 0 or max_unit_length <= 0:
        return []

    seq_upper = seq.upper()
    n = len(seq_upper)
    results: List[Tuple[str, int, int, int]] = []

    # Track which positions have already been claimed by a longer/earlier
    # repeat so we can report non-overlapping, longest-unit-first hits.
    covered = [False] * n

    # Search from longest unit to shortest so longer motifs take priority.
    for unit_len in range(max_unit_length, 0, -1):
        for start in range(n - unit_len * min_repeats + 1):
            if covered[start]:
                continue

            unit = seq_upper[start : start + unit_len]
            # Extend as far as the unit keeps repeating
            pos = start + unit_len
            repeat_count = 1
            while pos + unit_len <= n and seq_upper[pos : pos + unit_len] == unit:
                repeat_count += 1
                pos += unit_len

            if repeat_count >= min_repeats:
                end = start + unit_len * repeat_count
                # Verify this is a true repeat unit -- not a sub-multiple of
                # a shorter unit that would produce the same pattern.  E.g.
                # "AAAAAA" with unit "AA" is really unit "A" repeated 6x.
                is_reducible = False
                for sub_len in range(1, unit_len):
                    if unit_len % sub_len == 0:
                        sub_unit = unit[:sub_len]
                        if sub_unit * (unit_len // sub_len) == unit:
                            is_reducible = True
                            break

                if not is_reducible:
                    results.append((unit, start, end, repeat_count))
                    for idx in range(start, end):
                        covered[idx] = True

    results.sort(key=lambda x: x[1])
    return results


# ---------------------------------------------------------------------------
# Sequence complexity profiling
# ---------------------------------------------------------------------------


def sequence_complexity_profile(seq: str, window_size: int = 50, step: int = 10) -> List[Tuple[int, float]]:
    """Compute linguistic complexity in sliding windows.

    Linguistic complexity is defined as the ratio of distinct k-mers
    observed in the window to the maximum possible number of k-mers for
    that window, using a *k* of 3 (trinucleotides).  Values range from
    0.0 (completely repetitive) to 1.0 (maximum diversity).

    Args:
        seq: DNA sequence string.
        window_size: Number of bases per window.  Defaults to 50.
        step: Step size for the sliding window.  Defaults to 10.

    Returns:
        List of ``(start_position, complexity)`` tuples.

    Example:
        >>> profile = sequence_complexity_profile("ATCG" * 20, window_size=20, step=10)
        >>> len(profile) >= 1
        True
    """
    if not seq or window_size <= 0 or step <= 0:
        return []

    k = 3  # trinucleotide complexity
    seq_upper = seq.upper()
    n = len(seq_upper)

    if window_size < k:
        return []
    if window_size > n:
        return []

    results: List[Tuple[int, float]] = []

    for start in range(0, n - window_size + 1, step):
        window = seq_upper[start : start + window_size]
        num_kmers_in_window = window_size - k + 1
        observed = len(set(window[i : i + k] for i in range(num_kmers_in_window)))

        # Maximum possible distinct k-mers is min(4**k, num_kmers_in_window)
        max_possible = min(4**k, num_kmers_in_window)
        complexity = observed / max_possible if max_possible > 0 else 0.0
        results.append((start, complexity))

    return results


# ---------------------------------------------------------------------------
# Low-complexity masking
# ---------------------------------------------------------------------------


def mask_low_complexity(seq: str, window_size: int = 50, threshold: float = 0.5) -> str:
    """Mask low-complexity regions of a DNA sequence with 'N'.

    Scans the sequence in sliding windows of *window_size* bases.  Any
    window whose linguistic complexity (trinucleotide diversity ratio) falls
    below *threshold* has all its bases replaced with 'N'.

    Args:
        seq: DNA sequence string.
        window_size: Sliding window size.  Defaults to 50.
        threshold: Complexity threshold below which a window is masked.
            Defaults to 0.5.

    Returns:
        Masked sequence with low-complexity regions replaced by 'N'.

    Example:
        >>> masked = mask_low_complexity("A" * 100 + "ATCGATCGATCG" * 10, window_size=20, threshold=0.3)
        >>> masked[:100] == "N" * 100
        True
    """
    if not seq or window_size <= 0:
        return seq

    seq_upper = seq.upper()
    n = len(seq_upper)

    if window_size > n:
        # Single window: evaluate entire sequence
        profile = sequence_complexity_profile(seq_upper, window_size=n, step=1)
        if profile and profile[0][1] < threshold:
            return "N" * n
        return seq_upper

    k = 3  # trinucleotide complexity (matches sequence_complexity_profile)
    if window_size < k:
        return seq_upper

    masked = list(seq_upper)

    for start in range(0, n - window_size + 1):
        window = seq_upper[start : start + window_size]
        num_kmers = window_size - k + 1
        observed = len(set(window[i : i + k] for i in range(num_kmers)))
        max_possible = min(4**k, num_kmers)
        complexity = observed / max_possible if max_possible > 0 else 0.0

        if complexity < threshold:
            for idx in range(start, start + window_size):
                masked[idx] = "N"

    return "".join(masked)
