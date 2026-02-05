"""Quality metrics for long-read sequencing data.

Provides N50/Nx statistics, read length distributions, quality score analysis,
accuracy estimation, and throughput calculations. All computations are real
mathematical implementations with no placeholder data.

Optional dependencies:
    - numpy: For efficient numerical computation on large datasets
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Sequence

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

try:
    import numpy as np  # type: ignore[import-untyped]
except ImportError:
    np = None  # type: ignore[assignment]


@dataclass
class ReadLengthStatistics:
    """Comprehensive statistics for a read length distribution.

    Attributes:
        count: Total number of reads.
        total_bases: Sum of all read lengths.
        mean_length: Arithmetic mean of read lengths.
        median_length: Median read length.
        min_length: Shortest read.
        max_length: Longest read.
        std_dev: Standard deviation of read lengths.
        n50: N50 read length.
        n90: N90 read length.
        l50: L50 (minimum number of reads whose lengths sum to >= 50% total).
        percentiles: Dictionary of read length percentiles (p5, p25, p75, p95).
    """

    count: int = 0
    total_bases: int = 0
    mean_length: float = 0.0
    median_length: float = 0.0
    min_length: int = 0
    max_length: int = 0
    std_dev: float = 0.0
    n50: int = 0
    n90: int = 0
    l50: int = 0
    percentiles: dict[str, float] = field(default_factory=dict)


@dataclass
class QualityDistribution:
    """Distribution of Phred quality scores across reads.

    Attributes:
        mean_quality: Mean Phred quality score.
        median_quality: Median Phred quality score.
        min_quality: Minimum quality score.
        max_quality: Maximum quality score.
        q7_fraction: Fraction of reads with mean Q >= 7.
        q10_fraction: Fraction of reads with mean Q >= 10.
        q15_fraction: Fraction of reads with mean Q >= 15.
        q20_fraction: Fraction of reads with mean Q >= 20.
        per_read_means: List of per-read mean quality scores.
        histogram: Quality score histogram (score -> count).
    """

    mean_quality: float = 0.0
    median_quality: float = 0.0
    min_quality: float = 0.0
    max_quality: float = 0.0
    q7_fraction: float = 0.0
    q10_fraction: float = 0.0
    q15_fraction: float = 0.0
    q20_fraction: float = 0.0
    per_read_means: list[float] = field(default_factory=list)
    histogram: dict[int, int] = field(default_factory=dict)


def calculate_n50(read_lengths: Sequence[int]) -> int:
    """Calculate the N50 statistic for a set of read lengths.

    N50 is defined as the length N such that 50% of the total assembled
    bases are in sequences of length >= N. Equivalently, it is the length
    of the shortest sequence in the set of longest sequences that together
    cover at least 50% of the total length.

    Algorithm:
        1. Sort lengths in descending order
        2. Accumulate lengths until sum >= 50% of total
        3. Return the length that crosses the 50% threshold

    Args:
        read_lengths: Sequence of integer read lengths.

    Returns:
        N50 value. Returns 0 if input is empty.
    """
    return calculate_nx(read_lengths, x=50)


def calculate_nx(read_lengths: Sequence[int], x: int = 50) -> int:
    """Calculate the Nx statistic for a set of read lengths.

    Nx is the generalization of N50: the length N such that x% of the total
    bases are in sequences of length >= N.

    Args:
        read_lengths: Sequence of integer read lengths.
        x: Percentage threshold (0-100). Default is 50 (N50).

    Returns:
        Nx value. Returns 0 if input is empty.

    Raises:
        ValueError: If x is not between 0 and 100.
    """
    if x < 0 or x > 100:
        raise ValueError(f"x must be between 0 and 100, got {x}")

    lengths = [l for l in read_lengths if l > 0]
    if not lengths:
        return 0

    sorted_lengths = sorted(lengths, reverse=True)
    total = sum(sorted_lengths)
    threshold = total * x / 100.0

    cumulative = 0
    for length in sorted_lengths:
        cumulative += length
        if cumulative >= threshold:
            return length

    return sorted_lengths[-1]


def _calculate_lx(read_lengths: Sequence[int], x: int = 50) -> int:
    """Calculate Lx: minimum number of reads whose lengths sum to >= x% of total.

    Args:
        read_lengths: Sequence of integer read lengths.
        x: Percentage threshold (0-100).

    Returns:
        Lx value (count of reads). Returns 0 if input is empty.
    """
    lengths = [l for l in read_lengths if l > 0]
    if not lengths:
        return 0

    sorted_lengths = sorted(lengths, reverse=True)
    total = sum(sorted_lengths)
    threshold = total * x / 100.0

    cumulative = 0
    for i, length in enumerate(sorted_lengths, 1):
        cumulative += length
        if cumulative >= threshold:
            return i

    return len(sorted_lengths)


def read_length_stats(reads: Sequence[dict[str, Any] | int]) -> ReadLengthStatistics:
    """Calculate comprehensive read length statistics.

    Accepts either a sequence of integers (read lengths) or a sequence of
    dictionaries with a 'length' or 'sequence' key.

    Args:
        reads: Sequence of read lengths (int) or read dictionaries containing
            'length' (int) or 'sequence' (str) keys.

    Returns:
        ReadLengthStatistics with all computed metrics.
    """
    # Extract lengths
    lengths: list[int] = []
    for r in reads:
        if isinstance(r, int):
            lengths.append(r)
        elif isinstance(r, dict):
            if "length" in r:
                lengths.append(int(r["length"]))
            elif "sequence" in r and r["sequence"]:
                lengths.append(len(r["sequence"]))
        elif hasattr(r, "sequence") and r.sequence:
            lengths.append(len(r.sequence))
        elif hasattr(r, "query_length"):
            lengths.append(int(r.query_length))

    if not lengths:
        return ReadLengthStatistics()

    lengths = [l for l in lengths if l > 0]
    if not lengths:
        return ReadLengthStatistics()

    sorted_lengths = sorted(lengths)
    count = len(sorted_lengths)
    total = sum(sorted_lengths)

    # Basic statistics
    mean_len = total / count
    min_len = sorted_lengths[0]
    max_len = sorted_lengths[-1]

    # Median
    if count % 2 == 0:
        median_len = (sorted_lengths[count // 2 - 1] + sorted_lengths[count // 2]) / 2.0
    else:
        median_len = float(sorted_lengths[count // 2])

    # Standard deviation
    if count > 1:
        variance = sum((l - mean_len) ** 2 for l in sorted_lengths) / (count - 1)
        std = math.sqrt(variance)
    else:
        std = 0.0

    # Percentiles
    def _percentile(data: list[int], p: float) -> float:
        """Calculate percentile using linear interpolation."""
        n = len(data)
        k = (n - 1) * p / 100.0
        f = math.floor(k)
        c = math.ceil(k)
        if f == c:
            return float(data[int(k)])
        return data[f] * (c - k) + data[c] * (k - f)

    percentiles = {
        "p5": _percentile(sorted_lengths, 5),
        "p25": _percentile(sorted_lengths, 25),
        "p75": _percentile(sorted_lengths, 75),
        "p95": _percentile(sorted_lengths, 95),
    }

    n50 = calculate_nx(sorted_lengths, 50)
    n90 = calculate_nx(sorted_lengths, 90)
    l50 = _calculate_lx(sorted_lengths, 50)

    return ReadLengthStatistics(
        count=count,
        total_bases=total,
        mean_length=mean_len,
        median_length=median_len,
        min_length=min_len,
        max_length=max_len,
        std_dev=std,
        n50=n50,
        n90=n90,
        l50=l50,
        percentiles=percentiles,
    )


def quality_score_distribution(
    reads: Sequence[dict[str, Any] | str],
) -> QualityDistribution:
    """Compute Phred quality score distribution across reads.

    Accepts quality strings (ASCII-encoded Phred+33) or read dictionaries
    with a 'quality' or 'quality_string' key.

    Args:
        reads: Sequence of quality strings or read dictionaries.

    Returns:
        QualityDistribution with computed metrics.
    """
    per_read_means: list[float] = []
    histogram: dict[int, int] = {}

    for r in reads:
        qual_str: str | None = None
        if isinstance(r, str):
            qual_str = r
        elif isinstance(r, dict):
            qual_str = r.get("quality") or r.get("quality_string")
        elif hasattr(r, "quality_string"):
            qual_str = r.quality_string

        if not qual_str:
            continue

        scores = [ord(c) - 33 for c in qual_str]
        if not scores:
            continue

        mean_q = sum(scores) / len(scores)
        per_read_means.append(mean_q)

        for s in scores:
            histogram[s] = histogram.get(s, 0) + 1

    if not per_read_means:
        return QualityDistribution()

    sorted_means = sorted(per_read_means)
    n = len(sorted_means)

    mean_quality = sum(sorted_means) / n
    if n % 2 == 0:
        median_quality = (sorted_means[n // 2 - 1] + sorted_means[n // 2]) / 2.0
    else:
        median_quality = sorted_means[n // 2]

    q7_count = sum(1 for q in sorted_means if q >= 7)
    q10_count = sum(1 for q in sorted_means if q >= 10)
    q15_count = sum(1 for q in sorted_means if q >= 15)
    q20_count = sum(1 for q in sorted_means if q >= 20)

    return QualityDistribution(
        mean_quality=mean_quality,
        median_quality=median_quality,
        min_quality=sorted_means[0],
        max_quality=sorted_means[-1],
        q7_fraction=q7_count / n,
        q10_fraction=q10_count / n,
        q15_fraction=q15_count / n,
        q20_fraction=q20_count / n,
        per_read_means=per_read_means,
        histogram=histogram,
    )


def estimate_accuracy(quality_scores: Sequence[int] | str) -> float:
    """Estimate base-level accuracy from Phred quality scores.

    Converts Phred scores to error probabilities, then computes the
    mean accuracy across all bases.

    The relationship is:
        P(error) = 10^(-Q/10)
        Accuracy = 1 - P(error)

    For a set of bases, the mean accuracy is:
        mean_accuracy = 1 - mean(P(error))

    Args:
        quality_scores: Sequence of integer Phred scores, or an
            ASCII quality string (Phred+33 encoded).

    Returns:
        Estimated accuracy as a fraction (0.0 to 1.0).
        Returns 0.0 if input is empty.
    """
    scores: list[int]
    if isinstance(quality_scores, str):
        scores = [ord(c) - 33 for c in quality_scores]
    else:
        scores = list(quality_scores)

    if not scores:
        return 0.0

    # Compute mean error probability
    error_probs = [10.0 ** (-q / 10.0) for q in scores]
    mean_error = sum(error_probs) / len(error_probs)

    return 1.0 - mean_error


def calculate_throughput(
    reads: Sequence[dict[str, Any] | int],
    run_duration: float | None = None,
) -> dict[str, float]:
    """Calculate sequencing throughput metrics.

    Computes total bases, bases per hour, and reads per hour from a set
    of reads and an optional run duration.

    Args:
        reads: Sequence of read lengths (int) or read dictionaries containing
            'length' (int) or 'sequence' (str) keys.
        run_duration: Total run duration in hours. If None, only total metrics
            are returned without per-hour rates.

    Returns:
        Dictionary with throughput metrics:
            - total_bases: Total sequenced bases.
            - total_reads: Total number of reads.
            - bases_per_hour: Bases per hour (if run_duration provided).
            - reads_per_hour: Reads per hour (if run_duration provided).
            - gigabases_per_hour: Gigabases per hour (if run_duration provided).
            - mean_read_length: Average read length.
    """
    lengths: list[int] = []
    for r in reads:
        if isinstance(r, int):
            lengths.append(r)
        elif isinstance(r, dict):
            if "length" in r:
                lengths.append(int(r["length"]))
            elif "sequence" in r and r["sequence"]:
                lengths.append(len(r["sequence"]))
        elif hasattr(r, "sequence") and r.sequence:
            lengths.append(len(r.sequence))
        elif hasattr(r, "query_length"):
            lengths.append(int(r.query_length))

    total_bases = sum(lengths)
    total_reads = len(lengths)
    mean_length = total_bases / total_reads if total_reads > 0 else 0.0

    result: dict[str, float] = {
        "total_bases": float(total_bases),
        "total_reads": float(total_reads),
        "mean_read_length": mean_length,
    }

    if run_duration is not None and run_duration > 0:
        result["bases_per_hour"] = total_bases / run_duration
        result["reads_per_hour"] = total_reads / run_duration
        result["gigabases_per_hour"] = total_bases / run_duration / 1e9

    return result
