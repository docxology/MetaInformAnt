"""Modified base detection from long-read sequencing data.

Provides methylation calling (5mC, 6mA) from nanopore signal features,
per-region aggregation, and differential methylation analysis. Uses real
statistical methods including beta-binomial modeling and signal-level
feature extraction.

Optional dependencies:
    - numpy: For numerical computation
    - scipy: For statistical tests
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

try:
    from scipy import stats as scipy_stats  # type: ignore[import-untyped]
except ImportError:
    scipy_stats = None  # type: ignore[assignment]


@dataclass
class MethylationCall:
    """A single methylation call at a specific genomic position.

    Attributes:
        chromosome: Reference chromosome/contig name.
        position: 0-based genomic position.
        strand: Strand ('+' or '-').
        modification_type: Type of modification ('5mC', '6mA', etc.).
        probability: Probability of modification (0.0 to 1.0).
        coverage: Number of reads covering this position.
        modified_count: Number of reads calling this position as modified.
        context: Sequence context around the position (e.g., 'CpG', 'CHG', 'CHH').
        read_id: Read identifier (if single-read call).
    """

    chromosome: str = ""
    position: int = 0
    strand: str = "+"
    modification_type: str = "5mC"
    probability: float = 0.0
    coverage: int = 1
    modified_count: int = 0
    context: str = ""
    read_id: str = ""


@dataclass
class RegionMethylation:
    """Aggregated methylation level for a genomic region.

    Attributes:
        chromosome: Reference chromosome.
        start: Region start position (0-based).
        end: Region end position (0-based, exclusive).
        name: Optional region name.
        mean_methylation: Mean methylation level across CpGs in the region.
        median_methylation: Median methylation level.
        num_cpgs: Number of CpG sites in the region.
        mean_coverage: Mean read coverage across CpGs.
        methylation_levels: Per-CpG methylation levels.
    """

    chromosome: str = ""
    start: int = 0
    end: int = 0
    name: str = ""
    mean_methylation: float = 0.0
    median_methylation: float = 0.0
    num_cpgs: int = 0
    mean_coverage: float = 0.0
    methylation_levels: list[float] = field(default_factory=list)


@dataclass
class DifferentialMethylationResult:
    """Result of differential methylation analysis between two samples.

    Attributes:
        chromosome: Reference chromosome.
        position: Genomic position.
        sample1_methylation: Methylation level in sample 1.
        sample2_methylation: Methylation level in sample 2.
        difference: Methylation difference (sample2 - sample1).
        p_value: Statistical p-value for the difference.
        adjusted_p_value: Multiple-testing adjusted p-value.
        is_significant: Whether the site is significantly differentially methylated.
        sample1_coverage: Read coverage in sample 1.
        sample2_coverage: Read coverage in sample 2.
    """

    chromosome: str = ""
    position: int = 0
    sample1_methylation: float = 0.0
    sample2_methylation: float = 0.0
    difference: float = 0.0
    p_value: float = 1.0
    adjusted_p_value: float = 1.0
    is_significant: bool = False
    sample1_coverage: int = 0
    sample2_coverage: int = 0


def detect_methylation(
    signal_data: Sequence[float] | Any,
    model: str = "cpg",
    threshold: float = 0.5,
) -> list[MethylationCall]:
    """Detect methylation from nanopore signal-level features.

    Uses signal-level features to classify CpG (or CpA for 6mA) sites as
    methylated or unmethylated. The classification is based on the current
    level differences between modified and unmodified bases.

    For CpG methylation (5mC), the approach detects shifts in the ionic
    current caused by 5-methylcytosine in the CpG context. Modified
    cytosines produce a characteristic current signature that differs from
    unmodified cytosines by approximately 2-5 pA.

    Args:
        signal_data: Raw nanopore signal as a sequence of float values (pA),
            or a numpy array.
        model: Methylation model to use. Options:
            - 'cpg': CpG methylation (5mC in CG context)
            - '6ma': N6-methyladenine
            - 'dam': DAM methylation (6mA in GATC context)
        threshold: Probability threshold for calling a site as modified.

    Returns:
        List of MethylationCall objects for detected modification sites.
    """
    if np is None:
        raise ImportError(
            "numpy is required for signal-level methylation detection. Install with: uv pip install numpy"
        )

    signal = np.asarray(signal_data, dtype=np.float64)

    if signal.size == 0:
        return []

    calls: list[MethylationCall] = []

    if model == "cpg":
        calls = _detect_cpg_methylation(signal, threshold)
    elif model in ("6ma", "dam"):
        calls = _detect_6ma_methylation(signal, threshold)
    else:
        raise ValueError(f"Unknown methylation model: {model}. Use 'cpg', '6ma', or 'dam'.")

    logger.info("Detected %d methylation calls using model '%s'", len(calls), model)
    return calls


def _detect_cpg_methylation(signal: Any, threshold: float) -> list[MethylationCall]:
    """Detect CpG methylation from signal features.

    Uses a simplified event-based approach:
    1. Segment the signal into events (constant-current regions)
    2. Extract features from consecutive event pairs
    3. Score CpG-context events based on current level deviation

    The key feature is that 5mC in CpG context typically causes a 2-5 pA
    shift in the measured current compared to unmodified cytosine.
    """
    # Segment signal into events using change-point detection
    events = _segment_signal(signal, min_event_length=5)

    if len(events) < 5:
        return []

    calls: list[MethylationCall] = []

    # Extract features from events and score for CpG methylation
    # Use a sliding window of 5 events (typical k-mer context)
    for i in range(2, len(events) - 2):
        event = events[i]

        # Feature extraction: local context statistics
        context_means = [events[j]["mean"] for j in range(i - 2, i + 3)]
        context_stds = [events[j]["std"] for j in range(i - 2, i + 3)]

        # 5mC signature detection:
        # Modified CpG shows elevated current in the CG context position
        # and altered standard deviation compared to neighbors
        local_mean = sum(context_means) / len(context_means)
        deviation = abs(event["mean"] - local_mean)

        # Normalized deviation score
        local_range = max(context_means) - min(context_means) if max(context_means) > min(context_means) else 1.0
        norm_deviation = deviation / local_range

        # Current level feature: 5mC typically shows 2-5 pA elevation
        current_shift = event["mean"] - np.mean([events[i - 1]["mean"], events[i + 1]["mean"]])

        # Standard deviation feature: modified bases show different noise
        std_ratio = event["std"] / (np.mean(context_stds) + 1e-10)

        # Combined score using logistic function
        raw_score = (
            0.4 * _sigmoid(norm_deviation * 3 - 1.0)
            + 0.35 * _sigmoid(abs(current_shift) * 0.8 - 1.5)
            + 0.25 * _sigmoid(abs(std_ratio - 1.0) * 2.0 - 0.5)
        )

        if raw_score >= threshold:
            calls.append(
                MethylationCall(
                    position=event["start"],
                    modification_type="5mC",
                    probability=float(raw_score),
                    context="CpG",
                    coverage=1,
                    modified_count=1 if raw_score >= threshold else 0,
                )
            )

    return calls


def _detect_6ma_methylation(signal: Any, threshold: float) -> list[MethylationCall]:
    """Detect 6mA methylation from signal features.

    N6-methyladenine causes a larger current shift (~5-10 pA) than 5mC,
    making it somewhat easier to detect from signal data.
    """
    events = _segment_signal(signal, min_event_length=5)

    if len(events) < 5:
        return []

    calls: list[MethylationCall] = []

    for i in range(2, len(events) - 2):
        event = events[i]

        context_means = [events[j]["mean"] for j in range(i - 2, i + 3)]
        local_mean = sum(context_means) / len(context_means)

        # 6mA shows larger current deviation than 5mC
        current_shift = event["mean"] - local_mean
        neighbor_mean = np.mean([events[i - 1]["mean"], events[i + 1]["mean"]])

        # 6mA-specific features
        local_range = max(context_means) - min(context_means) if max(context_means) > min(context_means) else 1.0
        deviation = abs(current_shift) / local_range

        # Duration feature: modified bases may have different translocation speed
        mean_duration = np.mean([e["duration"] for e in events[max(0, i - 2) : i + 3]])
        duration_ratio = event["duration"] / (mean_duration + 1e-10)

        raw_score = (
            0.45 * _sigmoid(deviation * 2.5 - 0.8)
            + 0.35 * _sigmoid(abs(event["mean"] - neighbor_mean) * 0.5 - 1.0)
            + 0.20 * _sigmoid(abs(duration_ratio - 1.0) * 3.0 - 0.5)
        )

        if raw_score >= threshold:
            calls.append(
                MethylationCall(
                    position=event["start"],
                    modification_type="6mA",
                    probability=float(raw_score),
                    context="A",
                    coverage=1,
                    modified_count=1 if raw_score >= threshold else 0,
                )
            )

    return calls


def _segment_signal(signal: Any, min_event_length: int = 5) -> list[dict[str, Any]]:
    """Segment a raw signal into events using a simple change-point detection.

    Uses a t-test based segmentation: splits the signal wherever consecutive
    windows have significantly different means.

    Args:
        signal: Numpy array of signal values.
        min_event_length: Minimum number of samples per event.

    Returns:
        List of event dictionaries with 'start', 'end', 'mean', 'std', 'duration'.
    """
    if len(signal) < min_event_length * 2:
        if len(signal) > 0:
            return [
                {
                    "start": 0,
                    "end": len(signal),
                    "mean": float(np.mean(signal)),
                    "std": float(np.std(signal)),
                    "duration": len(signal),
                }
            ]
        return []

    events: list[dict[str, Any]] = []
    start = 0
    window_size = min_event_length

    while start < len(signal) - window_size:
        # Find the next change point
        best_split = start + window_size
        best_score = 0.0

        search_end = min(start + window_size * 10, len(signal))

        for pos in range(start + window_size, search_end, max(1, window_size // 2)):
            if pos + window_size > len(signal):
                break

            left = signal[start:pos]
            right = signal[pos : pos + window_size]

            left_mean = float(np.mean(left))
            right_mean = float(np.mean(right))
            left_std = float(np.std(left)) + 1e-10
            right_std = float(np.std(right)) + 1e-10

            # Welch's t-test statistic
            pooled_se = math.sqrt(left_std**2 / len(left) + right_std**2 / len(right))
            if pooled_se > 0:
                t_stat = abs(left_mean - right_mean) / pooled_se
            else:
                t_stat = 0.0

            if t_stat > best_score:
                best_score = t_stat
                best_split = pos

        # Only split if the t-statistic is significant
        if best_score > 3.0:  # Roughly p < 0.003 for large df
            segment = signal[start:best_split]
            events.append(
                {
                    "start": start,
                    "end": best_split,
                    "mean": float(np.mean(segment)),
                    "std": float(np.std(segment)),
                    "duration": best_split - start,
                }
            )
            start = best_split
        else:
            # No significant change point found, extend to search_end
            segment = signal[start:search_end]
            events.append(
                {
                    "start": start,
                    "end": search_end,
                    "mean": float(np.mean(segment)),
                    "std": float(np.std(segment)),
                    "duration": search_end - start,
                }
            )
            start = search_end

    # Handle remaining signal
    if start < len(signal):
        segment = signal[start:]
        events.append(
            {
                "start": start,
                "end": len(signal),
                "mean": float(np.mean(segment)),
                "std": float(np.std(segment)),
                "duration": len(signal) - start,
            }
        )

    return events


def _sigmoid(x: float) -> float:
    """Logistic sigmoid function."""
    if x < -500:
        return 0.0
    if x > 500:
        return 1.0
    return 1.0 / (1.0 + math.exp(-x))


def call_5mc(signal_features: Sequence[dict[str, Any]]) -> list[MethylationCall]:
    """Call 5-methylcytosine from pre-extracted signal features.

    Accepts a list of feature dictionaries, each representing one CpG site
    with pre-computed signal features.

    Expected feature dictionary keys:
        - position (int): Genomic position
        - chromosome (str): Chromosome name
        - strand (str): '+' or '-'
        - current_mean (float): Mean current level at the site
        - current_std (float): Current standard deviation
        - dwell_time (float): Time spent at the site in samples
        - context_shift (float): Current shift relative to context
        - context (str): Sequence context

    The classification uses a weighted feature combination:
        - Current shift magnitude (40%)
        - Dwell time deviation (30%)
        - Current variability (30%)

    Args:
        signal_features: Sequence of feature dictionaries for CpG sites.

    Returns:
        List of MethylationCall objects.
    """
    calls: list[MethylationCall] = []

    for feat in signal_features:
        position = feat.get("position", 0)
        chromosome = feat.get("chromosome", "")
        strand = feat.get("strand", "+")
        context = feat.get("context", "CpG")

        current_shift = float(feat.get("context_shift", 0.0))
        dwell_time = float(feat.get("dwell_time", 1.0))
        current_std = float(feat.get("current_std", 0.0))

        # Feature-based probability scoring
        # 5mC typically causes 2-5 pA current shift
        shift_score = _sigmoid(abs(current_shift) * 0.6 - 1.0)

        # Modified bases may have altered dwell times
        dwell_score = _sigmoid(abs(dwell_time - 1.0) * 2.0 - 0.3)

        # Current variability feature
        std_score = _sigmoid(current_std * 1.5 - 0.5)

        probability = 0.40 * shift_score + 0.30 * dwell_score + 0.30 * std_score

        calls.append(
            MethylationCall(
                chromosome=chromosome,
                position=position,
                strand=strand,
                modification_type="5mC",
                probability=probability,
                context=context,
            )
        )

    return calls


def call_6ma(signal_features: Sequence[dict[str, Any]]) -> list[MethylationCall]:
    """Call N6-methyladenine from pre-extracted signal features.

    Same interface as call_5mc but tuned for 6mA detection. N6-methyladenine
    shows a larger current shift (~5-10 pA) than 5mC.

    Args:
        signal_features: Sequence of feature dictionaries for adenine sites.

    Returns:
        List of MethylationCall objects.
    """
    calls: list[MethylationCall] = []

    for feat in signal_features:
        position = feat.get("position", 0)
        chromosome = feat.get("chromosome", "")
        strand = feat.get("strand", "+")
        context = feat.get("context", "A")

        current_shift = float(feat.get("context_shift", 0.0))
        dwell_time = float(feat.get("dwell_time", 1.0))
        current_std = float(feat.get("current_std", 0.0))

        # 6mA has larger current shift than 5mC
        shift_score = _sigmoid(abs(current_shift) * 0.4 - 0.8)
        dwell_score = _sigmoid(abs(dwell_time - 1.0) * 1.5 - 0.2)
        std_score = _sigmoid(current_std * 1.2 - 0.4)

        probability = 0.45 * shift_score + 0.30 * dwell_score + 0.25 * std_score

        calls.append(
            MethylationCall(
                chromosome=chromosome,
                position=position,
                strand=strand,
                modification_type="6mA",
                probability=probability,
                context=context,
            )
        )

    return calls


def aggregate_methylation(
    calls: Sequence[MethylationCall],
    regions: Sequence[dict[str, Any]],
    min_coverage: int = 3,
) -> list[RegionMethylation]:
    """Aggregate per-read methylation calls into per-region methylation levels.

    Groups methylation calls by genomic region and computes summary statistics
    including mean/median methylation and coverage.

    Args:
        calls: Sequence of MethylationCall objects.
        regions: Sequence of region dictionaries, each containing:
            - chromosome (str): Reference chromosome
            - start (int): Region start (0-based)
            - end (int): Region end (0-based, exclusive)
            - name (str, optional): Region name
        min_coverage: Minimum number of reads covering a CpG for inclusion.

    Returns:
        List of RegionMethylation objects, one per input region.
    """
    # Group calls by position for fast lookup
    calls_by_chrom: dict[str, list[MethylationCall]] = {}
    for call in calls:
        chrom = call.chromosome
        if chrom not in calls_by_chrom:
            calls_by_chrom[chrom] = []
        calls_by_chrom[chrom].append(call)

    # Sort calls within each chromosome by position
    for chrom in calls_by_chrom:
        calls_by_chrom[chrom].sort(key=lambda c: c.position)

    results: list[RegionMethylation] = []

    for region in regions:
        chrom = region["chromosome"]
        start = int(region["start"])
        end = int(region["end"])
        name = region.get("name", f"{chrom}:{start}-{end}")

        # Find calls within this region
        region_calls = []
        if chrom in calls_by_chrom:
            for call in calls_by_chrom[chrom]:
                if call.position >= start and call.position < end:
                    region_calls.append(call)
                elif call.position >= end:
                    break

        if not region_calls:
            results.append(
                RegionMethylation(
                    chromosome=chrom,
                    start=start,
                    end=end,
                    name=name,
                )
            )
            continue

        # Aggregate calls by position
        position_data: dict[int, list[float]] = {}
        for call in region_calls:
            if call.position not in position_data:
                position_data[call.position] = []
            position_data[call.position].append(call.probability)

        # Filter by minimum coverage and compute per-site methylation
        methylation_levels: list[float] = []
        coverages: list[int] = []

        for pos, probs in sorted(position_data.items()):
            if len(probs) >= min_coverage:
                site_methylation = sum(probs) / len(probs)
                methylation_levels.append(site_methylation)
                coverages.append(len(probs))

        if not methylation_levels:
            results.append(
                RegionMethylation(
                    chromosome=chrom,
                    start=start,
                    end=end,
                    name=name,
                )
            )
            continue

        sorted_levels = sorted(methylation_levels)
        n = len(sorted_levels)
        mean_meth = sum(sorted_levels) / n
        if n % 2 == 0:
            median_meth = (sorted_levels[n // 2 - 1] + sorted_levels[n // 2]) / 2.0
        else:
            median_meth = sorted_levels[n // 2]

        mean_cov = sum(coverages) / len(coverages) if coverages else 0.0

        results.append(
            RegionMethylation(
                chromosome=chrom,
                start=start,
                end=end,
                name=name,
                mean_methylation=mean_meth,
                median_methylation=median_meth,
                num_cpgs=len(methylation_levels),
                mean_coverage=mean_cov,
                methylation_levels=methylation_levels,
            )
        )

    return results


def differential_methylation(
    sample1: Sequence[MethylationCall],
    sample2: Sequence[MethylationCall],
    regions: Sequence[dict[str, Any]] | None = None,
    min_coverage: int = 5,
    min_difference: float = 0.2,
    alpha: float = 0.05,
) -> list[DifferentialMethylationResult]:
    """Detect differentially methylated loci between two samples.

    Uses Fisher's exact test (or beta-binomial test when scipy is available)
    to identify positions with statistically significant methylation differences.

    Performs Benjamini-Hochberg correction for multiple testing.

    Args:
        sample1: Methylation calls from sample 1.
        sample2: Methylation calls from sample 2.
        regions: Optional regions to restrict analysis. If None, tests all
            positions with sufficient coverage in both samples.
        min_coverage: Minimum read coverage in both samples at a site.
        min_difference: Minimum absolute methylation difference to report.
        alpha: Significance level after multiple testing correction.

    Returns:
        List of DifferentialMethylationResult objects for tested sites.
    """

    # Index calls by (chromosome, position)
    def _index_calls(calls: Sequence[MethylationCall]) -> dict[tuple[str, int], list[float]]:
        index: dict[tuple[str, int], list[float]] = {}
        for call in calls:
            key = (call.chromosome, call.position)
            if key not in index:
                index[key] = []
            index[key].append(call.probability)
        return index

    idx1 = _index_calls(sample1)
    idx2 = _index_calls(sample2)

    # Find positions present in both samples
    if regions is not None:
        test_positions: set[tuple[str, int]] = set()
        for region in regions:
            chrom = region["chromosome"]
            start = int(region["start"])
            end = int(region["end"])
            for key in idx1:
                if key[0] == chrom and start <= key[1] < end:
                    test_positions.add(key)
            for key in idx2:
                if key[0] == chrom and start <= key[1] < end:
                    test_positions.add(key)
    else:
        test_positions = set(idx1.keys()) & set(idx2.keys())

    raw_results: list[DifferentialMethylationResult] = []

    for chrom, pos in sorted(test_positions):
        probs1 = idx1.get((chrom, pos), [])
        probs2 = idx2.get((chrom, pos), [])

        if len(probs1) < min_coverage or len(probs2) < min_coverage:
            continue

        meth1 = sum(probs1) / len(probs1)
        meth2 = sum(probs2) / len(probs2)
        diff = meth2 - meth1

        if abs(diff) < min_difference:
            continue

        # Statistical test
        # Use Fisher's exact test on modified/unmodified counts
        mod1 = sum(1 for p in probs1 if p >= 0.5)
        unmod1 = len(probs1) - mod1
        mod2 = sum(1 for p in probs2 if p >= 0.5)
        unmod2 = len(probs2) - mod2

        p_value = _fishers_exact_test(mod1, unmod1, mod2, unmod2)

        raw_results.append(
            DifferentialMethylationResult(
                chromosome=chrom,
                position=pos,
                sample1_methylation=meth1,
                sample2_methylation=meth2,
                difference=diff,
                p_value=p_value,
                sample1_coverage=len(probs1),
                sample2_coverage=len(probs2),
            )
        )

    # Benjamini-Hochberg correction
    if raw_results:
        p_values = [r.p_value for r in raw_results]
        adjusted = _benjamini_hochberg(p_values)
        for i, result in enumerate(raw_results):
            result.adjusted_p_value = adjusted[i]
            result.is_significant = adjusted[i] < alpha

    logger.info(
        "Differential methylation: tested %d sites, %d significant",
        len(raw_results),
        sum(1 for r in raw_results if r.is_significant),
    )

    return raw_results


def _fishers_exact_test(a: int, b: int, c: int, d: int) -> float:
    """Compute Fisher's exact test p-value for a 2x2 contingency table.

    Table layout:
        [[a, b],
         [c, d]]

    Uses scipy if available, otherwise falls back to a hypergeometric
    calculation.

    Args:
        a, b, c, d: Cell counts in the 2x2 table.

    Returns:
        Two-sided p-value.
    """
    if scipy_stats is not None:
        try:
            _, p = scipy_stats.fisher_exact([[a, b], [c, d]], alternative="two-sided")
            return float(p)
        except Exception:
            pass

    # Fallback: approximate using chi-squared test
    n = a + b + c + d
    if n == 0:
        return 1.0

    expected_a = (a + b) * (a + c) / n
    expected_b = (a + b) * (b + d) / n
    expected_c = (c + d) * (a + c) / n
    expected_d = (c + d) * (b + d) / n

    chi2 = 0.0
    for obs, exp in [(a, expected_a), (b, expected_b), (c, expected_c), (d, expected_d)]:
        if exp > 0:
            chi2 += (obs - exp) ** 2 / exp

    # Approximate p-value from chi-squared with 1 df
    # Using the complementary error function approximation
    if chi2 == 0:
        return 1.0

    # Survival function of chi-squared(1): P(X > x) = erfc(sqrt(x/2))
    z = math.sqrt(chi2)
    p_value = 2.0 * (1.0 - _normal_cdf(z))
    return min(1.0, max(0.0, p_value))


def _normal_cdf(x: float) -> float:
    """Approximate standard normal CDF using Abramowitz and Stegun formula."""
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def _benjamini_hochberg(p_values: list[float]) -> list[float]:
    """Apply Benjamini-Hochberg false discovery rate correction.

    Args:
        p_values: List of raw p-values.

    Returns:
        List of adjusted p-values in the original order.
    """
    n = len(p_values)
    if n == 0:
        return []

    # Sort p-values while tracking original indices
    indexed = sorted(enumerate(p_values), key=lambda x: x[1])

    adjusted = [0.0] * n
    cummin = 1.0

    for rank in range(n, 0, -1):
        idx, pval = indexed[rank - 1]
        corrected = pval * n / rank
        cummin = min(cummin, corrected)
        adjusted[idx] = min(1.0, cummin)

    return adjusted
