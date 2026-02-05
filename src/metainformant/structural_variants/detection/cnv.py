"""Copy Number Variation (CNV) detection from read depth data.

Implements circular binary segmentation (CBS) for CNV detection from
read depth profiles, log2 ratio computation with GC correction,
copy number state assignment, and adjacent segment merging.

References:
    Olshen AB et al. (2004) Circular binary segmentation for the analysis
    of array-based DNA copy number data. Biostatistics 5(4):557-572.
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass, field
from typing import Any, Sequence

from metainformant.core.utils.logging import get_logger

try:
    import numpy as np
    from numpy.typing import NDArray
except ImportError:
    np = None  # type: ignore[assignment]
    NDArray = Any  # type: ignore[assignment, misc]

try:
    from scipy import stats as scipy_stats
except ImportError:
    scipy_stats = None  # type: ignore[assignment]

logger = get_logger(__name__)

# Config prefix for environment variable overrides
_ENV_PREFIX = "SV_"


@dataclass
class CNVSegment:
    """A segment of constant copy number.

    Attributes:
        chrom: Chromosome name.
        start: 0-based start position (inclusive).
        end: 0-based end position (exclusive).
        mean_log2ratio: Mean log2 ratio for this segment.
        n_bins: Number of bins (windows) in the segment.
        state: Copy number state ('DEL', 'DUP', 'NEUTRAL', 'AMP', 'HOMODEL').
        cn: Estimated integer copy number.
        confidence: Confidence score [0, 1].
    """

    chrom: str
    start: int
    end: int
    mean_log2ratio: float
    n_bins: int
    state: str = "NEUTRAL"
    cn: int = 2
    confidence: float = 0.0


@dataclass
class CNVResult:
    """Result of CNV detection.

    Attributes:
        segments: List of CNV segments detected.
        log2_ratios: The log2 ratio array used for detection.
        bin_size: Size of each bin/window in base pairs.
        method: Detection method used.
        parameters: Dictionary of parameters used.
    """

    segments: list[CNVSegment]
    log2_ratios: list[float] | None = None
    bin_size: int = 1000
    method: str = "segmentation"
    parameters: dict[str, Any] = field(default_factory=dict)


def detect_cnv_from_depth(
    depth_data: dict[str, list[float] | Any],
    window_size: int = 1000,
    method: str = "segmentation",
    significance: float = 0.01,
    ploidy: int = 2,
    min_segment_bins: int = 3,
    merge_distance: int = 1000,
) -> dict[str, CNVResult]:
    """Detect copy number variations from read depth data.

    Analyzes per-chromosome depth arrays to identify regions of copy number
    gain or loss using circular binary segmentation (CBS).

    Args:
        depth_data: Dictionary mapping chromosome names to arrays of read depth
            values. Each value represents the mean depth in a genomic window.
            Can be numpy arrays or plain lists of floats.
        window_size: Size of each depth bin in base pairs.
        method: Detection method. Currently supports 'segmentation' (CBS).
        significance: P-value threshold for CBS change-point detection.
        ploidy: Expected ploidy of the organism (default 2 for diploid).
        min_segment_bins: Minimum number of bins for a valid segment.
        merge_distance: Maximum gap (in bp) to merge adjacent same-state segments.

    Returns:
        Dictionary mapping chromosome names to CNVResult objects containing
        detected segments with copy number states.

    Raises:
        ValueError: If depth_data is empty or method is unsupported.
    """
    if np is None:
        raise ImportError("numpy is required for CNV detection. Install with: uv pip install numpy")

    if not depth_data:
        raise ValueError("depth_data must not be empty")

    if method not in ("segmentation",):
        raise ValueError(f"Unsupported method '{method}'. Supported: 'segmentation'")

    # Read environment overrides
    env_significance = os.environ.get(f"{_ENV_PREFIX}CNV_SIGNIFICANCE")
    if env_significance is not None:
        significance = float(env_significance)

    results: dict[str, CNVResult] = {}

    for chrom, depths in depth_data.items():
        depth_array = np.asarray(depths, dtype=np.float64)

        if len(depth_array) < min_segment_bins:
            logger.warning(f"Chromosome {chrom}: too few bins ({len(depth_array)}), skipping")
            continue

        # Compute log2 ratio relative to median (self-normalization)
        median_depth = np.median(depth_array[depth_array > 0]) if np.any(depth_array > 0) else 1.0
        if median_depth <= 0:
            median_depth = 1.0

        log2_ratios = np.log2(depth_array / median_depth + 1e-10)
        # Replace -inf/inf with bounded values
        log2_ratios = np.clip(log2_ratios, -5.0, 5.0)

        # Segment using CBS
        segments = segment_coverage(log2_ratios, significance=significance)

        # Filter small segments
        segments = [s for s in segments if s[1] - s[0] >= min_segment_bins]

        # Build CNVSegment objects with states
        cnv_segments: list[CNVSegment] = []
        for seg_start, seg_end, seg_mean in segments:
            state, cn = _assign_cnv_state(seg_mean, ploidy)
            confidence = _calculate_segment_confidence(
                log2_ratios[seg_start:seg_end], seg_mean, state
            )
            cnv_segments.append(
                CNVSegment(
                    chrom=chrom,
                    start=seg_start * window_size,
                    end=seg_end * window_size,
                    mean_log2ratio=float(seg_mean),
                    n_bins=seg_end - seg_start,
                    state=state,
                    cn=cn,
                    confidence=confidence,
                )
            )

        # Merge adjacent same-state segments
        cnv_segments = merge_adjacent_segments(cnv_segments, max_gap=merge_distance)

        results[chrom] = CNVResult(
            segments=cnv_segments,
            log2_ratios=log2_ratios.tolist(),
            bin_size=window_size,
            method=method,
            parameters={
                "significance": significance,
                "ploidy": ploidy,
                "min_segment_bins": min_segment_bins,
                "merge_distance": merge_distance,
            },
        )

    logger.info(
        f"CNV detection complete: {sum(len(r.segments) for r in results.values())} "
        f"segments across {len(results)} chromosomes"
    )
    return results


def segment_coverage(
    coverage_array: Any,
    significance: float = 0.01,
    max_iterations: int = 100,
    min_width: int = 2,
) -> list[tuple[int, int, float]]:
    """Circular Binary Segmentation (CBS) algorithm for coverage segmentation.

    Recursively finds change-points in a coverage array by testing for
    differences in means across circular permutations of the data, splitting
    at the most significant change-point, and recursing on each half.

    Args:
        coverage_array: Array of log2 ratio values to segment. Can be a numpy
            array or list of floats.
        significance: P-value threshold for accepting a change-point. Lower
            values produce fewer, more confident segments.
        max_iterations: Maximum recursion depth to prevent infinite loops.
        min_width: Minimum segment width in bins.

    Returns:
        List of (start, end, mean) tuples representing segments, where start
        is inclusive, end is exclusive, and mean is the segment mean.
    """
    if np is None:
        raise ImportError("numpy is required for CBS segmentation. Install with: uv pip install numpy")

    data = np.asarray(coverage_array, dtype=np.float64)
    n = len(data)

    if n < 2 * min_width:
        return [(0, n, float(np.mean(data)))]

    segments: list[tuple[int, int, float]] = []
    _cbs_recursive(data, 0, n, segments, significance, max_iterations, min_width, depth=0)

    # Sort by start position
    segments.sort(key=lambda x: x[0])
    return segments


def _cbs_recursive(
    data: Any,
    start: int,
    end: int,
    segments: list[tuple[int, int, float]],
    significance: float,
    max_iterations: int,
    min_width: int,
    depth: int,
) -> None:
    """Recursive CBS implementation.

    Tests for a change-point in data[start:end] using the maximum circular
    binary segmentation statistic. If significant, splits and recurses.

    Args:
        data: Full data array.
        start: Start index (inclusive).
        end: End index (exclusive).
        segments: Accumulator list for resulting segments.
        significance: P-value threshold.
        max_iterations: Maximum recursion depth.
        min_width: Minimum segment width.
        depth: Current recursion depth.
    """
    n = end - start

    if n < 2 * min_width or depth >= max_iterations:
        segments.append((start, end, float(np.mean(data[start:end]))))
        return

    # Compute the CBS test statistic
    segment_data = data[start:end]
    best_t_stat = 0.0
    best_i = -1
    best_j = -1

    cumsum = np.cumsum(segment_data - np.mean(segment_data))
    cumsum = np.insert(cumsum, 0, 0.0)  # prepend 0

    # Variance of the segment
    variance = np.var(segment_data, ddof=1) if n > 1 else 1.0
    if variance < 1e-12:
        # Constant segment, no change-point
        segments.append((start, end, float(np.mean(segment_data))))
        return

    std_dev = math.sqrt(variance)

    # Test all possible arc (i, j) pairs for circular binary segmentation
    # Optimization: use vectorized approach for the single change-point scan first
    # For efficiency, test single change-point model (equivalent to standard BS)
    # then extend to circular model

    # Standard binary segmentation: find single best split point
    best_single_stat = 0.0
    best_split = -1

    for k in range(min_width, n - min_width + 1):
        # T-statistic for splitting at position k
        left_mean = np.mean(segment_data[:k])
        right_mean = np.mean(segment_data[k:])
        n_left = k
        n_right = n - k

        # Weighted difference of means
        t_stat = abs(left_mean - right_mean) * math.sqrt(n_left * n_right / n) / std_dev

        if t_stat > best_single_stat:
            best_single_stat = t_stat
            best_split = k

    # Circular scan: test arc segments (i, j) where the arc mean differs from rest
    best_arc_stat = 0.0
    best_arc_i = -1
    best_arc_j = -1

    # Efficient circular scan using cumulative sums
    for i in range(0, n - min_width + 1):
        for j in range(i + min_width, min(i + n // 2 + 1, n + 1)):
            arc_length = j - i
            rest_length = n - arc_length
            if rest_length < min_width:
                continue

            arc_sum = cumsum[j] - cumsum[i]
            arc_mean_diff = arc_sum / arc_length - (cumsum[n] - arc_sum) / rest_length if rest_length > 0 else 0.0

            t_stat = abs(arc_mean_diff) * math.sqrt(arc_length * rest_length / n) / std_dev

            if t_stat > best_arc_stat:
                best_arc_stat = t_stat
                best_arc_i = i
                best_arc_j = j

    # Use the better of single split vs arc
    use_arc = best_arc_stat > best_single_stat

    if use_arc:
        best_t_stat = best_arc_stat
        best_i = best_arc_i
        best_j = best_arc_j
    else:
        best_t_stat = best_single_stat
        best_i = 0
        best_j = best_split

    # Assess significance via permutation-derived critical values
    # Approximation using the asymptotic distribution of the CBS statistic
    # For large n, the statistic follows approximately a Gumbel distribution
    p_value = _cbs_pvalue(best_t_stat, n)

    if p_value > significance:
        # No significant change-point
        segments.append((start, end, float(np.mean(segment_data))))
        return

    if use_arc:
        # Split into three parts: [start, start+i), [start+i, start+j), [start+j, end)
        abs_i = start + best_i
        abs_j = start + best_j

        if best_i > 0:
            _cbs_recursive(data, start, abs_i, segments, significance, max_iterations, min_width, depth + 1)
        _cbs_recursive(data, abs_i, abs_j, segments, significance, max_iterations, min_width, depth + 1)
        if abs_j < end:
            _cbs_recursive(data, abs_j, end, segments, significance, max_iterations, min_width, depth + 1)
    else:
        # Standard split
        abs_split = start + best_j
        _cbs_recursive(data, start, abs_split, segments, significance, max_iterations, min_width, depth + 1)
        _cbs_recursive(data, abs_split, end, segments, significance, max_iterations, min_width, depth + 1)


def _cbs_pvalue(statistic: float, n: int) -> float:
    """Approximate p-value for CBS test statistic.

    Uses the asymptotic approximation from Olshen et al. (2004) where the
    maximum of the standardized CBS statistic converges to a Gumbel
    distribution with parameters depending on the segment length.

    Args:
        statistic: CBS test statistic value.
        n: Number of data points in the segment.

    Returns:
        Approximate p-value.
    """
    if n < 4:
        return 1.0

    # Asymptotic approximation: the max of |T| over all arcs
    # follows approximately a distribution related to the maximum of
    # a Brownian bridge. Use Gumbel approximation.
    # Location and scale parameters derived from segment length
    a_n = math.sqrt(2 * math.log(n))
    b_n = 2 * math.log(n) + 0.5 * math.log(math.log(n)) - 0.5 * math.log(math.pi)

    if a_n < 1e-12:
        return 1.0

    # Standardize
    z = (statistic * a_n) - b_n

    # Gumbel survival function: P(X > x) = 1 - exp(-exp(-x))
    try:
        p_value = 1.0 - math.exp(-math.exp(-z))
    except OverflowError:
        p_value = 0.0 if z > 0 else 1.0

    return max(0.0, min(1.0, p_value))


def call_cnv_states(
    segments: list[tuple[int, int, float]] | list[CNVSegment],
    ploidy: int = 2,
    del_threshold: float = -0.3,
    dup_threshold: float = 0.3,
    amp_threshold: float = 1.0,
    homodel_threshold: float = -1.5,
) -> list[CNVSegment]:
    """Assign copy number states to segments based on log2 ratio thresholds.

    Classifies each segment as homozygous deletion, deletion, neutral,
    duplication, or amplification based on its mean log2 ratio value
    relative to configurable thresholds.

    Args:
        segments: List of segments as (start, end, mean_log2ratio) tuples
            or CNVSegment objects.
        ploidy: Expected ploidy (default 2 for diploid).
        del_threshold: Log2 ratio threshold below which a segment is called
            as a deletion (default -0.3, approximately CN=1 for diploid).
        dup_threshold: Log2 ratio threshold above which a segment is called
            as a duplication (default 0.3, approximately CN=3 for diploid).
        amp_threshold: Log2 ratio threshold for high-level amplification
            (default 1.0, approximately CN=4+ for diploid).
        homodel_threshold: Log2 ratio threshold for homozygous deletion
            (default -1.5, approximately CN=0 for diploid).

    Returns:
        List of CNVSegment objects with state and cn fields populated.
    """
    # Read environment overrides for thresholds
    env_del = os.environ.get(f"{_ENV_PREFIX}CNV_DEL_THRESHOLD")
    env_dup = os.environ.get(f"{_ENV_PREFIX}CNV_DUP_THRESHOLD")
    if env_del is not None:
        del_threshold = float(env_del)
    if env_dup is not None:
        dup_threshold = float(env_dup)

    result: list[CNVSegment] = []

    for seg in segments:
        if isinstance(seg, CNVSegment):
            mean_lr = seg.mean_log2ratio
            chrom = seg.chrom
            start = seg.start
            end = seg.end
            n_bins = seg.n_bins
        else:
            start_idx, end_idx, mean_lr = seg
            chrom = "unknown"
            start = start_idx
            end = end_idx
            n_bins = end_idx - start_idx

        state, cn = _assign_cnv_state(
            mean_lr, ploidy, del_threshold, dup_threshold, amp_threshold, homodel_threshold
        )

        result.append(
            CNVSegment(
                chrom=chrom,
                start=start,
                end=end,
                mean_log2ratio=float(mean_lr),
                n_bins=n_bins,
                state=state,
                cn=cn,
                confidence=0.0,
            )
        )

    return result


def _assign_cnv_state(
    mean_log2ratio: float,
    ploidy: int = 2,
    del_threshold: float = -0.3,
    dup_threshold: float = 0.3,
    amp_threshold: float = 1.0,
    homodel_threshold: float = -1.5,
) -> tuple[str, int]:
    """Assign a copy number state and integer CN from a log2 ratio.

    Args:
        mean_log2ratio: Mean log2 ratio value for the segment.
        ploidy: Expected ploidy.
        del_threshold: Deletion threshold.
        dup_threshold: Duplication threshold.
        amp_threshold: Amplification threshold.
        homodel_threshold: Homozygous deletion threshold.

    Returns:
        Tuple of (state_string, integer_copy_number).
    """
    # Convert log2 ratio to copy number: CN = ploidy * 2^(log2ratio)
    cn_float = ploidy * (2 ** mean_log2ratio)
    cn = max(0, round(cn_float))

    if mean_log2ratio <= homodel_threshold:
        return "HOMODEL", 0
    elif mean_log2ratio <= del_threshold:
        return "DEL", max(0, min(cn, ploidy - 1))
    elif mean_log2ratio >= amp_threshold:
        return "AMP", cn
    elif mean_log2ratio >= dup_threshold:
        return "DUP", max(ploidy + 1, cn)
    else:
        return "NEUTRAL", ploidy


def _calculate_segment_confidence(
    log2_ratios: Any,
    segment_mean: float,
    state: str,
) -> float:
    """Calculate confidence score for a CNV segment.

    Confidence is based on: (1) deviation from neutral, (2) variance within
    the segment (lower is better), and (3) number of supporting bins.

    Args:
        log2_ratios: Array of log2 ratio values within the segment.
        segment_mean: Mean log2 ratio for the segment.
        state: Assigned CNV state.

    Returns:
        Confidence score between 0 and 1.
    """
    if np is None:
        return 0.0

    data = np.asarray(log2_ratios)
    n = len(data)

    if n == 0 or state == "NEUTRAL":
        return 0.0

    # Component 1: deviation from neutral (higher deviation = higher confidence)
    deviation_score = min(1.0, abs(segment_mean) / 1.5)

    # Component 2: low variance (tighter distribution = more confident)
    if n > 1:
        segment_std = float(np.std(data, ddof=1))
        variance_score = max(0.0, 1.0 - segment_std / (abs(segment_mean) + 0.01))
    else:
        variance_score = 0.5

    # Component 3: length support (more bins = more confident, diminishing returns)
    length_score = min(1.0, n / 20.0)

    # Weighted combination
    confidence = 0.4 * deviation_score + 0.35 * variance_score + 0.25 * length_score
    return max(0.0, min(1.0, confidence))


def merge_adjacent_segments(
    segments: list[CNVSegment],
    max_gap: int = 1000,
) -> list[CNVSegment]:
    """Merge adjacent segments with the same copy number state.

    Combines neighboring segments on the same chromosome that share
    the same CNV state and are separated by no more than max_gap base pairs.

    Args:
        segments: List of CNVSegment objects, ideally sorted by position.
        max_gap: Maximum gap in base pairs between segments to allow merging.

    Returns:
        List of merged CNVSegment objects.
    """
    if not segments:
        return []

    # Sort by chromosome and start position
    sorted_segs = sorted(segments, key=lambda s: (s.chrom, s.start))

    merged: list[CNVSegment] = [sorted_segs[0]]

    for seg in sorted_segs[1:]:
        prev = merged[-1]

        # Check if we can merge
        if (
            seg.chrom == prev.chrom
            and seg.state == prev.state
            and (seg.start - prev.end) <= max_gap
        ):
            # Merge: extend the previous segment
            total_bins = prev.n_bins + seg.n_bins
            weighted_mean = (
                prev.mean_log2ratio * prev.n_bins + seg.mean_log2ratio * seg.n_bins
            ) / total_bins

            merged[-1] = CNVSegment(
                chrom=prev.chrom,
                start=prev.start,
                end=seg.end,
                mean_log2ratio=weighted_mean,
                n_bins=total_bins,
                state=prev.state,
                cn=prev.cn,
                confidence=max(prev.confidence, seg.confidence),
            )
        else:
            merged.append(seg)

    logger.debug(f"Merged {len(sorted_segs)} segments into {len(merged)}")
    return merged


def calculate_log2_ratio(
    tumor_depth: Any,
    normal_depth: Any,
    gc_content: Any | None = None,
    pseudocount: float = 1.0,
) -> Any:
    """Compute log2 ratio between tumor and normal read depths.

    Calculates the log2(tumor/normal) ratio with optional GC-content bias
    correction using LOESS-like smoothing.

    Args:
        tumor_depth: Array of tumor read depth values per window.
        normal_depth: Array of normal (reference) read depth values per window.
        gc_content: Optional array of GC content fractions per window (0-1)
            for GC bias correction. If provided, the ratio is corrected
            using median normalization within GC bins.
        pseudocount: Small value added to avoid division by zero.

    Returns:
        Numpy array of log2 ratio values. If numpy is not available,
        returns a list of floats.

    Raises:
        ValueError: If arrays have different lengths.
    """
    if np is None:
        # Pure Python fallback
        if len(tumor_depth) != len(normal_depth):
            raise ValueError(
                f"Array length mismatch: tumor={len(tumor_depth)}, normal={len(normal_depth)}"
            )
        ratios: list[float] = []
        for t, n in zip(tumor_depth, normal_depth):
            ratio = (t + pseudocount) / (n + pseudocount)
            ratios.append(math.log2(ratio))
        return ratios

    tumor = np.asarray(tumor_depth, dtype=np.float64)
    normal = np.asarray(normal_depth, dtype=np.float64)

    if tumor.shape != normal.shape:
        raise ValueError(
            f"Array shape mismatch: tumor={tumor.shape}, normal={normal.shape}"
        )

    # Compute raw log2 ratio
    ratio = (tumor + pseudocount) / (normal + pseudocount)
    log2_ratio = np.log2(ratio)

    # GC correction if GC content provided
    if gc_content is not None:
        gc = np.asarray(gc_content, dtype=np.float64)
        if gc.shape != tumor.shape:
            raise ValueError(
                f"GC content shape mismatch: gc={gc.shape}, tumor={tumor.shape}"
            )
        log2_ratio = _gc_correct(log2_ratio, gc)

    # Center the log2 ratios around 0 (median normalization)
    median_ratio = np.median(log2_ratio[np.isfinite(log2_ratio)])
    log2_ratio -= median_ratio

    return log2_ratio


def _gc_correct(
    log2_ratios: Any,
    gc_content: Any,
    n_bins: int = 20,
) -> Any:
    """GC-content bias correction using binned median normalization.

    Divides the GC content range into bins, computes the median log2 ratio
    within each GC bin, and subtracts the bin-specific median to remove
    GC-dependent bias.

    Args:
        log2_ratios: Array of log2 ratio values.
        gc_content: Array of GC content fractions (0-1).
        n_bins: Number of GC bins for correction.

    Returns:
        GC-corrected log2 ratio array.
    """
    corrected = log2_ratios.copy()
    gc_bins = np.linspace(0, 1, n_bins + 1)

    for i in range(n_bins):
        mask = (gc_content >= gc_bins[i]) & (gc_content < gc_bins[i + 1])
        if np.sum(mask) > 0:
            bin_median = np.median(log2_ratios[mask])
            corrected[mask] -= bin_median

    # Recenter
    global_median = np.median(corrected[np.isfinite(corrected)])
    corrected -= global_median

    return corrected
