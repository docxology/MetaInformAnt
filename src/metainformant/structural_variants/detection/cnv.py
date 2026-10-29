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
import statistics
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

from metainformant.core.utils.logging import get_logger

np: Any  # None when numpy is unavailable (optional dependency)
try:
    import numpy as np
    from numpy.typing import NDArray
except ImportError:
    np = None
    if not TYPE_CHECKING:
        NDArray = Any  # runtime-only stand-in; type checkers keep the real numpy alias

scipy_stats: Any  # None when scipy is unavailable (optional dependency)
try:
    from scipy import stats as scipy_stats
except ImportError:
    scipy_stats = None

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
            values. Each value represents the mean depth in a genomic window and
            must be finite and non-negative (zero is a bin with no coverage).
            Can be numpy arrays or plain lists of floats. Log2 ratios are not
            depths: pass them to :func:`segment_coverage` instead.
        window_size: Size of each depth bin in base pairs.
        method: Detection method. Currently supports 'segmentation' (CBS).
        significance: P-value threshold for CBS change-point detection. The
            environment variable ``SV_CNV_SIGNIFICANCE``, when set, overrides it.
        ploidy: Expected ploidy of the organism (default 2 for diploid).
        min_segment_bins: Minimum number of bins for a valid segment.
        merge_distance: Maximum gap (in bp) to merge adjacent same-state segments.

    Returns:
        Dictionary mapping chromosome names to CNVResult objects containing
        detected segments with copy number states.

    Raises:
        ValueError: If depth_data is empty, method is unsupported, or a chromosome
            that would be analysed holds negative or non-finite depth values.
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

    # Check every chromosome that will be analysed before analysing any of them.
    depth_arrays: dict[str, Any] = {}
    for chrom, depths in depth_data.items():
        depth_array = np.asarray(depths, dtype=np.float64)

        if len(depth_array) < min_segment_bins:
            logger.warning(f"Chromosome {chrom}: too few bins ({len(depth_array)}), skipping")
            continue

        _check_depths(chrom, depth_array)
        depth_arrays[chrom] = depth_array

    for chrom, depth_array in depth_arrays.items():
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
            confidence = _calculate_segment_confidence(log2_ratios[seg_start:seg_end], seg_mean, state)
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


def _check_depths(chrom: str, depth_array: Any) -> None:
    """Raise a ValueError naming the problem if the values cannot be read-depth counts.

    Args:
        chrom: Chromosome name (for the message).
        depth_array: Depth values as a float array.

    Raises:
        ValueError: If any value is negative or not finite.
    """
    n_negative = int(np.count_nonzero(depth_array < 0))
    n_nonfinite = int(np.count_nonzero(~np.isfinite(depth_array)))
    if n_negative == 0 and n_nonfinite == 0:
        return
    raise ValueError(
        f"depth_data[{chrom!r}] has {n_negative} negative and {n_nonfinite} non-finite depth value(s) "
        f"out of {len(depth_array)}: read depths must be finite and >= 0 (0 is an uncovered bin). "
        "Negative values usually mean log2 ratios were passed instead of depths; give "
        "detect_cnv_from_depth the depths, or run segment_coverage() and call_cnv_states() on the log2 ratios."
    )


def segment_coverage(
    coverage_array: Any,
    significance: float = 0.01,
    max_iterations: int = 100,
    min_width: int = 2,
) -> list[tuple[int, int, float]]:
    """Circular Binary Segmentation (CBS) algorithm for coverage segmentation.

    Recursively finds change-points in a coverage array by scoring the
    difference in means between every circular arc of the data and the rest,
    splitting at the best-scoring arc, and recursing on the pieces. The
    p-value of the best arc is an asymptotic approximation (see
    :func:`_cbs_pvalue`); no permutation test is run.

    Every arc of at most half the data is scored, so the work grows with the
    square of the segment length; the scan is vectorised over all arcs of one
    length, which makes it about two orders of magnitude faster than scoring
    the arcs one at a time (and gives exactly the same segmentation).

    Args:
        coverage_array: Array of log2 ratio values to segment. Can be a numpy
            array or list of floats.
        significance: P-value threshold for accepting a change-point. Lower
            values produce fewer, more confident segments.
        max_iterations: Maximum recursion depth to prevent infinite loops.
        min_width: Minimum segment width in bins (at least 1).

    Returns:
        List of (start, end, mean) tuples representing segments, where start
        is inclusive, end is exclusive, and mean is the segment mean.

    Raises:
        ValueError: If min_width is less than 1.
    """
    if np is None:
        raise ImportError("numpy is required for CBS segmentation. Install with: uv pip install numpy")
    if min_width < 1:
        raise ValueError(f"min_width must be at least 1, got {min_width}")

    data = np.asarray(coverage_array, dtype=np.float64)
    n = len(data)

    if n < 2 * min_width:
        return [(0, n, float(np.mean(data)))]

    segments: list[tuple[int, int, float]] = []
    _cbs_recursive(data, 0, n, segments, significance, max_iterations, min_width, depth=0)

    # Sort by start position
    segments.sort(key=lambda x: x[0])
    return segments


def _cbs_best_arc(segment_data: Any, std_dev: float, min_width: int) -> tuple[float, int, int]:
    """Find the arc whose mean differs most from the rest of the segment.

    Scores every arc ``[i, j)`` with ``min_width <= j - i <= n // 2`` (and at least
    ``min_width`` bins outside it) by the standardised difference between the mean
    inside and outside the arc,
    ``|mean_in - mean_out| * sqrt(len_in * len_out / n) / std_dev``. A change point
    is the arc that touches an end of the segment, so the single-split scan of plain
    binary segmentation is the special case of this one.

    The arcs of one length are scored together with array arithmetic that performs
    the same floating-point operations in the same order as scoring them one at a
    time, so the statistics are bit-for-bit the same. The first arc in
    ``(i, j)`` order wins a tie.

    Args:
        segment_data: The segment's values (at least ``2 * min_width`` of them).
        std_dev: Sample standard deviation of the segment (positive).
        min_width: Minimum arc length, and minimum length outside the arc.

    Returns:
        ``(statistic, i, j)`` of the best arc, or ``(0.0, -1, -1)`` if no arc scores above zero.
    """
    n = len(segment_data)
    cumsum = np.cumsum(segment_data - np.mean(segment_data))
    cumsum = np.insert(cumsum, 0, 0.0)  # prepend 0
    total = cumsum[n]

    best_stat = 0.0
    best_i = -1
    best_j = -1
    for arc_length in range(min_width, n // 2 + 1):
        rest_length = n - arc_length

        # Arc (i, i + arc_length) for every i: its sum, and the standardised mean difference.
        arc_sum = cumsum[arc_length:] - cumsum[: n + 1 - arc_length]
        arc_mean_diff = arc_sum / arc_length - (total - arc_sum) / rest_length
        stats = np.abs(arc_mean_diff) * math.sqrt(arc_length * rest_length / n) / std_dev

        i = int(np.argmax(stats))  # first arc of this length with the largest statistic
        stat = float(stats[i])
        if stat > best_stat or (stat == best_stat and 0 <= best_i and i < best_i):
            best_stat = stat
            best_i = i
            best_j = i + arc_length

    return best_stat, best_i, best_j


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

    segment_data = data[start:end]

    # Variance of the segment
    variance = np.var(segment_data, ddof=1) if n > 1 else 1.0
    if not math.isfinite(variance) or variance < 1e-12:
        # Constant (or non-finite) segment, no change-point
        segments.append((start, end, float(np.mean(segment_data))))
        return

    std_dev = math.sqrt(variance)

    # Circular scan: the arc (i, j) whose mean differs most from the rest. Splitting at a
    # single change-point is the case where the arc touches an end of the segment.
    best_t_stat, best_i, best_j = _cbs_best_arc(segment_data, std_dev, min_width)

    # Assess significance with an asymptotic approximation of the null distribution of the maximum statistic
    # (a Gumbel law for large n); this is not a permutation test, so no random numbers are involved.
    p_value = _cbs_pvalue(best_t_stat, n)

    if p_value > significance:
        # No significant change-point
        segments.append((start, end, float(np.mean(segment_data))))
        return

    # Split into up to three parts: [start, start+i), [start+i, start+j), [start+j, end)
    abs_i = start + best_i
    abs_j = start + best_j

    if best_i > 0:
        _cbs_recursive(data, start, abs_i, segments, significance, max_iterations, min_width, depth + 1)
    _cbs_recursive(data, abs_i, abs_j, segments, significance, max_iterations, min_width, depth + 1)
    if abs_j < end:
        _cbs_recursive(data, abs_j, end, segments, significance, max_iterations, min_width, depth + 1)


def _cbs_pvalue(statistic: float, n: int) -> float:
    """Approximate p-value for CBS test statistic.

    Uses an asymptotic approximation in the spirit of Olshen et al. (2004):
    the maximum of the standardized CBS statistic is taken to follow a Gumbel
    distribution with parameters depending on the segment length. It is an
    approximation, not the permutation p-value of the original method.

    Args:
        statistic: CBS test statistic value.
        n: Number of data points in the segment.

    Returns:
        Approximate p-value. A segment of fewer than 4 points is too short for
        the approximation: 1.0 is returned, which means "no change point can be
        claimed here", not a measured probability.
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
            as a deletion (default -0.3, approximately CN=1 for diploid). The
            environment variable ``SV_CNV_DEL_THRESHOLD``, when set, overrides it.
        dup_threshold: Log2 ratio threshold above which a segment is called
            as a duplication (default 0.3, approximately CN=3 for diploid). The
            environment variable ``SV_CNV_DUP_THRESHOLD``, when set, overrides it.
        amp_threshold: Log2 ratio threshold for high-level amplification
            (default 1.0, approximately CN=4+ for diploid).
        homodel_threshold: Log2 ratio threshold for homozygous deletion
            (default -1.5, approximately CN=0 for diploid).

    Returns:
        List of CNVSegment objects with state and cn fields populated. Tuples carry
        no chromosome and positions in bins, so their segments have
        ``chrom == "unknown"`` and ``start``/``end`` equal to the bin indices;
        ``confidence`` is 0.0 for every segment (it needs the per-bin values, which
        are not given here).
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

        state, cn = _assign_cnv_state(mean_lr, ploidy, del_threshold, dup_threshold, amp_threshold, homodel_threshold)

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
    cn_float = ploidy * (2**mean_log2ratio)
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
        if seg.chrom == prev.chrom and seg.state == prev.state and (seg.start - prev.end) <= max_gap:
            # Merge: extend the previous segment
            total_bins = prev.n_bins + seg.n_bins
            weighted_mean = (prev.mean_log2ratio * prev.n_bins + seg.mean_log2ratio * seg.n_bins) / total_bins

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
    correction by binned median normalization (the median of each GC bin is
    subtracted; there is no LOESS or other smooth fit), and centers the ratios
    on their median.

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
        # Pure Python fallback: same algorithm as the numpy path below
        if len(tumor_depth) != len(normal_depth):
            raise ValueError(f"Array length mismatch: tumor={len(tumor_depth)}, normal={len(normal_depth)}")
        ratios: list[float] = []
        for t, n in zip(tumor_depth, normal_depth):
            ratio = (t + pseudocount) / (n + pseudocount)
            ratios.append(math.log2(ratio))
        if gc_content is not None:
            if len(gc_content) != len(ratios):
                raise ValueError(f"GC content length mismatch: gc={len(gc_content)}, tumor={len(ratios)}")
            ratios = _gc_correct_py(ratios, list(gc_content))
        return _center_on_median(ratios)

    tumor = np.asarray(tumor_depth, dtype=np.float64)
    normal = np.asarray(normal_depth, dtype=np.float64)

    if tumor.shape != normal.shape:
        raise ValueError(f"Array shape mismatch: tumor={tumor.shape}, normal={normal.shape}")

    # Compute raw log2 ratio
    ratio = (tumor + pseudocount) / (normal + pseudocount)
    log2_ratio = np.log2(ratio)

    # GC correction if GC content provided
    if gc_content is not None:
        gc = np.asarray(gc_content, dtype=np.float64)
        if gc.shape != tumor.shape:
            raise ValueError(f"GC content shape mismatch: gc={gc.shape}, tumor={tumor.shape}")
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


def _center_on_median(values: list[float]) -> list[float]:
    """Subtract the median of the finite values (pure-Python median centering).

    Args:
        values: Log2 ratio values.

    Returns:
        Values shifted so that the median of the finite entries is zero.
    """
    finite = [v for v in values if math.isfinite(v)]
    center = statistics.median(finite) if finite else math.nan
    return [v - center for v in values]


def _gc_correct_py(log2_ratios: list[float], gc_content: list[float], n_bins: int = 20) -> list[float]:
    """Pure-Python twin of :func:`_gc_correct` (binned median normalization).

    Args:
        log2_ratios: Log2 ratio values.
        gc_content: GC content fractions (0-1), one per log2 ratio.
        n_bins: Number of GC bins for correction.

    Returns:
        GC-corrected, median-centered log2 ratios.
    """
    corrected = list(log2_ratios)
    # Bin edges identical to numpy.linspace(0, 1, n_bins + 1)
    step = 1.0 / n_bins
    edges = [i * step for i in range(n_bins)] + [1.0]

    for i in range(n_bins):
        members = [j for j, gc in enumerate(gc_content) if edges[i] <= gc < edges[i + 1]]
        if members:
            bin_median = statistics.median(log2_ratios[j] for j in members)
            for j in members:
                corrected[j] -= bin_median

    return _center_on_median(corrected)
