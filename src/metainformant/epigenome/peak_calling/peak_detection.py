"""Peak detection and calling for ChIP-seq and ATAC-seq data.

This module implements signal-based peak calling algorithms including
simple narrow peak calling with Poisson statistics, broad domain
detection for diffuse histone marks, summit refinement, peak merging
and filtering, FRiP quality metrics, and differential peak analysis.
"""

from __future__ import annotations

import math
import statistics
from typing import Any, Optional

from metainformant.core import logging

logger = logging.get_logger(__name__)

# Optional numpy for vectorized operations
try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None  # type: ignore[assignment]

# Optional scipy for statistical tests
try:
    from scipy import stats as scipy_stats

    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    scipy_stats = None  # type: ignore[assignment]


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _poisson_pvalue(observed: float, expected: float) -> float:
    """Compute one-sided Poisson p-value P(X >= observed | lambda=expected).

    Uses scipy when available, otherwise falls back to a pure-Python
    cumulative Poisson CDF computation.

    Args:
        observed: Observed count (signal value in a window).
        expected: Expected count (lambda from background estimate).

    Returns:
        One-sided p-value.
    """
    if expected <= 0:
        return 1.0 if observed <= 0 else 0.0

    if HAS_SCIPY:
        return float(1.0 - scipy_stats.poisson.cdf(int(observed) - 1, expected))

    # Pure-Python fallback: sum CDF up to observed - 1
    k_max = max(0, int(observed) - 1)
    cdf = 0.0
    log_lambda = math.log(expected)
    log_fact = 0.0  # log(0!)
    for k in range(k_max + 1):
        if k > 0:
            log_fact += math.log(k)
        cdf += math.exp(k * log_lambda - expected - log_fact)
    return max(0.0, 1.0 - cdf)


def _benjamini_hochberg(pvalues: list[float]) -> list[float]:
    """Apply Benjamini-Hochberg FDR correction to a list of p-values.

    Args:
        pvalues: Raw p-values.

    Returns:
        List of q-values (FDR-adjusted p-values) in the original order.
    """
    n = len(pvalues)
    if n == 0:
        return []

    # Sort indices by p-value
    indexed = sorted(enumerate(pvalues), key=lambda x: x[1])
    qvalues = [0.0] * n

    cummin = 1.0
    for rank_idx in range(n - 1, -1, -1):
        original_idx, pval = indexed[rank_idx]
        rank = rank_idx + 1  # 1-based rank
        adjusted = pval * n / rank
        cummin = min(cummin, adjusted)
        qvalues[original_idx] = min(cummin, 1.0)

    return qvalues


def _fold_enrichment(signal_val: float, background_val: float) -> float:
    """Compute fold enrichment of signal over background.

    Args:
        signal_val: Observed signal value.
        background_val: Background (control or local lambda) value.

    Returns:
        Fold enrichment, clamped to a minimum of 0.0.
    """
    if background_val <= 0:
        return float(signal_val) if signal_val > 0 else 0.0
    return max(0.0, signal_val / background_val)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def compute_local_lambda(
    signal: list[float],
    position: int,
    window_sizes: list[int] | None = None,
) -> float:
    """Estimate local background lambda using multiple window sizes.

    Computes the mean signal in windows of varying sizes centred on
    *position* and returns the maximum.  This mirrors MACS2's approach
    of taking the maximum local lambda across 1 kb, 5 kb, and 10 kb
    windows to avoid under-estimating background in regions with local
    enrichment.

    Args:
        signal: Signal values (one value per genomic position or bin).
        position: Centre position for the lambda estimate.
        window_sizes: Half-window sizes to evaluate.  Defaults to
            ``[500, 2500, 5000]`` (corresponding to 1 kb, 5 kb, 10 kb
            full windows).

    Returns:
        Maximum local lambda across all window sizes.

    Raises:
        ValueError: If signal is empty or position is out of range.
    """
    if not signal:
        raise ValueError("Signal list must not be empty")
    if position < 0 or position >= len(signal):
        raise ValueError(f"Position {position} is out of range [0, {len(signal) - 1}]")

    if window_sizes is None:
        window_sizes = [500, 2500, 5000]

    max_lambda = 0.0
    n = len(signal)

    for half_w in window_sizes:
        start = max(0, position - half_w)
        end = min(n, position + half_w + 1)
        window = signal[start:end]
        if window:
            window_mean = sum(window) / len(window)
            max_lambda = max(max_lambda, window_mean)

    # Floor: genome-wide mean as absolute minimum
    genome_mean = sum(signal) / len(signal)
    max_lambda = max(max_lambda, genome_mean)

    return max_lambda


def call_peaks_simple(
    signal: list[float],
    control: list[float] | None = None,
    threshold: float = 5.0,
    min_length: int = 100,
    merge_distance: int = 50,
) -> list[dict]:
    """Call narrow peaks from a signal track using Poisson enrichment.

    Scans the signal for positions that are significantly enriched over
    the local background (estimated from *control* or from the signal
    itself when no control is provided).  Enriched positions are grouped
    into contiguous peaks, filtered by minimum length, and optionally
    merged when they are within *merge_distance* of each other.

    Args:
        signal: Signal values (e.g. ChIP-seq read counts per bin).
        control: Optional control/input signal of the same length.
            When ``None`` the local lambda is estimated from *signal*.
        threshold: Minimum -log10(p-value) to call a position as
            enriched.  Default ``5.0`` corresponds to p < 1e-5.
        min_length: Minimum peak length (number of consecutive
            enriched positions) to retain.
        merge_distance: Maximum gap between enriched runs to merge
            into a single peak.

    Returns:
        List of peak dictionaries with keys: ``chrom``, ``start``,
        ``end``, ``summit``, ``score``, ``p_value``, ``fold_enrichment``,
        ``signal_value``.

    Raises:
        ValueError: If signal is empty or control length mismatches.
    """
    if not signal:
        raise ValueError("Signal list must not be empty")
    if control is not None and len(control) != len(signal):
        raise ValueError(f"Control length ({len(control)}) must match signal length " f"({len(signal)})")

    logger.info(
        "Calling narrow peaks: threshold=%.1f, min_length=%d, " "merge_distance=%d, signal_length=%d",
        threshold,
        min_length,
        merge_distance,
        len(signal),
    )

    n = len(signal)
    # Genome-wide background lambda from control or signal
    if control is not None:
        global_lambda = sum(control) / n if n > 0 else 1.0
    else:
        global_lambda = sum(signal) / n if n > 0 else 1.0

    # Phase 1: identify enriched positions
    enriched: list[bool] = []
    pvalues: list[float] = []

    for i in range(n):
        if control is not None:
            local_bg = compute_local_lambda(control, i, window_sizes=[500, 2500, 5000])
        else:
            local_bg = compute_local_lambda(signal, i, window_sizes=[500, 2500, 5000])

        bg = max(local_bg, global_lambda, 1e-10)
        pval = _poisson_pvalue(signal[i], bg)
        pvalues.append(pval)

        neg_log_p = -math.log10(pval) if pval > 0 else 300.0
        enriched.append(neg_log_p >= threshold)

    # Phase 2: group contiguous enriched positions into candidate peaks
    candidates: list[dict] = []
    in_peak = False
    peak_start = 0

    for i in range(n):
        if enriched[i] and not in_peak:
            peak_start = i
            in_peak = True
        elif not enriched[i] and in_peak:
            candidates.append({"start": peak_start, "end": i})
            in_peak = False
    if in_peak:
        candidates.append({"start": peak_start, "end": n})

    # Phase 3: merge nearby candidates
    if merge_distance > 0 and len(candidates) > 1:
        merged_candidates: list[dict] = [candidates[0]]
        for cand in candidates[1:]:
            prev = merged_candidates[-1]
            if cand["start"] - prev["end"] <= merge_distance:
                prev["end"] = cand["end"]
            else:
                merged_candidates.append(cand)
        candidates = merged_candidates

    # Phase 4: filter by minimum length and compute peak statistics
    peaks: list[dict] = []
    raw_pvals: list[float] = []

    for cand in candidates:
        length = cand["end"] - cand["start"]
        if length < min_length:
            continue

        region_signal = signal[cand["start"] : cand["end"]]
        max_signal = max(region_signal)
        summit_offset = region_signal.index(max_signal)
        summit = cand["start"] + summit_offset

        # Use minimum p-value in the peak region
        region_pvals = pvalues[cand["start"] : cand["end"]]
        best_pval = min(region_pvals)

        # Fold enrichment at summit
        if control is not None:
            bg_at_summit = max(control[summit], global_lambda, 1e-10)
        else:
            bg_at_summit = max(global_lambda, 1e-10)
        fold = _fold_enrichment(signal[summit], bg_at_summit)

        score = -math.log10(best_pval) if best_pval > 0 else 300.0

        peak = {
            "chrom": "chr_unknown",
            "start": cand["start"],
            "end": cand["end"],
            "summit": summit,
            "score": round(score, 4),
            "p_value": best_pval,
            "fold_enrichment": round(fold, 4),
            "signal_value": round(max_signal, 4),
        }
        peaks.append(peak)
        raw_pvals.append(best_pval)

    # Phase 5: BH FDR correction
    if raw_pvals:
        qvalues = _benjamini_hochberg(raw_pvals)
        for i, peak in enumerate(peaks):
            peak["q_value"] = qvalues[i]

    logger.info("Called %d narrow peaks from %d positions", len(peaks), n)
    return peaks


def call_peaks_broad(
    signal: list[float],
    control: list[float] | None = None,
    p_threshold: float = 0.1,
    broad_cutoff: float = 0.1,
    min_length: int = 200,
) -> list[dict]:
    """Call broad peaks / domains for diffuse histone marks.

    Uses a two-pass approach: first identifies seed peaks at a stringent
    threshold (``p_threshold / 10``), then extends each seed at the
    relaxed ``broad_cutoff`` threshold to capture the full domain.
    Suitable for marks like H3K27me3 and H3K36me3.

    Args:
        signal: Signal values (e.g. read counts per bin).
        control: Optional control signal of the same length.
        p_threshold: Relaxed p-value threshold for domain extension.
        broad_cutoff: Secondary p-value cutoff for broad region
            boundary determination.
        min_length: Minimum broad peak length to retain.

    Returns:
        List of broad peak dictionaries with keys: ``chrom``, ``start``,
        ``end``, ``summit``, ``score``, ``p_value``, ``fold_enrichment``,
        ``signal_value``, ``peak_type``.

    Raises:
        ValueError: If signal is empty or control length mismatches.
    """
    if not signal:
        raise ValueError("Signal list must not be empty")
    if control is not None and len(control) != len(signal):
        raise ValueError(f"Control length ({len(control)}) must match signal length " f"({len(signal)})")

    logger.info(
        "Calling broad peaks: p_threshold=%.3f, broad_cutoff=%.3f, " "min_length=%d, signal_length=%d",
        p_threshold,
        broad_cutoff,
        min_length,
        len(signal),
    )

    n = len(signal)
    if control is not None:
        global_lambda = sum(control) / n if n > 0 else 1.0
    else:
        global_lambda = sum(signal) / n if n > 0 else 1.0

    # Compute p-values at every position
    pvalues: list[float] = []
    for i in range(n):
        if control is not None:
            local_bg = compute_local_lambda(control, i, window_sizes=[500, 2500, 5000])
        else:
            local_bg = compute_local_lambda(signal, i, window_sizes=[500, 2500, 5000])
        bg = max(local_bg, global_lambda, 1e-10)
        pval = _poisson_pvalue(signal[i], bg)
        pvalues.append(pval)

    # Stringent threshold for seed peaks
    stringent_threshold = p_threshold / 10.0

    # Pass 1: find seed peaks at stringent threshold
    seeds: list[dict] = []
    in_seed = False
    seed_start = 0

    for i in range(n):
        if pvalues[i] <= stringent_threshold and not in_seed:
            seed_start = i
            in_seed = True
        elif (pvalues[i] > stringent_threshold or i == n - 1) and in_seed:
            seed_end = i + 1 if (i == n - 1 and pvalues[i] <= stringent_threshold) else i
            seeds.append({"start": seed_start, "end": seed_end})
            in_seed = False

    logger.debug("Found %d seed regions at stringent threshold", len(seeds))

    # Pass 2: extend each seed at the relaxed threshold
    broad_peaks: list[dict] = []
    raw_pvals: list[float] = []

    for seed in seeds:
        # Extend left
        ext_start = seed["start"]
        while ext_start > 0 and pvalues[ext_start - 1] <= broad_cutoff:
            ext_start -= 1

        # Extend right
        ext_end = seed["end"]
        while ext_end < n and pvalues[ext_end] <= broad_cutoff:
            ext_end += 1

        length = ext_end - ext_start
        if length < min_length:
            continue

        region_signal = signal[ext_start:ext_end]
        max_signal = max(region_signal)
        summit_offset = region_signal.index(max_signal)
        summit = ext_start + summit_offset

        region_pvals = pvalues[ext_start:ext_end]
        best_pval = min(region_pvals)

        if control is not None:
            bg_at_summit = max(control[summit], global_lambda, 1e-10)
        else:
            bg_at_summit = max(global_lambda, 1e-10)
        fold = _fold_enrichment(signal[summit], bg_at_summit)

        score = -math.log10(best_pval) if best_pval > 0 else 300.0

        peak = {
            "chrom": "chr_unknown",
            "start": ext_start,
            "end": ext_end,
            "summit": summit,
            "score": round(score, 4),
            "p_value": best_pval,
            "fold_enrichment": round(fold, 4),
            "signal_value": round(max_signal, 4),
            "peak_type": "broad",
        }
        broad_peaks.append(peak)
        raw_pvals.append(best_pval)

    # Merge overlapping broad peaks
    if len(broad_peaks) > 1:
        broad_peaks = merge_peaks(broad_peaks, distance=0)
        # Recompute raw_pvals after merge
        raw_pvals = [p["p_value"] for p in broad_peaks]

    # BH FDR correction
    if raw_pvals:
        qvalues = _benjamini_hochberg(raw_pvals)
        for i, peak in enumerate(broad_peaks):
            peak["q_value"] = qvalues[i]

    logger.info("Called %d broad peaks from %d positions", len(broad_peaks), n)
    return broad_peaks


def peak_summit_refinement(
    signal: list[float],
    peak: dict,
    window: int = 50,
) -> dict:
    """Refine peak summit position using parabolic interpolation.

    Fits a parabola to the signal values around the current summit
    and adjusts the summit to the vertex of the parabola.  This yields
    sub-bin resolution for the summit location.

    Args:
        signal: Full signal array.
        peak: Peak dictionary with at least ``start``, ``end``, and
            ``summit`` keys.
        window: Half-window around the summit to use for fitting.

    Returns:
        A new peak dictionary with updated ``summit`` and
        ``refined_summit_position`` (float, sub-bin resolution).

    Raises:
        ValueError: If signal is empty or peak summit is out of range.
    """
    if not signal:
        raise ValueError("Signal list must not be empty")

    summit = peak.get("summit", (peak["start"] + peak["end"]) // 2)
    if summit < 0 or summit >= len(signal):
        raise ValueError(f"Summit position {summit} is out of signal range " f"[0, {len(signal) - 1}]")

    # Restrict window to peak boundaries and signal bounds
    win_start = max(peak["start"], summit - window, 0)
    win_end = min(peak["end"], summit + window + 1, len(signal))

    region = signal[win_start:win_end]
    if len(region) < 3:
        # Not enough points for parabolic fit; return unchanged
        refined = dict(peak)
        refined["refined_summit_position"] = float(summit)
        return refined

    # Find the maximum in the window
    local_max_val = max(region)
    local_max_idx = region.index(local_max_val)
    abs_max_idx = win_start + local_max_idx

    # Parabolic interpolation requires the point and its two neighbours
    refined_position = float(abs_max_idx)

    if 0 < local_max_idx < len(region) - 1:
        y_left = region[local_max_idx - 1]
        y_centre = region[local_max_idx]
        y_right = region[local_max_idx + 1]

        denominator = 2.0 * (2.0 * y_centre - y_left - y_right)
        if abs(denominator) > 1e-12:
            offset = (y_left - y_right) / denominator
            refined_position = abs_max_idx + offset

    refined = dict(peak)
    refined["summit"] = int(round(refined_position))
    refined["refined_summit_position"] = round(refined_position, 4)
    # Update signal value at refined summit
    clamped_summit = max(0, min(len(signal) - 1, refined["summit"]))
    refined["signal_value"] = round(signal[clamped_summit], 4)

    logger.debug(
        "Refined summit from %d to %.2f (integer %d)",
        summit,
        refined_position,
        refined["summit"],
    )
    return refined


def merge_peaks(
    peaks: list[dict],
    distance: int = 0,
) -> list[dict]:
    """Merge overlapping or nearby peaks, keeping the best summit.

    Peaks are sorted by start position, then consecutive peaks whose
    gap is at most *distance* are merged.  The merged peak retains the
    summit with the highest signal value and the best (lowest) p-value.

    Args:
        peaks: List of peak dictionaries.
        distance: Maximum gap between peaks to merge (0 = overlapping
            only).

    Returns:
        List of merged peak dictionaries sorted by start position.
    """
    if not peaks:
        return []

    logger.debug("Merging %d peaks with distance=%d", len(peaks), distance)

    # Sort by start position
    sorted_peaks = sorted(peaks, key=lambda p: p["start"])

    merged: list[dict] = [dict(sorted_peaks[0])]

    for peak in sorted_peaks[1:]:
        prev = merged[-1]
        if peak["start"] <= prev["end"] + distance:
            # Merge: extend end
            prev["end"] = max(prev["end"], peak["end"])

            # Keep summit with highest signal
            if peak.get("signal_value", 0) > prev.get("signal_value", 0):
                prev["summit"] = peak.get("summit", prev.get("summit"))
                prev["signal_value"] = peak.get("signal_value", prev.get("signal_value", 0))

            # Keep best (lowest) p-value
            prev_pval = prev.get("p_value", 1.0)
            peak_pval = peak.get("p_value", 1.0)
            if peak_pval < prev_pval:
                prev["p_value"] = peak_pval

            # Keep best score (highest)
            prev["score"] = max(prev.get("score", 0), peak.get("score", 0))

            # Keep best fold enrichment (highest)
            prev["fold_enrichment"] = max(prev.get("fold_enrichment", 0), peak.get("fold_enrichment", 0))
        else:
            merged.append(dict(peak))

    logger.debug("Merged %d peaks into %d", len(peaks), len(merged))
    return merged


def filter_peaks(
    peaks: list[dict],
    min_fold: float = 2.0,
    max_qvalue: float = 0.05,
    blacklist_regions: list[dict] | None = None,
) -> list[dict]:
    """Filter peaks by fold enrichment, q-value, and blacklist exclusion.

    Applies sequential filters: first by fold enrichment, then by
    q-value (or p-value when q-value is absent), and finally excludes
    peaks that overlap blacklisted genomic regions.

    Args:
        peaks: List of peak dictionaries.
        min_fold: Minimum fold enrichment to retain.
        max_qvalue: Maximum q-value (FDR) to retain.
        blacklist_regions: Optional list of region dictionaries with
            ``chrom``, ``start``, ``end`` keys.  Peaks overlapping any
            blacklist region are removed.

    Returns:
        Filtered list of peak dictionaries.
    """
    if not peaks:
        return []

    logger.info(
        "Filtering %d peaks: min_fold=%.2f, max_qvalue=%.4f, " "blacklist_regions=%d",
        len(peaks),
        min_fold,
        max_qvalue,
        len(blacklist_regions) if blacklist_regions else 0,
    )

    filtered = list(peaks)

    # Filter by fold enrichment
    filtered = [p for p in filtered if p.get("fold_enrichment", 0) >= min_fold]
    logger.debug("%d peaks pass fold enrichment filter", len(filtered))

    # Filter by q-value (fall back to p-value if q-value absent)
    filtered = [p for p in filtered if p.get("q_value", p.get("p_value", 1.0)) <= max_qvalue]
    logger.debug("%d peaks pass q-value filter", len(filtered))

    # Filter by blacklist
    if blacklist_regions:

        def _overlaps_blacklist(peak: dict) -> bool:
            peak_chrom = peak.get("chrom", "")
            for bl in blacklist_regions:  # type: ignore[union-attr]
                if bl.get("chrom", "") != peak_chrom:
                    continue
                # Check overlap
                if peak["start"] < bl["end"] and peak["end"] > bl["start"]:
                    return True
            return False

        before = len(filtered)
        filtered = [p for p in filtered if not _overlaps_blacklist(p)]
        logger.debug(
            "%d peaks removed by blacklist, %d remain",
            before - len(filtered),
            len(filtered),
        )

    logger.info("Filtered to %d peaks", len(filtered))
    return filtered


def compute_frip(reads_in_peaks: int, total_reads: int) -> dict:
    """Compute Fraction of Reads in Peaks (FRiP) quality metric.

    FRiP is a standard quality metric for ChIP-seq experiments.
    ENCODE guidelines suggest FRiP >= 0.01 for a usable experiment,
    and FRiP >= 0.05 for a good experiment.

    Args:
        reads_in_peaks: Number of reads falling within called peaks.
        total_reads: Total number of mapped reads.

    Returns:
        Dictionary with ``frip``, ``reads_in_peaks``, ``total_reads``,
        ``reads_outside_peaks``, and ``quality_assessment``.

    Raises:
        ValueError: If total_reads is zero or negative, or
            reads_in_peaks exceeds total_reads.
    """
    if total_reads <= 0:
        raise ValueError(f"total_reads must be positive, got {total_reads}")
    if reads_in_peaks < 0:
        raise ValueError(f"reads_in_peaks must be non-negative, got {reads_in_peaks}")
    if reads_in_peaks > total_reads:
        raise ValueError(f"reads_in_peaks ({reads_in_peaks}) cannot exceed " f"total_reads ({total_reads})")

    frip = reads_in_peaks / total_reads
    reads_outside = total_reads - reads_in_peaks

    # Quality assessment per ENCODE guidelines
    if frip >= 0.05:
        quality = "good"
    elif frip >= 0.01:
        quality = "acceptable"
    else:
        quality = "low"

    result = {
        "frip": round(frip, 6),
        "reads_in_peaks": reads_in_peaks,
        "total_reads": total_reads,
        "reads_outside_peaks": reads_outside,
        "percent_in_peaks": round(frip * 100, 4),
        "quality_assessment": quality,
    }

    logger.info(
        "FRiP = %.4f (%s): %d / %d reads in peaks",
        frip,
        quality,
        reads_in_peaks,
        total_reads,
    )
    return result


def differential_peaks(
    peaks_a: list[dict],
    peaks_b: list[dict],
    signals_a: list[float],
    signals_b: list[float],
) -> list[dict]:
    """Identify differential peaks between two conditions.

    Compares peak sets from two conditions, identifies peaks unique to
    each condition, and for shared peaks computes fold change and a
    statistical test (Welch's t-test or simple fold-change threshold)
    on the signal values.

    Args:
        peaks_a: Peaks from condition A.
        peaks_b: Peaks from condition B.
        signals_a: Full signal track for condition A.
        signals_b: Full signal track for condition B.

    Returns:
        List of differential peak dictionaries with keys: ``chrom``,
        ``start``, ``end``, ``log2_fold_change``, ``p_value``,
        ``direction`` ("up_in_A", "up_in_B", or "unchanged"),
        ``status`` ("shared", "unique_A", or "unique_B"),
        ``mean_signal_a``, ``mean_signal_b``.

    Raises:
        ValueError: If signal lengths differ.
    """
    if len(signals_a) != len(signals_b):
        raise ValueError(f"Signal lengths must match: A={len(signals_a)}, B={len(signals_b)}")

    logger.info(
        "Computing differential peaks: %d peaks in A, %d peaks in B",
        len(peaks_a),
        len(peaks_b),
    )

    results: list[dict] = []

    # Build interval index for peaks_b for overlap detection
    def _peaks_by_region(peaks: list[dict]) -> dict[str, list[dict]]:
        """Index peaks by chromosome."""
        index: dict[str, list[dict]] = {}
        for p in peaks:
            chrom = p.get("chrom", "chr_unknown")
            if chrom not in index:
                index[chrom] = []
            index[chrom].append(p)
        return index

    index_b = _peaks_by_region(peaks_b)
    matched_b: set[int] = set()  # Track which B peaks are matched

    # For each peak in A, find overlapping peaks in B
    for peak_a in peaks_a:
        chrom_a = peak_a.get("chrom", "chr_unknown")
        a_start = peak_a["start"]
        a_end = peak_a["end"]

        # Extract mean signal for the peak region
        region_start = max(0, a_start)
        region_end = min(len(signals_a), a_end)
        if region_end > region_start:
            mean_a = sum(signals_a[region_start:region_end]) / (region_end - region_start)
            mean_b = sum(signals_b[region_start:region_end]) / (region_end - region_start)
        else:
            mean_a = 0.0
            mean_b = 0.0

        # Find overlapping peak in B
        overlap_found = False
        candidates_b = index_b.get(chrom_a, [])

        for idx_b, peak_b in enumerate(candidates_b):
            b_start = peak_b["start"]
            b_end = peak_b["end"]
            # Check overlap
            if a_start < b_end and a_end > b_start:
                overlap_found = True
                # Compute using the union region
                union_start = min(a_start, b_start)
                union_end = max(a_end, b_end)
                u_start = max(0, union_start)
                u_end = min(len(signals_a), union_end)

                if u_end > u_start:
                    mean_a = sum(signals_a[u_start:u_end]) / (u_end - u_start)
                    mean_b = sum(signals_b[u_start:u_end]) / (u_end - u_start)

                # Log2 fold change
                pseudo = 1e-10
                log2fc = math.log2((mean_a + pseudo) / (mean_b + pseudo))

                # Simple statistical test on region values
                vals_a = signals_a[u_start:u_end]
                vals_b = signals_b[u_start:u_end]
                pval = _compute_differential_pvalue(vals_a, vals_b)

                if abs(log2fc) >= 1.0 and pval < 0.05:
                    direction = "up_in_A" if log2fc > 0 else "up_in_B"
                else:
                    direction = "unchanged"

                results.append(
                    {
                        "chrom": chrom_a,
                        "start": union_start,
                        "end": union_end,
                        "log2_fold_change": round(log2fc, 4),
                        "p_value": pval,
                        "direction": direction,
                        "status": "shared",
                        "mean_signal_a": round(mean_a, 4),
                        "mean_signal_b": round(mean_b, 4),
                    }
                )

                # Mark this B peak as matched
                # Use id() for tracking since peaks are dicts
                matched_b.add(id(peak_b))
                break

        if not overlap_found:
            # Peak unique to A
            results.append(
                {
                    "chrom": chrom_a,
                    "start": a_start,
                    "end": a_end,
                    "log2_fold_change": float("inf") if mean_a > 0 else 0.0,
                    "p_value": 0.0,
                    "direction": "up_in_A",
                    "status": "unique_A",
                    "mean_signal_a": round(mean_a, 4),
                    "mean_signal_b": round(mean_b, 4),
                }
            )

    # Find B peaks not matched to any A peak
    for chrom_b, peaks_in_chrom in index_b.items():
        for peak_b in peaks_in_chrom:
            if id(peak_b) not in matched_b:
                b_start = peak_b["start"]
                b_end = peak_b["end"]
                region_start = max(0, b_start)
                region_end = min(len(signals_b), b_end)

                if region_end > region_start:
                    mean_a = sum(signals_a[region_start:region_end]) / (region_end - region_start)
                    mean_b = sum(signals_b[region_start:region_end]) / (region_end - region_start)
                else:
                    mean_a = 0.0
                    mean_b = 0.0

                results.append(
                    {
                        "chrom": chrom_b,
                        "start": b_start,
                        "end": b_end,
                        "log2_fold_change": float("-inf") if mean_b > 0 else 0.0,
                        "p_value": 0.0,
                        "direction": "up_in_B",
                        "status": "unique_B",
                        "mean_signal_a": round(mean_a, 4),
                        "mean_signal_b": round(mean_b, 4),
                    }
                )

    logger.info(
        "Identified %d differential peaks: %d shared, %d unique_A, %d unique_B",
        len(results),
        sum(1 for r in results if r["status"] == "shared"),
        sum(1 for r in results if r["status"] == "unique_A"),
        sum(1 for r in results if r["status"] == "unique_B"),
    )
    return results


def _compute_differential_pvalue(
    vals_a: list[float],
    vals_b: list[float],
) -> float:
    """Compute p-value for differential signal between two conditions.

    Uses Welch's t-test via scipy when available, otherwise falls back
    to a simple z-test approximation.

    Args:
        vals_a: Signal values from condition A.
        vals_b: Signal values from condition B.

    Returns:
        Two-sided p-value.
    """
    if len(vals_a) < 2 or len(vals_b) < 2:
        return 1.0

    if HAS_SCIPY:
        stat, pval = scipy_stats.ttest_ind(vals_a, vals_b, equal_var=False)
        return float(pval) if not math.isnan(pval) else 1.0

    # Pure-Python fallback: z-test approximation
    mean_a = sum(vals_a) / len(vals_a)
    mean_b = sum(vals_b) / len(vals_b)

    if len(vals_a) > 1:
        var_a = sum((x - mean_a) ** 2 for x in vals_a) / (len(vals_a) - 1)
    else:
        var_a = 0.0

    if len(vals_b) > 1:
        var_b = sum((x - mean_b) ** 2 for x in vals_b) / (len(vals_b) - 1)
    else:
        var_b = 0.0

    se = math.sqrt(var_a / len(vals_a) + var_b / len(vals_b))
    if se < 1e-12:
        return 1.0

    z = abs(mean_a - mean_b) / se
    # Approximate two-sided p-value from z-score
    pval = 2.0 * (1.0 - 0.5 * (1.0 + math.erf(z / math.sqrt(2.0))))
    return max(0.0, min(1.0, pval))
