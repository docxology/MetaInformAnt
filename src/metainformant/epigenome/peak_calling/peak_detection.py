"""Peak detection and calling for ChIP-seq and ATAC-seq data.

This module implements signal-based peak calling algorithms including
simple narrow peak calling with Poisson statistics, broad domain
detection for diffuse histone marks, summit refinement, peak merging
and filtering, FRiP quality metrics, and differential peak analysis.
"""

from __future__ import annotations

import math
import operator
from bisect import bisect_left
from collections.abc import Mapping, Sequence
from typing import Any, cast

from metainformant.core.utils import logging
from metainformant.epigenome.assays.chipseq import _poisson_upper_tail

logger = logging.get_logger(__name__)


# Optional scipy for statistical tests
try:
    from scipy import stats as scipy_stats

    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    scipy_stats = cast(Any, None)  # optional dependency; callers guard on HAS_SCIPY


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _poisson_pvalue(observed: float, expected: float) -> float:
    """Compute one-sided Poisson p-value P(X >= observed | lambda=expected).

    The upper tail is summed directly in log space (the same routine the ChIP-seq enrichment uses), with or without
    SciPy installed. The former ``1 - cdf`` could not represent a p-value below about 1e-16: a count of 20 over a
    background of 1 (p = 1.6e-19) came out as exactly 0.0, so every strong peak tied at p = 0 and a score of 300,
    and their q-values were all 0 whatever the evidence. ``observed`` is truncated to an integer count.

    Args:
        observed: Observed count (signal value in a window).
        expected: Expected count (lambda from background estimate).

    Returns:
        One-sided p-value (0.0 only where it underflows, or for a background of 0 with a positive count).
    """
    if expected <= 0:
        return 1.0 if observed <= 0 else 0.0

    return _poisson_upper_tail(int(observed), expected)[0]


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
    chrom: str = "chr_unknown",
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
        chrom: Chromosome the signal belongs to; it is written to the
            ``chrom`` key of every peak.  Peak coordinates are 0-based
            indices into *signal*, so they are genomic positions when the
            track starts at position 0 of *chrom*.  The default
            ``"chr_unknown"`` is only a placeholder for signals with no
            chromosome.

    Returns:
        List of peak dictionaries with keys: ``chrom``, ``start``,
        ``end``, ``summit``, ``score``, ``p_value``, ``fold_enrichment``,
        ``signal_value``, ``q_value``.

    Raises:
        ValueError: If signal is empty, control length mismatches or
            *chrom* is empty.
    """
    if not signal:
        raise ValueError("Signal list must not be empty")
    if not chrom:
        raise ValueError("chrom must not be empty")
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
            "chrom": chrom,
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
    chrom: str = "chr_unknown",
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
        chrom: Chromosome the signal belongs to; it is written to the
            ``chrom`` key of every peak.  Peak coordinates are 0-based
            indices into *signal*, so they are genomic positions when the
            track starts at position 0 of *chrom*.  The default
            ``"chr_unknown"`` is only a placeholder for signals with no
            chromosome.

    Returns:
        List of broad peak dictionaries with keys: ``chrom``, ``start``,
        ``end``, ``summit``, ``score``, ``p_value``, ``fold_enrichment``,
        ``signal_value``, ``peak_type``, ``q_value``.

    Raises:
        ValueError: If signal is empty, control length mismatches or
            *chrom* is empty.
    """
    if not signal:
        raise ValueError("Signal list must not be empty")
    if not chrom:
        raise ValueError("chrom must not be empty")
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

    if in_seed:
        # Seed opened at the final position and never closed by the loop
        seeds.append({"start": seed_start, "end": n})

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
            "chrom": chrom,
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

    Finds the highest bin within *window* of the current summit (inside the peak and the signal), fits the parabola
    through that bin and its two neighbours, and moves the summit to the vertex of the parabola.  This yields
    sub-bin resolution for the summit location, never more than half a bin from the highest bin.  When the highest
    bin sits at the edge of the searched region, or the three bins are collinear, the highest bin itself is the
    summit.

    Args:
        signal: Full signal array.
        peak: Peak dictionary with at least ``start``, ``end``, and
            ``summit`` keys.
        window: Half-window around the summit in which to look for the highest bin (the parabola itself always
            uses just that bin and its two neighbours).

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

        # Vertex of the parabola through (-1, y_left), (0, y_centre) and (1, y_right). The curvature is negative at
        # a maximum; the vertex lies towards the higher neighbour (the offset had the opposite sign, so a summit
        # was moved away from the stronger side)
        curvature = y_left - 2.0 * y_centre + y_right
        if abs(curvature) > 1e-12:
            offset = 0.5 * (y_left - y_right) / curvature
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
    """Merge overlapping or nearby peaks on the same chromosome, keeping the best summit.

    Peaks are grouped by their ``chrom`` key (``"chr_unknown"`` for a peak
    without one, as in the rest of this module) and sorted by start position
    within each chromosome; consecutive peaks whose gap is at most *distance*
    are then merged.  Peaks on different chromosomes are never merged.  The
    merged peak retains the summit with the highest signal value, the best
    (lowest) ``p_value`` and ``q_value`` of its members, and the highest
    ``score`` and ``fold_enrichment``; every other key comes from the first
    peak.  The q-value is the lowest of the members' q-values, so it stays
    consistent with the retained p-value; it is not re-estimated for the
    merged set.

    Args:
        peaks: List of peak dictionaries.
        distance: Maximum gap between peaks to merge (0 = overlapping
            only).

    Returns:
        List of merged peak dictionaries: chromosomes in order of first
        appearance, each sorted by start position.
    """
    if not peaks:
        return []

    logger.debug("Merging %d peaks with distance=%d", len(peaks), distance)

    # Group by chromosome; a merge never crosses chromosomes
    by_chrom: dict[str, list[dict]] = {}
    for peak in peaks:
        by_chrom.setdefault(peak.get("chrom", "chr_unknown"), []).append(peak)

    merged: list[dict] = []

    for chrom_peaks in by_chrom.values():
        # Sort by start position
        sorted_peaks = sorted(chrom_peaks, key=lambda p: p["start"])

        merged.append(dict(sorted_peaks[0]))

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

                # Keep best (lowest) q-value, so it stays in step with the retained p-value
                peak_qval = peak.get("q_value")
                if peak_qval is not None:
                    prev_qval = prev.get("q_value")
                    if prev_qval is None or peak_qval < prev_qval:
                        prev["q_value"] = peak_qval

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
    peaks that overlap blacklisted genomic regions.  A peak without
    ``fold_enrichment``, or with neither ``q_value`` nor ``p_value``,
    cannot pass the corresponding filter and is removed; a warning
    counts such peaks instead of dropping them silently.

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
    without_fold = sum(1 for p in filtered if "fold_enrichment" not in p)
    if without_fold:
        logger.warning(
            "%d of %d peaks have no fold_enrichment and are removed by the fold filter", without_fold, len(peaks)
        )
    filtered = [p for p in filtered if p.get("fold_enrichment", 0) >= min_fold]
    logger.debug("%d peaks pass fold enrichment filter", len(filtered))

    # Filter by q-value (fall back to p-value if q-value absent)
    without_significance = sum(1 for p in filtered if "q_value" not in p and "p_value" not in p)
    if without_significance:
        logger.warning(
            "%d of %d peaks have neither q_value nor p_value and are removed by the significance filter",
            without_significance,
            len(filtered),
        )
    filtered = [p for p in filtered if p.get("q_value", p.get("p_value", 1.0)) <= max_qvalue]
    logger.debug("%d peaks pass q-value filter", len(filtered))

    # Filter by blacklist
    if blacklist_regions:

        def _overlaps_blacklist(peak: dict) -> bool:
            peak_chrom = peak.get("chrom", "")
            for bl in blacklist_regions:
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


def _peak_extents(peaks: list[dict], label: str) -> list[tuple[str, int, int]]:
    """The ``(chrom, start, end)`` of every peak, checked.

    Args:
        peaks: Peak dictionaries with integer ``start`` and ``end`` (``chrom`` defaults to ``"chr_unknown"``).
        label: Name of the peak set, for error messages.

    Raises:
        ValueError: If a peak has no integer ``start`` and ``end``, starts below 0 or has no extent (``end <= start``).
            Such a peak used to give zero signal means and a significant-looking row instead of an error.
    """
    extents: list[tuple[str, int, int]] = []
    for position, peak in enumerate(peaks):
        try:
            start = operator.index(peak["start"])
            end = operator.index(peak["end"])
        except (KeyError, TypeError):
            raise ValueError(
                f"Peak {position} of {label} needs integer 'start' and 'end' values, got {peak!r}"
            ) from None
        if start < 0 or end <= start:
            raise ValueError(f"Peak {position} of {label} must satisfy 0 <= start < end, got start={start}, end={end}")
        extents.append((peak.get("chrom", "chr_unknown"), start, end))
    return extents


def _best_overlaps(
    queries: list[tuple[str, int, int]],
    targets: list[tuple[str, int, int]],
) -> list[int | None]:
    """For every query peak, the index of the target peak it shares the most bases with.

    Peaks overlap when they are on the same chromosome and share at least one base (intervals are half-open). A tie
    in the number of shared bases goes to the target that starts first, then ends first, then comes first in
    *targets*, so the answer does not depend on the order of the targets.

    Args:
        queries: ``(chrom, start, end)`` of the peaks to look up.
        targets: ``(chrom, start, end)`` of the peaks to look in.

    Returns:
        One entry per query: the index into *targets*, or ``None`` when the query overlaps no target.
    """
    grouped: dict[str, list[tuple[int, int, int]]] = {}
    for position, (chrom, start, end) in enumerate(targets):
        grouped.setdefault(chrom, []).append((start, end, position))

    # Per chromosome: the targets by start, their starts, and the running maximum of their ends. A target can only
    # overlap a query that starts before it ends, and once the running maximum end is at or before the query's
    # start no earlier target can overlap it, so the scan below stops there
    index: dict[str, tuple[list[tuple[int, int, int]], list[int], list[int]]] = {}
    for chrom, items in grouped.items():
        items.sort()
        running_end: list[int] = []
        highest = items[0][1]
        for _, end, _ in items:
            highest = max(highest, end)
            running_end.append(highest)
        index[chrom] = (items, [item[0] for item in items], running_end)

    best: list[int | None] = []
    for chrom, query_start, query_end in queries:
        choice: int | None = None
        entry = index.get(chrom)
        if entry is not None:
            items, starts, running_end = entry
            best_key: tuple[int, int, int, int] | None = None
            for k in range(bisect_left(starts, query_end) - 1, -1, -1):
                if running_end[k] <= query_start:
                    break
                start, end, position = items[k]
                shared = min(query_end, end) - max(query_start, start)
                if shared <= 0:
                    continue
                key = (-shared, start, end, position)
                if best_key is None or key < best_key:
                    best_key = key
                    choice = position
        best.append(choice)
    return best


def _signal_tracks(
    signals_a: Sequence[float] | Mapping[str, Sequence[float]],
    signals_b: Sequence[float] | Mapping[str, Sequence[float]],
    chromosomes: set[str],
) -> dict[str, tuple[Sequence[float], Sequence[float]]]:
    """The two signal tracks to read for each chromosome that has peaks.

    Raises:
        ValueError: If the tracks differ in length, only one of the two arguments is a mapping, a chromosome with
            peaks has no track, or one plain track would have to describe peaks on several chromosomes.
    """
    if isinstance(signals_a, Mapping) or isinstance(signals_b, Mapping):
        if not (isinstance(signals_a, Mapping) and isinstance(signals_b, Mapping)):
            raise ValueError("signals_a and signals_b must both be lists or both be dicts of chromosome -> list")
        if set(signals_a) != set(signals_b):
            raise ValueError(f"Signal chromosomes must match: A={sorted(signals_a)}, B={sorted(signals_b)}")
        for chrom in signals_a:
            if len(signals_a[chrom]) != len(signals_b[chrom]):
                raise ValueError(
                    f"Signal lengths must match for {chrom}: A={len(signals_a[chrom])}, B={len(signals_b[chrom])}"
                )
        missing = sorted(chromosomes - set(signals_a))
        if missing:
            raise ValueError(f"No signal track for the chromosomes of peaks: {missing}")
        return {chrom: (signals_a[chrom], signals_b[chrom]) for chrom in chromosomes}

    if len(signals_a) != len(signals_b):
        raise ValueError(f"Signal lengths must match: A={len(signals_a)}, B={len(signals_b)}")
    if len(chromosomes) > 1:
        raise ValueError(
            f"One signal track cannot describe peaks on several chromosomes ({sorted(chromosomes)}): "
            "pass signals_a and signals_b as dicts of chromosome -> signal list"
        )
    return {chrom: (signals_a, signals_b) for chrom in chromosomes}


def _differential_record(
    chrom: str,
    start: int,
    end: int,
    status: str,
    track_a: Sequence[float],
    track_b: Sequence[float],
) -> dict:
    """Compare the two conditions' signal over ``[start, end)`` of one chromosome.

    The part of the region beyond the end of the tracks is ignored.

    Raises:
        ValueError: If the region lies wholly outside the tracks, or a mean signal is negative (no fold change).
    """
    low = max(0, start)
    high = min(len(track_a), end)
    if high <= low:
        raise ValueError(f"Peak {chrom}:{start}-{end} lies outside the signal tracks (length {len(track_a)})")

    values_a = track_a[low:high]
    values_b = track_b[low:high]
    mean_a = math.fsum(values_a) / len(values_a)
    mean_b = math.fsum(values_b) / len(values_b)
    if mean_a < 0 or mean_b < 0:
        raise ValueError(
            f"Mean signal over {chrom}:{start}-{end} is negative (A={mean_a:.4g}, B={mean_b:.4g}): "
            "a fold change needs non-negative signal"
        )

    pseudo = 1e-10
    log2fc = math.log2((mean_a + pseudo) / (mean_b + pseudo))
    p_value = _compute_differential_pvalue(values_a, values_b)

    if abs(log2fc) >= 1.0 and p_value < 0.05:
        direction = "up_in_A" if log2fc > 0 else "up_in_B"
    else:
        direction = "unchanged"

    return {
        "chrom": chrom,
        "start": start,
        "end": end,
        "log2_fold_change": round(log2fc, 4),
        "p_value": p_value,
        "direction": direction,
        "status": status,
        "mean_signal_a": round(mean_a, 4),
        "mean_signal_b": round(mean_b, 4),
    }


def differential_peaks(
    peaks_a: list[dict],
    peaks_b: list[dict],
    signals_a: Sequence[float] | Mapping[str, Sequence[float]],
    signals_b: Sequence[float] | Mapping[str, Sequence[float]],
) -> list[dict]:
    """Identify differential peaks between two conditions.

    Every peak is paired with the peak of the other condition it overlaps most, and the signal of the two
    conditions is compared over each pair.

    **Pairing.** Two peaks overlap when they are on the same chromosome and share at least one base. Every A peak is
    paired with the B peak it shares the most bases with, and every B peak with the A peak it shares the most bases
    with (ties go to the peak that starts first, then ends first, then comes first in its list); a pair chosen from
    both sides is reported once. A peak that overlaps several peaks of the other condition is therefore compared with
    the best-matching one, not with whichever comes first in the list, and a peak is ``"shared"`` as soon as it
    overlaps anything: ``"unique_A"`` and ``"unique_B"`` are exactly the peaks that overlap no peak of the other
    condition. A peak that is the best partner of several peaks, or whose best partner prefers another peak,
    appears in several records. The records do not depend on the order of the peaks in either list.

    **Statistics.** The region of a record is the union of its pair (the peak itself for a unique peak).
    ``mean_signal_a`` and ``mean_signal_b`` are the mean signal of each condition over the region and
    ``log2_fold_change`` is ``log2((mean_a + 1e-10) / (mean_b + 1e-10))``, finite even when one mean is 0.
    ``p_value`` is the two-sided Welch's t-test of the per-base signal of A against B over the region (bases are
    autocorrelated, so it ranks regions but is not a calibrated probability); with SciPy missing, the t distribution
    is replaced by a normal one. It is ``nan`` for a one-base region, which cannot be tested; ``0.0`` when neither
    condition varies within the region and the two levels differ (an infinite t statistic); and ``1.0`` when
    neither varies and the levels are equal. ``direction`` is ``"up_in_A"`` or ``"up_in_B"`` when
    ``|log2_fold_change| >= 1`` and ``p_value < 0.05``, else ``"unchanged"``; unique peaks are tested the same way
    (their other-condition signal is read from the track, not assumed to be zero), so a unique peak is only ``up``
    when its signal really is higher. There is no fold-change threshold argument.

    Args:
        peaks_a: Peaks from condition A (dicts with integer ``start`` and ``end``, end exclusive, and ``chrom``,
            ``"chr_unknown"`` when absent).
        peaks_b: Peaks from condition B.
        signals_a: Signal track of condition A, indexed by position: one list for peaks on a single chromosome, or
            a dict of chromosome -> list when the peaks span several.
        signals_b: Signal track(s) of condition B, of the same layout and lengths.

    Returns:
        List of differential peak dictionaries with keys: ``chrom``,
        ``start``, ``end``, ``log2_fold_change``, ``p_value``,
        ``direction`` ("up_in_A", "up_in_B", or "unchanged"),
        ``status`` ("shared", "unique_A", or "unique_B"),
        ``mean_signal_a``, ``mean_signal_b``. The records of the A peaks come first, in the order of ``peaks_a``,
        then the records of B peaks whose best partner is an A peak that chose another B peak, and the unique B
        peaks, in the order of ``peaks_b``.

    Raises:
        ValueError: If the signal tracks differ in length, one plain track would have to describe several
            chromosomes, a peak lacks integer ``start`` and ``end``, has no extent or lies wholly outside the
            tracks, or a region's mean signal is negative.
    """
    # Validate the tracks before anything else, so a length mismatch is reported even without peaks
    extents_a = _peak_extents(peaks_a, "peaks_a")
    extents_b = _peak_extents(peaks_b, "peaks_b")
    tracks = _signal_tracks(
        signals_a, signals_b, {chrom for chrom, _, _ in extents_a} | {chrom for chrom, _, _ in extents_b}
    )

    logger.info(
        "Computing differential peaks: %d peaks in A, %d peaks in B",
        len(peaks_a),
        len(peaks_b),
    )

    best_b_for_a = _best_overlaps(extents_a, extents_b)
    best_a_for_b = _best_overlaps(extents_b, extents_a)

    results: list[dict] = []

    # The A peaks first, each with its best B partner
    listed: set[tuple[int, int]] = set()
    for index_a, (chrom, start_a, end_a) in enumerate(extents_a):
        track_a, track_b = tracks[chrom]
        partner_b = best_b_for_a[index_a]
        if partner_b is None:
            results.append(_differential_record(chrom, start_a, end_a, "unique_A", track_a, track_b))
            continue
        listed.add((index_a, partner_b))
        _, start_b, end_b = extents_b[partner_b]
        results.append(
            _differential_record(chrom, min(start_a, start_b), max(end_a, end_b), "shared", track_a, track_b)
        )

    # Then the B peaks: those with no overlap, and those whose best partner is an A peak that preferred another B peak
    for index_b, (chrom, start_b, end_b) in enumerate(extents_b):
        track_a, track_b = tracks[chrom]
        partner_a = best_a_for_b[index_b]
        if partner_a is None:
            results.append(_differential_record(chrom, start_b, end_b, "unique_B", track_a, track_b))
        elif (partner_a, index_b) not in listed:
            _, start_a, end_a = extents_a[partner_a]
            results.append(
                _differential_record(chrom, min(start_a, start_b), max(end_a, end_b), "shared", track_a, track_b)
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
    vals_a: Sequence[float],
    vals_b: Sequence[float],
) -> float:
    """Compute p-value for differential signal between two conditions.

    Welch's t-test (unequal variances) of the two samples; the t distribution comes from SciPy when it is
    available, otherwise it is approximated by a normal distribution. The degenerate cases get their own answers
    instead of a silent 1.0: a sample of fewer than two values cannot be tested (``nan``), and when neither sample
    varies the t statistic is infinite if their values differ (p = 0.0) and undefined but null if they are equal
    (p = 1.0).

    Args:
        vals_a: Signal values from condition A.
        vals_b: Signal values from condition B.

    Returns:
        Two-sided p-value, or ``nan`` if fewer than two values are available in either sample.
    """
    n_a, n_b = len(vals_a), len(vals_b)
    if n_a < 2 or n_b < 2:
        return float("nan")

    def mean_and_variance(values: Sequence[float]) -> tuple[float, float]:
        # A sample whose values are all equal has exactly zero variance (summing floats would leave rounding noise)
        if min(values) == max(values):
            return float(values[0]), 0.0
        mean = math.fsum(values) / len(values)
        return mean, math.fsum((x - mean) ** 2 for x in values) / (len(values) - 1)

    mean_a, var_a = mean_and_variance(vals_a)
    mean_b, var_b = mean_and_variance(vals_b)

    term_a = var_a / n_a
    term_b = var_b / n_b
    squared_error = term_a + term_b
    if squared_error == 0.0:
        return 1.0 if mean_a == mean_b else 0.0

    t_statistic = abs(mean_a - mean_b) / math.sqrt(squared_error)

    if HAS_SCIPY:
        # Welch-Satterthwaite degrees of freedom, written with the two terms as shares of the squared error so that
        # tiny variances cannot underflow the squares; a sample with zero variance contributes no term
        share_a = term_a / squared_error
        share_b = term_b / squared_error
        degrees_of_freedom = 1.0 / (share_a**2 / (n_a - 1) + share_b**2 / (n_b - 1))
        return float(2.0 * scipy_stats.t.sf(t_statistic, degrees_of_freedom))

    # Pure-Python fallback: normal approximation of the two-sided tail
    return math.erfc(t_statistic / math.sqrt(2.0))
