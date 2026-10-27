"""Regression tests for CNV detection: depth validation and the vectorised CBS.

* Negative depths (log2 ratios passed by mistake) became NaN under ``log2`` and then crashed in
  ``_assign_cnv_state`` with ``ValueError: cannot convert float NaN to integer``. Depths are now validated
  and the error says what is wrong. Zero depth stays valid: it is a bin with no coverage.
* CBS scored every arc of the data one at a time in pure Python (about 40-90 s for 8,000 bins). The scan is now
  vectorised over the arcs of one length. It performs the same floating-point operations in the same order,
  so the segmentation is identical; ``_reference_*`` below is the original implementation, kept here verbatim
  as the oracle that the new one is compared with.
"""

from __future__ import annotations

import math
import time
from typing import Any

import numpy as np
import pytest

from metainformant.structural_variants.detection import cnv
from metainformant.structural_variants.detection.cnv import detect_cnv_from_depth, segment_coverage

# ---------------------------------------------------------------------------
# The original (pure-Python) CBS, as the oracle
# ---------------------------------------------------------------------------


def _reference_cbs_recursive(
    data: Any,
    start: int,
    end: int,
    segments: list[tuple[int, int, float]],
    significance: float,
    max_iterations: int,
    min_width: int,
    depth: int,
) -> None:
    n = end - start

    if n < 2 * min_width or depth >= max_iterations:
        segments.append((start, end, float(np.mean(data[start:end]))))
        return

    segment_data = data[start:end]
    best_t_stat = 0.0
    best_i = -1
    best_j = -1

    cumsum = np.cumsum(segment_data - np.mean(segment_data))
    cumsum = np.insert(cumsum, 0, 0.0)  # prepend 0

    variance = np.var(segment_data, ddof=1) if n > 1 else 1.0
    if variance < 1e-12:
        segments.append((start, end, float(np.mean(segment_data))))
        return

    std_dev = math.sqrt(variance)

    # Standard binary segmentation: find single best split point
    best_single_stat = 0.0
    best_split = -1

    for k in range(min_width, n - min_width + 1):
        left_mean = np.mean(segment_data[:k])
        right_mean = np.mean(segment_data[k:])
        n_left = k
        n_right = n - k

        t_stat = abs(left_mean - right_mean) * math.sqrt(n_left * n_right / n) / std_dev

        if t_stat > best_single_stat:
            best_single_stat = t_stat
            best_split = k

    # Circular scan: test arc segments (i, j) where the arc mean differs from rest
    best_arc_stat = 0.0
    best_arc_i = -1
    best_arc_j = -1

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

    use_arc = best_arc_stat > best_single_stat

    if use_arc:
        best_t_stat = best_arc_stat
        best_i = best_arc_i
        best_j = best_arc_j
    else:
        best_t_stat = best_single_stat
        best_i = 0
        best_j = best_split

    p_value = cnv._cbs_pvalue(best_t_stat, n)

    if p_value > significance:
        segments.append((start, end, float(np.mean(segment_data))))
        return

    if use_arc:
        abs_i = start + best_i
        abs_j = start + best_j

        if best_i > 0:
            _reference_cbs_recursive(data, start, abs_i, segments, significance, max_iterations, min_width, depth + 1)
        _reference_cbs_recursive(data, abs_i, abs_j, segments, significance, max_iterations, min_width, depth + 1)
        if abs_j < end:
            _reference_cbs_recursive(data, abs_j, end, segments, significance, max_iterations, min_width, depth + 1)
    else:
        abs_split = start + best_j
        _reference_cbs_recursive(data, start, abs_split, segments, significance, max_iterations, min_width, depth + 1)
        _reference_cbs_recursive(data, abs_split, end, segments, significance, max_iterations, min_width, depth + 1)


def _reference_segment_coverage(
    coverage: Any, significance: float = 0.01, max_iterations: int = 100, min_width: int = 2
) -> list[tuple[int, int, float]]:
    data = np.asarray(coverage, dtype=np.float64)
    n = len(data)
    if n < 2 * min_width:
        return [(0, n, float(np.mean(data)))]
    segments: list[tuple[int, int, float]] = []
    _reference_cbs_recursive(data, 0, n, segments, significance, max_iterations, min_width, depth=0)
    segments.sort(key=lambda x: x[0])
    return segments


def _profile(rng: np.random.Generator, n: int, noise: float, changes: int) -> np.ndarray:
    values = rng.normal(0.0, noise, n)
    for _ in range(changes):
        begin = int(rng.integers(0, max(1, n - 2)))
        finish = int(min(n, begin + rng.integers(2, max(3, n // 2))))
        values[begin:finish] += rng.choice([-1.0, -0.6, 0.5, 1.0])
    return values


# Integer-valued and repeating data: these produce exact ties between candidate change points, where only
# the order in which the arcs are scanned decides, so they are where a different scan order would show.
TIE_PRONE = {
    "step": [0] * 50 + [2] * 50,
    "plateau": [0] * 50 + [2] * 50 + [0] * 50,
    "two-plateaus": [0] * 50 + [2] * 50 + [0] * 50 + [2] * 50,
    "square-wave": [0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1],
    "spikes": [1, 1, 1, 5, 5, 5, 1, 1, 1, 5, 5, 5, 1, 1, 1],
    "decimal-plateau": [0.1] * 30 + [0.3] * 30 + [0.1] * 30,
    "five-blocks": [0.1] * 20 + [0.7] * 20 + [0.1] * 20 + [0.7] * 20 + [0.1] * 20,
    "ramp": list(range(40)),
    "short-bump": [0] * 10 + [3] * 3 + [0] * 10,
    "half-and-half": [0] * 6 + [1] * 6,
    "alternating": [0, 1] * 20,
    "last-bin": [0] * 99 + [5],
    "first-bin": [5] + [0] * 99,
    "constant": [0] * 30,
}


class TestVectorisedSegmentationMatchesTheOriginal:
    @pytest.mark.parametrize("name", sorted(TIE_PRONE))
    @pytest.mark.parametrize("min_width", [1, 2, 5])
    def test_tie_prone_data(self, name: str, min_width: int) -> None:
        data = np.array(TIE_PRONE[name], dtype=float)
        assert segment_coverage(data, min_width=min_width) == _reference_segment_coverage(data, min_width=min_width)

    @pytest.mark.parametrize("seed", range(12))
    def test_random_profiles_with_change_points(self, seed: int) -> None:
        rng = np.random.default_rng(20261027 + seed)
        n = int(rng.choice([6, 9, 14, 25, 40, 70, 110, 160]))
        data = _profile(rng, n, noise=float(rng.choice([0.05, 0.2, 0.5])), changes=int(rng.integers(0, 4)))
        for significance in (0.05, 0.01, 0.001):
            for min_width in (1, 2, 3):
                new = segment_coverage(data, significance=significance, min_width=min_width)
                assert new == _reference_segment_coverage(data, significance=significance, min_width=min_width)

    def test_the_fixture_profile_of_the_cnv_tests(self) -> None:
        rng = np.random.default_rng(42)
        depths = rng.normal(100, 10, 500)
        depths[100:150] = np.maximum(0, rng.normal(40, 5, 50))
        depths[300:350] = rng.normal(160, 10, 50)
        log2_ratios = np.clip(np.log2(depths / np.median(depths) + 1e-10), -5.0, 5.0)

        new = segment_coverage(log2_ratios)
        assert new == _reference_segment_coverage(log2_ratios)
        assert [(start, end) for start, end, _ in new] != [(0, 500)]  # the comparison is not of a single segment

    def test_the_segments_cover_the_data_exactly_once(self) -> None:
        rng = np.random.default_rng(5)
        data = _profile(rng, 300, noise=0.1, changes=3)
        segments = segment_coverage(data)
        assert segments[0][0] == 0 and segments[-1][1] == 300
        assert all(a[1] == b[0] for a, b in zip(segments, segments[1:]))
        for start, end, mean in segments:
            assert mean == float(np.mean(data[start:end]))

    def test_non_finite_input_is_one_segment_as_before(self) -> None:
        data = np.array([0.0, 1.0, np.nan, 2.0, 3.0, 1.0, 0.0, 2.0])
        with np.errstate(all="ignore"):
            segments = segment_coverage(data)
        assert [(s, e) for s, e, _ in segments] == [(0, 8)]

    @pytest.mark.parametrize("min_width", [0, -1])
    def test_a_minimum_width_below_one_is_rejected(self, min_width: int) -> None:
        with pytest.raises(ValueError, match="min_width"):
            segment_coverage([0.0, 1.0, 2.0, 3.0], min_width=min_width)


def test_the_vectorised_scan_is_an_order_of_magnitude_faster_than_the_scalar_one() -> None:
    """On the same data, in the same process, the scalar scan takes 10 to 100 times as long as the vectorised one.

    The cost still grows with the square of the segment length (every arc is scored), so what the vectorised scan
    changes is the constant: 8,000 bins took 37-86 s, and now take a fraction of a second. Timing the original
    next to the new code keeps the comparison independent of the machine's speed.
    """
    rng = np.random.default_rng(11)
    data = rng.normal(0.0, 0.2, 800)
    data[100:160] -= 0.9
    data[400:520] += 0.6
    data[640:650] -= 1.5

    started = time.perf_counter()
    expected = _reference_segment_coverage(data)
    scalar_time = time.perf_counter() - started
    started = time.perf_counter()
    segments = segment_coverage(data)
    vectorised_time = time.perf_counter() - started

    assert segments == expected
    boundaries = {start for start, _, _ in segments} | {end for _, end, _ in segments}
    for planted in (100, 160, 400, 520, 640, 650):
        assert min(abs(boundary - planted) for boundary in boundaries) <= 3, planted
    assert scalar_time > 10 * vectorised_time, f"scalar {scalar_time:.3f}s, vectorised {vectorised_time:.3f}s"


# ---------------------------------------------------------------------------
# End to end
# ---------------------------------------------------------------------------


def _depth_profile(n: int = 500, seed: int = 3) -> np.ndarray:
    """Depth around 100x; profiles of 350 bins or more carry a lost and a gained stretch."""
    rng = np.random.default_rng(seed)
    depth = rng.normal(100.0, 2.0, n)
    if n >= 350:
        depth[100:150] = rng.normal(50.0, 2.0, 50)  # one copy lost
        depth[300:350] = rng.normal(150.0, 2.0, 50)  # one copy gained
    return depth


def test_detect_cnv_recovers_the_planted_copy_number_changes() -> None:
    result = detect_cnv_from_depth({"chr1": _depth_profile()}, window_size=1_000)["chr1"]
    called = [(s.start, s.end, s.state, s.cn) for s in result.segments if s.state != "NEUTRAL"]
    assert called == [(100_000, 150_000, "DEL", 1), (300_000, 350_000, "DUP", 3)]


class TestDepthValidation:
    def test_negative_depths_get_a_clear_error(self) -> None:
        log2_ratios = np.random.default_rng(1).normal(0.0, 0.3, 100)  # what someone might pass by mistake
        with pytest.raises(ValueError, match=r"depth_data\['chr1'\] has \d+ negative") as error:
            detect_cnv_from_depth({"chr1": log2_ratios.tolist()})
        message = str(error.value)
        assert "log2 ratios" in message
        assert "segment_coverage" in message
        assert "cannot convert" not in message

    def test_the_error_counts_the_offending_values(self) -> None:
        depth = _depth_profile(60)
        depth[[3, 9]] = -1.0
        depth[17] = np.nan
        with pytest.raises(ValueError, match=r"2 negative and 1 non-finite depth value\(s\) out of 60"):
            detect_cnv_from_depth({"chr1": depth})

    @pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf])
    def test_non_finite_depths_are_rejected(self, bad: float) -> None:
        depth = _depth_profile(60)
        depth[5] = bad
        with pytest.raises(ValueError, match="non-finite"):
            detect_cnv_from_depth({"chr1": depth})

    def test_zero_depth_is_valid(self) -> None:
        depth = np.concatenate([np.zeros(20), np.full(80, 100.0)])
        result = detect_cnv_from_depth({"chr1": depth}, window_size=1_000)["chr1"]
        assert result.segments[0].start == 0
        assert result.segments[0].state == "HOMODEL"  # the uncovered stretch
        assert min(result.log2_ratios or []) == -5.0

    def test_every_analysed_chromosome_is_checked_before_any_is_segmented(self) -> None:
        bad = _depth_profile(60)
        bad[0] = -3.0
        with pytest.raises(ValueError, match=r"depth_data\['chr2'\]"):
            detect_cnv_from_depth({"chr1": _depth_profile(), "chr2": bad})

    def test_a_chromosome_that_is_skipped_for_having_too_few_bins_is_not_checked(self) -> None:
        results = detect_cnv_from_depth({"chr1": _depth_profile(), "chrM": [-1.0, -2.0]}, min_segment_bins=3)
        assert set(results) == {"chr1"}

    def test_integer_depths_are_accepted(self) -> None:
        depth = [int(v) for v in _depth_profile(200)]
        assert "chr1" in detect_cnv_from_depth({"chr1": depth}, window_size=1_000)
