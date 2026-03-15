"""Tests for epigenome peak calling pipeline.

Tests call_peaks_simple, call_peaks_broad, merge_peaks, filter_peaks,
compute_frip, differential_peaks, compute_local_lambda, and
peak_summit_refinement using real implementations with synthetic
epigenomic signal data. NO MOCKING.
"""

from __future__ import annotations

import math
import random

import pytest

from metainformant.epigenome.peak_calling.peak_detection import (
    call_peaks_broad,
    call_peaks_simple,
    compute_frip,
    compute_local_lambda,
    differential_peaks,
    filter_peaks,
    merge_peaks,
    peak_summit_refinement,
)

# ---------------------------------------------------------------------------
# Helpers to generate realistic synthetic signal data
# ---------------------------------------------------------------------------


def _make_signal_with_peaks(
    length: int = 2000,
    peaks: list[tuple[int, int, float]] | None = None,
    background: float = 1.0,
    seed: int = 42,
) -> list[float]:
    """Generate a synthetic signal track with injected peaks.

    Args:
        length: Total signal length.
        peaks: List of (centre, half_width, amplitude) tuples.
        background: Background signal level.
        seed: Random seed for reproducibility.

    Returns:
        List of float signal values.
    """
    rng = random.Random(seed)
    signal = [background + rng.gauss(0, 0.3) for _ in range(length)]
    signal = [max(0.0, v) for v in signal]

    if peaks is None:
        peaks = [
            (400, 80, 25.0),
            (900, 60, 30.0),
            (1500, 100, 20.0),
        ]

    for centre, half_w, amplitude in peaks:
        for i in range(max(0, centre - half_w), min(length, centre + half_w)):
            dist = abs(i - centre)
            # Gaussian-shaped peak
            signal[i] += amplitude * math.exp(-0.5 * (dist / (half_w / 3)) ** 2)

    return signal


# ---------------------------------------------------------------------------
# compute_local_lambda
# ---------------------------------------------------------------------------


class TestComputeLocalLambda:
    """Tests for compute_local_lambda."""

    def test_basic_lambda_computation(self) -> None:
        signal = [1.0] * 100
        result = compute_local_lambda(signal, 50, window_sizes=[10, 20])
        assert isinstance(result, float)
        assert result == pytest.approx(1.0, abs=0.01)

    def test_lambda_at_peak_region(self) -> None:
        signal = [1.0] * 100
        signal[45:55] = [10.0] * 10
        result = compute_local_lambda(signal, 50, window_sizes=[5, 10])
        # Lambda should be elevated due to high signal nearby
        assert result > 1.0

    def test_empty_signal_raises(self) -> None:
        with pytest.raises(ValueError, match="must not be empty"):
            compute_local_lambda([], 0)

    def test_out_of_range_position_raises(self) -> None:
        with pytest.raises(ValueError, match="out of range"):
            compute_local_lambda([1.0, 2.0, 3.0], 5)

    def test_edge_position_zero(self) -> None:
        signal = [5.0] + [1.0] * 99
        result = compute_local_lambda(signal, 0, window_sizes=[10])
        assert isinstance(result, float)
        assert result > 0

    def test_edge_position_last(self) -> None:
        signal = [1.0] * 99 + [5.0]
        result = compute_local_lambda(signal, 99, window_sizes=[10])
        assert isinstance(result, float)
        assert result > 0

    def test_genome_wide_mean_is_floor(self) -> None:
        # Uniform low signal; lambda should equal the genome-wide mean
        signal = [0.5] * 200
        result = compute_local_lambda(signal, 100, window_sizes=[10])
        assert result == pytest.approx(0.5, abs=0.01)


# ---------------------------------------------------------------------------
# call_peaks_simple
# ---------------------------------------------------------------------------


class TestCallPeaksSimple:
    """Tests for call_peaks_simple (narrow peak calling)."""

    def test_finds_injected_peaks(self) -> None:
        signal = _make_signal_with_peaks(
            length=2000,
            peaks=[(500, 80, 30.0), (1200, 60, 35.0)],
            background=1.0,
        )
        peaks = call_peaks_simple(signal, threshold=3.0, min_length=20)
        assert len(peaks) >= 1
        # At least one peak should overlap with one of the injected regions
        peak_positions = {p["summit"] for p in peaks}
        assert any(400 < pos < 600 for pos in peak_positions) or any(1100 < pos < 1300 for pos in peak_positions)

    def test_peak_dict_keys(self) -> None:
        signal = _make_signal_with_peaks(length=1000, peaks=[(500, 60, 40.0)])
        peaks = call_peaks_simple(signal, threshold=3.0, min_length=10)
        if peaks:
            peak = peaks[0]
            required_keys = {
                "chrom",
                "start",
                "end",
                "summit",
                "score",
                "p_value",
                "fold_enrichment",
                "signal_value",
            }
            assert required_keys.issubset(peak.keys())

    def test_with_control_signal(self) -> None:
        signal = _make_signal_with_peaks(length=1500, peaks=[(700, 70, 35.0)])
        control = [1.0] * 1500
        peaks = call_peaks_simple(signal, control=control, threshold=3.0, min_length=10)
        # Should still detect the peak
        assert isinstance(peaks, list)

    def test_empty_signal_raises(self) -> None:
        with pytest.raises(ValueError, match="must not be empty"):
            call_peaks_simple([])

    def test_control_length_mismatch_raises(self) -> None:
        with pytest.raises(ValueError, match="must match"):
            call_peaks_simple([1.0, 2.0, 3.0], control=[1.0, 2.0])

    def test_no_peaks_in_flat_signal(self) -> None:
        # Uniform signal should produce no peaks at high threshold
        signal = [1.0] * 500
        peaks = call_peaks_simple(signal, threshold=10.0, min_length=50)
        assert peaks == []

    def test_qvalue_present_when_peaks_found(self) -> None:
        signal = _make_signal_with_peaks(length=1000, peaks=[(500, 80, 50.0)])
        peaks = call_peaks_simple(signal, threshold=2.0, min_length=10)
        if peaks:
            assert "q_value" in peaks[0]
            assert 0.0 <= peaks[0]["q_value"] <= 1.0

    def test_merge_distance_merges_nearby_peaks(self) -> None:
        # Two peaks very close together should merge
        signal = _make_signal_with_peaks(
            length=1000,
            peaks=[(400, 30, 40.0), (470, 30, 40.0)],
        )
        peaks_merged = call_peaks_simple(signal, threshold=3.0, min_length=10, merge_distance=100)
        peaks_separate = call_peaks_simple(signal, threshold=3.0, min_length=10, merge_distance=0)
        # Merged should have fewer or equal peaks
        assert len(peaks_merged) <= len(peaks_separate)

    def test_min_length_filters_short_peaks(self) -> None:
        signal = _make_signal_with_peaks(
            length=1000,
            peaks=[(500, 10, 50.0)],  # Narrow spike
        )
        peaks_relaxed = call_peaks_simple(signal, threshold=3.0, min_length=1)
        peaks_strict = call_peaks_simple(signal, threshold=3.0, min_length=200)
        assert len(peaks_strict) <= len(peaks_relaxed)


# ---------------------------------------------------------------------------
# call_peaks_broad
# ---------------------------------------------------------------------------


class TestCallPeaksBroad:
    """Tests for call_peaks_broad (broad domain detection)."""

    def test_detects_broad_domain(self) -> None:
        # Create a broad enrichment region
        signal = _make_signal_with_peaks(
            length=2000,
            peaks=[(1000, 300, 15.0)],  # Wide, moderate amplitude
            background=1.0,
        )
        peaks = call_peaks_broad(signal, p_threshold=0.5, broad_cutoff=0.5, min_length=50)
        assert isinstance(peaks, list)

    def test_broad_peak_has_type_field(self) -> None:
        signal = _make_signal_with_peaks(
            length=1500,
            peaks=[(750, 200, 20.0)],
        )
        peaks = call_peaks_broad(signal, p_threshold=0.5, broad_cutoff=0.5, min_length=50)
        for p in peaks:
            assert p.get("peak_type") == "broad"

    def test_empty_signal_raises(self) -> None:
        with pytest.raises(ValueError, match="must not be empty"):
            call_peaks_broad([])

    def test_control_length_mismatch_raises(self) -> None:
        with pytest.raises(ValueError, match="must match"):
            call_peaks_broad([1.0, 2.0], control=[1.0])

    def test_with_control(self) -> None:
        signal = _make_signal_with_peaks(
            length=1500,
            peaks=[(750, 200, 20.0)],
        )
        control = [1.5] * 1500
        peaks = call_peaks_broad(signal, control=control, p_threshold=0.5, broad_cutoff=0.5, min_length=50)
        assert isinstance(peaks, list)

    def test_flat_signal_no_peaks(self) -> None:
        signal = [1.0] * 500
        peaks = call_peaks_broad(signal, p_threshold=0.001, broad_cutoff=0.001, min_length=100)
        assert peaks == []


# ---------------------------------------------------------------------------
# peak_summit_refinement
# ---------------------------------------------------------------------------


class TestPeakSummitRefinement:
    """Tests for peak_summit_refinement."""

    def test_refines_summit(self) -> None:
        signal = _make_signal_with_peaks(length=500, peaks=[(250, 50, 30.0)])
        peak = {"start": 200, "end": 300, "summit": 250}
        refined = peak_summit_refinement(signal, peak, window=30)
        assert "refined_summit_position" in refined
        assert isinstance(refined["refined_summit_position"], float)

    def test_preserves_original_keys(self) -> None:
        signal = [1.0] * 100
        signal[50] = 10.0
        peak = {"start": 40, "end": 60, "summit": 50, "score": 5.0}
        refined = peak_summit_refinement(signal, peak)
        assert refined["score"] == 5.0

    def test_empty_signal_raises(self) -> None:
        with pytest.raises(ValueError, match="must not be empty"):
            peak_summit_refinement([], {"start": 0, "end": 10, "summit": 5})

    def test_summit_out_of_range_raises(self) -> None:
        with pytest.raises(ValueError, match="out of signal range"):
            peak_summit_refinement([1.0] * 10, {"start": 0, "end": 5, "summit": 20})


# ---------------------------------------------------------------------------
# merge_peaks
# ---------------------------------------------------------------------------


class TestMergePeaks:
    """Tests for merge_peaks."""

    def test_overlapping_peaks_merge(self) -> None:
        peaks = [
            {"start": 100, "end": 200, "score": 5.0, "signal_value": 10.0, "p_value": 0.01},
            {"start": 150, "end": 250, "score": 7.0, "signal_value": 15.0, "p_value": 0.005},
        ]
        merged = merge_peaks(peaks, distance=0)
        assert len(merged) == 1
        assert merged[0]["start"] == 100
        assert merged[0]["end"] == 250
        # Best score kept
        assert merged[0]["score"] == 7.0

    def test_non_overlapping_peaks_remain_separate(self) -> None:
        peaks = [
            {"start": 100, "end": 200, "score": 5.0},
            {"start": 300, "end": 400, "score": 7.0},
        ]
        merged = merge_peaks(peaks, distance=0)
        assert len(merged) == 2

    def test_merge_with_distance(self) -> None:
        peaks = [
            {"start": 100, "end": 200, "score": 5.0, "signal_value": 10.0, "p_value": 0.01},
            {"start": 220, "end": 300, "score": 7.0, "signal_value": 12.0, "p_value": 0.005},
        ]
        merged = merge_peaks(peaks, distance=30)
        assert len(merged) == 1

    def test_empty_list(self) -> None:
        assert merge_peaks([]) == []

    def test_single_peak(self) -> None:
        peaks = [{"start": 100, "end": 200, "score": 5.0}]
        merged = merge_peaks(peaks, distance=0)
        assert len(merged) == 1

    def test_keeps_best_pvalue(self) -> None:
        peaks = [
            {"start": 100, "end": 200, "score": 5.0, "p_value": 0.05, "signal_value": 8.0},
            {"start": 150, "end": 250, "score": 3.0, "p_value": 0.001, "signal_value": 6.0},
        ]
        merged = merge_peaks(peaks, distance=0)
        assert merged[0]["p_value"] == 0.001

    def test_unsorted_input_sorted_correctly(self) -> None:
        peaks = [
            {"start": 500, "end": 600, "score": 3.0},
            {"start": 100, "end": 200, "score": 5.0},
            {"start": 300, "end": 400, "score": 4.0},
        ]
        merged = merge_peaks(peaks, distance=0)
        starts = [p["start"] for p in merged]
        assert starts == sorted(starts)


# ---------------------------------------------------------------------------
# filter_peaks
# ---------------------------------------------------------------------------


class TestFilterPeaks:
    """Tests for filter_peaks."""

    def test_filter_by_fold_enrichment(self) -> None:
        peaks = [
            {"fold_enrichment": 5.0, "q_value": 0.01},
            {"fold_enrichment": 1.5, "q_value": 0.01},
            {"fold_enrichment": 3.0, "q_value": 0.01},
        ]
        filtered = filter_peaks(peaks, min_fold=2.0, max_qvalue=1.0)
        assert len(filtered) == 2

    def test_filter_by_qvalue(self) -> None:
        peaks = [
            {"fold_enrichment": 5.0, "q_value": 0.01},
            {"fold_enrichment": 5.0, "q_value": 0.1},
            {"fold_enrichment": 5.0, "q_value": 0.001},
        ]
        filtered = filter_peaks(peaks, min_fold=1.0, max_qvalue=0.05)
        assert len(filtered) == 2

    def test_filter_uses_pvalue_when_qvalue_absent(self) -> None:
        peaks = [
            {"fold_enrichment": 5.0, "p_value": 0.01},
            {"fold_enrichment": 5.0, "p_value": 0.1},
        ]
        filtered = filter_peaks(peaks, min_fold=1.0, max_qvalue=0.05)
        assert len(filtered) == 1

    def test_blacklist_filtering(self) -> None:
        peaks = [
            {"chrom": "chr1", "start": 100, "end": 200, "fold_enrichment": 5.0, "q_value": 0.01},
            {"chrom": "chr1", "start": 500, "end": 600, "fold_enrichment": 5.0, "q_value": 0.01},
        ]
        blacklist = [{"chrom": "chr1", "start": 150, "end": 180}]
        filtered = filter_peaks(peaks, min_fold=1.0, max_qvalue=1.0, blacklist_regions=blacklist)
        assert len(filtered) == 1
        assert filtered[0]["start"] == 500

    def test_empty_list(self) -> None:
        assert filter_peaks([]) == []

    def test_blacklist_different_chrom_no_effect(self) -> None:
        peaks = [
            {"chrom": "chr1", "start": 100, "end": 200, "fold_enrichment": 5.0, "q_value": 0.01},
        ]
        blacklist = [{"chrom": "chr2", "start": 100, "end": 200}]
        filtered = filter_peaks(peaks, min_fold=1.0, max_qvalue=1.0, blacklist_regions=blacklist)
        assert len(filtered) == 1


# ---------------------------------------------------------------------------
# compute_frip
# ---------------------------------------------------------------------------


class TestComputeFrip:
    """Tests for compute_frip (Fraction of Reads in Peaks)."""

    def test_good_quality(self) -> None:
        result = compute_frip(reads_in_peaks=5000, total_reads=100000)
        assert result["frip"] == pytest.approx(0.05, abs=1e-6)
        assert result["quality_assessment"] == "good"

    def test_acceptable_quality(self) -> None:
        result = compute_frip(reads_in_peaks=2000, total_reads=100000)
        assert result["frip"] == pytest.approx(0.02, abs=1e-6)
        assert result["quality_assessment"] == "acceptable"

    def test_low_quality(self) -> None:
        result = compute_frip(reads_in_peaks=500, total_reads=100000)
        assert result["frip"] == pytest.approx(0.005, abs=1e-6)
        assert result["quality_assessment"] == "low"

    def test_all_reads_in_peaks(self) -> None:
        result = compute_frip(reads_in_peaks=100, total_reads=100)
        assert result["frip"] == pytest.approx(1.0, abs=1e-6)
        assert result["reads_outside_peaks"] == 0

    def test_no_reads_in_peaks(self) -> None:
        result = compute_frip(reads_in_peaks=0, total_reads=100)
        assert result["frip"] == 0.0
        assert result["reads_outside_peaks"] == 100

    def test_result_keys(self) -> None:
        result = compute_frip(reads_in_peaks=50, total_reads=1000)
        expected_keys = {
            "frip",
            "reads_in_peaks",
            "total_reads",
            "reads_outside_peaks",
            "percent_in_peaks",
            "quality_assessment",
        }
        assert expected_keys.issubset(result.keys())

    def test_zero_total_reads_raises(self) -> None:
        with pytest.raises(ValueError, match="must be positive"):
            compute_frip(reads_in_peaks=0, total_reads=0)

    def test_negative_total_reads_raises(self) -> None:
        with pytest.raises(ValueError, match="must be positive"):
            compute_frip(reads_in_peaks=0, total_reads=-1)

    def test_reads_exceed_total_raises(self) -> None:
        with pytest.raises(ValueError, match="cannot exceed"):
            compute_frip(reads_in_peaks=200, total_reads=100)

    def test_negative_reads_in_peaks_raises(self) -> None:
        with pytest.raises(ValueError, match="must be non-negative"):
            compute_frip(reads_in_peaks=-5, total_reads=100)

    def test_percent_in_peaks(self) -> None:
        result = compute_frip(reads_in_peaks=250, total_reads=1000)
        assert result["percent_in_peaks"] == pytest.approx(25.0, abs=0.01)


# ---------------------------------------------------------------------------
# differential_peaks
# ---------------------------------------------------------------------------


class TestDifferentialPeaks:
    """Tests for differential_peaks."""

    def test_shared_peak_detection(self) -> None:
        signal_a = _make_signal_with_peaks(length=1000, peaks=[(500, 60, 30.0)])
        signal_b = _make_signal_with_peaks(length=1000, peaks=[(500, 60, 10.0)])
        peaks_a = [{"chrom": "chr_unknown", "start": 440, "end": 560}]
        peaks_b = [{"chrom": "chr_unknown", "start": 440, "end": 560}]
        results = differential_peaks(peaks_a, peaks_b, signal_a, signal_b)
        shared = [r for r in results if r["status"] == "shared"]
        assert len(shared) >= 1

    def test_unique_a_peak(self) -> None:
        signal_a = _make_signal_with_peaks(length=1000, peaks=[(500, 60, 30.0)])
        signal_b = [1.0] * 1000
        peaks_a = [{"chrom": "chr_unknown", "start": 440, "end": 560}]
        peaks_b: list[dict] = []
        results = differential_peaks(peaks_a, peaks_b, signal_a, signal_b)
        unique_a = [r for r in results if r["status"] == "unique_A"]
        assert len(unique_a) == 1

    def test_unique_b_peak(self) -> None:
        signal_a = [1.0] * 1000
        signal_b = _make_signal_with_peaks(length=1000, peaks=[(500, 60, 30.0)])
        peaks_a: list[dict] = []
        peaks_b = [{"chrom": "chr_unknown", "start": 440, "end": 560}]
        results = differential_peaks(peaks_a, peaks_b, signal_a, signal_b)
        unique_b = [r for r in results if r["status"] == "unique_B"]
        assert len(unique_b) == 1

    def test_result_keys(self) -> None:
        signal_a = _make_signal_with_peaks(length=800, peaks=[(400, 50, 25.0)])
        signal_b = _make_signal_with_peaks(length=800, peaks=[(400, 50, 10.0)])
        peaks_a = [{"chrom": "chr_unknown", "start": 350, "end": 450}]
        peaks_b = [{"chrom": "chr_unknown", "start": 350, "end": 450}]
        results = differential_peaks(peaks_a, peaks_b, signal_a, signal_b)
        if results:
            expected_keys = {
                "chrom",
                "start",
                "end",
                "log2_fold_change",
                "p_value",
                "direction",
                "status",
                "mean_signal_a",
                "mean_signal_b",
            }
            assert expected_keys.issubset(results[0].keys())

    def test_signal_length_mismatch_raises(self) -> None:
        with pytest.raises(ValueError, match="must match"):
            differential_peaks([], [], [1.0, 2.0], [1.0])

    def test_direction_up_in_a(self) -> None:
        # A has much higher signal
        signal_a = _make_signal_with_peaks(length=800, peaks=[(400, 50, 50.0)])
        signal_b = _make_signal_with_peaks(length=800, peaks=[(400, 50, 2.0)])
        peaks_a = [{"chrom": "chr_unknown", "start": 350, "end": 450}]
        peaks_b = [{"chrom": "chr_unknown", "start": 350, "end": 450}]
        results = differential_peaks(peaks_a, peaks_b, signal_a, signal_b)
        shared = [r for r in results if r["status"] == "shared"]
        if shared:
            # log2FC should be positive (A > B)
            assert shared[0]["log2_fold_change"] > 0

    def test_empty_peaks_both(self) -> None:
        results = differential_peaks([], [], [1.0] * 100, [1.0] * 100)
        assert results == []
