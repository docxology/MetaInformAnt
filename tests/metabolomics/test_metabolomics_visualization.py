"""Tests for metabolomics visualization utilities."""
from __future__ import annotations

import numpy as np
import pytest

from metainformant.metabolomics.visualization.plots import (
    ChromatogramPeak,
    MassSpectrumPlotData,
    VolcanoPoint,
    detect_chromatographic_peaks,
    intensity_heatmap_data,
    pca_metabolomics,
    prepare_spectrum_plot,
    retention_time_alignment_data,
    volcano_plot_data,
)


class TestVolcanoPlotData:
    """Tests for volcano plot data preparation."""

    def test_significant_classification(self) -> None:
        names = ["A", "B", "C"]
        fc = np.array([2.0, 0.1, -2.0])
        pvals = np.array([0.001, 0.5, 0.01])
        points = volcano_plot_data(names, fc, pvals, fc_threshold=1.0, p_threshold=0.05)
        assert len(points) == 3
        assert bool(points[0].significant) is True  # fc=2, p=0.001
        assert bool(points[1].significant) is False  # fc=0.1
        assert bool(points[2].significant) is True  # fc=-2, p=0.01

    def test_point_fields(self) -> None:
        points = volcano_plot_data(["X"], np.array([1.5]), np.array([0.01]))
        p = points[0]
        assert isinstance(p, VolcanoPoint)
        assert p.name == "X"
        assert p.log2_fc == pytest.approx(1.5)
        assert p.neg_log10_p == pytest.approx(2.0)

    def test_tiny_pvalue(self) -> None:
        """Very small p-values should not produce inf."""
        points = volcano_plot_data(["A"], np.array([1.0]), np.array([1e-300]))
        assert np.isfinite(points[0].neg_log10_p)

    def test_empty_input(self) -> None:
        points = volcano_plot_data([], np.array([]), np.array([]))
        assert points == []


class TestPCAMetabolomics:
    """Tests for PCA on metabolomics data."""

    def test_output_shapes(self) -> None:
        """PCA should return correct shapes."""
        data = np.random.default_rng(42).random((50, 10))  # 50 metabolites, 10 samples
        scores, var_ratio = pca_metabolomics(data, n_components=2)
        assert scores.shape == (10, 2)
        assert var_ratio.shape == (2,)

    def test_variance_ratio_sums(self) -> None:
        """Variance ratios should be between 0 and 1."""
        data = np.random.default_rng(42).random((20, 8))
        _, var_ratio = pca_metabolomics(data, n_components=3)
        assert all(0 <= v <= 1 for v in var_ratio)
        assert sum(var_ratio) <= 1.0 + 1e-10

    def test_single_component(self) -> None:
        data = np.random.default_rng(42).random((10, 5))
        scores, var_ratio = pca_metabolomics(data, n_components=1)
        assert scores.shape == (5, 1)
        assert var_ratio.shape == (1,)

    def test_constant_data(self) -> None:
        """All-constant data should yield zero variance."""
        data = np.ones((10, 5))
        scores, var_ratio = pca_metabolomics(data, n_components=2)
        np.testing.assert_allclose(var_ratio, [0.0, 0.0], atol=1e-10)


class TestIntensityHeatmapData:
    """Tests for heatmap data preparation."""

    def test_z_score_normalization(self) -> None:
        data = np.array([[10.0, 20.0, 30.0], [100.0, 100.0, 100.0]])
        result = intensity_heatmap_data(data, ["A", "B"], ["S1", "S2", "S3"], normalize=True)
        matrix = result["matrix"]
        # Row means should be ~0 after z-scoring
        assert abs(matrix[0].mean()) < 1e-10
        # Constant row (B) should be all zeros
        np.testing.assert_allclose(matrix[1], [0.0, 0.0, 0.0], atol=1e-10)

    def test_no_normalization(self) -> None:
        data = np.array([[10.0, 20.0]])
        result = intensity_heatmap_data(data, ["A"], ["S1", "S2"], normalize=False)
        np.testing.assert_array_equal(result["matrix"], data)

    def test_labels_preserved(self) -> None:
        data = np.array([[1.0, 2.0]])
        result = intensity_heatmap_data(data, ["met1"], ["sampleA", "sampleB"])
        assert result["row_labels"] == ["met1"]
        assert result["col_labels"] == ["sampleA", "sampleB"]


class TestDetectChromatographicPeaks:
    """Tests for peak detection in chromatograms."""

    def test_single_gaussian_peak(self) -> None:
        """Should detect a clear Gaussian peak."""
        x = np.linspace(0, 100, 200)
        y = 1000.0 * np.exp(-((x - 50) ** 2) / (2 * 5**2))
        peaks = detect_chromatographic_peaks(x, y, min_intensity=10.0)
        assert len(peaks) >= 1
        best = max(peaks, key=lambda p: p.intensity)
        assert 45 < best.retention_time < 55
        assert best.area > 0

    def test_multiple_peaks(self) -> None:
        x = np.linspace(0, 200, 400)
        y = (500.0 * np.exp(-((x - 50) ** 2) / 50)
             + 800.0 * np.exp(-((x - 150) ** 2) / 50))
        peaks = detect_chromatographic_peaks(x, y, min_intensity=50.0)
        assert len(peaks) >= 2

    def test_no_peaks_in_noise(self) -> None:
        x = np.linspace(0, 100, 50)
        y = np.ones(50) * 5.0
        peaks = detect_chromatographic_peaks(x, y, min_intensity=10.0)
        assert len(peaks) == 0

    def test_short_array(self) -> None:
        peaks = detect_chromatographic_peaks(np.array([1.0, 2.0]), np.array([10.0, 20.0]))
        assert peaks == []

    def test_peak_fields(self) -> None:
        x = np.array([0, 1, 2, 3, 4, 5, 6], dtype=float)
        y = np.array([0, 10, 50, 100, 50, 10, 0], dtype=float)
        peaks = detect_chromatographic_peaks(x, y)
        if peaks:
            p = peaks[0]
            assert isinstance(p, ChromatogramPeak)
            assert p.left_idx <= p.peak_idx <= p.right_idx
            assert p.area > 0


class TestPrepareSpectrumPlot:
    """Tests for mass spectrum plot preparation."""

    def test_base_peak_normalization(self) -> None:
        mz = np.array([100.0, 200.0, 300.0])
        intensity = np.array([500.0, 1000.0, 250.0])
        result = prepare_spectrum_plot(mz, intensity)
        assert isinstance(result, MassSpectrumPlotData)
        assert result.base_peak_mz == 200.0
        assert result.base_peak_intensity == 1000.0
        # Intensities should be on 0-100 scale
        assert max(result.intensities) == pytest.approx(100.0)

    def test_top_n_filtering(self) -> None:
        mz = np.arange(100.0, 200.0, 1.0)  # 100 peaks
        intensity = np.random.default_rng(42).random(100) * 1000
        result = prepare_spectrum_plot(mz, intensity, top_n=10)
        assert len(result.mz_values) == 10

    def test_annotation_matching(self) -> None:
        mz = np.array([100.0, 200.0])
        intensity = np.array([500.0, 1000.0])
        annotations = {100.0: "fragment_A", 200.0: "fragment_B"}
        result = prepare_spectrum_plot(mz, intensity, annotation_db=annotations)
        assert 100.0 in result.annotations
        assert result.annotations[100.0] == "fragment_A"

    def test_empty_spectrum(self) -> None:
        result = prepare_spectrum_plot(np.array([]), np.array([]))
        assert len(result.mz_values) == 0
        assert result.base_peak_mz == 0.0


class TestRetentionTimeAlignment:
    """Tests for RT alignment assessment."""

    def test_perfect_alignment(self) -> None:
        ref = np.array([10.0, 20.0, 30.0, 40.0, 50.0])
        sample = np.array([10.0, 20.0, 30.0, 40.0, 50.0])
        deviations, corrected, rmse = retention_time_alignment_data(ref, sample)
        np.testing.assert_allclose(deviations, 0.0, atol=1e-10)
        assert rmse < 0.01

    def test_constant_shift(self) -> None:
        ref = np.array([10.0, 20.0, 30.0, 40.0])
        sample = ref + 5.0
        deviations, corrected, rmse = retention_time_alignment_data(ref, sample)
        assert all(d == pytest.approx(5.0) for d in deviations)
        # After correction, RMSE should be small
        assert rmse < 1.0

    def test_linear_drift(self) -> None:
        ref = np.array([10.0, 20.0, 30.0, 40.0, 50.0])
        sample = ref + np.linspace(0, 10, 5)  # growing RT shift
        deviations, corrected, rmse = retention_time_alignment_data(ref, sample)
        assert deviations[0] < deviations[-1]  # drift increases
        # Polynomial correction should reduce error
        assert rmse < 5.0
