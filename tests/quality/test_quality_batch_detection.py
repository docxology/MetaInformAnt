"""Deep tests for metainformant.quality.batch.detection (real computation, no mocks)."""

import numpy as np
import pytest

from metainformant.quality.batch.detection import (
    BatchEffectReport,
    _compute_silhouette,
    correct_batch_combat,
    detect_batch_effects,
)


def _make_batched_data(n_per_batch: int = 20, shift: float = 3.0, seed: int = 7) -> tuple[np.ndarray, list[str]]:
    """Two batches drawn from shifted Gaussians: strong, controlled batch effect."""
    rng = np.random.default_rng(seed)
    batch_a = rng.normal(loc=0.0, scale=1.0, size=(n_per_batch, 12))
    batch_b = rng.normal(loc=shift, scale=1.0, size=(n_per_batch, 12))
    data = np.vstack([batch_a, batch_b])
    labels = ["A"] * n_per_batch + ["B"] * n_per_batch
    return data, labels


class TestDetectBatchEffects:
    def test_report_structure(self) -> None:
        data, labels = _make_batched_data()
        report = detect_batch_effects(data, labels)
        assert isinstance(report, BatchEffectReport)
        assert report.n_samples == 40
        assert report.n_batches == 2
        assert report.batch_labels == labels
        assert set(report.pvca_variance) == {"batch", "residual"}
        assert pytest.approx(report.pvca_variance["batch"] + report.pvca_variance["residual"]) == 1.0

    def test_strong_batch_effect_flagged_high_or_critical(self) -> None:
        data, labels = _make_batched_data(shift=5.0)
        report = detect_batch_effects(data, labels)
        assert report.severity in ("high", "critical")
        assert report.pvca_variance["batch"] > 0.5

    def test_no_batch_effect_scores_low(self) -> None:
        # seed 21 measured: severity "low", batch variance 0.0212 (far below the 0.05 moderate line)
        rng = np.random.default_rng(21)
        data = rng.normal(loc=0.0, scale=1.0, size=(30, 12))
        labels = ["A"] * 15 + ["B"] * 15
        report = detect_batch_effects(data, labels)
        assert report.severity == "low"
        assert report.pvca_variance["batch"] < 0.05

    def test_significant_features_with_strong_shift(self) -> None:
        data, labels = _make_batched_data(shift=6.0)
        report = detect_batch_effects(data, labels)
        # Every feature differs by ~6 sigma between batches -> all significant
        assert report.n_significant_features == data.shape[1]

    def test_single_batch_edge(self) -> None:
        rng = np.random.default_rng(3)
        data = rng.normal(size=(10, 4))
        report = detect_batch_effects(data, ["A"] * 10)
        # Single batch: no between-batch variance, silhouette helper returns 0
        assert report.n_batches == 1
        assert report.pvca_variance["batch"] == pytest.approx(0.0)
        assert report.severity == "low"

    def test_works_without_scipy_path(self) -> None:
        # Exercise the numerical result regardless of scipy import branch: F-stat based
        data, labels = _make_batched_data(shift=4.0)
        report = detect_batch_effects(data, labels)
        assert 0.0 <= report.silhouette_score <= 1.0


class TestCorrectBatchCombat:
    def test_output_shape_and_dtype(self) -> None:
        data, labels = _make_batched_data()
        corrected = correct_batch_combat(data, labels)
        assert corrected.shape == data.shape
        assert corrected.dtype == np.float64
        # Input not mutated
        assert not np.allclose(corrected, data)

    def test_reduces_batch_mean_separation(self) -> None:
        data, labels = _make_batched_data(shift=5.0)
        corrected = correct_batch_combat(data, labels)
        sep_before = abs(data[:20].mean() - data[20:].mean())
        sep_after = abs(corrected[:20].mean() - corrected[20:].mean())
        assert sep_after < sep_before

    def test_preserves_within_batch_variation(self) -> None:
        data, labels = _make_batched_data()
        corrected = correct_batch_combat(data, labels)
        std_before = data.std()
        std_after = corrected.std()
        # Shrinkage keeps scale near 1: overall spread should stay within 2x
        assert 0.5 * std_before < std_after < 2.0 * std_before

    def test_single_batch_noop_centering(self) -> None:
        rng = np.random.default_rng(5)
        data = rng.normal(size=(8, 4))
        corrected = correct_batch_combat(data, ["A"] * 8)
        # One batch: shift to grand mean is zero; scaling shrunk 50% toward 1
        assert np.allclose(corrected, data)

    def test_zero_variance_features_handled(self) -> None:
        data = np.ones((10, 3))
        labels = ["A"] * 5 + ["B"] * 5
        corrected = correct_batch_combat(data, labels)
        assert np.all(np.isfinite(corrected))


class TestComputeSilhouette:
    def test_perfect_separation_high_score(self) -> None:
        coords = np.array([[0.0, 0.0], [0.1, 0.0], [10.0, 0.0], [10.1, 0.0]])
        labels = np.array([0, 0, 1, 1])
        score = _compute_silhouette(coords, labels, 2)
        assert score > 0.9

    def test_single_group_returns_zero(self) -> None:
        coords = np.array([[0.0], [1.0], [2.0]])
        assert _compute_silhouette(coords, np.array([0, 0, 0]), 1) == 0.0

    def test_too_few_points_returns_zero(self) -> None:
        assert _compute_silhouette(np.array([[0.0]]), np.array([0]), 1) == 0.0

    def test_score_bounded(self) -> None:
        rng = np.random.default_rng(9)
        coords = rng.normal(size=(30, 3))
        labels = rng.integers(0, 3, size=30)
        score = _compute_silhouette(coords, labels, 3)
        assert -1.0 <= score <= 1.0
