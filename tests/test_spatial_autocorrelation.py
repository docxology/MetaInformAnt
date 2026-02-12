"""Tests for spatial autocorrelation: Moran's I, Geary's C, LISA, Getis-Ord, variogram.

NO MOCKING POLICY: All tests use real implementations.
"""
from __future__ import annotations

import numpy as np
import pytest

from metainformant.spatial.analysis.autocorrelation import (
    GearyCResult,
    GetisOrdResult,
    LocalMoransResult,
    MoransIResult,
    VariogramResult,
    gearys_c,
    getis_ord_g,
    local_morans_i,
    morans_i,
    spatial_variogram,
    spatial_weights_matrix,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_grid_data(n=25, seed=42):
    """Create a 5x5 grid with spatially correlated values."""
    rng = np.random.RandomState(seed)
    side = int(np.sqrt(n))
    coords = np.array([[i, j] for i in range(side) for j in range(side)], dtype=np.float64)
    # Create spatially autocorrelated values: smooth gradient + noise
    values = np.array([float(i + j) for i in range(side) for j in range(side)]) + rng.randn(n) * 0.5
    return values, coords


def _make_random_data(n=30, seed=42):
    """Create random spatial data with no autocorrelation."""
    rng = np.random.RandomState(seed)
    coords = rng.rand(n, 2) * 100
    values = rng.randn(n)
    return values, coords


# ---------------------------------------------------------------------------
# spatial_weights_matrix
# ---------------------------------------------------------------------------


class TestSpatialWeightsMatrix:
    def test_knn_method(self):
        _, coords = _make_grid_data()
        W = spatial_weights_matrix(coords, method="knn", k=4)
        assert W.shape == (25, 25)
        assert W.nnz > 0

    def test_distance_method(self):
        _, coords = _make_grid_data()
        W = spatial_weights_matrix(coords, method="distance", bandwidth=2.0)
        assert W.shape == (25, 25)

    def test_binary_method(self):
        _, coords = _make_grid_data()
        W = spatial_weights_matrix(coords, method="binary", bandwidth=2.0)
        assert W.shape == (25, 25)

    def test_row_standardized(self):
        _, coords = _make_grid_data()
        W = spatial_weights_matrix(coords, method="knn", k=4, row_standardize=True)
        row_sums = np.array(W.sum(axis=1)).flatten()
        # Non-zero rows should sum to ~1
        nonzero_rows = row_sums[row_sums > 0]
        np.testing.assert_allclose(nonzero_rows, 1.0, atol=1e-10)

    def test_not_row_standardized(self):
        _, coords = _make_grid_data()
        W = spatial_weights_matrix(coords, method="knn", k=4, row_standardize=False)
        row_sums = np.array(W.sum(axis=1)).flatten()
        # Without standardization, sums should be integers (neighbor counts)
        assert np.all(row_sums >= 1)

    def test_unknown_method_raises(self):
        _, coords = _make_grid_data()
        with pytest.raises(ValueError, match="Unknown method"):
            spatial_weights_matrix(coords, method="invalid")

    def test_auto_bandwidth(self):
        _, coords = _make_random_data()
        W = spatial_weights_matrix(coords, method="distance")
        assert W.shape[0] == 30


# ---------------------------------------------------------------------------
# morans_i
# ---------------------------------------------------------------------------


class TestMoransI:
    def test_positive_autocorrelation(self):
        values, coords = _make_grid_data()
        W = spatial_weights_matrix(coords, method="knn", k=4)
        result = morans_i(values, W)
        assert isinstance(result, MoransIResult)
        # Gradient data should show positive autocorrelation
        assert result.I > 0

    def test_expected_I_negative(self):
        values, coords = _make_grid_data()
        W = spatial_weights_matrix(coords, method="knn", k=4)
        result = morans_i(values, W)
        assert result.expected_I < 0  # -1/(n-1)

    def test_p_value_range(self):
        values, coords = _make_grid_data()
        W = spatial_weights_matrix(coords, method="knn", k=4)
        result = morans_i(values, W)
        assert 0.0 <= result.p_value <= 1.0

    def test_n_correct(self):
        values, coords = _make_grid_data()
        W = spatial_weights_matrix(coords, method="knn", k=4)
        result = morans_i(values, W)
        assert result.n == 25

    def test_constant_values_zero_I(self):
        coords = np.random.RandomState(42).rand(20, 2)
        values = np.ones(20)
        W = spatial_weights_matrix(coords, method="knn", k=4)
        result = morans_i(values, W)
        assert result.I == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# gearys_c
# ---------------------------------------------------------------------------


class TestGearysC:
    def test_positive_autocorrelation(self):
        values, coords = _make_grid_data()
        W = spatial_weights_matrix(coords, method="knn", k=4)
        result = gearys_c(values, W)
        assert isinstance(result, GearyCResult)
        # Positive autocorrelation => C < 1
        assert result.C < 1.0

    def test_expected_C_is_one(self):
        values, coords = _make_grid_data()
        W = spatial_weights_matrix(coords, method="knn", k=4)
        result = gearys_c(values, W)
        assert result.expected_C == 1.0

    def test_p_value_range(self):
        values, coords = _make_grid_data()
        W = spatial_weights_matrix(coords, method="knn", k=4)
        result = gearys_c(values, W)
        assert 0.0 <= result.p_value <= 1.0

    def test_constant_values(self):
        coords = np.random.RandomState(42).rand(20, 2)
        values = np.ones(20) * 5.0
        W = spatial_weights_matrix(coords, method="knn", k=4)
        result = gearys_c(values, W)
        assert result.C == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# local_morans_i (LISA)
# ---------------------------------------------------------------------------


class TestLocalMoransI:
    def test_basic_lisa(self):
        values, coords = _make_grid_data()
        W = spatial_weights_matrix(coords, method="knn", k=4)
        result = local_morans_i(values, W)
        assert isinstance(result, LocalMoransResult)
        assert len(result.local_I) == 25
        assert len(result.cluster_labels) == 25

    def test_cluster_labels_valid(self):
        values, coords = _make_grid_data()
        W = spatial_weights_matrix(coords, method="knn", k=4)
        result = local_morans_i(values, W)
        valid_labels = {"HH", "LL", "HL", "LH", "NS"}
        for label in result.cluster_labels:
            assert label in valid_labels

    def test_p_values_range(self):
        values, coords = _make_grid_data()
        W = spatial_weights_matrix(coords, method="knn", k=4)
        result = local_morans_i(values, W)
        assert np.all(result.p_values >= 0.0)
        assert np.all(result.p_values <= 1.0)

    def test_custom_significance(self):
        values, coords = _make_grid_data()
        W = spatial_weights_matrix(coords, method="knn", k=4)
        result = local_morans_i(values, W, significance=0.01)
        assert result.significance_level == 0.01

    def test_constant_values_all_ns(self):
        coords = np.random.RandomState(42).rand(20, 2)
        values = np.ones(20)
        W = spatial_weights_matrix(coords, method="knn", k=4)
        result = local_morans_i(values, W)
        assert all(label == "NS" for label in result.cluster_labels)


# ---------------------------------------------------------------------------
# getis_ord_g
# ---------------------------------------------------------------------------


class TestGetisOrdG:
    def test_basic_getis(self):
        values, coords = _make_grid_data()
        W = spatial_weights_matrix(coords, method="knn", k=4, row_standardize=False)
        result = getis_ord_g(values, W)
        assert isinstance(result, GetisOrdResult)
        assert len(result.g_star) == 25

    def test_hot_cold_spots(self):
        values, coords = _make_grid_data()
        W = spatial_weights_matrix(coords, method="knn", k=4, row_standardize=False)
        result = getis_ord_g(values, W)
        assert isinstance(result.hot_spots, np.ndarray)
        assert isinstance(result.cold_spots, np.ndarray)

    def test_p_values_range(self):
        values, coords = _make_grid_data()
        W = spatial_weights_matrix(coords, method="knn", k=4, row_standardize=False)
        result = getis_ord_g(values, W)
        assert np.all(result.p_values >= 0.0)
        assert np.all(result.p_values <= 1.0)

    def test_constant_values_no_hotspots(self):
        coords = np.random.RandomState(42).rand(20, 2)
        values = np.ones(20)
        W = spatial_weights_matrix(coords, method="knn", k=4, row_standardize=False)
        result = getis_ord_g(values, W)
        assert not np.any(result.hot_spots)
        assert not np.any(result.cold_spots)


# ---------------------------------------------------------------------------
# spatial_variogram
# ---------------------------------------------------------------------------


class TestSpatialVariogram:
    def test_basic_variogram(self):
        values, coords = _make_grid_data()
        result = spatial_variogram(values, coords, n_bins=10)
        assert isinstance(result, VariogramResult)
        assert len(result.bin_centers) == 10
        assert len(result.semivariance) == 10

    def test_sill_and_nugget(self):
        values, coords = _make_grid_data()
        result = spatial_variogram(values, coords, n_bins=10)
        assert result.sill >= 0
        assert result.nugget >= 0

    def test_n_pairs_non_negative(self):
        values, coords = _make_grid_data()
        result = spatial_variogram(values, coords, n_bins=10)
        assert np.all(result.n_pairs >= 0)

    def test_custom_max_distance(self):
        values, coords = _make_grid_data()
        result = spatial_variogram(values, coords, n_bins=5, max_distance=3.0)
        assert result.bin_centers[-1] < 3.5

    def test_single_point(self):
        values = np.array([1.0])
        coords = np.array([[0.0, 0.0]])
        result = spatial_variogram(values, coords, n_bins=5)
        assert len(result.bin_centers) == 0

    def test_random_data_variogram(self):
        values, coords = _make_random_data(n=40)
        result = spatial_variogram(values, coords, n_bins=10)
        assert isinstance(result, VariogramResult)
        assert result.sill >= 0
