"""Tests for spatial autocorrelation: Moran's I, Geary's C, LISA, Getis-Ord, variogram.

REAL IMPLEMENTATION POLICY: All tests use real implementations.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest
from scipy import sparse as sp

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

    def test_constant_values_have_no_morans_i(self):
        """0 / 0 is undefined: NaN and a status, not an I of 0 and a p-value of 1 that read as a measurement."""
        coords = np.random.RandomState(42).rand(20, 2)
        values = np.ones(20)
        W = spatial_weights_matrix(coords, method="knn", k=4)
        result = morans_i(values, W)
        assert np.isnan(result.I) and np.isnan(result.z_score) and np.isnan(result.p_value)
        assert result.expected_I == pytest.approx(-1.0 / 19)
        assert "constant" in result.status


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
        assert np.isnan(result.C) and np.isnan(result.z_score) and np.isnan(result.p_value)
        assert result.expected_C == 1.0  # the null expectation is still defined
        assert "constant" in result.status


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
# Sparse weight formats
# ---------------------------------------------------------------------------

SPARSE_FORMATS = ["csr", "csc", "coo", "lil", "dok", "dia", "bsr"]
SPARSE_FAMILIES = ["matrix", "array"]


def _sparse_weights(weights, fmt: str, family: str):
    """Convert CSR weights to the requested scipy sparse format and container family."""
    return getattr(sp, f"{fmt}_{family}")(weights)


def _gearys_c_from_definition(values, dense_weights) -> float:
    """Geary's C computed from its definition on a dense matrix (independent of the library code path)."""
    x = np.asarray(values, dtype=np.float64)
    n = len(x)
    squared_differences = (x[:, None] - x[None, :]) ** 2
    numerator = float((dense_weights * squared_differences).sum())
    denominator = float(((x - x.mean()) ** 2).sum())
    return (n - 1) * numerator / (2.0 * float(dense_weights.sum()) * denominator)


def _assert_results_close(actual, expected) -> None:
    """Compare two result dataclasses field by field (labels exactly, numbers to 1e-10 relative)."""
    for field in dataclasses.fields(expected):
        got, want = getattr(actual, field.name), getattr(expected, field.name)
        if isinstance(want, (list, str)):
            assert got == want, field.name
        else:
            np.testing.assert_allclose(np.asarray(got), np.asarray(want), rtol=1e-10, atol=1e-12, err_msg=field.name)


class TestSparseWeightFormats:
    """Every statistic accepts any scipy sparse format and container family, not only CSR."""

    @pytest.mark.parametrize("family", SPARSE_FAMILIES)
    @pytest.mark.parametrize("fmt", SPARSE_FORMATS)
    def test_gearys_c_matches_definition(self, fmt: str, family: str) -> None:
        # Regression: gearys_c indexed W[sources, targets], which raises TypeError or
        # NotImplementedError for COO, DIA, BSR, LIL and DOK weights (matrix and array families).
        values, coords = _make_grid_data()
        csr = spatial_weights_matrix(coords, method="knn", k=4, row_standardize=False)
        result = gearys_c(values, _sparse_weights(csr, fmt, family))
        assert result.C == pytest.approx(_gearys_c_from_definition(values, csr.toarray()), rel=1e-12)

    @pytest.mark.parametrize("family", SPARSE_FAMILIES)
    @pytest.mark.parametrize("fmt", SPARSE_FORMATS)
    def test_all_statistics_independent_of_sparse_format(self, fmt: str, family: str) -> None:
        values, coords = _make_grid_data()
        csr = spatial_weights_matrix(coords, method="knn", k=4, row_standardize=False)
        weights = _sparse_weights(csr, fmt, family)

        for statistic in (morans_i, gearys_c, local_morans_i, getis_ord_g):
            _assert_results_close(statistic(values, weights), statistic(values, csr))

    def test_gearys_c_dense_weights_match_sparse(self) -> None:
        values, coords = _make_grid_data()
        csr = spatial_weights_matrix(coords, method="knn", k=4)
        _assert_results_close(gearys_c(values, csr.toarray()), gearys_c(values, csr))

    def test_gearys_c_sums_duplicate_coo_entries(self) -> None:
        # A COO matrix may store one weight as several entries; the matrix value is their sum.
        values, coords = _make_grid_data()
        csr = spatial_weights_matrix(coords, method="knn", k=4, row_standardize=False)
        coo = csr.tocoo()
        split = sp.coo_matrix(
            (
                np.concatenate([coo.data / 2.0, coo.data / 2.0]),
                (np.concatenate([coo.row, coo.row]), np.concatenate([coo.col, coo.col])),
            ),
            shape=csr.shape,
        )
        assert split.nnz == 2 * csr.nnz
        _assert_results_close(gearys_c(values, split), gearys_c(values, csr))


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

    def test_kdtree_path_large_dataset(self):
        # n > 500 exercises the KDTree branch of spatial_variogram
        rng = np.random.RandomState(0)
        coords = rng.uniform(0, 100, (510, 2))
        values = 0.5 * (coords[:, 0] + coords[:, 1]) + rng.normal(0, 0.1, 510)
        result = spatial_variogram(values, coords, n_bins=5)
        assert isinstance(result, VariogramResult)
        assert len(result.bin_centers) == 5
        assert int(result.n_pairs.sum()) > 0
        assert result.sill > 0
