"""Spatial statistics that cannot be computed are NaN with a status, never a made-up value.

Moran's I, Geary's C, local Moran's I and Getis-Ord G* used to return I = 0 (or C = 1, G* = 0) and p = 1 for
constant values or an empty weights matrix, local Moran's I clamped a non-positive variance to 1e-20 (so an
observation without neighbors got a z-score near 1e9 and p = 0), and the variogram reported a semivariance of 0
for bins that hold no pair. A made-up p of 1 also entered the Benjamini-Hochberg family of the spatially variable
genes. Real data and real calls only; the definitions are re-derived here with numpy for the regular cases.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy import sparse, stats as scipy_stats

from metainformant.spatial.analysis.autocorrelation import (
    gearys_c,
    getis_ord_g,
    local_morans_i,
    morans_i,
    spatial_variogram,
    spatial_weights_matrix,
)
from metainformant.spatial.analysis.spatially_variable_genes import detect_spatially_variable_genes


def _grid_weights(side: int = 5, row_standardize: bool = False) -> tuple[np.ndarray, sparse.csr_matrix]:
    xs, ys = np.meshgrid(np.arange(side), np.arange(side))
    coords = np.column_stack([xs.ravel(), ys.ravel()]).astype(float)
    return coords, spatial_weights_matrix(coords, method="knn", k=4, row_standardize=row_standardize)


def _gradient(coords: np.ndarray, seed: int = 0) -> np.ndarray:
    return coords[:, 0] + coords[:, 1] + 0.1 * np.random.default_rng(seed).standard_normal(len(coords))


class TestMoransI:
    def test_regular_input_is_unchanged_and_ok(self) -> None:
        coords, weights = _grid_weights()
        values = _gradient(coords)

        result = morans_i(values, weights)

        z = values - values.mean()
        dense = weights.toarray()
        expected = len(values) / dense.sum() * (z @ dense @ z) / (z @ z)
        assert result.I == pytest.approx(expected, rel=1e-12)
        assert result.status == "ok"
        assert 0.0 <= result.p_value < 0.01

    def test_constant_values_have_no_statistic_not_zero(self) -> None:
        coords, weights = _grid_weights()

        result = morans_i(np.full(len(coords), 3.0), weights)

        assert math.isnan(result.I) and math.isnan(result.z_score) and math.isnan(result.p_value)
        assert math.isnan(result.variance_I)
        assert result.expected_I == pytest.approx(-1.0 / (len(coords) - 1))  # the null expectation is still defined
        assert "constant" in result.status

    def test_a_weights_matrix_without_edges_has_no_statistic(self) -> None:
        coords, _ = _grid_weights()

        result = morans_i(_gradient(coords), sparse.csr_matrix((len(coords), len(coords))))

        assert math.isnan(result.I) and math.isnan(result.p_value)
        assert "no edges" in result.status

    def test_three_observations_have_a_statistic_but_no_variance(self) -> None:
        weights = np.array([[0.0, 1.0, 0.0], [1.0, 0.0, 1.0], [0.0, 1.0, 0.0]])

        result = morans_i(np.array([1.0, 2.0, 4.0]), weights)

        assert math.isfinite(result.I)
        assert math.isnan(result.variance_I) and math.isnan(result.z_score) and math.isnan(result.p_value)
        assert "fewer than 4" in result.status

    def test_a_single_observation(self) -> None:
        result = morans_i(np.array([1.0]), np.zeros((1, 1)))

        assert math.isnan(result.I) and math.isnan(result.expected_I) and math.isnan(result.p_value)
        assert "fewer than 2" in result.status


class TestGearysC:
    def test_regular_input_is_ok(self) -> None:
        coords, weights = _grid_weights()
        values = _gradient(coords)

        result = gearys_c(values, weights)

        dense = weights.toarray()
        z = values - values.mean()
        squared = (values[:, None] - values[None, :]) ** 2
        expected = (len(values) - 1) * (dense * squared).sum() / (2.0 * dense.sum() * (z @ z))
        assert result.C == pytest.approx(expected, rel=1e-12)
        assert result.status == "ok" and 0.0 <= result.p_value < 0.01

    def test_constant_values_have_no_statistic_not_c_equal_one(self) -> None:
        coords, weights = _grid_weights()

        result = gearys_c(np.full(len(coords), 5.0), weights)

        assert math.isnan(result.C) and math.isnan(result.z_score) and math.isnan(result.p_value)
        assert result.expected_C == 1.0
        assert "constant" in result.status

    def test_a_weights_matrix_without_edges_has_no_statistic(self) -> None:
        coords, _ = _grid_weights()

        result = gearys_c(_gradient(coords), np.zeros((len(coords), len(coords))))

        assert math.isnan(result.C) and math.isnan(result.p_value)
        assert "no edges" in result.status

    def test_three_observations_have_a_statistic_but_no_variance(self) -> None:
        weights = np.array([[0.0, 1.0, 0.0], [1.0, 0.0, 1.0], [0.0, 1.0, 0.0]])

        result = gearys_c(np.array([1.0, 2.0, 4.0]), weights)

        assert math.isfinite(result.C)
        assert math.isnan(result.variance_C) and math.isnan(result.p_value)
        assert "fewer than 4" in result.status


class TestLocalMoransI:
    @staticmethod
    def _ring_with_an_isolated_spot() -> np.ndarray:
        """Spots 0-6 form a ring (row-standardized); spot 7 has no neighbor at all."""
        weights = np.zeros((8, 8))
        for i in range(7):
            weights[i, (i + 1) % 7] = weights[i, (i - 1) % 7] = 0.5
        return weights

    def test_a_spot_without_neighbors_has_no_z_score_and_no_p_value(self) -> None:
        values = np.random.default_rng(1).standard_normal(8)

        result = local_morans_i(values, self._ring_with_an_isolated_spot())

        assert math.isnan(result.z_scores[7]) and math.isnan(result.p_values[7])  # was z ~ 1e9 and p = 0
        assert result.cluster_labels[7] == "NS"
        assert np.isfinite(result.z_scores[:7]).all() and np.isfinite(result.p_values[:7]).all()
        assert result.status == "ok"

    def test_a_nan_p_value_never_becomes_a_cluster(self) -> None:
        values = np.zeros(8)
        values[7] = 50.0  # a huge value at the isolated spot would be "HH" if its NaN p-value counted as small
        values[:7] = np.random.default_rng(2).standard_normal(7)

        result = local_morans_i(values, self._ring_with_an_isolated_spot())

        assert result.cluster_labels[7] == "NS"

    def test_constant_values_leave_every_statistic_undefined(self) -> None:
        coords, weights = _grid_weights(row_standardize=True)

        result = local_morans_i(np.ones(len(coords)), weights)

        assert np.isnan(result.local_I).all() and np.isnan(result.z_scores).all() and np.isnan(result.p_values).all()
        assert set(result.cluster_labels) == {"NS"}
        assert "constant" in result.status

    def test_regular_input_has_finite_statistics(self) -> None:
        coords, weights = _grid_weights(row_standardize=True)

        result = local_morans_i(_gradient(coords), weights)

        assert np.isfinite(result.z_scores).all() and np.isfinite(result.p_values).all()
        assert ((0.0 <= result.p_values) & (result.p_values <= 1.0)).all()
        assert result.status == "ok"

    def test_fewer_than_three_observations(self) -> None:
        result = local_morans_i(np.array([1.0, 2.0]), np.array([[0.0, 1.0], [1.0, 0.0]]))

        assert np.isnan(result.p_values).all()
        assert "fewer than 3" in result.status


class TestGetisOrdG:
    def test_a_complete_graph_has_no_spread_to_standardize_by(self) -> None:
        """Every row of W* is all ones, so n sum w^2 - (sum w)^2 = 0: G* is undefined, not 0 with p = 1."""
        n = 10
        weights = np.ones((n, n)) - np.eye(n)
        values = np.random.default_rng(3).standard_normal(n)

        result = getis_ord_g(values, weights)

        assert np.isnan(result.g_star).all() and np.isnan(result.p_values).all()
        assert not result.hot_spots.any() and not result.cold_spots.any()

    def test_constant_values_leave_every_statistic_undefined(self) -> None:
        coords, weights = _grid_weights()

        result = getis_ord_g(np.full(len(coords), 2.0), weights)

        assert np.isnan(result.g_star).all() and np.isnan(result.p_values).all()
        assert "constant" in result.status
        assert not result.hot_spots.any()

    def test_regular_input_finds_the_hot_corner(self) -> None:
        coords, weights = _grid_weights(side=6)
        values = np.where((coords[:, 0] >= 4) & (coords[:, 1] >= 4), 10.0, 0.0)
        values = values + 0.01 * np.random.default_rng(4).standard_normal(len(values))

        result = getis_ord_g(values, weights)

        assert np.isfinite(result.g_star).all()
        assert result.hot_spots[(coords[:, 0] == 5) & (coords[:, 1] == 5)].all()
        assert result.status == "ok"


class TestVariogram:
    def test_bins_without_pairs_have_no_semivariance_not_zero(self) -> None:
        # Two tight groups 100 apart; with the default maximum distance (half the largest distance) only the
        # within-group pairs are binned, so every bin beyond the first few holds no pair
        coords = np.array([[0.0, 0.0], [0.5, 0.0], [1.0, 0.0], [100.0, 0.0], [100.5, 0.0], [101.0, 0.0]])
        values = np.array([1.0, 2.0, 3.0, 10.0, 11.0, 12.0])

        result = spatial_variogram(values, coords, n_bins=10)

        empty = result.n_pairs == 0
        assert empty.any() and (~empty).any()
        assert np.isnan(result.semivariance[empty]).all()
        assert np.isfinite(result.semivariance[~empty]).all() and (result.semivariance[~empty] > 0).all()

    def test_constant_values_have_a_sill_of_zero_and_no_range(self) -> None:
        coords = np.random.default_rng(5).uniform(0, 10, (30, 2))

        result = spatial_variogram(np.full(30, 4.0), coords, n_bins=6)

        assert result.sill == 0.0
        assert math.isnan(result.range_param)  # no range can be read off a flat variogram; it was 0.0

    def test_fewer_than_two_points_have_no_parameters(self) -> None:
        result = spatial_variogram(np.array([1.0]), np.array([[0.0, 0.0]]), n_bins=4)

        assert len(result.bin_centers) == 0
        assert math.isnan(result.nugget) and math.isnan(result.sill) and math.isnan(result.range_param)


class TestSpatiallyVariableGenes:
    @staticmethod
    def _slide(n_genes: int = 6, side: int = 6) -> tuple[np.ndarray, np.ndarray]:
        xs, ys = np.meshgrid(np.arange(side), np.arange(side))
        coords = np.column_stack([xs.ravel(), ys.ravel()]).astype(float)
        rng = np.random.default_rng(6)
        expression = rng.poisson(5.0, (n_genes, len(coords))).astype(float)
        expression[0] += 3.0 * coords[:, 0]  # one spatially structured gene
        return expression, coords

    def test_a_weights_matrix_without_edges_gives_nan_p_and_q_values(self) -> None:
        expression, coords = self._slide()

        result = detect_spatially_variable_genes(expression, coords, k=0)

        assert len(result.genes) == 6
        assert np.isnan(result.p_values).all()
        assert result.q_values is not None and np.isnan(result.q_values).all()

    def test_a_gene_without_a_p_value_is_outside_the_benjamini_hochberg_family(self) -> None:
        expression, coords = self._slide()
        expression[3] = 4.0  # constant; min_variance=-1 lets it through the variance filter

        with_constant = detect_spatially_variable_genes(expression, coords, genes=list("abcdef"), min_variance=-1.0)
        without = detect_spatially_variable_genes(
            np.delete(expression, 3, axis=0), coords, genes=list("abcef"), min_variance=-1.0
        )

        index = with_constant.genes.index("d")
        assert math.isnan(with_constant.p_values[index]) and math.isnan(with_constant.q_values[index])
        assert index == len(with_constant.genes) - 1  # ranked last: its Moran's I is NaN
        finite = ~np.isnan(with_constant.p_values)
        expected = scipy_stats.false_discovery_control(with_constant.p_values[finite], method="bh")
        np.testing.assert_allclose(with_constant.q_values[finite], expected, rtol=1e-12)
        # ... and the other genes' q-values are those of the data without the constant gene
        by_gene = dict(zip(with_constant.genes, with_constant.q_values))
        for gene, q in zip(without.genes, without.q_values):
            assert by_gene[gene] == pytest.approx(q, rel=1e-12)

    def test_regular_genes_have_q_values_and_the_structured_one_is_top(self) -> None:
        expression, coords = self._slide()

        result = detect_spatially_variable_genes(expression, coords, genes=list("abcdef"))

        assert result.top_genes(1) == ["a"]
        assert result.q_values is not None and np.isfinite(result.q_values).all()
        assert result.p_values[0] < 1e-6 and result.q_values[0] >= result.p_values[0]
