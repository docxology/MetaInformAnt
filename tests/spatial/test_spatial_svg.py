"""Tests for batch spatially variable gene detection (spatial analysis).

Real synthetic spatial expression matrices — a gradient gene (high Moran's I),
a random gene (near-expected I), and a constant gene (filtered) — verified
against the per-vector morans_i implementation.
"""

from __future__ import annotations

import numpy as np
import pytest

from metainformant.spatial.analysis.autocorrelation import morans_i, spatial_weights_matrix
from metainformant.spatial.analysis.spatially_variable_genes import (
    SVGResult,
    detect_spatially_variable_genes,
)


def _grid_coords(n_per_side: int = 8) -> np.ndarray:
    step = 1.0 / n_per_side
    xs = np.arange(n_per_side) * step
    ys = np.arange(n_per_side) * step
    xx, yy = np.meshgrid(xs, ys)
    return np.column_stack([xx.ravel(), yy.ravel()])


def _gradient_expression(coords: np.ndarray) -> np.ndarray:
    """Smooth spatial gradient — strongly spatially autocorrelated."""
    return (coords[:, 0] + coords[:, 1]).reshape(1, -1)


def test_gradient_gene_ranks_above_random_gene():
    rng = np.random.default_rng(42)
    coords = _grid_coords()
    gradient = _gradient_expression(coords)
    noise = rng.normal(size=(1, coords.shape[0]))

    expression = np.vstack([gradient, noise])
    result = detect_spatially_variable_genes(expression, coords, genes=["GRAD", "NOISE"], k=4)

    assert isinstance(result, SVGResult)
    assert result.n_spots == coords.shape[0]
    assert result.genes == ["GRAD", "NOISE"]  # sorted by Moran's I descending
    assert result.morans_i[0] > result.morans_i[1]
    assert result.morans_i[0] > 0.8  # near-perfect positive autocorrelation
    assert result.p_values[0] < 1e-6
    assert result.q_values is not None
    assert result.q_values[0] >= result.p_values[0]  # BH adjusts upward
    assert result.q_values[0] < 0.05  # still significant after adjustment
    # z-score consistent with the statistic
    assert result.z_scores[0] > 5


def test_matches_per_vector_morans_i():
    """Batch result must equal the single-vector morans_i on the same weights."""
    rng = np.random.default_rng(7)
    coords = _grid_coords(6)
    values = rng.normal(size=coords.shape[0])
    W = spatial_weights_matrix(coords, method="knn", k=4)
    expected = morans_i(values, W)

    result = detect_spatially_variable_genes(values.reshape(1, -1), coords, genes=["X"], k=4)
    assert result.genes == ["X"]
    assert result.morans_i[0] == pytest.approx(expected.I)
    assert result.z_scores[0] == pytest.approx(expected.z_score)
    assert result.p_values[0] == pytest.approx(expected.p_value)


def test_spot_by_gene_orientation_and_constant_filter():
    coords = _grid_coords(6)
    gradient = _gradient_expression(coords).ravel()
    constant = np.full(coords.shape[0], 5.0)
    expression = np.column_stack([gradient, constant])  # spots x genes

    result = detect_spatially_variable_genes(expression, coords, genes=["GRAD", "CONST"], k=4)
    # constant gene dropped (zero variance)
    assert result.genes == ["GRAD"]
    assert result.morans_i[0] > 0.8


def test_default_gene_names_and_top_genes():
    rng = np.random.default_rng(3)
    coords = _grid_coords(5)
    expression = rng.normal(size=(4, coords.shape[0]))
    result = detect_spatially_variable_genes(expression, coords, k=3)
    # Results are sorted by Moran's I descending, so default names come reordered.
    assert sorted(result.genes) == [f"Gene{i}" for i in range(4)]
    top = result.top_genes(2)
    assert len(top) == 2
    assert top[0] == result.genes[0]


def test_min_filters_can_drop_all_genes():
    coords = _grid_coords(4)
    expression = rng_less_simple(coords.shape[0])
    result = detect_spatially_variable_genes(expression, coords, k=3, min_mean_expression=1e9)
    assert result.genes == []
    assert result.morans_i.size == 0


def rng_less_simple(n: int) -> np.ndarray:
    return np.random.default_rng(11).normal(size=(1, n))


def test_shape_mismatch_raises():
    coords = _grid_coords(4)
    with pytest.raises(ValueError):
        detect_spatially_variable_genes(np.ones((3, 7)), coords, k=3)
    with pytest.raises(ValueError):
        detect_spatially_variable_genes(np.ones((2, coords.shape[0])), coords, genes=["a"], k=3)
