"""Randomness audit of the single-cell analysis paths.

Library code must not read or write numpy's global random state, every random draw must be controllable
through a seed, and a numerical failure must not be papered over with random numbers. Real data and real
library calls only.
"""

from __future__ import annotations

from typing import Any, Callable

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from metainformant.singlecell.analysis.clustering import compute_cluster_silhouette, hierarchical_clustering
from metainformant.singlecell.analysis.nonlinear_methods import compute_dimensionality_metrics, diffusion_map_reduction
from metainformant.singlecell.data.preprocessing import SingleCellData
from metainformant.singlecell.visualization.visualization import plot_marker_expression

matplotlib.use("Agg")  # Non-interactive backend for tests; selected before any figure is created


@pytest.fixture(autouse=True)
def _close_plots():
    yield
    plt.close("all")


def _global_state() -> tuple[Any, ...]:
    return np.random.get_state()


def _same_state(a: tuple[Any, ...], b: tuple[Any, ...]) -> bool:
    return a[0] == b[0] and np.array_equal(a[1], b[1]) and a[2:] == b[2:]


def _structureless(n_cells: int, n_genes: int = 5, seed: int = 0) -> SingleCellData:
    """Uniform noise: any change in which cells are drawn changes the result."""
    X = np.random.default_rng(seed).random((n_cells, n_genes))
    obs = pd.DataFrame({"cluster": [f"c{i % 3}" for i in range(n_cells)]}, index=[f"cell{i}" for i in range(n_cells)])
    return SingleCellData(X, obs=obs, var=pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)]))


def _leaves_global_state_alone(call: Callable[[], Any]) -> bool:
    np.random.seed(12345)
    before = _global_state()
    call()
    return _same_state(before, _global_state())


class TestHierarchicalClusteringSubsample:
    """Above 1000 cells the linkage is built on a random subset of cells."""

    def test_the_same_seed_gives_the_same_clusters(self) -> None:
        data = _structureless(1200)

        first = hierarchical_clustering(data, n_clusters=4, random_state=7)
        again = hierarchical_clustering(data, n_clusters=4, random_state=7)

        np.testing.assert_array_equal(first.obs["hierarchical_cluster"], again.obs["hierarchical_cluster"])
        assert first.uns["hierarchical_clustering"]["random_state"] == 7

    def test_a_different_seed_draws_a_different_subset(self) -> None:
        data = _structureless(1200)

        first = hierarchical_clustering(data, n_clusters=4, random_state=7)
        other = hierarchical_clustering(data, n_clusters=4, random_state=8)

        assert not np.array_equal(first.obs["hierarchical_cluster"], other.obs["hierarchical_cluster"])

    def test_the_global_random_state_is_neither_read_nor_advanced(self) -> None:
        data = _structureless(1200)

        assert _leaves_global_state_alone(lambda: hierarchical_clustering(data, n_clusters=4, random_state=7))
        assert _leaves_global_state_alone(lambda: hierarchical_clustering(data, n_clusters=4))


class TestSilhouetteSubsample:
    @staticmethod
    def _labelled(n_cells: int = 200) -> SingleCellData:
        data = _structureless(n_cells)
        data.obs["cluster"] = [i % 4 for i in range(n_cells)]
        return data

    def test_the_same_seed_gives_the_same_score(self) -> None:
        data = self._labelled()

        first = compute_cluster_silhouette(data, "cluster", sample_size=50, random_state=3)
        again = compute_cluster_silhouette(data, "cluster", sample_size=50, random_state=3)

        assert first["overall_silhouette_score"] == again["overall_silhouette_score"]

    def test_a_different_seed_draws_a_different_subsample(self) -> None:
        data = self._labelled()

        scores = {
            compute_cluster_silhouette(data, "cluster", sample_size=50, random_state=seed)["overall_silhouette_score"]
            for seed in range(5)
        }

        assert len(scores) > 1

    def test_the_global_random_state_is_neither_read_nor_advanced(self) -> None:
        data = self._labelled()

        assert _leaves_global_state_alone(lambda: compute_cluster_silhouette(data, "cluster", sample_size=50))


class TestDimensionalityMetricsSubsample:
    @staticmethod
    def _with_embedding(n_cells: int = 1100) -> SingleCellData:
        data = _structureless(n_cells, n_genes=6)
        embedding = np.random.default_rng(1).random((n_cells, 2))
        data.obs["Emb1"], data.obs["Emb2"] = embedding[:, 0], embedding[:, 1]
        return data

    def test_the_same_seed_gives_the_same_distance_correlation(self) -> None:
        data = self._with_embedding()

        first = compute_dimensionality_metrics(data, ["Emb1", "Emb2"], random_state=5)
        again = compute_dimensionality_metrics(data, ["Emb1", "Emb2"], random_state=5)
        other = compute_dimensionality_metrics(data, ["Emb1", "Emb2"], random_state=6)

        assert first["distance_correlation"] == again["distance_correlation"]
        assert first["distance_correlation"] != other["distance_correlation"]

    def test_the_global_random_state_is_neither_read_nor_advanced(self) -> None:
        data = self._with_embedding(n_cells=80)

        assert _leaves_global_state_alone(lambda: compute_dimensionality_metrics(data, ["Emb1", "Emb2"]))


class TestMarkerHeatmapSubsample:
    @staticmethod
    def _shown(data: SingleCellData, **kwargs: Any) -> np.ndarray:
        fig = plot_marker_expression(data, marker_genes=["gene_0", "gene_1"], method="heatmap", max_cells=20, **kwargs)
        return np.asarray(fig.axes[0].images[0].get_array())

    def test_the_same_seed_shows_the_same_cells(self) -> None:
        data = _structureless(100)

        first = self._shown(data, random_state=4)
        again = self._shown(data, random_state=4)

        assert first.shape == (2, 20)
        np.testing.assert_array_equal(first, again)

    def test_a_different_seed_shows_other_cells(self) -> None:
        data = _structureless(100)

        assert not np.array_equal(self._shown(data, random_state=4), self._shown(data, random_state=5))

    def test_the_global_random_state_is_neither_read_nor_advanced(self) -> None:
        data = _structureless(100)

        assert _leaves_global_state_alone(lambda: self._shown(data))


class TestDiffusionMapReduction:
    """A failed eigendecomposition used to be replaced by random vectors and a re-seeded global RNG.

    Whether LAPACK reports failure for a non-finite matrix depends on the library build and the matrix
    size, so the tests first check that it does for the exact all-NaN matrix the function will build
    (a NaN anywhere makes the kernel width, hence every kernel entry, NaN) and skip where it does not.
    """

    N_CELLS = 6

    @classmethod
    def _non_finite_data(cls) -> SingleCellData:
        data = _structureless(cls.N_CELLS, n_genes=3)
        data.X[2, 1] = np.nan
        try:
            np.linalg.eigh(np.full((cls.N_CELLS, cls.N_CELLS), np.nan))
        except np.linalg.LinAlgError:
            return data
        pytest.skip("this LAPACK build returns NaN instead of failing on a non-finite matrix")

    def test_a_failed_eigendecomposition_raises_instead_of_returning_random_vectors(self) -> None:
        data = self._non_finite_data()

        with pytest.raises(np.linalg.LinAlgError):
            diffusion_map_reduction(data, n_components=2, n_neighbors=5)

    def test_the_global_random_state_is_not_reseeded(self) -> None:
        data = self._non_finite_data()

        def run() -> None:
            with pytest.raises(np.linalg.LinAlgError):
                diffusion_map_reduction(data, n_components=2, n_neighbors=5)

        assert _leaves_global_state_alone(run)
