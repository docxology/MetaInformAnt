"""Spatial Leiden and Louvain never substitute the greedy modularity optimization silently.

``leiden_clustering`` and ``louvain_clustering`` used to run an in-package greedy optimization when
python-igraph/leidenalg or python-louvain/networkx could not be imported, logging one INFO line, and
``spatial_cluster`` then labelled the result ``method="leiden"`` or ``"louvain"`` anyway. Without the packages
the functions now raise unless ``allow_greedy_fallback=True``, which warns, and the result records the method
that ran. The environment decides which half of the contract can be observed, so each half is gated on what is
installed and the other half is skipped. Real data and real library calls only.
"""

from __future__ import annotations

import importlib
import warnings

import numpy as np
import pytest
from scipy import sparse

from metainformant.spatial.analysis import clustering
from metainformant.spatial.analysis.clustering import (
    build_spatial_graph,
    leiden_clustering,
    louvain_clustering,
    spatial_cluster,
)


def _importable(*names: str) -> bool:
    """Whether every module imports (an upstream docstring can make some leidenalg releases raise SyntaxError)."""
    for name in names:
        try:
            importlib.import_module(name)
        except (ImportError, SyntaxError):
            return False
    return True


LEIDEN_AVAILABLE = _importable("igraph", "leidenalg")
LOUVAIN_AVAILABLE = _importable("community", "networkx")


def _grid() -> tuple[np.ndarray, np.ndarray]:
    xs, ys = np.meshgrid(np.arange(8), np.arange(5))
    coordinates = np.column_stack([xs.ravel(), ys.ravel()]).astype(float)
    expression = np.random.default_rng(3).random((coordinates.shape[0], 10))
    return expression, coordinates


def _adjacency() -> sparse.csr_matrix:
    return build_spatial_graph(np.random.default_rng(2).uniform(0, 50, (30, 2)), method="knn", n_neighbors=4)


class TestWithoutTheLeidenPackages:
    pytestmark = pytest.mark.skipif(LEIDEN_AVAILABLE, reason="needs an environment without python-igraph/leidenalg")

    def test_leiden_clustering_raises_instead_of_running_another_algorithm(self) -> None:
        with pytest.raises(ImportError, match="allow_greedy_fallback"):
            leiden_clustering(_adjacency())

    def test_spatial_cluster_raises_too(self) -> None:
        expression, coordinates = _grid()

        with pytest.raises(ImportError, match="python-igraph and leidenalg"):
            spatial_cluster(expression, coordinates, method="leiden")

    def test_the_explicit_fallback_runs_the_greedy_optimization_and_says_so(self) -> None:
        with pytest.warns(UserWarning, match="greedy modularity"):
            labels, modularity = leiden_clustering(_adjacency(), allow_greedy_fallback=True)

        assert len(labels) == 30 and isinstance(modularity, float)
        greedy_labels, greedy_modularity = clustering._greedy_modularity_clustering(_adjacency(), seed=42)
        np.testing.assert_array_equal(labels, greedy_labels)
        assert modularity == greedy_modularity

    def test_the_result_records_that_the_greedy_optimization_ran(self) -> None:
        expression, coordinates = _grid()

        with pytest.warns(UserWarning, match="Leiden clustering is unavailable"):
            result = spatial_cluster(expression, coordinates, method="leiden", seed=1, allow_greedy_fallback=True)

        assert result.method == "greedy_modularity"
        assert result.metadata["requested_method"] == "leiden"
        assert result.metadata["fallback_used"] is True
        assert len(result.labels) == 40

    def test_the_warning_points_at_the_callers_line(self) -> None:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            leiden_clustering(_adjacency(), allow_greedy_fallback=True)

        relevant = [w for w in caught if "greedy modularity" in str(w.message)]
        assert len(relevant) == 1
        assert relevant[0].filename == __file__

    def test_zero_iterations_are_still_rejected_first(self) -> None:
        with pytest.raises(ValueError, match="n_iterations"):
            leiden_clustering(_adjacency(), n_iterations=0)


class TestWithoutTheLouvainPackages:
    pytestmark = pytest.mark.skipif(LOUVAIN_AVAILABLE, reason="needs an environment without python-louvain/networkx")

    def test_louvain_clustering_raises_instead_of_running_another_algorithm(self) -> None:
        with pytest.raises(ImportError, match="allow_greedy_fallback"):
            louvain_clustering(_adjacency())

    def test_spatial_cluster_raises_too(self) -> None:
        expression, coordinates = _grid()

        with pytest.raises(ImportError, match="python-louvain and networkx"):
            spatial_cluster(expression, coordinates, method="louvain")

    def test_the_explicit_fallback_runs_the_greedy_optimization_and_says_so(self) -> None:
        with pytest.warns(UserWarning, match="greedy modularity"):
            labels, _ = louvain_clustering(_adjacency(), allow_greedy_fallback=True)

        assert len(labels) == 30

    def test_the_result_records_that_the_greedy_optimization_ran(self) -> None:
        expression, coordinates = _grid()

        with pytest.warns(UserWarning, match="Louvain clustering is unavailable"):
            result = spatial_cluster(expression, coordinates, method="louvain", seed=1, allow_greedy_fallback=True)

        assert result.method == "greedy_modularity"
        assert result.metadata["requested_method"] == "louvain"
        assert result.metadata["fallback_used"] is True


class TestWithTheGraphPackages:
    @pytest.mark.skipif(not LEIDEN_AVAILABLE, reason="needs python-igraph and leidenalg")
    def test_the_flag_is_inert_and_the_result_says_leiden(self) -> None:
        expression, coordinates = _grid()

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = spatial_cluster(expression, coordinates, method="leiden", seed=1, allow_greedy_fallback=True)

        assert not [w for w in caught if "greedy" in str(w.message)]
        assert result.method == "leiden"
        assert result.metadata["requested_method"] == "leiden"
        assert result.metadata["fallback_used"] is False

    @pytest.mark.skipif(not LOUVAIN_AVAILABLE, reason="needs python-louvain and networkx")
    def test_louvain_flag_is_inert_and_the_result_says_louvain(self) -> None:
        expression, coordinates = _grid()

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = spatial_cluster(expression, coordinates, method="louvain", seed=1, allow_greedy_fallback=True)

        assert not [w for w in caught if "greedy" in str(w.message)]
        assert result.method == "louvain"
        assert result.metadata["fallback_used"] is False


class TestEveryEnvironment:
    def test_kmeans_needs_no_graph_package_and_is_not_a_fallback(self) -> None:
        expression, coordinates = _grid()

        result = spatial_cluster(expression, coordinates, n_clusters=3, method="kmeans", seed=0)

        assert result.method == "kmeans"
        assert "fallback_used" not in result.metadata

    def test_the_greedy_optimization_separates_two_cliques_joined_by_one_edge(self) -> None:
        """The fallback itself is a working optimizer; it is only never chosen without being asked for."""
        n = 10
        block = np.ones((n, n)) - np.eye(n)
        dense = np.zeros((2 * n, 2 * n))
        dense[:n, :n] = block
        dense[n:, n:] = block
        dense[n - 1, n] = dense[n, n - 1] = 1.0

        labels, modularity = clustering._greedy_modularity_clustering(sparse.csr_matrix(dense), seed=0)

        assert len(set(labels[:n])) == 1 and len(set(labels[n:])) == 1
        assert labels[0] != labels[-1]
        assert modularity > 0.4

    def test_the_default_is_to_refuse(self) -> None:
        import inspect

        for function in (leiden_clustering, louvain_clustering, spatial_cluster):
            parameter = inspect.signature(function).parameters["allow_greedy_fallback"]
            assert parameter.default is False
            assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
