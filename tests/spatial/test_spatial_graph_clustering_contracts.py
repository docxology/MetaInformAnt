"""Contracts of the spatial graph and clustering helpers.

- kNN graph weights depend on the distance of a pair only, not on whether the neighbor relation is mutual;
- ``spatial_cluster(graph_method="radius")`` takes a ``radius`` and runs;
- ``leiden_clustering`` hands leidenalg the iteration count it documents (-1 = until convergence).

Real data and real library calls only; nothing is replaced by a test double.
"""

from __future__ import annotations

import numpy as np
import pytest

from metainformant.spatial.analysis.clustering import build_spatial_graph, leiden_clustering, spatial_cluster

try:
    import igraph  # noqa: F401
    import leidenalg  # noqa: F401

    LEIDENALG_AVAILABLE = True
except (ImportError, SyntaxError):
    LEIDENALG_AVAILABLE = False


def _leidenalg_iterations(n_iterations: int) -> int:
    """The count handed to leidenalg for a public ``n_iterations`` (resolved lazily: a private helper)."""
    from metainformant.spatial.analysis import clustering

    return int(clustering._leidenalg_iterations(n_iterations))


class TestKnnWeightsAreNotDoubled:
    def test_mutual_and_one_way_pairs_get_their_own_distance_weight(self) -> None:
        # Spots 0 and 1 are each other's nearest neighbor (mutual); spot 2 lists spot 1 but not vice versa.
        coordinates = np.array([[0.0, 0.0], [10.0, 0.0], [100.0, 0.0]])

        adjacency = build_spatial_graph(coordinates, method="knn", n_neighbors=1).toarray()

        assert adjacency[0, 1] == adjacency[1, 0] == pytest.approx(1.0 / 11.0)  # distance 10, mutual
        assert adjacency[1, 2] == adjacency[2, 1] == pytest.approx(1.0 / 91.0)  # distance 90, one-way
        assert adjacency[0, 2] == 0.0

    def test_every_edge_weight_is_a_function_of_its_endpoint_distance(self) -> None:
        coordinates = np.random.default_rng(0).uniform(0.0, 100.0, size=(60, 2))

        edges = build_spatial_graph(coordinates, method="knn", n_neighbors=5).tocoo()

        distance = np.linalg.norm(coordinates[edges.row] - coordinates[edges.col], axis=1)
        np.testing.assert_allclose(edges.data, 1.0 / (1.0 + distance))

    def test_the_graph_is_symmetric_and_keeps_every_one_way_edge(self) -> None:
        coordinates = np.random.default_rng(1).uniform(0.0, 100.0, size=(40, 2))

        adjacency = build_spatial_graph(coordinates, method="knn", n_neighbors=4).toarray()

        np.testing.assert_array_equal(adjacency, adjacency.T)
        assert np.all((adjacency > 0).sum(axis=1) >= 4)  # each spot keeps at least its own 4 neighbors

    def test_close_mutual_neighbors_are_not_capped_at_one(self) -> None:
        # With distances below 1 the doubled weights used to be clipped to exactly 1.0
        coordinates = np.array([[0.0, 0.0], [0.2, 0.0], [5.0, 0.0]])

        adjacency = build_spatial_graph(coordinates, method="knn", n_neighbors=1).toarray()

        assert adjacency[0, 1] == pytest.approx(1.0 / 1.2)
        assert adjacency[0, 1] < 1.0


class TestSpatialClusterRadiusGraph:
    @staticmethod
    def _grid() -> tuple[np.ndarray, np.ndarray]:
        xs, ys = np.meshgrid(np.arange(8), np.arange(5))
        coordinates = np.column_stack([xs.ravel(), ys.ravel()]).astype(float)
        expression = np.random.default_rng(3).random((coordinates.shape[0], 10))
        return expression, coordinates

    @pytest.mark.parametrize("method", ["leiden", "louvain"])
    def test_a_radius_graph_can_be_used(self, method: str) -> None:
        expression, coordinates = self._grid()

        result = spatial_cluster(
            expression, coordinates, method=method, graph_method="radius", radius=1.5, seed=1  # type: ignore[arg-type]
        )

        assert len(result.labels) == 40
        assert result.metadata["graph_method"] == "radius"
        assert result.metadata["radius"] == 1.5

    def test_the_radius_reaches_the_spatial_graph(self) -> None:
        expression, coordinates = self._grid()

        result = spatial_cluster(
            expression, coordinates, method="louvain", graph_method="radius", radius=1.0, spatial_weight=1.0, seed=1
        )

        expected = build_spatial_graph(coordinates, method="radius", radius=1.0).toarray()
        np.testing.assert_allclose(result.spatial_graph.toarray(), expected)
        assert np.count_nonzero(expected) == 2 * (7 * 5 + 8 * 4)  # the 4-neighborhood of the 8 x 5 grid

    def test_a_larger_radius_connects_more_spots(self) -> None:
        expression, coordinates = self._grid()

        small = spatial_cluster(expression, coordinates, graph_method="radius", radius=1.0, spatial_weight=1.0)
        large = spatial_cluster(expression, coordinates, graph_method="radius", radius=2.0, spatial_weight=1.0)

        assert np.count_nonzero(large.spatial_graph.toarray()) > np.count_nonzero(small.spatial_graph.toarray())

    @pytest.mark.parametrize("method", ["leiden", "louvain"])
    def test_a_radius_graph_without_a_radius_is_refused_up_front(self, method: str) -> None:
        expression, coordinates = self._grid()

        with pytest.raises(ValueError, match="radius"):
            spatial_cluster(expression, coordinates, method=method, graph_method="radius")  # type: ignore[arg-type]

    def test_kmeans_does_not_use_the_graph_so_it_needs_no_radius(self) -> None:
        expression, coordinates = self._grid()

        result = spatial_cluster(expression, coordinates, n_clusters=3, method="kmeans", graph_method="radius")

        assert result.n_clusters == 3


class TestLeidenIterations:
    @pytest.mark.parametrize("requested, passed", [(-1, -1), (-5, -1), (1, 1), (2, 2), (7, 7)])
    def test_the_documented_counts_reach_leidenalg_unchanged(self, requested: int, passed: int) -> None:
        # -1 means "until convergence" and used to be turned into 2 iterations
        assert _leidenalg_iterations(requested) == passed

    def test_zero_iterations_are_rejected_not_rewritten(self) -> None:
        with pytest.raises(ValueError, match="n_iterations"):
            _leidenalg_iterations(0)

    def test_leiden_clustering_rejects_zero_iterations_whatever_is_installed(self) -> None:
        adjacency = build_spatial_graph(np.random.default_rng(2).uniform(0, 50, (30, 2)), method="knn", n_neighbors=4)

        with pytest.raises(ValueError, match="n_iterations"):
            leiden_clustering(adjacency, n_iterations=0)

    @pytest.mark.skipif(not LEIDENALG_AVAILABLE, reason="needs python-igraph and leidenalg")
    def test_running_to_convergence_is_at_least_as_good_as_two_iterations(self) -> None:
        adjacency = build_spatial_graph(np.random.default_rng(4).uniform(0, 100, (200, 2)), method="knn", n_neighbors=6)

        _, converged = leiden_clustering(adjacency, n_iterations=-1, seed=7)
        _, two = leiden_clustering(adjacency, n_iterations=2, seed=7)

        assert converged >= two - 1e-12
