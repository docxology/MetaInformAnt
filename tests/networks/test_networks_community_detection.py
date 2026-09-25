"""Direct tests for the community detection wrappers in networks.analysis.community.

Covers the public per-algorithm functions that the ``detect_communities``
dispatcher delegates to:

- ``louvain_communities`` and ``leiden_communities`` (optional dependencies:
  tests skip cleanly when ``python-louvain`` / ``leidenalg`` are absent and
  exercise the real dependency guard otherwise),
- ``girvan_newman_communities`` (deterministic NetworkX algorithm),
- ``asyn_lpa_communities`` (seeded, deterministic),
- ``fluid_communities`` (seeded, deterministic).

All tests use real NetworkX graphs; no test doubles.
"""

from __future__ import annotations

import networkx as nx
import pytest

from metainformant.networks.analysis.community import (
    asyn_lpa_communities,
    fluid_communities,
    girvan_newman_communities,
    leiden_communities,
    louvain_communities,
)

try:
    import community as community_louvain  # type: ignore[import-not-found]  # noqa: F401

    HAS_LOUVAIN = True
except ImportError:
    HAS_LOUVAIN = False

try:
    import leidenalg  # type: ignore[import-not-found]  # noqa: F401

    HAS_LEIDEN = True
except ImportError:
    HAS_LEIDEN = False


def _two_cluster_graph() -> nx.Graph:
    """Build two 5-cliques joined by a single bridge edge (deterministic)."""
    g: nx.Graph = nx.Graph()
    cluster_a, cluster_b = list(range(5)), list(range(5, 10))
    g.add_edges_from((i, j) for i in cluster_a for j in cluster_a if i < j)
    g.add_edges_from((i, j) for i in cluster_b for j in cluster_b if i < j)
    g.add_edge(4, 5)
    return g


def _assert_two_cluster_split(community_lists: list[list[str]]) -> None:
    """Assert the partition is exactly the two dense clusters."""
    assert len(community_lists) == 2
    assert sorted(map(sorted, community_lists)) == [[0, 1, 2, 3, 4], [5, 6, 7, 8, 9]]


class TestGirvanNewmanCommunities:
    """Girvan-Newman is fully deterministic (edge-betweenness peeling)."""

    def test_two_cluster_split(self) -> None:
        communities = girvan_newman_communities(_two_cluster_graph(), n_communities=2)
        _assert_two_cluster_split(communities)

    def test_edgeless_graph_falls_back_to_singletons(self) -> None:
        graph = nx.Graph()
        graph.add_nodes_from(["x", "y", "z"])
        assert girvan_newman_communities(graph) == [["x"], ["y"], ["z"]]

    def test_empty_graph_returns_empty_partition(self) -> None:
        assert girvan_newman_communities(nx.Graph()) == []


class TestAsynLpaCommunities:
    """Asynchronous label propagation with a fixed seed is deterministic."""

    def test_two_cluster_split_with_seed(self) -> None:
        communities = asyn_lpa_communities(_two_cluster_graph(), seed=42)
        _assert_two_cluster_split(communities)

    def test_edgeless_graph_falls_back_to_singletons(self) -> None:
        graph = nx.Graph()
        graph.add_nodes_from(["x", "y", "z"])
        assert asyn_lpa_communities(graph) == [["x"], ["y"], ["z"]]

    def test_empty_graph_returns_empty_partition(self) -> None:
        assert asyn_lpa_communities(nx.Graph()) == []


class TestFluidCommunities:
    """Fluid communities with a fixed seed is deterministic."""

    def test_two_cluster_split_with_seed(self) -> None:
        communities = fluid_communities(_two_cluster_graph(), k=2, seed=42)
        _assert_two_cluster_split(communities)

    def test_empty_graph_returns_empty_partition(self) -> None:
        assert fluid_communities(nx.Graph(), k=2) == []

    def test_disconnected_graph_raises_clean_networkx_error(self) -> None:
        """asyn_fluidc requires a connected graph; the wrapper must surface it."""
        graph = nx.Graph()
        graph.add_nodes_from(["x", "y", "z"])
        with pytest.raises(nx.NetworkXError, match="connected"):
            fluid_communities(graph, k=2)


class TestLouvainCommunities:
    """Louvain requires the optional ``python-louvain`` package."""

    def test_dependency_guard_is_clean_when_missing(self) -> None:
        if HAS_LOUVAIN:
            pytest.skip("python-louvain installed; dependency guard path not reachable")
        with pytest.raises(ImportError, match="python-louvain"):
            louvain_communities(_two_cluster_graph())

    def test_empty_graph_returns_empty_partition(self) -> None:
        if not HAS_LOUVAIN:
            pytest.skip("python-louvain not available")
        assert louvain_communities(nx.Graph()) == []

    def test_edgeless_graph_falls_back_to_singletons(self) -> None:
        if not HAS_LOUVAIN:
            pytest.skip("python-louvain not available")
        graph = nx.Graph()
        graph.add_nodes_from(["x", "y", "z"])
        assert louvain_communities(graph) == [["x"], ["y"], ["z"]]

    def test_two_cluster_split_with_random_state(self) -> None:
        if not HAS_LOUVAIN:
            pytest.skip("python-louvain not available")
        communities = louvain_communities(_two_cluster_graph(), random_state=42)
        _assert_two_cluster_split(communities)


class TestLeidenCommunities:
    """Leiden requires the optional ``leidenalg``/``igraph`` packages."""

    def test_empty_graph_returns_empty_partition(self) -> None:
        """The trivial fallback runs before the optional dependency import."""
        assert leiden_communities(nx.Graph()) == []

    def test_edgeless_graph_falls_back_to_singletons(self) -> None:
        graph = nx.Graph()
        graph.add_nodes_from(["x", "y", "z"])
        assert leiden_communities(graph) == [["x"], ["y"], ["z"]]

    def test_dependency_guard_is_clean_when_missing(self) -> None:
        if HAS_LEIDEN:
            pytest.skip("leidenalg installed; dependency guard path not reachable")
        with pytest.raises(ImportError, match="leidenalg"):
            leiden_communities(_two_cluster_graph())

    def test_two_cluster_split_with_random_state(self) -> None:
        if not HAS_LEIDEN:
            pytest.skip("leidenalg/igraph not available")
        communities = leiden_communities(_two_cluster_graph(), random_state=42)
        _assert_two_cluster_split(communities)
