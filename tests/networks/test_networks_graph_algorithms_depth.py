"""Depth tests for metainformant.networks.analysis.graph_algorithms (Round-4 T3).

Zero-mock, value-pinned: real networkx graphs with known topology, exact
centrality/similarity values, export/import roundtrips, and invariant checks.
"""

from __future__ import annotations

import json

import networkx as nx
import pytest

from metainformant.networks.analysis.graph_algorithms import (
    centrality_measures,
    export_network,
    filter_network,
    get_connected_components,
    import_network,
    network_intersection,
    network_similarity,
    network_union,
    shortest_paths,
)


def _star_graph() -> nx.Graph:
    """5-node star: center 0 connected to 1-4. Known exact centralities."""
    g = nx.Graph()
    g.add_edges_from([(0, 1), (0, 2), (0, 3), (0, 4)])
    return g


class TestCentralityMeasures:
    """Exact values on the star graph."""

    def test_star_degree_centrality_exact(self) -> None:
        result = centrality_measures(_star_graph())
        assert result["degree"][0] == 1.0  # connected to all 4 others
        assert all(result["degree"][i] == pytest.approx(0.25) for i in range(1, 5))

    def test_star_betweenness_exact(self) -> None:
        result = centrality_measures(_star_graph())
        # Center lies on all leaf-leaf shortest paths: normalized to 1.0
        assert result["betweenness"][0] == pytest.approx(1.0)
        assert all(result["betweenness"][i] == pytest.approx(0.0) for i in range(1, 5))

    def test_star_closeness_exact(self) -> None:
        result = centrality_measures(_star_graph())
        # Center: (n-1)/sum(d) = 4/4 = 1.0; leaves: 4/7
        assert result["closeness"][0] == pytest.approx(1.0)
        assert result["closeness"][1] == pytest.approx(4.0 / 7.0)

    def test_pagerank_sums_to_one(self) -> None:
        result = centrality_measures(_star_graph())
        assert sum(result["pagerank"].values()) == pytest.approx(1.0)

    def test_empty_graph_shape(self) -> None:
        assert centrality_measures(nx.Graph()) == {}

    def test_disconnected_eigenvector_falls_back(self) -> None:
        # Eigenvector centrality may fail to converge on disconnected graphs;
        # the API must still return all keys.
        g = nx.Graph()
        g.add_edges_from([(0, 1)])
        g.add_edge(2, 3)
        result = centrality_measures(g)
        assert set(result.keys()) == {"degree", "betweenness", "closeness", "eigenvector", "pagerank"}


class TestShortestPaths:
    def test_single_pair(self) -> None:
        g = _star_graph()
        assert shortest_paths(g, source=1, target=2) == {1: {2: 2}}

    def test_single_source_leaf_distances_exact(self) -> None:
        g = _star_graph()
        result = shortest_paths(g, source=1)
        assert result[1] == {1: 0, 0: 1, 2: 2, 3: 2, 4: 2}

    def test_source_equals_target(self) -> None:
        g = _star_graph()
        assert shortest_paths(g, source=0, target=0) == {0: {0: 0}}

    def test_unreachable_pair_returns_empty_dict(self) -> None:
        # Regression: NetworkXNoPath is not a NetworkXError subclass, so an
        # unreachable pair previously raised instead of returning the
        # documented empty dict.
        g = nx.Graph()
        g.add_edge(0, 1)
        g.add_edge(2, 3)
        assert shortest_paths(g, source=0, target=2) == {}
        # Reachable queries are unaffected by the fix.
        assert shortest_paths(g, source=0, target=1) == {0: {1: 1}}


class TestSimilarity:
    def test_summary_jaccards_exact(self) -> None:
        g1 = nx.Graph([(0, 1), (1, 2)])
        g2 = nx.Graph([(1, 2), (2, 3)])
        sim = network_similarity(g1, g2, method="summary")
        # Nodes {0,1,2} vs {1,2,3}: intersection 2, union 4 => 0.5
        assert sim["node_jaccard"] == pytest.approx(0.5)
        # Edges {(0,1),(1,2)} vs {(1,2),(2,3)}: intersection 1, union 3 => 1/3
        assert sim["edge_jaccard"] == pytest.approx(1.0 / 3.0)

    def test_identical_graphs_are_one(self) -> None:
        g = _star_graph()
        assert network_similarity(g, g, method="jaccard") == pytest.approx(1.0)
        assert network_similarity(g, g, method="dice") == pytest.approx(1.0)
        assert network_similarity(g, g, method="overlap") == pytest.approx(1.0)

    def test_dice_edge_symmetry(self) -> None:
        g1 = nx.Graph([(0, 1), (1, 2), (2, 3)])
        g2 = nx.Graph([(1, 2)])
        # 2*1/(3+1) = 0.5
        assert network_similarity(g1, g2, method="dice") == pytest.approx(0.5)

    def test_disjoint_graphs_zero(self) -> None:
        g1 = nx.Graph([(0, 1)])
        g2 = nx.Graph([(2, 3)])
        assert network_similarity(g1, g2, method="jaccard") == 0.0

    def test_unknown_method_raises(self) -> None:
        with pytest.raises(ValueError, match="Unsupported"):
            network_similarity(nx.Graph([(0, 1)]), nx.Graph([(0, 1)]), method="nope")


class TestUnionIntersection:
    def test_union_sums_duplicate_edge_weights(self) -> None:
        g1 = nx.Graph()
        g1.add_edge(0, 1, weight=2.0)
        g2 = nx.Graph()
        g2.add_edge(0, 1, weight=3.0)
        u = network_union(g1, g2)
        assert u.graph[0][1]["weight"] == 5.0

    def test_union_keeps_disjoint_nodes(self) -> None:
        g1 = nx.Graph([(0, 1)])
        g2 = nx.Graph([(2, 3)])
        u = network_union(g1, g2)
        assert set(u.graph.nodes()) == {0, 1, 2, 3}
        assert u.graph.number_of_edges() == 2

    def test_intersection_takes_common_edges_only(self) -> None:
        g1 = nx.Graph([(0, 1), (1, 2)])
        g2 = nx.Graph([(0, 1), (2, 3)])
        inter = network_intersection(g1, g2)
        assert set(inter.graph.nodes()) == {0, 1, 2}  # node 3 not in g1
        assert set(inter.graph.edges()) == {(0, 1)}

    def test_intersection_disjoint_returns_empty(self) -> None:
        g1 = nx.Graph([(0, 1)])
        g2 = nx.Graph([(2, 3)])
        inter = network_intersection(g1, g2)
        assert inter.graph.number_of_nodes() == 0

    def test_union_preserves_directedness(self) -> None:
        d1 = nx.DiGraph([(0, 1)])
        d2 = nx.DiGraph([(1, 2)])
        u = network_union(d1, d2)
        assert u.graph.is_directed()


class TestConnectedComponents:
    def test_two_components_sorted_content(self) -> None:
        g = nx.Graph()
        g.add_edges_from([(0, 1), (2, 3), (3, 4)])
        comps = get_connected_components(g)
        assert sorted(sorted(c) for c in comps) == [[0, 1], [2, 3, 4]]

    def test_directed_uses_weakly_connected(self) -> None:
        g = nx.DiGraph()
        g.add_edges_from([(0, 1), (2, 1)])
        comps = get_connected_components(g)
        assert sorted(comps, key=lambda c: sorted(c)) == [[0, 1, 2]]


class TestExportImportRoundtrip:
    """Export then import must reproduce topology and weights."""

    @pytest.mark.parametrize("fmt", ["json", "graphml", "gml", "edgelist", "csv"])
    def test_roundtrip_preserves_graph(self, tmp_path, fmt: str) -> None:
        g = nx.Graph()
        g.add_edges_from([(0, 1), (1, 2), (2, 3)], weight=2.5)
        path = tmp_path / f"net.{fmt}"
        export_network(g, path, format=fmt)
        restored = import_network(path, format=fmt)
        assert set(map(str, restored.graph.nodes())) == {"0", "1", "2", "3"}
        assert restored.graph.number_of_edges() == 3

    def test_json_roundtrip_preserves_weights(self, tmp_path) -> None:
        g = nx.Graph()
        g.add_edge("A", "B", weight=3.5)
        path = tmp_path / "w.json"
        export_network(g, path, format="json")
        data = json.loads(path.read_text())
        assert data["edges"][0]["weight"] == 3.5
        restored = import_network(path, format="json")
        assert restored.graph["A"]["B"]["weight"] == 3.5

    def test_unsupported_export_format_raises(self, tmp_path) -> None:
        with pytest.raises(ValueError, match="Unsupported"):
            export_network(nx.Graph([(0, 1)]), tmp_path / "x.xml", format="xml")

    def test_unsupported_import_format_raises_for_existing_file(self, tmp_path) -> None:
        # The format check happens after existence; give a real file so the
        # ValueError (not FileNotFoundError) is what surfaces.
        path = tmp_path / "x.xml"
        path.write_text("<net/>")
        with pytest.raises(ValueError, match="Unsupported"):
            import_network(path, format="xml")

    def test_missing_file_raises(self, tmp_path) -> None:
        with pytest.raises(FileNotFoundError):
            import_network(tmp_path / "absent.json", format="json")


class TestFilterNetwork:
    def test_min_degree_removes_leaves(self) -> None:
        g = _star_graph()
        filtered = filter_network(g, min_degree=4)
        assert set(filtered.graph.nodes()) == {0}

    def test_max_degree_keeps_leaves_only(self) -> None:
        g = _star_graph()
        filtered = filter_network(g, min_degree=1, max_degree=1)
        assert set(filtered.graph.nodes()) == {1, 2, 3, 4}

    def test_weight_filter_drops_light_edges(self) -> None:
        g = nx.Graph()
        g.add_edge(0, 1, weight=0.1)
        g.add_edge(1, 2, weight=0.9)
        filtered = filter_network(g, min_weight=0.5)
        assert set(filtered.graph.edges()) == {(1, 2)}

    def test_min_edge_weight_alias(self) -> None:
        g = nx.Graph()
        g.add_edge(0, 1, weight=0.2)
        g.add_edge(1, 2, weight=0.8)
        filtered = filter_network(g, min_edge_weight=0.5)
        assert set(filtered.graph.edges()) == {(1, 2)}

    def test_filter_is_not_in_place(self) -> None:
        g = _star_graph()
        filter_network(g, min_degree=4)
        assert g.number_of_nodes() == 5  # original untouched
