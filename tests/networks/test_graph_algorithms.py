"""Contract tests for metainformant.networks.analysis.graph_algorithms."""

from __future__ import annotations

import networkx as nx
import pytest

from metainformant.networks.analysis.graph_algorithms import shortest_paths
from metainformant.networks.analysis.graph_core import BiologicalNetwork


class TestShortestPathsAllPairs:
    def test_path_graph_distances_hand_enumerated(self) -> None:
        g = nx.Graph()
        g.add_edges_from([("A", "B"), ("B", "C"), ("C", "D")])
        sp = shortest_paths(g)
        assert sp["A"] == {"A": 0, "B": 1, "C": 2, "D": 3}
        assert sp["D"] == {"D": 0, "C": 1, "B": 2, "A": 3}

    def test_weights_ignored_hop_counts_only(self) -> None:
        # The wrapper passes weight=None through to networkx: lengths are
        # hop counts, not weighted distances.
        g = nx.Graph()
        g.add_edge("A", "B", weight=4.0)
        g.add_edge("B", "C", weight=1.0)
        g.add_edge("A", "C", weight=2.0)
        sp = shortest_paths(g)
        assert sp["A"]["C"] == 1
        assert sp["A"]["B"] == 1

    def test_disconnected_pairs_omitted(self) -> None:
        g = nx.Graph()
        g.add_edge("A", "B")
        g.add_edge("C", "D")
        g.add_node("E")
        sp = shortest_paths(g)
        assert "C" not in sp["A"] and "E" not in sp["A"]
        assert sp["C"] == {"C": 0, "D": 1}
        assert sp["E"] == {"E": 0}

    def test_biological_network_uses_same_sparse_contract(self) -> None:
        net = BiologicalNetwork()
        net.add_edge(0, 1)
        net.add_node(2)  # isolated: unreachable from everything
        distances = shortest_paths(net)
        assert distances[0][1] == 1
        assert distances[1][0] == 1
        # Aligned contract: unreachable pairs are omitted rather than
        # inf-filled; the diagonal is the 0-length path to itself.
        assert 2 not in distances[0]
        assert 0 not in distances[2]
        assert distances[2] == {2: 0}

    def test_directed_graph_respects_direction(self) -> None:
        g = nx.DiGraph()
        g.add_edges_from([("A", "B"), ("B", "C")])
        sp = shortest_paths(g)
        assert sp["A"]["C"] == 2
        assert "A" not in sp["C"]  # C cannot reach A


class TestShortestPathsQueries:
    def test_single_source_reachable_only(self) -> None:
        g = nx.Graph()
        g.add_edges_from([("A", "B"), ("B", "C")])
        g.add_node("D")
        assert shortest_paths(g, source="A") == {"A": {"A": 0, "B": 1, "C": 2}}

    def test_single_pair(self) -> None:
        g = nx.Graph()
        g.add_edges_from([("A", "B"), ("B", "C")])
        assert shortest_paths(g, source="A", target="C") == {"A": {"C": 2}}
        assert shortest_paths(g, source="C", target="A") == {
            "C": {"A": 2}
        }  # undirected

    def test_source_equals_target(self) -> None:
        g = nx.Graph()
        g.add_edge("A", "B")
        assert shortest_paths(g, source="A", target="A") == {"A": {"A": 0}}

    def test_unreachable_pair_returns_empty_dict(self) -> None:
        g = nx.Graph()
        g.add_node("A")
        g.add_edge("B", "C")
        assert shortest_paths(g, source="A", target="C") == {}
        assert shortest_paths(g, source="A", target="A") == {"A": {"A": 0}}

    def test_empty_graph_returns_empty(self) -> None:
        assert shortest_paths(nx.Graph()) == {}

    def test_missing_nodes_raise_with_detail(self) -> None:
        g = nx.Graph()
        g.add_edge("A", "B")
        with pytest.raises(ValueError, match="Z"):
            shortest_paths(g, source="Z")
        with pytest.raises(ValueError, match="Q"):
            shortest_paths(g, source="A", target="Q")
