"""Comprehensive tests for the networks module.

Tests cover BiologicalNetwork, graph creation and manipulation, metrics,
community detection, pathway analysis, protein-protein interactions,
gene regulatory networks, configuration, workflow orchestration, and edge cases.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from metainformant.networks.analysis.graph import (
    BiologicalNetwork,
    add_edges_from_correlation,
    add_edges_from_dataframe,
    add_edges_from_interactions,
    add_nodes_from_dataframe,
    centrality_measures,
    create_network,
    export_network,
    extract_subgraph,
    filter_network,
    get_connected_components,
    get_network_summary,
    import_network,
    load_network,
    network_intersection,
    network_metrics,
    network_similarity,
    network_union,
    save_network,
    shortest_paths,
    validate_network,
)
from metainformant.networks.analysis.community import (
    community_metrics,
    detect_communities,
    evaluate_communities,
    greedy_modularity_communities,
    label_propagation_communities,
    modularity,
)
from metainformant.networks.analysis.pathway import (
    PathwayNetwork,
    load_pathway_database,
    pathway_enrichment,
    pathway_enrichment_analysis,
)
from metainformant.networks.interaction.ppi import (
    ProteinNetwork,
    ppi_network_analysis,
    predict_interactions,
)
from metainformant.networks.interaction.regulatory import (
    GeneRegulatoryNetwork,
    infer_grn,
    regulatory_motifs,
)
from metainformant.networks.config import (
    CommunityDetectionConfig,
    GRNConfig,
    NetworkConfig,
    NetworkWorkflowConfig,
    PPIConfig,
    PathwayEnrichmentConfig,
)
from metainformant.networks.workflow import NetworkWorkflow

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_triangle_network() -> BiologicalNetwork:
    """Build a simple undirected triangle A-B-C."""
    net = BiologicalNetwork(directed=False)
    net.add_edge("A", "B", weight=1.0)
    net.add_edge("B", "C", weight=0.8)
    net.add_edge("A", "C", weight=0.6)
    return net


def _make_two_cliques_graph():
    """Create an nx.Graph with two dense cliques joined by a single bridge edge."""
    edges = [
        ("A", "B"),
        ("A", "C"),
        ("B", "C"),
        ("D", "E"),
        ("D", "F"),
        ("E", "F"),
        ("C", "D"),  # bridge
    ]
    return create_network(edges, directed=False)


# ---------------------------------------------------------------------------
# 1. TestBiologicalNetwork
# ---------------------------------------------------------------------------


class TestBiologicalNetwork:
    """Tests for the BiologicalNetwork class."""

    def test_init_undirected(self) -> None:
        net = BiologicalNetwork(directed=False)
        assert net.number_of_nodes() == 0
        assert net.number_of_edges() == 0
        assert net.is_directed() is False

    def test_init_directed(self) -> None:
        net = BiologicalNetwork(directed=True)
        assert net.is_directed() is True

    def test_add_node(self) -> None:
        net = BiologicalNetwork()
        net.add_node("X", color="red")
        assert "X" in net
        assert net.number_of_nodes() == 1
        assert net.graph.nodes["X"]["color"] == "red"

    def test_add_edge_auto_creates_nodes(self) -> None:
        net = BiologicalNetwork()
        net.add_edge("A", "B", weight=0.5)
        assert net.number_of_nodes() == 2
        assert net.number_of_edges() == 1

    def test_remove_node(self) -> None:
        net = _make_triangle_network()
        net.remove_node("C")
        assert "C" not in net
        assert net.number_of_nodes() == 2
        # Removing C also removes edges A-C and B-C
        assert net.number_of_edges() == 1

    def test_remove_edge(self) -> None:
        net = _make_triangle_network()
        net.remove_edge("A", "B")
        assert net.number_of_edges() == 2

    def test_get_nodes(self) -> None:
        net = _make_triangle_network()
        nodes = net.get_nodes()
        assert set(nodes) == {"A", "B", "C"}

    def test_get_edges(self) -> None:
        net = _make_triangle_network()
        edges = net.get_edges()
        assert len(edges) == 3

    def test_size_returns_edge_count(self) -> None:
        net = _make_triangle_network()
        assert net.size() == 3

    def test_size_with_weight(self) -> None:
        net = _make_triangle_network()
        total_weight = net.size(weight="weight")
        assert total_weight == pytest.approx(1.0 + 0.8 + 0.6, rel=1e-6)

    def test_nodes_data(self) -> None:
        net = BiologicalNetwork()
        net.add_node("G1", gene_type="TF")
        data_list = list(net.nodes(data=True))
        assert data_list[0] == ("G1", {"gene_type": "TF"})

    def test_edges_data(self) -> None:
        net = BiologicalNetwork()
        net.add_edge("A", "B", weight=0.9)
        data_list = list(net.edges(data=True))
        assert data_list[0][2]["weight"] == pytest.approx(0.9)

    def test_neighbors(self) -> None:
        net = _make_triangle_network()
        neighbors = set(net.neighbors("A"))
        assert neighbors == {"B", "C"}

    def test_degree_single_node(self) -> None:
        net = _make_triangle_network()
        assert net.degree("A") == 2

    def test_degree_all_nodes(self) -> None:
        net = _make_triangle_network()
        deg = dict(net.degree())
        assert deg["A"] == 2
        assert deg["B"] == 2
        assert deg["C"] == 2

    def test_len(self) -> None:
        net = _make_triangle_network()
        assert len(net) == 3

    def test_contains(self) -> None:
        net = _make_triangle_network()
        assert "A" in net
        assert "Z" not in net

    def test_iter(self) -> None:
        net = _make_triangle_network()
        assert set(net) == {"A", "B", "C"}

    def test_directed_edge_semantics(self) -> None:
        net = BiologicalNetwork(directed=True)
        net.add_edge("A", "B")
        assert net.graph.has_edge("A", "B")
        assert not net.graph.has_edge("B", "A")


# ---------------------------------------------------------------------------
# 2. TestGraphCreation
# ---------------------------------------------------------------------------


class TestGraphCreation:
    """Tests for create_network and edge-adding utilities."""

    def test_create_network_from_edges(self) -> None:
        G = create_network([("A", "B"), ("B", "C")], directed=False)
        assert G.number_of_nodes() == 3
        assert G.number_of_edges() == 2

    def test_create_directed_network(self) -> None:
        G = create_network([("A", "B")], directed=True)
        assert G.is_directed() is True

    def test_create_network_empty(self) -> None:
        G = create_network([], directed=False)
        assert G.number_of_nodes() == 0

    def test_add_edges_from_correlation_dataframe(self) -> None:
        import networkx as nx

        labels = ["G1", "G2", "G3", "G4"]
        corr = pd.DataFrame(
            [
                [1.0, 0.9, 0.2, 0.1],
                [0.9, 1.0, 0.3, 0.05],
                [0.2, 0.3, 1.0, 0.85],
                [0.1, 0.05, 0.85, 1.0],
            ],
            index=labels,
            columns=labels,
        )
        G = nx.Graph()
        add_edges_from_correlation(G, corr, threshold=0.7)
        # Edges above 0.7: G1-G2 (0.9), G3-G4 (0.85)
        assert G.number_of_edges() == 2
        assert G.has_edge("G1", "G2")
        assert G.has_edge("G3", "G4")

    def test_add_edges_from_correlation_numpy(self) -> None:
        import networkx as nx

        corr_arr = np.array([[1.0, 0.8, 0.1], [0.8, 1.0, 0.3], [0.1, 0.3, 1.0]])
        G = nx.Graph()
        add_edges_from_correlation(G, corr_arr, threshold=0.5)
        assert G.number_of_edges() == 1  # Only Node_0-Node_1 at 0.8

    def test_add_edges_from_interactions(self) -> None:
        import networkx as nx

        G = nx.Graph()
        interactions = [("P1", "P2", 0.9), ("P2", "P3", 0.4), ("P3", "P4", 0.1)]
        add_edges_from_interactions(G, interactions, weight_threshold=0.3)
        # P1-P2 (0.9) and P2-P3 (0.4) pass; P3-P4 (0.1) does not
        assert G.has_edge("P1", "P2")
        assert G.has_edge("P2", "P3")
        assert not G.has_edge("P3", "P4")

    def test_add_edges_from_interactions_zero_threshold(self) -> None:
        import networkx as nx

        G = nx.Graph()
        interactions = [("X", "Y", 0.01)]
        add_edges_from_interactions(G, interactions, weight_threshold=0.0)
        assert G.has_edge("X", "Y")

    def test_add_nodes_from_dataframe(self) -> None:
        import networkx as nx

        df = pd.DataFrame({"gene": ["G1", "G2", "G3"], "score": [1.0, 2.0, 3.0]})
        G = nx.Graph()
        add_nodes_from_dataframe(G, df, node_column="gene", attribute_columns=["score"])
        assert set(G.nodes()) == {"G1", "G2", "G3"}
        assert G.nodes["G1"]["score"] == 1.0

    def test_add_edges_from_dataframe(self) -> None:
        import networkx as nx

        df = pd.DataFrame({"src": ["A", "B"], "tgt": ["B", "C"], "w": [0.5, 0.9]})
        G = nx.Graph()
        add_edges_from_dataframe(G, df, source_column="src", target_column="tgt", attribute_columns=["w"])
        assert G.number_of_edges() == 2
        assert G.edges["A", "B"]["w"] == 0.5


# ---------------------------------------------------------------------------
# 3. TestNetworkMetrics
# ---------------------------------------------------------------------------


class TestNetworkMetrics:
    """Tests for network_metrics, centrality, shortest_paths, summary, validate."""

    def _star_graph(self):
        """Hub A connected to B, C, D, E."""
        return create_network([("A", "B"), ("A", "C"), ("A", "D"), ("A", "E")], directed=False)

    def test_network_metrics_basic(self) -> None:
        G = self._star_graph()
        m = network_metrics(G)
        assert m["num_nodes"] == 5
        assert m["num_edges"] == 4
        assert m["avg_degree"] == pytest.approx(1.6)
        assert m["max_degree"] == 4

    def test_network_metrics_density(self) -> None:
        G = create_network([("A", "B"), ("B", "C"), ("A", "C")], directed=False)
        m = network_metrics(G)
        # 3 nodes, 3 edges, undirected max = 3*(3-1)/2 = 3, density = 3/3 = 1.0
        assert m["density"] == pytest.approx(1.0)

    def test_centrality_measures_hub(self) -> None:
        G = self._star_graph()
        c = centrality_measures(G)
        assert "degree" in c
        assert "betweenness" in c
        assert "closeness" in c
        # A should have highest degree centrality
        assert c["degree"]["A"] > c["degree"]["B"]

    def test_centrality_empty_graph(self) -> None:
        import networkx as nx

        G = nx.Graph()
        c = centrality_measures(G)
        assert c == {}

    def test_shortest_paths_all_pairs(self) -> None:
        G = create_network([("A", "B"), ("B", "C"), ("C", "D")], directed=False)
        sp = shortest_paths(G)
        assert sp["A"]["D"] == 3
        assert sp["A"]["A"] == 0

    def test_shortest_paths_single_pair(self) -> None:
        G = create_network([("A", "B"), ("B", "C")], directed=False)
        sp = shortest_paths(G, source="A", target="C")
        assert sp["A"]["C"] == 2

    def test_shortest_paths_single_source(self) -> None:
        G = create_network([("A", "B"), ("B", "C")], directed=False)
        sp = shortest_paths(G, source="A")
        assert "A" in sp
        assert sp["A"]["B"] == 1

    def test_get_network_summary_undirected(self) -> None:
        G = create_network([("A", "B"), ("B", "C"), ("A", "C")], directed=False)
        s = get_network_summary(G)
        assert s["n_nodes"] == 3
        assert s["n_edges"] == 3
        assert s["directed"] is False
        assert "degree_stats" in s
        assert s["connected_components"]["n_components"] == 1

    def test_validate_network_clean(self) -> None:
        G = create_network([("A", "B"), ("B", "C")], directed=False)
        is_valid, errors = validate_network(G)
        assert is_valid is True
        assert errors == []

    def test_validate_network_with_isolated(self) -> None:
        import networkx as nx

        G = nx.Graph()
        G.add_nodes_from(["A", "B", "C"])
        G.add_edge("A", "B")
        is_valid, errors = validate_network(G)
        assert is_valid is False
        assert any("isolated" in e.lower() for e in errors)

    def test_validate_network_with_self_loop(self) -> None:
        import networkx as nx

        G = nx.Graph()
        G.add_edge("A", "A")
        is_valid, errors = validate_network(G)
        assert is_valid is False
        assert any("self-loop" in e.lower() for e in errors)


# ---------------------------------------------------------------------------
# 4. TestNetworkIO
# ---------------------------------------------------------------------------


class TestNetworkIO:
    """Tests for export_network / import_network and save_network / load_network."""

    def test_export_import_json_roundtrip(self, tmp_path: Path) -> None:
        net = _make_triangle_network()
        filepath = tmp_path / "net.json"
        export_network(net, str(filepath), format="json")
        assert filepath.exists()

        loaded = import_network(str(filepath), format="json")
        assert isinstance(loaded, BiologicalNetwork)
        assert loaded.number_of_nodes() == 3
        assert loaded.number_of_edges() == 3

    def test_export_import_preserves_weight(self, tmp_path: Path) -> None:
        net = BiologicalNetwork()
        net.add_edge("X", "Y", weight=0.42)
        filepath = tmp_path / "weighted.json"
        export_network(net, str(filepath), format="json")
        loaded = import_network(str(filepath), format="json")
        # Check edge weight is preserved
        edge_data = loaded.graph.edges["X", "Y"]
        assert edge_data["weight"] == pytest.approx(0.42)

    def test_save_load_json(self, tmp_path: Path) -> None:
        G = create_network([("A", "B"), ("B", "C")], directed=False)
        path = tmp_path / "saved.json"
        save_network(G, str(path), format="json")
        assert path.exists()
        loaded = load_network(str(path), format="json")
        assert loaded.number_of_nodes() == 3

    def test_save_load_edgelist(self, tmp_path: Path) -> None:
        G = create_network([("A", "B"), ("B", "C")], directed=False)
        path = tmp_path / "saved.edgelist"
        save_network(G, str(path), format="edgelist")
        loaded = load_network(str(path), format="edgelist")
        assert loaded.number_of_edges() == 2

    def test_import_nonexistent_file_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            import_network(str(tmp_path / "no_such_file.json"), format="json")

    def test_export_unsupported_format_raises(self, tmp_path: Path) -> None:
        net = BiologicalNetwork()
        with pytest.raises(ValueError):
            export_network(net, str(tmp_path / "bad.xyz"), format="xyz")

    def test_export_nx_graph_directly(self, tmp_path: Path) -> None:
        """export_network should accept a raw nx.Graph too."""
        G = create_network([("A", "B")], directed=False)
        filepath = tmp_path / "raw.json"
        export_network(G, str(filepath), format="json")
        loaded = import_network(str(filepath), format="json")
        assert loaded.number_of_edges() == 1


# ---------------------------------------------------------------------------
# 5. TestNetworkOperations
# ---------------------------------------------------------------------------


class TestNetworkOperations:
    """Tests for subgraph, filter, components, union, intersection, similarity."""

    def test_extract_subgraph(self) -> None:
        net = _make_triangle_network()
        sub = extract_subgraph(net, ["A", "B"])
        assert isinstance(sub, BiologicalNetwork)
        assert sub.number_of_nodes() == 2
        assert sub.number_of_edges() == 1

    def test_filter_network_min_degree(self) -> None:
        G = create_network([("A", "B"), ("A", "C"), ("A", "D"), ("B", "C")], directed=False)
        filtered = filter_network(G, min_degree=2)
        assert isinstance(filtered, BiologicalNetwork)
        # D has degree 1 so is removed
        assert "D" not in filtered

    def test_filter_network_max_degree(self) -> None:
        G = create_network([("hub", "A"), ("hub", "B"), ("hub", "C"), ("A", "B")], directed=False)
        filtered = filter_network(G, max_degree=2)
        # hub has degree 3 so is removed
        assert "hub" not in filtered

    def test_filter_network_min_weight(self) -> None:
        import networkx as nx

        G = nx.Graph()
        G.add_edge("A", "B", weight=0.9)
        G.add_edge("B", "C", weight=0.2)
        filtered = filter_network(G, min_weight=0.5)
        assert filtered.graph.has_edge("A", "B")
        assert not filtered.graph.has_edge("B", "C")

    def test_get_connected_components(self) -> None:
        import networkx as nx

        G = nx.Graph()
        G.add_edges_from([("A", "B"), ("C", "D")])
        comps = get_connected_components(G)
        assert len(comps) == 2

    def test_get_connected_components_biological_network(self) -> None:
        net = BiologicalNetwork()
        net.add_edge("A", "B")
        net.add_edge("C", "D")
        comps = get_connected_components(net)
        assert len(comps) == 2

    def test_network_union(self) -> None:
        n1 = BiologicalNetwork()
        n1.add_edge("A", "B")
        n2 = BiologicalNetwork()
        n2.add_edge("C", "D")
        union = network_union(n1, n2)
        assert union.number_of_nodes() == 4
        assert union.number_of_edges() == 2

    def test_network_union_overlapping(self) -> None:
        n1 = BiologicalNetwork()
        n1.add_edge("A", "B")
        n2 = BiologicalNetwork()
        n2.add_edge("A", "C")
        union = network_union(n1, n2)
        assert union.number_of_nodes() == 3
        assert union.number_of_edges() == 2

    def test_network_intersection_common(self) -> None:
        n1 = BiologicalNetwork()
        n1.add_edge("A", "B")
        n1.add_edge("B", "C")
        n2 = BiologicalNetwork()
        n2.add_edge("A", "B")
        n2.add_edge("B", "D")
        inter = network_intersection(n1, n2)
        assert "A" in inter
        assert "B" in inter
        assert inter.graph.has_edge("A", "B")
        # C and D are unique
        assert "C" not in inter
        assert "D" not in inter

    def test_network_intersection_disjoint(self) -> None:
        n1 = BiologicalNetwork()
        n1.add_edge("A", "B")
        n2 = BiologicalNetwork()
        n2.add_edge("C", "D")
        inter = network_intersection(n1, n2)
        assert inter.number_of_nodes() == 0

    def test_network_similarity_jaccard(self) -> None:
        n1 = BiologicalNetwork()
        n1.add_edge("A", "B")
        n1.add_edge("B", "C")
        n2 = BiologicalNetwork()
        n2.add_edge("A", "B")
        n2.add_edge("B", "D")
        sim = network_similarity(n1, n2, method="jaccard")
        # nodes: n1={A,B,C}, n2={A,B,D}, inter={A,B}, union={A,B,C,D}
        assert sim == pytest.approx(2 / 4)

    def test_network_similarity_overlap(self) -> None:
        n1 = BiologicalNetwork()
        n1.add_edge("A", "B")
        n2 = BiologicalNetwork()
        n2.add_edge("A", "B")
        n2.add_edge("B", "C")
        sim = network_similarity(n1, n2, method="overlap")
        # nodes: inter=2, min_size=2 => 1.0
        assert sim == pytest.approx(1.0)

    @pytest.mark.parametrize("method", ["jaccard", "dice", "overlap"])
    def test_network_similarity_methods_exist(self, method: str) -> None:
        n1 = BiologicalNetwork()
        n1.add_edge("A", "B")
        n2 = BiologicalNetwork()
        n2.add_edge("A", "B")
        sim = network_similarity(n1, n2, method=method)
        assert 0.0 <= sim <= 1.0

    def test_network_similarity_invalid_method_raises(self) -> None:
        n1 = BiologicalNetwork()
        n2 = BiologicalNetwork()
        with pytest.raises(ValueError):
            network_similarity(n1, n2, method="bogus")


# ---------------------------------------------------------------------------
# 6. TestCommunityDetection
# ---------------------------------------------------------------------------


class TestCommunityDetection:
    """Tests for community detection algorithms."""

    def _two_cliques_graph(self):
        return _make_two_cliques_graph()

    def test_greedy_modularity_communities(self) -> None:
        G = self._two_cliques_graph()
        comms = greedy_modularity_communities(G)
        assert isinstance(comms, list)
        assert len(comms) >= 2  # Should identify at least 2 clusters

    def test_label_propagation_communities(self) -> None:
        G = self._two_cliques_graph()
        comms = label_propagation_communities(G)
        assert isinstance(comms, list)
        assert len(comms) >= 1

    def test_detect_communities_greedy(self) -> None:
        G = self._two_cliques_graph()
        mapping = detect_communities(G, method="greedy")
        assert isinstance(mapping, dict)
        assert len(mapping) == 6  # all 6 nodes assigned

    def test_detect_communities_label_propagation(self) -> None:
        G = self._two_cliques_graph()
        mapping = detect_communities(G, method="label_propagation")
        assert isinstance(mapping, dict)
        assert set(mapping.keys()) == {"A", "B", "C", "D", "E", "F"}

    def test_detect_communities_invalid_method_raises(self) -> None:
        G = create_network([("A", "B")], directed=False)
        with pytest.raises(ValueError):
            detect_communities(G, method="nonexistent")

    def test_evaluate_communities(self) -> None:
        G = self._two_cliques_graph()
        comms = [["A", "B", "C"], ["D", "E", "F"]]
        evaluation = evaluate_communities(G, comms)
        assert "n_communities" in evaluation
        assert evaluation["n_communities"] == 2
        assert "modularity" in evaluation
        assert "coverage" in evaluation

    def test_modularity_two_clusters(self) -> None:
        G = self._two_cliques_graph()
        comms = [["A", "B", "C"], ["D", "E", "F"]]
        m = modularity(G, comms)
        # Good partition should have positive modularity
        assert m > 0.0

    def test_modularity_single_community(self) -> None:
        G = self._two_cliques_graph()
        comms = [["A", "B", "C", "D", "E", "F"]]
        m = modularity(G, comms)
        # Single community should have zero or near-zero modularity
        assert m <= 0.01

    def test_community_metrics(self) -> None:
        G = self._two_cliques_graph()
        comms = [["A", "B", "C"], ["D", "E", "F"]]
        m = community_metrics(G, comms)
        assert m["n_communities"] == 2
        assert "modularity" in m
        assert "coverage" in m
        assert m["coverage"] == pytest.approx(1.0)
        assert "community_sizes" in m
        assert m["community_sizes"]["mean"] == pytest.approx(3.0)


# ---------------------------------------------------------------------------
# 7. TestPathwayNetwork
# ---------------------------------------------------------------------------


class TestPathwayNetwork:
    """Tests for PathwayNetwork class, enrichment, and database loading."""

    def test_init_empty(self) -> None:
        pn = PathwayNetwork(name="test")
        assert pn.name == "test"
        assert len(pn) == 0

    def test_init_with_pathways(self) -> None:
        pn = PathwayNetwork(
            name="my_pw",
            pathways={"pw1": ["G1", "G2"], "pw2": ["G3"]},
        )
        assert len(pn) == 2

    def test_add_and_get_pathway(self) -> None:
        pn = PathwayNetwork()
        pn.add_pathway("apoptosis", ["TP53", "BAX", "BCL2"])
        assert "apoptosis" in pn
        assert pn.get_pathway("apoptosis") == ["TP53", "BAX", "BCL2"]

    def test_get_pathway_missing_returns_empty(self) -> None:
        pn = PathwayNetwork()
        assert pn.get_pathway("nope") == []

    def test_remove_pathway(self) -> None:
        pn = PathwayNetwork(pathways={"pw1": ["A"], "pw2": ["B"]})
        assert pn.remove_pathway("pw1") is True
        assert len(pn) == 1
        assert pn.remove_pathway("nonexistent") is False

    def test_find_pathways_containing_gene(self) -> None:
        pn = PathwayNetwork(pathways={"pw1": ["G1", "G2"], "pw2": ["G2", "G3"], "pw3": ["G4"]})
        found = pn.find_pathways_containing_gene("G2")
        assert set(found) == {"pw1", "pw2"}

    def test_filter_pathways_by_size(self) -> None:
        pn = PathwayNetwork(
            pathways={
                "small": ["G1"],
                "medium": ["G1", "G2", "G3"],
                "large": ["G1", "G2", "G3", "G4", "G5"],
            }
        )
        filtered = pn.filter_pathways_by_size(min_size=2, max_size=4)
        assert "small" not in filtered.pathways
        assert "medium" in filtered.pathways
        assert "large" not in filtered.pathways

    def test_getitem(self) -> None:
        pn = PathwayNetwork(pathways={"pw1": ["A", "B"]})
        assert pn["pw1"] == ["A", "B"]

    def test_pathway_overlap_matrix(self) -> None:
        pn = PathwayNetwork(pathways={"pw1": ["A", "B", "C"], "pw2": ["B", "C", "D"]})
        matrix = pn.pathway_overlap_matrix()
        assert matrix["pw1"]["pw1"] == 1.0
        assert matrix["pw1"]["pw2"] == matrix["pw2"]["pw1"]
        # Jaccard: |{B,C}| / |{A,B,C,D}| = 2/4 = 0.5
        assert matrix["pw1"]["pw2"] == pytest.approx(0.5)

    def test_add_gene_to_pathway(self) -> None:
        pn = PathwayNetwork(pathways={"pw1": ["A"]})
        pn.add_gene_to_pathway("pw1", "B")
        assert "B" in pn.get_pathway("pw1")
        # Adding duplicate should not create duplicate
        pn.add_gene_to_pathway("pw1", "B")
        assert pn.get_pathway("pw1").count("B") == 1

    def test_load_pathway_database_direct_format(self) -> None:
        data = {
            "pw1": {"genes": ["G1", "G2"]},
            "pw2": {"genes": ["G3", "G4", "G5"]},
        }
        pn = load_pathway_database(data, name="test_db")
        assert len(pn) == 2
        assert pn.get_pathway("pw1") == ["G1", "G2"]

    def test_load_pathway_database_nested_format(self) -> None:
        data = {
            "pathways": {
                "pathway_A": ["X", "Y"],
                "pathway_B": ["Y", "Z"],
            }
        }
        pn = load_pathway_database(data, name="nested")
        assert len(pn) == 2

    def test_load_from_database_json(self, tmp_path: Path) -> None:
        import json

        db_file = tmp_path / "pathways.json"
        payload = {
            "name": "test_db",
            "pathways": {"glycolysis": ["HK1", "PFK", "PK"]},
        }
        with open(db_file, "w") as f:
            json.dump(payload, f)

        pn = PathwayNetwork.load_from_database(str(db_file), format="json")
        assert "glycolysis" in pn
        assert pn.get_pathway("glycolysis") == ["HK1", "PFK", "PK"]

    def test_pathway_enrichment_basic(self) -> None:
        pn = PathwayNetwork(
            pathways={
                "response": ["G1", "G2", "G3", "G4"],
                "metabolism": ["G5", "G6", "G7"],
            }
        )
        background = ["G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9", "G10"]
        query = ["G1", "G2", "G3"]

        results = pathway_enrichment(
            gene_list=query,
            pathway_network=pn,
            background_genes=background,
            method="fisher",
            min_overlap=1,
        )
        assert "response" in results
        assert results["response"]["overlap_size"] >= 3

    def test_pathway_enrichment_analysis_function(self) -> None:
        genes = ["G1", "G2", "G3"]
        background = ["G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8"]
        pathways = {"pw1": ["G1", "G2", "G3", "G4"], "pw2": ["G5", "G6"]}
        results = pathway_enrichment_analysis(genes=genes, background_genes=background, pathways=pathways)
        assert isinstance(results, list)
        assert len(results) == 2
        # results sorted by p-value
        pw_names = [r["pathway"] for r in results]
        assert "pw1" in pw_names


# ---------------------------------------------------------------------------
# 8. TestProteinNetwork
# ---------------------------------------------------------------------------


class TestProteinNetwork:
    """Tests for ProteinNetwork class and related functions."""

    def test_init_empty(self) -> None:
        pn = ProteinNetwork(name="test_ppi")
        assert len(pn) == 0
        assert pn.get_proteins() == []
        assert pn.get_interactions() == []

    def test_from_interactions(self) -> None:
        interactions = [("P1", "P2"), ("P2", "P3"), ("P3", "P4")]
        pn = ProteinNetwork.from_interactions(interactions, name="small_ppi")
        assert len(pn) == 4
        assert pn.name == "small_ppi"
        assert len(pn.get_interactions()) == 3

    def test_add_interaction(self) -> None:
        pn = ProteinNetwork()
        pn.add_interaction("A", "B", weight=0.9)
        assert "A" in pn
        assert "B" in pn
        assert pn.graph.edges["A", "B"]["weight"] == 0.9

    def test_find_neighbors(self) -> None:
        pn = ProteinNetwork.from_interactions([("A", "B"), ("A", "C")])
        neighbors = pn.find_neighbors("A")
        assert set(neighbors) == {"B", "C"}

    def test_find_neighbors_missing_protein(self) -> None:
        pn = ProteinNetwork()
        assert pn.find_neighbors("Z") == []

    def test_shortest_path(self) -> None:
        pn = ProteinNetwork.from_interactions([("A", "B"), ("B", "C"), ("C", "D")])
        path = pn.shortest_path("A", "D")
        assert path == ["A", "B", "C", "D"]

    def test_shortest_path_no_path(self) -> None:
        pn = ProteinNetwork.from_interactions([("A", "B")])
        pn.graph.add_node("Z")  # isolated node
        assert pn.shortest_path("A", "Z") is None

    def test_degree_distribution(self) -> None:
        pn = ProteinNetwork.from_interactions([("A", "B"), ("A", "C")])
        dd = pn.degree_distribution()
        assert dd["A"] == 2
        assert dd["B"] == 1

    def test_connected_components(self) -> None:
        pn = ProteinNetwork.from_interactions([("A", "B"), ("C", "D")])
        comps = pn.connected_components()
        assert len(comps) == 2

    def test_network_summary(self) -> None:
        pn = ProteinNetwork.from_interactions([("P1", "P2"), ("P2", "P3"), ("P1", "P3")])
        summary = pn.network_summary()
        assert summary["n_proteins"] == 3
        assert summary["n_interactions"] == 3
        assert summary["n_components"] == 1

    def test_ppi_network_analysis(self) -> None:
        import networkx as nx

        G = nx.Graph()
        G.add_edges_from([("P1", "P2"), ("P2", "P3"), ("P3", "P1"), ("P3", "P4")])
        analysis = ppi_network_analysis(G)
        assert analysis["basic_stats"]["n_proteins"] == 4
        assert analysis["basic_stats"]["n_interactions"] == 4
        assert "degree_distribution" in analysis

    def test_predict_interactions_correlation(self) -> None:
        np.random.seed(42)
        features = {
            "P1": list(np.random.randn(20)),
            "P2": list(np.random.randn(20)),
            "P3": list(np.random.randn(20)),
        }
        # Make P1 and P2 correlated
        features["P2"] = [x + 0.1 * np.random.randn() for x in features["P1"]]

        predictions = predict_interactions(
            target_proteins=["P1"],
            features=features,
            method="correlation",
            threshold=0.5,
        )
        assert isinstance(predictions, dict)
        assert "P1" in predictions

    def test_predict_interactions_similarity(self) -> None:
        known = ProteinNetwork.from_interactions([("P1", "P2"), ("P1", "P3"), ("P2", "P3"), ("P3", "P4")])
        predictions = predict_interactions(
            target_proteins=["P1"],
            known_network=known,
            method="similarity",
            threshold=0.0,
        )
        assert isinstance(predictions, dict)

    def test_protein_metadata(self) -> None:
        pn = ProteinNetwork()
        pn.add_protein_metadata("P1", function="kinase", species="human")
        meta = pn.get_protein_metadata("P1")
        assert meta["function"] == "kinase"
        assert meta["species"] == "human"


# ---------------------------------------------------------------------------
# 9. TestGeneRegulatoryNetwork
# ---------------------------------------------------------------------------


class TestGeneRegulatoryNetwork:
    """Tests for GeneRegulatoryNetwork class, inference, and motifs."""

    def test_init_empty(self) -> None:
        grn = GeneRegulatoryNetwork(name="test_grn")
        assert len(grn) == 0
        assert grn.name == "test_grn"

    def test_add_regulation(self) -> None:
        grn = GeneRegulatoryNetwork()
        grn.add_regulation("TF1", "GENE1", weight=0.9)
        assert "TF1" in grn
        assert "GENE1" in grn
        assert len(grn) == 1

    def test_from_interactions(self) -> None:
        interactions = [("TF1", "G1", 0.9), ("TF1", "G2", 0.7), ("TF2", "G1", 0.8)]
        grn = GeneRegulatoryNetwork.from_interactions(interactions, name="inferred")
        assert len(grn) == 3
        assert grn.name == "inferred"

    def test_get_transcription_factors(self) -> None:
        grn = GeneRegulatoryNetwork()
        grn.add_regulation("TF1", "G1")
        grn.add_regulation("TF2", "G2")
        tfs = grn.get_transcription_factors()
        assert set(tfs) == {"TF1", "TF2"}

    def test_get_targets(self) -> None:
        grn = GeneRegulatoryNetwork()
        grn.add_regulation("TF1", "G1")
        grn.add_regulation("TF1", "G2")
        targets = grn.get_targets()
        assert set(targets) == {"G1", "G2"}

    def test_get_tf_targets(self) -> None:
        grn = GeneRegulatoryNetwork()
        grn.add_regulation("TF1", "G1")
        grn.add_regulation("TF1", "G2")
        grn.add_regulation("TF2", "G3")
        assert set(grn.get_tf_targets("TF1")) == {"G1", "G2"}

    def test_get_target_regulators(self) -> None:
        grn = GeneRegulatoryNetwork()
        grn.add_regulation("TF1", "G1")
        grn.add_regulation("TF2", "G1")
        regs = grn.get_target_regulators("G1")
        assert set(regs) == {"TF1", "TF2"}

    def test_find_common_targets(self) -> None:
        grn = GeneRegulatoryNetwork()
        grn.add_regulation("TF1", "G1")
        grn.add_regulation("TF1", "G2")
        grn.add_regulation("TF2", "G2")
        grn.add_regulation("TF2", "G3")
        common = grn.find_common_targets("TF1", "TF2")
        assert "G2" in common

    def test_find_shared_regulators(self) -> None:
        grn = GeneRegulatoryNetwork()
        grn.add_regulation("TF1", "G1")
        grn.add_regulation("TF1", "G2")
        grn.add_regulation("TF2", "G1")
        shared = grn.find_shared_regulators("G1", "G2")
        assert "TF1" in shared

    def test_regulatory_motifs_method(self) -> None:
        grn = GeneRegulatoryNetwork()
        # Feed-forward: TF1 -> TF2 -> G1, and TF1 -> G1
        grn.add_regulation("TF1", "TF2")
        grn.add_regulation("TF2", "G1")
        grn.add_regulation("TF1", "G1")
        # Also make TF2 a TF by ensuring it's in tf_targets
        motifs = grn.regulatory_motifs()
        assert "feed_forward" in motifs
        assert len(motifs["feed_forward"]) >= 1

    def test_network_summary(self) -> None:
        grn = GeneRegulatoryNetwork(name="summary_test")
        grn.add_regulation("TF1", "G1")
        grn.add_regulation("TF1", "G2")
        summary = grn.network_summary()
        assert summary["name"] == "summary_test"
        assert summary["n_tfs"] == 1
        assert summary["n_targets"] == 2
        assert summary["is_directed"] is True

    def test_infer_grn_correlation(self) -> None:
        np.random.seed(42)
        n_samples = 50
        n_genes = 10
        expression = np.random.randn(n_samples, n_genes)
        # Make gene 0 and gene 1 correlated
        expression[:, 1] = 0.9 * expression[:, 0] + 0.1 * np.random.randn(n_samples)

        gene_names = [f"G{i}" for i in range(n_genes)]
        tf_genes = ["G0", "G1"]

        grn = infer_grn(expression, gene_names, method="correlation", threshold=0.5, tf_genes=tf_genes)
        assert isinstance(grn, GeneRegulatoryNetwork)
        assert len(grn) > 0

    def test_regulatory_motifs_standalone(self) -> None:
        grn = GeneRegulatoryNetwork()
        # Build feed-forward: A -> B -> C, A -> C
        grn.add_regulation("A", "B")
        grn.add_regulation("B", "C")
        grn.add_regulation("A", "C")
        motifs_list = regulatory_motifs(grn)
        assert isinstance(motifs_list, list)
        ffl = [m for m in motifs_list if m["motif_type"] == "feed_forward_loop"]
        assert len(ffl) >= 1

    def test_auto_regulation_motif(self) -> None:
        grn = GeneRegulatoryNetwork()
        grn.add_regulation("TF1", "TF1")  # self-regulation
        motifs_list = regulatory_motifs(grn)
        auto = [m for m in motifs_list if m["motif_type"] == "auto_regulation"]
        assert len(auto) == 1


# ---------------------------------------------------------------------------
# 10. TestConfig
# ---------------------------------------------------------------------------


class TestConfig:
    """Tests for all configuration dataclasses."""

    def test_network_config_defaults(self) -> None:
        cfg = NetworkConfig()
        assert cfg.directed is False
        assert cfg.weighted is True
        assert cfg.min_edge_weight == 0.0
        assert cfg.max_nodes is None
        assert cfg.remove_isolates is False

    def test_community_detection_config_defaults(self) -> None:
        cfg = CommunityDetectionConfig()
        assert cfg.method == "greedy"
        assert cfg.resolution == 1.0
        assert cfg.n_communities is None
        assert cfg.compare_methods is False

    def test_pathway_enrichment_config_defaults(self) -> None:
        cfg = PathwayEnrichmentConfig()
        assert cfg.method == "fisher"
        assert cfg.correction == "bonferroni"
        assert cfg.min_overlap == 1
        assert cfg.max_p_value == 0.05

    def test_ppi_config_defaults(self) -> None:
        cfg = PPIConfig()
        assert cfg.confidence_threshold == 0.4
        assert cfg.prediction_method == "similarity"

    def test_grn_config_defaults(self) -> None:
        cfg = GRNConfig()
        assert cfg.method == "correlation"
        assert cfg.threshold == 0.5
        assert "feed_forward" in cfg.motif_types

    def test_workflow_config_defaults(self) -> None:
        cfg = NetworkWorkflowConfig()
        assert isinstance(cfg.network, NetworkConfig)
        assert isinstance(cfg.community, CommunityDetectionConfig)
        assert isinstance(cfg.pathway, PathwayEnrichmentConfig)
        assert isinstance(cfg.ppi, PPIConfig)
        assert isinstance(cfg.grn, GRNConfig)
        assert cfg.export_format == "json"
        assert cfg.verbose is False

    def test_workflow_config_from_dict(self) -> None:
        d = {
            "network": {"directed": True, "weighted": False},
            "community": {"method": "label_propagation", "resolution": 2.0},
            "pathway": {"method": "hypergeometric", "correction": "fdr"},
            "ppi": {"confidence_threshold": 0.8},
            "grn": {"method": "mutual_info", "threshold": 0.3},
            "output_dir": "/tmp/results",
            "export_format": "graphml",
            "verbose": True,
        }
        cfg = NetworkWorkflowConfig.from_dict(d)
        assert cfg.network.directed is True
        assert cfg.community.method == "label_propagation"
        assert cfg.pathway.method == "hypergeometric"
        assert cfg.ppi.confidence_threshold == 0.8
        assert cfg.grn.method == "mutual_info"
        assert cfg.output_dir == "/tmp/results"
        assert cfg.verbose is True

    def test_workflow_config_from_dict_empty(self) -> None:
        cfg = NetworkWorkflowConfig.from_dict({})
        assert cfg.network.directed is False
        assert cfg.community.method == "greedy"


# ---------------------------------------------------------------------------
# 11. TestNetworkWorkflow
# ---------------------------------------------------------------------------


class TestNetworkWorkflow:
    """Tests for NetworkWorkflow chaining and summary."""

    def test_build_network_from_edges(self) -> None:
        wf = NetworkWorkflow()
        wf.build_network(edges=[("A", "B"), ("B", "C"), ("A", "C")])
        assert "network" in wf.results
        net = wf.results["network"]
        assert net.number_of_nodes() == 3
        assert net.number_of_edges() == 3

    def test_build_network_from_correlation(self) -> None:
        corr = np.array([[1.0, 0.9, 0.1], [0.9, 1.0, 0.2], [0.1, 0.2, 1.0]])
        wf = NetworkWorkflow()
        wf.build_network(
            correlation_matrix=corr,
            node_names=["A", "B", "C"],
            threshold=0.5,
        )
        net = wf.results["network"]
        assert net.number_of_edges() >= 1  # At least A-B at 0.9

    def test_chaining_build_and_metrics(self) -> None:
        wf = NetworkWorkflow()
        wf.build_network(edges=[("A", "B"), ("B", "C"), ("A", "C")]).analyze_metrics()

        assert "metrics" in wf.results
        assert "centrality" in wf.results
        assert wf.results["metrics"]["num_nodes"] == 3

    def test_chaining_build_detect_communities(self) -> None:
        edges = [
            ("A", "B"),
            ("A", "C"),
            ("B", "C"),
            ("D", "E"),
            ("D", "F"),
            ("E", "F"),
            ("C", "D"),
        ]
        wf = NetworkWorkflow()
        wf.build_network(edges=edges).detect_communities()

        assert "communities" in wf.results
        assert "community_evaluation" in wf.results

    def test_detect_communities_before_build_raises(self) -> None:
        wf = NetworkWorkflow()
        with pytest.raises(RuntimeError):
            wf.detect_communities()

    def test_analyze_metrics_before_build_raises(self) -> None:
        wf = NetworkWorkflow()
        with pytest.raises(RuntimeError):
            wf.analyze_metrics()

    def test_summary(self) -> None:
        wf = NetworkWorkflow()
        wf.build_network(edges=[("A", "B"), ("B", "C")]).analyze_metrics()

        s = wf.summary()
        assert "build_network" in s["steps_completed"]
        assert "analyze_metrics" in s["steps_completed"]
        assert s["network"]["n_nodes"] == 3

    def test_summary_empty_workflow(self) -> None:
        wf = NetworkWorkflow()
        s = wf.summary()
        assert s["steps_completed"] == []
        assert s["n_steps"] == 0

    def test_export_results(self, tmp_path: Path) -> None:
        wf = NetworkWorkflow()
        wf.build_network(edges=[("A", "B"), ("B", "C"), ("A", "C")]).analyze_metrics()

        exported = wf.export_results(output_dir=str(tmp_path / "net_out"))
        assert "network" in exported
        assert Path(exported["network"]).exists()
        assert "metrics" in exported

    def test_workflow_with_config(self) -> None:
        cfg = NetworkWorkflowConfig.from_dict({"community": {"method": "greedy"}})
        wf = NetworkWorkflow(config=cfg)
        wf.build_network(edges=[("A", "B"), ("B", "C"), ("A", "C")]).detect_communities()
        assert "communities" in wf.results

    def test_pathway_enrichment_step(self) -> None:
        pn = PathwayNetwork(
            pathways={
                "pw1": ["G1", "G2", "G3"],
                "pw2": ["G4", "G5"],
            }
        )
        wf = NetworkWorkflow()
        wf.build_network(edges=[("G1", "G2")])
        wf.run_pathway_enrichment(
            gene_list=["G1", "G2"],
            pathway_network=pn,
        )
        assert "pathway_enrichment" in wf.results


# ---------------------------------------------------------------------------
# 12. TestEdgeCases
# ---------------------------------------------------------------------------


class TestEdgeCases:
    """Tests for empty graphs, single nodes, disconnected components."""

    def test_empty_biological_network(self) -> None:
        net = BiologicalNetwork()
        assert net.number_of_nodes() == 0
        assert net.number_of_edges() == 0
        assert net.get_nodes() == []
        assert net.get_edges() == []
        assert net.size() == 0
        assert len(net) == 0

    def test_single_node_network(self) -> None:
        net = BiologicalNetwork()
        net.add_node("solo")
        assert net.number_of_nodes() == 1
        assert net.number_of_edges() == 0
        assert net.degree("solo") == 0

    def test_single_node_metrics(self) -> None:
        import networkx as nx

        G = nx.Graph()
        G.add_node("X")
        m = network_metrics(G)
        assert m["num_nodes"] == 1
        assert m["num_edges"] == 0
        assert m["avg_degree"] == 0

    def test_single_node_centrality(self) -> None:
        import networkx as nx

        G = nx.Graph()
        G.add_node("X")
        c = centrality_measures(G)
        # networkx 3.2+ returns 1.0 for single-node degree centrality
        assert c["degree"]["X"] == pytest.approx(1.0)

    def test_disconnected_components_metrics(self) -> None:
        import networkx as nx

        G = nx.Graph()
        G.add_edges_from([("A", "B"), ("C", "D"), ("E", "F")])
        m = network_metrics(G)
        assert m["num_components"] == 3

    def test_disconnected_shortest_paths(self) -> None:
        import networkx as nx

        G = nx.Graph()
        G.add_edge("A", "B")
        G.add_edge("C", "D")
        sp = shortest_paths(G)
        # A and C are disconnected, so no path A->C in results
        assert "C" not in sp.get("A", {})

    def test_empty_graph_validate(self) -> None:
        import networkx as nx

        G = nx.Graph()
        is_valid, errors = validate_network(G)
        assert is_valid is True

    def test_empty_graph_summary(self) -> None:
        import networkx as nx

        G = nx.Graph()
        s = get_network_summary(G)
        assert s["n_nodes"] == 0
        assert s["n_edges"] == 0

    def test_empty_pathway_network_enrichment(self) -> None:
        pn = PathwayNetwork()
        results = pathway_enrichment(gene_list=["G1"], pathway_network=pn, background_genes=["G1"])
        assert results == {}

    def test_filter_network_removes_all_nodes(self) -> None:
        import networkx as nx

        G = nx.Graph()
        G.add_edge("A", "B")  # both degree 1
        filtered = filter_network(G, min_degree=5)
        assert filtered.number_of_nodes() == 0

    def test_extract_subgraph_empty_node_list(self) -> None:
        net = _make_triangle_network()
        sub = extract_subgraph(net, [])
        assert sub.number_of_nodes() == 0

    def test_network_union_with_empty(self) -> None:
        n1 = BiologicalNetwork()
        n1.add_edge("A", "B")
        n2 = BiologicalNetwork()
        union = network_union(n1, n2)
        assert union.number_of_nodes() == 2
        assert union.number_of_edges() == 1

    def test_network_similarity_empty_graphs(self) -> None:
        n1 = BiologicalNetwork()
        n2 = BiologicalNetwork()
        sim = network_similarity(n1, n2, method="jaccard")
        assert sim == 0.0

    def test_community_single_node(self) -> None:
        import networkx as nx

        G = nx.Graph()
        G.add_node("A")
        comms = greedy_modularity_communities(G)
        assert len(comms) == 1

    def test_protein_network_empty_summary(self) -> None:
        pn = ProteinNetwork()
        summary = pn.network_summary()
        assert summary["n_proteins"] == 0
        assert summary["n_interactions"] == 0

    def test_grn_empty(self) -> None:
        grn = GeneRegulatoryNetwork()
        assert grn.get_transcription_factors() == []
        assert grn.get_targets() == []
        motifs = grn.regulatory_motifs()
        for motif_type, instances in motifs.items():
            assert instances == []

    def test_self_loop_graph(self) -> None:
        import networkx as nx

        G = nx.Graph()
        G.add_edge("A", "A")
        is_valid, errors = validate_network(G)
        assert is_valid is False

    @pytest.mark.parametrize("n_nodes", [2, 5, 10])
    def test_complete_graph_density(self, n_nodes: int) -> None:
        import networkx as nx

        G = nx.complete_graph(n_nodes)
        m = network_metrics(G)
        # Complete graph density should be 1.0
        expected_density = 1.0
        assert m["density"] == pytest.approx(expected_density, rel=1e-3)

    def test_directed_network_components(self) -> None:
        import networkx as nx

        G = nx.DiGraph()
        G.add_edges_from([("A", "B"), ("C", "D")])
        comps = get_connected_components(G)
        assert len(comps) == 2

    def test_export_import_empty_network(self, tmp_path: Path) -> None:
        net = BiologicalNetwork()
        filepath = tmp_path / "empty.json"
        export_network(net, str(filepath), format="json")
        loaded = import_network(str(filepath), format="json")
        assert loaded.number_of_nodes() == 0
        assert loaded.number_of_edges() == 0
