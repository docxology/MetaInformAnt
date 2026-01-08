"""Tests for biological network graph functionality.

Real implementation testing for network analysis methods.
No mocking used - all tests use real computational methods.
"""

from __future__ import annotations

import numpy as np
import pytest

from metainformant.networks.analysis.graph import (
    BiologicalNetwork,
    add_edges_from_correlation,
    add_edges_from_interactions,
    centrality_measures,
    create_network,
    network_metrics,
    shortest_paths,
)


class TestBiologicalNetwork:
    """Test basic BiologicalNetwork functionality."""

    def test_network_initialization(self):
        """Test network creation and basic properties."""
        # Undirected network
        network = BiologicalNetwork(directed=False)
        assert network.num_nodes() == 0
        assert network.num_edges() == 0
        assert not network.directed

        # Directed network
        network_dir = BiologicalNetwork(directed=True)
        assert network_dir.directed

    def test_add_nodes_and_edges(self):
        """Test adding nodes and edges to network."""
        network = BiologicalNetwork()

        # Add nodes
        network.add_node("gene1", expression=10.5)
        network.add_node("gene2", expression=8.2)

        assert network.num_nodes() == 2
        assert "gene1" in network.nodes
        assert "gene2" in network.nodes
        assert network.node_attrs["gene1"]["expression"] == 10.5

        # Add edge
        network.add_edge("gene1", "gene2", weight=0.8)
        assert network.num_edges() == 1
        assert network.get_edge_weight("gene1", "gene2") == 0.8
        assert network.get_edge_weight("gene2", "gene1") == 0.8  # Undirected

    def test_directed_vs_undirected_edges(self):
        """Test edge behavior in directed vs undirected networks."""
        # Undirected network
        undirected = BiologicalNetwork(directed=False)
        undirected.add_edge("A", "B", weight=1.0)

        assert undirected.get_edge_weight("A", "B") == 1.0
        assert undirected.get_edge_weight("B", "A") == 1.0  # Same edge

        # Directed network
        directed = BiologicalNetwork(directed=True)
        directed.add_edge("A", "B", weight=1.0)

        assert directed.get_edge_weight("A", "B") == 1.0
        assert directed.get_edge_weight("B", "A") is None  # Different edge

    def test_get_neighbors(self):
        """Test neighbor finding functionality."""
        network = BiologicalNetwork()
        network.add_edge("A", "B")
        network.add_edge("A", "C")
        network.add_edge("B", "D")

        neighbors_A = network.get_neighbors("A")
        assert set(neighbors_A) == {"B", "C"}

        neighbors_B = network.get_neighbors("B")
        assert set(neighbors_B) == {"A", "D"}

        neighbors_D = network.get_neighbors("D")
        assert neighbors_D == ["B"]

    def test_network_density(self):
        """Test network density calculation."""
        # Empty network
        network = BiologicalNetwork()
        assert network.density() == 0.0

        # Single node
        network.add_node("A")
        assert network.density() == 0.0

        # Complete graph with 3 nodes
        network.add_edge("A", "B")
        network.add_edge("B", "C")
        network.add_edge("A", "C")

        # 3 nodes can have maximum 3 edges in undirected graph
        # We have all 3 edges, so density = 1.0
        assert network.density() == 1.0


class TestNetworkCreation:
    """Test network creation utilities."""

    def test_create_network(self):
        """Test creating network from node list."""
        nodes = ["gene1", "gene2", "gene3"]
        network = create_network(nodes, directed=False)

        assert network.num_nodes() == 3
        assert network.num_edges() == 0
        assert set(network.nodes) == set(nodes)
        assert not network.directed

        # Directed network
        directed_net = create_network(nodes, directed=True)
        assert directed_net.directed

    def test_add_edges_from_correlation(self):
        """Test adding edges based on correlation matrix."""
        network = BiologicalNetwork()

        # Create correlation matrix
        correlation_matrix = np.array([[1.0, 0.8, 0.2], [0.8, 1.0, -0.6], [0.2, -0.6, 1.0]])
        node_names = ["gene1", "gene2", "gene3"]

        add_edges_from_correlation(network, correlation_matrix, node_names, threshold=0.5)

        # Should have edges for correlations > 0.5
        # gene1-gene2: 0.8 > 0.5 ✓
        # gene2-gene3: |-0.6| = 0.6 > 0.5 ✓
        # gene1-gene3: 0.2 < 0.5 ✗

        assert network.num_edges() == 2
        assert network.get_edge_weight("gene1", "gene2") == 0.8
        assert network.get_edge_weight("gene2", "gene3") == 0.6  # Absolute value
        assert network.get_edge_weight("gene1", "gene3") is None

    def test_add_edges_from_correlation_dimension_mismatch(self):
        """Test error handling for dimension mismatch."""
        network = BiologicalNetwork()
        correlation_matrix = np.array([[1.0, 0.5], [0.5, 1.0]])
        node_names = ["gene1", "gene2", "gene3"]  # Too many names

        with pytest.raises(ValueError, match="dimensions don't match"):
            add_edges_from_correlation(network, correlation_matrix, node_names)

    def test_add_edges_from_interactions(self):
        """Test adding edges from interaction list."""
        network = BiologicalNetwork()

        interactions = [("protein1", "protein2", 0.9), ("protein2", "protein3", 0.7), ("protein1", "protein3", 0.4)]

        add_edges_from_interactions(network, interactions)

        assert network.num_edges() == 3
        assert network.get_edge_weight("protein1", "protein2") == 0.9
        assert network.get_edge_weight("protein2", "protein3") == 0.7
        assert network.get_edge_weight("protein1", "protein3") == 0.4


class TestNetworkMetrics:
    """Test network metrics calculations."""

    def setup_method(self):
        """Set up test network."""
        self.network = BiologicalNetwork()
        # Create a small test network
        #   A --- B
        #   |     |
        #   C --- D
        edges = [("A", "B"), ("B", "D"), ("D", "C"), ("C", "A")]
        for node1, node2 in edges:
            self.network.add_edge(node1, node2)

    def test_basic_network_metrics(self):
        """Test basic network metric calculations."""
        metrics = network_metrics(self.network)

        assert metrics["num_nodes"] == 4
        assert metrics["num_edges"] == 4
        assert metrics["avg_degree"] == 2.0  # Each node has degree 2
        assert metrics["max_degree"] == 2
        assert metrics["min_degree"] == 2

        # Density: 4 edges out of 6 possible = 2/3
        expected_density = 4 / 6
        assert abs(metrics["density"] - expected_density) < 1e-6

    def test_metrics_empty_network(self):
        """Test metrics on empty network."""
        empty_network = BiologicalNetwork()
        metrics = network_metrics(empty_network)

        assert metrics["num_nodes"] == 0
        assert metrics["num_edges"] == 0
        assert metrics["avg_degree"] == 0.0
        assert metrics["density"] == 0.0

    def test_metrics_single_node(self):
        """Test metrics on single node network."""
        single_network = BiologicalNetwork()
        single_network.add_node("A")
        metrics = network_metrics(single_network)

        assert metrics["num_nodes"] == 1
        assert metrics["num_edges"] == 0
        assert metrics["avg_degree"] == 0.0
        assert metrics["density"] == 0.0


class TestCentralityMeasures:
    """Test centrality measure calculations."""

    def setup_method(self):
        """Set up test network with known centrality properties."""
        self.network = BiologicalNetwork()
        # Star network: B is central, connected to A, C, D
        edges = [("B", "A"), ("B", "C"), ("B", "D")]
        for node1, node2 in edges:
            self.network.add_edge(node1, node2)

    def test_degree_centrality(self):
        """Test degree centrality calculation."""
        centralities = centrality_measures(self.network)

        degree_cent = centralities["degree"]

        # B has degree 3 out of max 3 = 1.0
        assert degree_cent["B"] == 1.0

        # A, C, D each have degree 1 out of max 3 = 1/3
        for node in ["A", "C", "D"]:
            assert abs(degree_cent[node] - 1 / 3) < 1e-6

    def test_closeness_centrality(self):
        """Test closeness centrality calculation."""
        centralities = centrality_measures(self.network)

        closeness_cent = centralities["closeness"]

        # B is closest to all other nodes (distance 1)
        # A, C, D are farther from others (distance 1 to B, 2 to others)
        assert closeness_cent["B"] > closeness_cent["A"]
        assert closeness_cent["B"] > closeness_cent["C"]
        assert closeness_cent["B"] > closeness_cent["D"]

    def test_centrality_empty_network(self):
        """Test centrality measures on empty network."""
        empty_network = BiologicalNetwork()
        centralities = centrality_measures(empty_network)

        assert len(centralities["degree"]) == 0
        assert len(centralities["closeness"]) == 0
        assert len(centralities["betweenness"]) == 0


class TestShortestPaths:
    """Test shortest path calculations."""

    def setup_method(self):
        """Set up test network."""
        self.network = BiologicalNetwork()
        # Linear network: A - B - C - D
        edges = [("A", "B"), ("B", "C"), ("C", "D")]
        for node1, node2 in edges:
            self.network.add_edge(node1, node2)

    def test_shortest_paths_linear_network(self):
        """Test shortest paths in linear network."""
        distances = shortest_paths(self.network)

        # Check specific distances
        assert distances["A"]["A"] == 0.0
        assert distances["A"]["B"] == 1.0
        assert distances["A"]["C"] == 2.0
        assert distances["A"]["D"] == 3.0

        assert distances["B"]["C"] == 1.0
        assert distances["D"]["A"] == 3.0

    def test_shortest_paths_disconnected(self):
        """Test shortest paths with disconnected components."""
        network = BiologicalNetwork()
        # Two disconnected components: A-B and C-D
        network.add_edge("A", "B")
        network.add_edge("C", "D")

        distances = shortest_paths(network)

        # Within components
        assert distances["A"]["B"] == 1.0
        assert distances["C"]["D"] == 1.0

        # Between components (should be infinite)
        assert distances["A"]["C"] == float("inf")
        assert distances["B"]["D"] == float("inf")

    def test_shortest_paths_single_node(self):
        """Test shortest paths with single node."""
        network = BiologicalNetwork()
        network.add_node("A")

        distances = shortest_paths(network)
        assert distances["A"]["A"] == 0.0


class TestNetworkIntegration:
    """Integration tests for network functionality."""

    def test_full_network_analysis_workflow(self):
        """Test complete network analysis workflow."""
        # Create correlation-based gene network
        gene_names = ["BRCA1", "TP53", "MYC", "EGFR", "PTEN"]

        # Simulate expression correlation matrix
        np.random.seed(42)
        correlation_matrix = np.random.rand(5, 5)
        correlation_matrix = (correlation_matrix + correlation_matrix.T) / 2  # Symmetric
        np.fill_diagonal(correlation_matrix, 1.0)  # Perfect self-correlation

        # Create network
        network = create_network(gene_names)
        add_edges_from_correlation(network, correlation_matrix, gene_names, threshold=0.6)

        # Analyze network
        metrics = network_metrics(network)
        centralities = centrality_measures(network)
        distances = shortest_paths(network)

        # Verify results make sense
        assert metrics["num_nodes"] == 5
        assert metrics["num_edges"] >= 0
        assert metrics["density"] >= 0.0
        assert metrics["density"] <= 1.0

        # Check centrality measures exist for all nodes
        for measure in ["degree", "closeness", "betweenness", "eigenvector"]:
            assert len(centralities[measure]) == 5
            for node in gene_names:
                assert node in centralities[measure]
                # Degree centrality is normalized to [0, 1]
                # Other measures may exceed 1.0 depending on normalization
                if measure == "degree":
                    assert 0.0 <= centralities[measure][node] <= 1.0
                else:
                    assert centralities[measure][node] >= 0.0

        # Check distances
        assert len(distances) == 5
        for source in gene_names:
            assert len(distances[source]) == 5
            assert distances[source][source] == 0.0  # Self-distance is 0


class TestEdgeCases:
    """Test edge cases and error conditions."""

    def test_nonexistent_nodes(self):
        """Test operations on nonexistent nodes."""
        network = BiologicalNetwork()
        network.add_edge("A", "B")

        # Getting neighbors of nonexistent node
        neighbors = network.get_neighbors("NONEXISTENT")
        assert neighbors == []

        # Getting edge weight for nonexistent nodes
        weight = network.get_edge_weight("NONEXISTENT", "A")
        assert weight is None

    def test_zero_weight_edges(self):
        """Test handling of zero-weight edges."""
        network = BiologicalNetwork()
        network.add_edge("A", "B", weight=0.0)

        assert network.get_edge_weight("A", "B") == 0.0
        assert network.num_edges() == 1

        neighbors_A = network.get_neighbors("A")
        assert "B" in neighbors_A

    def test_large_weights(self):
        """Test handling of very large edge weights."""
        network = BiologicalNetwork()
        large_weight = 1e10
        network.add_edge("A", "B", weight=large_weight)

        assert network.get_edge_weight("A", "B") == large_weight

        # Should not break metrics calculations
        metrics = network_metrics(network)
        assert np.isfinite(metrics["density"])


class TestNetworkExportImport:
    """Test network export and import functionality."""

    def test_export_import_json(self, tmp_path):
        """Test JSON export and import."""
        from metainformant.networks.analysis.graph import export_network, import_network

        network = create_network(["A", "B", "C"], directed=False)
        network.add_edge("A", "B", weight=0.8)
        network.add_node("A", function="test")

        json_file = tmp_path / "test_network.json"
        export_network(network, str(json_file), format="json")

        loaded = import_network(str(json_file), format="json")
        assert loaded.num_nodes() == 3
        assert loaded.num_edges() == 1
        assert loaded.get_edge_weight("A", "B") == 0.8
        assert loaded.node_attrs["A"]["function"] == "test"

    def test_export_import_csv(self, tmp_path):
        """Test CSV export and import."""
        from metainformant.networks.analysis.graph import export_network, import_network

        network = create_network(["A", "B"], directed=False)
        network.add_edge("A", "B", weight=0.9)

        csv_file = tmp_path / "test_network.csv"
        export_network(network, str(csv_file), format="csv")

        loaded = import_network(str(csv_file), format="csv")
        assert loaded.num_nodes() == 2
        assert loaded.num_edges() == 1

    def test_export_invalid_format(self):
        """Test error handling for invalid export format."""
        from metainformant.networks.analysis.graph import export_network

        network = create_network(["A", "B"])
        with pytest.raises(ValueError, match="Unsupported export format"):
            export_network(network, "test.xyz", format="xyz")


class TestNetworkComparison:
    """Test network comparison utilities."""

    def test_network_similarity(self):
        """Test network similarity calculation."""
        from metainformant.networks.analysis.graph import network_similarity

        net1 = create_network(["A", "B", "C"], directed=False)
        net1.add_edge("A", "B", weight=0.8)

        net2 = create_network(["A", "B", "D"], directed=False)
        net2.add_edge("A", "B", weight=0.9)

        sim = network_similarity(net1, net2)
        assert "node_jaccard" in sim
        assert "edge_jaccard" in sim
        assert sim["edge_jaccard"] > 0.0

    def test_extract_subgraph(self):
        """Test subgraph extraction."""
        from metainformant.networks.analysis.graph import extract_subgraph

        network = create_network(["A", "B", "C", "D"], directed=False)
        network.add_edge("A", "B")
        network.add_edge("B", "C")
        network.add_edge("C", "D")

        subgraph = extract_subgraph(network, ["A", "B", "C"])
        assert subgraph.num_nodes() == 3
        assert subgraph.num_edges() == 2

    def test_filter_network(self):
        """Test network filtering."""
        from metainformant.networks.analysis.graph import filter_network

        network = create_network(["A", "B", "C", "D"], directed=False)
        network.add_edge("A", "B", weight=0.9)
        network.add_edge("C", "D", weight=0.3)

        filtered = filter_network(network, min_edge_weight=0.5)
        assert filtered.num_edges() == 1

    def test_get_connected_components(self):
        """Test connected component detection."""
        from metainformant.networks.analysis.graph import get_connected_components

        network = create_network(["A", "B", "C", "D"], directed=False)
        network.add_edge("A", "B")
        network.add_edge("C", "D")

        components = get_connected_components(network)
        assert len(components) == 2

    def test_network_union_intersection(self):
        """Test network union and intersection."""
        from metainformant.networks.analysis.graph import network_union, network_intersection

        net1 = create_network(["A", "B", "C"], directed=False)
        net1.add_edge("A", "B", weight=0.5)

        net2 = create_network(["A", "B", "D"], directed=False)
        net2.add_edge("A", "B", weight=0.3)

        union = network_union(net1, net2)
        assert union.num_nodes() == 4
        assert union.get_edge_weight("A", "B") == 0.8  # Sum

        intersection = network_intersection(net1, net2)
        assert intersection.num_nodes() == 2
        assert intersection.num_edges() == 1

    def test_remove_node_edge(self):
        """Test node and edge removal."""
        from metainformant.networks.analysis.graph import remove_node, remove_edge

        network = create_network(["A", "B", "C"], directed=False)
        network.add_edge("A", "B")
        network.add_edge("B", "C")

        remove_edge(network, "A", "B")
        assert network.num_edges() == 1

        remove_node(network, "B")
        assert network.num_nodes() == 2
        assert network.num_edges() == 0

    def test_duplicate_edges(self):
        """Test behavior with duplicate edge additions."""
        network = BiologicalNetwork()

        # Add same edge twice with different weights
        network.add_edge("A", "B", weight=1.0)
        network.add_edge("A", "B", weight=2.0)  # Should overwrite

        assert network.get_edge_weight("A", "B") == 2.0
        assert network.num_edges() == 1  # Still just one edge
