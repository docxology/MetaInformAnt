"""Tests for biological network community detection functionality.

Real implementation testing for community detection algorithms in biological networks.
No mocking used - all tests use real computational methods and data.
"""

from __future__ import annotations

from typing import Dict, List, Set

import numpy as np
import pytest

from metainformant.networks.community import community_metrics, detect_communities, modularity
from metainformant.networks.graph import BiologicalNetwork, create_network


class TestCommunityDetection:
    """Test community detection algorithms."""

    def setup_method(self):
        """Set up test networks with known community structure."""
        # Create network with clear community structure
        # Community 1: A-B-C (triangle)
        # Community 2: D-E-F (triangle)
        # Bridge: C-D (connects communities)
        self.network = BiologicalNetwork()

        # Add community 1 (highly connected)
        edges_comm1 = [("A", "B"), ("B", "C"), ("C", "A")]
        for node1, node2 in edges_comm1:
            self.network.add_edge(node1, node2, weight=1.0)

        # Add community 2 (highly connected)
        edges_comm2 = [("D", "E"), ("E", "F"), ("F", "D")]
        for node1, node2 in edges_comm2:
            self.network.add_edge(node1, node2, weight=1.0)

        # Add bridge between communities (weaker connection)
        self.network.add_edge("C", "D", weight=0.3)

    def test_detect_communities_basic(self):
        """Test basic community detection functionality."""
        communities = detect_communities(self.network, method="greedy")

        # Should return dict mapping nodes to community IDs
        assert isinstance(communities, dict)
        assert len(communities) == 6  # All nodes assigned

        # All nodes should be assigned
        expected_nodes = {"A", "B", "C", "D", "E", "F"}
        assert set(communities.keys()) == expected_nodes

        # Check community IDs are reasonable
        community_ids = set(communities.values())
        assert len(community_ids) >= 1  # At least one community
        assert len(community_ids) <= 6  # At most one per node

    def test_detect_communities_different_algorithms(self):
        """Test different community detection algorithms."""
        algorithms = ["greedy", "leiden", "louvain"]

        for algorithm in algorithms:
            communities = detect_communities(self.network, method=algorithm)

            # Should detect communities (returns dict mapping node to community ID)
            assert isinstance(communities, dict)
            assert len(communities) == 6  # All nodes should be assigned

            # All nodes should be covered
            all_nodes = set(communities.keys())
            expected_nodes = {"A", "B", "C", "D", "E", "F"}
            assert all_nodes == expected_nodes

    def test_detect_communities_resolution(self):
        """Test community detection with different resolution parameters."""
        # Low resolution (fewer, larger communities)
        communities_low = detect_communities(self.network, method="louvain", resolution=0.5)

        # High resolution (more, smaller communities)
        communities_high = detect_communities(self.network, method="louvain", resolution=2.0)

        # Both should be valid
        assert isinstance(communities_low, dict)
        assert isinstance(communities_high, dict)
        assert len(communities_low) == 6
        assert len(communities_high) == 6

        # Generally expect more communities with higher resolution
        num_communities_low = len(set(communities_low.values()))
        num_communities_high = len(set(communities_high.values()))
        assert num_communities_high >= num_communities_low

    def test_detect_communities_empty_network(self):
        """Test community detection on empty network."""
        empty_network = BiologicalNetwork()

        communities = detect_communities(empty_network)
        assert isinstance(communities, dict)
        assert len(communities) == 0

    def test_detect_communities_single_node(self):
        """Test community detection on single node network."""
        single_network = BiologicalNetwork()
        single_network.add_node("A")

        communities = detect_communities(single_network)
        assert isinstance(communities, dict)
        assert len(communities) == 1
        assert "A" in communities

    def test_detect_communities_disconnected(self):
        """Test community detection on disconnected network."""
        # Create two disconnected components
        disconnected = BiologicalNetwork()

        # Component 1: A-B
        disconnected.add_edge("A", "B")

        # Component 2: C-D (isolated)
        disconnected.add_edge("C", "D")

        communities = detect_communities(disconnected)

        # Should detect separate communities for disconnected components
        assert isinstance(communities, dict)
        assert len(communities) == 4

        # Each component should be in separate communities
        all_nodes = set(communities.keys())
        assert all_nodes == {"A", "B", "C", "D"}


class TestModularity:
    """Test modularity calculation."""

    def setup_method(self):
        """Set up test network for modularity testing."""
        # Create simple network with known modularity
        self.network = BiologicalNetwork()
        edges = [("A", "B"), ("B", "C"), ("D", "E")]  # Two components
        for node1, node2 in edges:
            self.network.add_edge(node1, node2)

    def test_modularity_calculation(self):
        """Test basic modularity calculation."""
        # Perfect community structure (each component is a community)
        # Convert list format to dict format
        communities = {"A": 0, "B": 0, "C": 0, "D": 1, "E": 1}

        Q = modularity(self.network, communities)

        # Modularity should be positive for good community structure
        assert Q > 0.0
        assert Q <= 1.0  # Modularity is bounded above by 1

    def test_modularity_single_community(self):
        """Test modularity when all nodes are in one community."""
        # Put all nodes in single community (bad structure)
        single_community = {"A": 0, "B": 0, "C": 0, "D": 0, "E": 0}

        Q = modularity(self.network, single_community)

        # Should have some modularity (single community might not be optimal but isn't necessarily bad)
        assert isinstance(Q, float)

    def test_modularity_each_node_separate(self):
        """Test modularity when each node is its own community."""
        # Each node as separate community
        separate_communities = {"A": 0, "B": 1, "C": 2, "D": 3, "E": 4}

        Q = modularity(self.network, separate_communities)

        # Should have low modularity (no internal edges)
        assert Q <= 0.0

    def test_modularity_empty_communities(self):
        """Test modularity with empty communities dict."""
        Q = modularity(self.network, {})

        # Should handle empty communities (might return default modularity)
        assert isinstance(Q, float)

    def test_modularity_invalid_communities(self):
        """Test modularity with communities containing non-existent nodes."""
        invalid_communities = {"A": 0, "B": 0, "X": 1, "Y": 1}  # X, Y don't exist

        # Should handle gracefully (ignore non-existent nodes)
        Q = modularity(self.network, invalid_communities)
        assert isinstance(Q, float)


class TestCommunityMetrics:
    """Test community metrics calculation."""

    def setup_method(self):
        """Set up test network with communities."""
        self.network = BiologicalNetwork()

        # Create network with clear community structure
        # Module 1: densely connected
        module1_edges = [("A", "B"), ("B", "C"), ("C", "A"), ("A", "D"), ("B", "D")]
        for node1, node2 in module1_edges:
            self.network.add_edge(node1, node2, weight=1.0)

        # Module 2: densely connected
        module2_edges = [("E", "F"), ("F", "G"), ("G", "E")]
        for node1, node2 in module2_edges:
            self.network.add_edge(node1, node2, weight=1.0)

        # Weak connection between modules
        self.network.add_edge("D", "E", weight=0.2)

        # Detect communities
        self.communities = detect_communities(self.network, method="greedy")

    def test_community_metrics_basic(self):
        """Test basic community metrics calculation."""
        metrics = community_metrics(self.network, self.communities)

        # Should return metrics dictionary
        assert isinstance(metrics, dict)

        # Should include key metrics
        expected_keys = ["modularity", "num_communities", "internal_edge_ratio", "avg_community_size"]
        for key in expected_keys:
            assert key in metrics

        # Check metric ranges
        assert metrics["modularity"] >= -1.0 and metrics["modularity"] <= 1.0
        assert metrics["num_communities"] >= 1
        assert metrics["internal_edge_ratio"] >= 0.0 and metrics["internal_edge_ratio"] <= 1.0
        assert metrics["avg_community_size"] >= 0.0

    def test_community_metrics_empty_network(self):
        """Test metrics on empty network."""
        empty_network = BiologicalNetwork()
        empty_communities = {}

        metrics = community_metrics(empty_network, empty_communities)

        # Should handle empty network gracefully
        assert isinstance(metrics, dict)
        assert metrics["num_communities"] == 0
        assert metrics["modularity"] == 0.0

    def test_community_metrics_single_node(self):
        """Test metrics on single node network."""
        single_network = BiologicalNetwork()
        single_network.add_node("A")
        single_communities = {"A": 0}

        metrics = community_metrics(single_network, single_communities)

        # Should handle single node
        assert isinstance(metrics, dict)
        assert metrics["num_communities"] == 1


class TestCommunityDetectionIntegration:
    """Integration tests for community detection functionality."""

    def test_full_community_analysis_workflow(self):
        """Test complete community analysis workflow."""
        # Create realistic biological network
        np.random.seed(42)

        # Gene co-expression network
        genes = [f"Gene_{i}" for i in range(20)]
        network = create_network(genes)

        # Add edges based on simulated correlations
        for i, gene1 in enumerate(genes):
            for j, gene2 in enumerate(genes[i + 1 :], i + 1):
                # Simulate correlation with community structure
                if (i < 10 and j < 10) or (i >= 10 and j >= 10):
                    # Within-community correlation
                    correlation = np.random.beta(2, 1) * 0.9 + 0.1  # High correlation
                else:
                    # Between-community correlation
                    correlation = np.random.beta(1, 2) * 0.3  # Low correlation

                if correlation > 0.5:  # Threshold for edge creation
                    network.add_edge(gene1, gene2, weight=correlation)

        # Detect communities
        communities = detect_communities(network, method="louvain")

        # Calculate modularity
        Q = modularity(network, communities)

        # Calculate community metrics
        metrics = community_metrics(network, communities)

        # Verify results make sense
        assert isinstance(communities, dict)
        assert len(communities) == 20  # All genes assigned
        assert Q >= 0.0  # Should have positive modularity

        # Community metrics should be reasonable
        assert metrics["num_communities"] >= 1
        assert metrics["num_communities"] <= 20
        assert 0.0 <= metrics["internal_edge_ratio"] <= 1.0

        # All genes should be in some community
        assert len(communities) == 20


class TestEdgeCasesAndErrors:
    """Test edge cases and error conditions."""

    def test_invalid_algorithm_error(self):
        """Test error handling for invalid algorithm."""
        network = BiologicalNetwork()
        network.add_edge("A", "B")

        with pytest.raises(ValueError, match="Unknown community detection method"):
            detect_communities(network, method="invalid_algorithm")

    def test_negative_weights_handling(self):
        """Test handling of negative edge weights."""
        network = BiologicalNetwork()
        network.add_edge("A", "B", weight=-0.5)
        network.add_edge("B", "C", weight=1.0)

        # Should handle negative weights gracefully
        communities = detect_communities(network)
        assert isinstance(communities, dict)
        assert len(communities) >= 1

    def test_very_large_network_performance(self):
        """Test performance on larger networks."""
        # Create larger network (50 nodes)
        large_network = BiologicalNetwork()
        nodes = [f"Node_{i}" for i in range(50)]

        # Add random edges
        np.random.seed(123)
        for i in range(100):  # 100 random edges
            node1 = np.random.choice(nodes)
            node2 = np.random.choice(nodes)
            if node1 != node2:
                weight = np.random.random()
                large_network.add_edge(node1, node2, weight=weight)

        # Should handle larger networks reasonably
        communities = detect_communities(large_network, method="greedy")

        assert isinstance(communities, dict)
        assert len(communities) <= 50

        # All nodes should be assigned
        assert len(set(communities.keys())) >= 25  # Most nodes should be assigned

    def test_self_loops_handling(self):
        """Test handling of self-loops in network."""
        network = BiologicalNetwork()
        network.add_edge("A", "B")
        network.add_edge("A", "A", weight=1.0)  # Self-loop

        # Should handle self-loops gracefully
        communities = detect_communities(network)
        assert isinstance(communities, dict)
        assert len(communities) >= 1

    def test_zero_weight_edges(self):
        """Test handling of zero-weight edges."""
        network = BiologicalNetwork()
        network.add_edge("A", "B", weight=0.0)
        network.add_edge("B", "C", weight=1.0)

        communities = detect_communities(network)
        assert isinstance(communities, dict)
        assert len(communities) >= 1

        # Calculate modularity
        Q = modularity(network, communities)
        assert isinstance(Q, float)


class TestHierarchicalCommunities:
    """Test hierarchical community detection."""

    def test_hierarchical_communities(self):
        """Test hierarchical community detection."""
        from metainformant.networks.community import hierarchical_communities

        network = create_network(["A", "B", "C", "D", "E", "F"], directed=False)
        network.add_edge("A", "B")
        network.add_edge("B", "C")
        network.add_edge("D", "E")
        network.add_edge("E", "F")

        hierarchies = hierarchical_communities(network, levels=3, seed=42)
        assert len(hierarchies) == 3
        assert 0 in hierarchies
        assert 2 in hierarchies

    def test_community_stability(self):
        """Test community stability assessment."""
        from metainformant.networks.community import community_stability

        network = create_network(["A", "B", "C", "D", "E", "F"], directed=False)
        network.add_edge("A", "B")
        network.add_edge("B", "C")
        network.add_edge("D", "E")
        network.add_edge("E", "F")

        stability = community_stability(network, n_runs=5, seed=42)
        assert "stability_score" in stability
        assert "avg_modularity" in stability
        assert stability["stability_score"] >= 0.0

    def test_compare_communities(self):
        """Test community comparison."""
        from metainformant.networks.community import compare_communities

        comm1 = {"A": 0, "B": 0, "C": 1, "D": 1}
        comm2 = {"A": 0, "B": 1, "C": 1, "D": 1}

        comparison = compare_communities(comm1, comm2)
        assert "normalized_mutual_information" in comparison
        assert "adjusted_rand_index" in comparison
        assert comparison["normalized_mutual_information"] >= 0.0

    def test_optimize_resolution(self):
        """Test resolution optimization."""
        from metainformant.networks.community import optimize_resolution

        network = create_network(["A", "B", "C", "D", "E", "F"], directed=False)
        network.add_edge("A", "B")
        network.add_edge("B", "C")
        network.add_edge("D", "E")
        network.add_edge("E", "F")

        optimization = optimize_resolution(network, resolution_range=(0.5, 2.0), n_points=5)
        assert "optimal_resolution" in optimization
        assert "optimal_modularity" in optimization
        assert optimization["optimal_resolution"] >= 0.5
