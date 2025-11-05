"""Comprehensive tests for networks module.

Tests cover graph operations, community detection, pathway analysis,
protein-protein interactions, and regulatory networks.
"""

from __future__ import annotations

import numpy as np
import pytest

from metainformant.networks.community import community_metrics, detect_communities, modularity
from metainformant.networks.graph import (
    BiologicalNetwork,
    add_edges_from_correlation,
    add_edges_from_interactions,
    centrality_measures,
    create_network,
    network_metrics,
    shortest_paths,
)
from metainformant.networks.pathway import PathwayNetwork, pathway_enrichment
from metainformant.networks.ppi import ProteinNetwork, predict_interactions
from metainformant.networks.regulatory import GeneRegulatoryNetwork, infer_grn, regulatory_motifs


class TestBiologicalNetwork:
    """Tests for BiologicalNetwork class."""

    def test_network_creation(self):
        """Test creating empty network."""
        network = BiologicalNetwork(directed=False)
        assert network.num_nodes() == 0
        assert network.num_edges() == 0

    def test_add_nodes_and_edges(self):
        """Test adding nodes and edges."""
        network = BiologicalNetwork(directed=False)
        network.add_node("A")
        network.add_node("B", function="transcription")
        network.add_edge("A", "B", weight=0.8)

        assert network.num_nodes() == 2
        assert network.num_edges() == 1
        assert network.get_edge_weight("A", "B") == 0.8

    def test_directed_vs_undirected(self):
        """Test directed vs undirected network behavior."""
        # Undirected
        undir = BiologicalNetwork(directed=False)
        undir.add_edge("A", "B")
        assert undir.get_edge_weight("B", "A") == 1.0

        # Directed
        dir_net = BiologicalNetwork(directed=True)
        dir_net.add_edge("A", "B")
        assert dir_net.get_edge_weight("A", "B") == 1.0
        # Should not have reverse edge automatically
        assert dir_net.get_edge_weight("B", "A") is None

    def test_network_density(self):
        """Test network density calculation."""
        network = create_network(["A", "B", "C", "D"], directed=False)
        network.add_edge("A", "B")
        network.add_edge("A", "C")
        network.add_edge("B", "C")

        # 3 edges out of 6 possible (for 4 nodes, undirected)
        assert abs(network.density() - (3 / 6)) < 1e-10

    def test_get_neighbors(self):
        """Test neighbor retrieval."""
        network = create_network(["A", "B", "C"], directed=False)
        network.add_edge("A", "B")
        network.add_edge("A", "C")

        neighbors = network.get_neighbors("A")
        assert set(neighbors) == {"B", "C"}


class TestNetworkMetrics:
    """Tests for network metrics calculations."""

    def test_network_metrics_basic(self):
        """Test basic network metrics."""
        network = create_network(["A", "B", "C", "D"], directed=False)
        network.add_edge("A", "B")
        network.add_edge("A", "C")
        network.add_edge("B", "C")
        network.add_edge("C", "D")

        metrics = network_metrics(network)

        assert metrics["num_nodes"] == 4
        assert metrics["num_edges"] == 4
        assert metrics["avg_degree"] > 0
        assert metrics["max_degree"] >= metrics["avg_degree"]

    def test_centrality_measures(self):
        """Test centrality calculations."""
        network = create_network(["A", "B", "C", "D"], directed=False)
        # A is hub (connected to all others)
        network.add_edge("A", "B")
        network.add_edge("A", "C")
        network.add_edge("A", "D")
        network.add_edge("B", "C")

        centralities = centrality_measures(network)

        assert "degree" in centralities
        assert "closeness" in centralities
        assert "betweenness" in centralities
        assert "eigenvector" in centralities

        # A should have highest degree centrality
        assert centralities["degree"]["A"] > centralities["degree"]["B"]

    def test_shortest_paths(self):
        """Test shortest path calculations."""
        network = create_network(["A", "B", "C", "D"], directed=False)
        network.add_edge("A", "B")
        network.add_edge("B", "C")
        network.add_edge("C", "D")

        paths = shortest_paths(network)

        assert paths["A"]["D"] == 3.0  # A -> B -> C -> D
        assert paths["A"]["A"] == 0.0  # Distance to self is 0


class TestCommunityDetection:
    """Tests for community detection algorithms."""

    def test_detect_communities_louvain(self):
        """Test Louvain community detection."""
        network = create_network(["A", "B", "C", "D", "E", "F"], directed=False)
        # Two clear communities: {A, B, C} and {D, E, F}
        network.add_edge("A", "B")
        network.add_edge("A", "C")
        network.add_edge("B", "C")
        network.add_edge("D", "E")
        network.add_edge("D", "F")
        network.add_edge("E", "F")

        communities = detect_communities(network, method="louvain", seed=42)

        assert len(communities) == 6
        # Nodes in same community should have same ID
        assert communities["A"] == communities["B"] == communities["C"]
        assert communities["D"] == communities["E"] == communities["F"]

    def test_modularity_calculation(self):
        """Test modularity calculation."""
        network = create_network(["A", "B", "C", "D"], directed=False)
        network.add_edge("A", "B")
        network.add_edge("C", "D")

        # Partition into two communities
        communities = {"A": 0, "B": 0, "C": 1, "D": 1}

        mod = modularity(network, communities)

        assert mod > 0.0  # Should have positive modularity for good partition

    def test_community_metrics(self):
        """Test community metrics calculation."""
        network = create_network(["A", "B", "C", "D", "E"], directed=False)
        network.add_edge("A", "B")
        network.add_edge("B", "C")
        network.add_edge("D", "E")

        communities = {"A": 0, "B": 0, "C": 0, "D": 1, "E": 1}

        metrics = community_metrics(network, communities)

        assert metrics["num_communities"] == 2
        assert metrics["avg_community_size"] == 2.5
        assert metrics["internal_edges"] >= 0


class TestPathwayNetwork:
    """Tests for pathway network functionality."""

    def test_pathway_network_creation(self):
        """Test creating pathway network."""
        pn = PathwayNetwork(name="test_pathways")
        pn.add_pathway("path1", ["GENE1", "GENE2", "GENE3"], metadata={"name": "Pathway 1"})

        assert len(pn.pathways) == 1
        assert pn.get_pathway_genes("path1") == {"GENE1", "GENE2", "GENE3"}

    def test_pathway_overlap(self):
        """Test pathway overlap calculation."""
        pn = PathwayNetwork()
        pn.add_pathway("path1", ["GENE1", "GENE2", "GENE3"])
        pn.add_pathway("path2", ["GENE2", "GENE3", "GENE4"])

        overlap, jaccard = pn.pathway_overlap("path1", "path2")

        assert overlap == {"GENE2", "GENE3"}
        assert jaccard > 0.0

    def test_pathway_enrichment(self):
        """Test pathway enrichment analysis."""
        pn = PathwayNetwork()
        pn.add_pathway("path1", ["GENE1", "GENE2", "GENE3", "GENE4"])
        pn.add_pathway("path2", ["GENE5", "GENE6", "GENE7"])

        query_genes = ["GENE1", "GENE2", "GENE5"]

        enrichment = pathway_enrichment(query_genes, pn)

        assert len(enrichment) > 0
        # path1 should be enriched (2 out of 3 query genes)
        enriched_pathways = [e["pathway_id"] for e in enrichment]
        assert "path1" in enriched_pathways


class TestProteinNetwork:
    """Tests for protein-protein interaction networks."""

    def test_protein_network_creation(self):
        """Test creating PPI network."""
        ppi = ProteinNetwork(name="test_ppi")
        ppi.add_interaction("P1", "P2", confidence=0.8, evidence_types=["experimental"])

        assert len(ppi.proteins) == 2
        assert len(ppi.interactions) == 1

    def test_create_network_with_confidence_filter(self):
        """Test creating network with confidence threshold."""
        ppi = ProteinNetwork()
        ppi.add_interaction("P1", "P2", confidence=0.9)
        ppi.add_interaction("P2", "P3", confidence=0.5)
        ppi.add_interaction("P3", "P4", confidence=0.3)

        # High confidence only
        network = ppi.create_network(min_confidence=0.7)

        assert network.num_edges() == 1  # Only P1-P2 above threshold

    def test_predict_interactions_correlation(self):
        """Test interaction prediction from features."""
        np.random.seed(42)
        # Create correlated feature vectors for some proteins
        n_proteins = 20
        n_features = 50
        features = np.random.randn(n_proteins, n_features)

        # Make some proteins highly correlated
        features[0] = features[1] + 0.1 * np.random.randn(n_features)

        protein_ids = [f"P{i:03d}" for i in range(n_proteins)]

        ppi = predict_interactions(features, protein_ids, method="correlation", threshold=0.8)

        assert len(ppi.proteins) > 0
        assert len(ppi.interactions) > 0


class TestGeneRegulatoryNetwork:
    """Tests for gene regulatory networks."""

    def test_grn_creation(self):
        """Test creating GRN."""
        grn = GeneRegulatoryNetwork(name="test_grn")
        grn.add_regulation("TF1", "GENE1", regulation_type="activation", strength=0.85, confidence=0.9)

        assert len(grn.genes) == 2
        assert "TF1" in grn.transcription_factors

    def test_get_regulators_and_targets(self):
        """Test retrieving regulators and targets."""
        grn = GeneRegulatoryNetwork()
        grn.add_regulation("TF1", "GENE1", regulation_type="activation")
        grn.add_regulation("TF2", "GENE1", regulation_type="repression")
        grn.add_regulation("TF1", "GENE2", regulation_type="activation")

        regulators = grn.get_regulators("GENE1")
        assert len(regulators) == 2

        targets = grn.get_targets("TF1")
        assert len(targets) == 2

    def test_infer_grn_correlation(self):
        """Test GRN inference from expression data."""
        np.random.seed(42)
        n_samples = 50
        n_genes = 20

        # Create expression data with some correlations
        expression = np.random.randn(n_samples, n_genes)
        # Make first gene correlated with others
        expression[:, 1] = 0.8 * expression[:, 0] + 0.2 * np.random.randn(n_samples)

        gene_names = [f"GENE_{i}" for i in range(n_genes)]
        tf_list = [f"GENE_{i}" for i in range(5)]  # First 5 are TFs

        grn = infer_grn(expression, gene_names, method="correlation", tf_list=tf_list, threshold=0.7)

        assert len(grn.genes) > 0
        assert len(grn.regulations) > 0

    def test_regulatory_motifs(self):
        """Test regulatory motif detection."""
        grn = GeneRegulatoryNetwork()
        # Create feed-forward loop: A -> B -> C, A -> C
        grn.add_regulation("A", "B", regulation_type="activation")
        grn.add_regulation("B", "C", regulation_type="activation")
        grn.add_regulation("A", "C", regulation_type="activation")

        motifs = regulatory_motifs(grn, motif_types=["feed_forward_loop"])

        assert "feed_forward_loop" in motifs
        assert len(motifs["feed_forward_loop"]) > 0


class TestGraphUtilities:
    """Tests for graph utility functions."""

    def test_add_edges_from_correlation(self):
        """Test adding edges based on correlation matrix."""
        network = create_network(["A", "B", "C", "D"], directed=False)

        # Create correlation matrix
        corr_matrix = np.array([[1.0, 0.8, 0.3, 0.2], [0.8, 1.0, 0.4, 0.1], [0.3, 0.4, 1.0, 0.9], [0.2, 0.1, 0.9, 1.0]])
        node_names = ["A", "B", "C", "D"]

        add_edges_from_correlation(network, corr_matrix, node_names, threshold=0.7)

        assert network.num_edges() > 0
        assert network.get_edge_weight("A", "B") is not None

    def test_add_edges_from_interactions(self):
        """Test adding edges from interaction list."""
        network = create_network(["A", "B", "C"], directed=False)

        interactions = [("A", "B", 0.8), ("B", "C", 0.6), ("A", "C", 0.9)]

        add_edges_from_interactions(network, interactions)

        assert network.num_edges() == 3
        assert network.get_edge_weight("A", "C") == 0.9

    def test_export_import_roundtrip(self, tmp_path):
        """Test network export and import round-trip."""
        from metainformant.networks.graph import export_network, import_network

        original = create_network(["A", "B", "C"], directed=False)
        original.add_edge("A", "B", weight=0.8)
        original.add_node("A", test_attr="test_value")

        json_file = tmp_path / "test_network.json"
        export_network(original, str(json_file), format="json")
        loaded = import_network(str(json_file), format="json")

        assert loaded.num_nodes() == original.num_nodes()
        assert loaded.num_edges() == original.num_edges()
        assert loaded.get_edge_weight("A", "B") == original.get_edge_weight("A", "B")
        assert loaded.node_attrs["A"]["test_attr"] == "test_value"

    def test_hierarchical_community_detection(self):
        """Test hierarchical community detection."""
        from metainformant.networks.community import hierarchical_communities

        network = create_network(["A", "B", "C", "D", "E", "F"], directed=False)
        network.add_edge("A", "B")
        network.add_edge("B", "C")
        network.add_edge("D", "E")
        network.add_edge("E", "F")

        hierarchies = hierarchical_communities(network, levels=3, seed=42)
        assert len(hierarchies) == 3

    def test_protein_complex_detection(self):
        """Test protein complex detection."""
        from metainformant.networks.ppi import detect_complexes

        ppi = ProteinNetwork()
        # Create dense subgraph
        for i in range(4):
            for j in range(i + 1, 4):
                ppi.add_interaction(f"P{i}", f"P{j}", confidence=0.8)

        complexes = detect_complexes(ppi, min_confidence=0.7, min_size=3)
        assert len(complexes) > 0

    def test_regulatory_cascade_detection(self):
        """Test regulatory cascade detection."""
        from metainformant.networks.regulatory import detect_regulatory_cascades

        grn = GeneRegulatoryNetwork()
        grn.add_regulation("TF1", "TF2", confidence=0.8)
        grn.add_regulation("TF2", "GENE1", confidence=0.9)

        cascades = detect_regulatory_cascades(grn, max_length=5)
        assert isinstance(cascades, list)

    def test_pathway_similarity_and_activity(self):
        """Test pathway similarity and activity scoring."""
        from metainformant.networks.pathway import pathway_similarity, pathway_activity_score

        pn = PathwayNetwork()
        pn.add_pathway("path1", ["GENE1", "GENE2", "GENE3"])
        pn.add_pathway("path2", ["GENE2", "GENE3", "GENE4"])

        path1_genes = pn.get_pathway_genes("path1")
        path2_genes = pn.get_pathway_genes("path2")

        similarity = pathway_similarity(path1_genes, path2_genes)
        assert 0.0 <= similarity <= 1.0

        expression = {"GENE1": 10.0, "GENE2": 8.0, "GENE3": 12.0}
        activity = pathway_activity_score(pn, "path1", expression)
        assert activity > 0.0





