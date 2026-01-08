"""Tests for biological pathway network functionality.

Real implementation testing for pathway analysis methods.
No mocking used - all tests use real computational methods and data.
"""

from __future__ import annotations

from typing import Dict, List, Set

import numpy as np
import pytest

from metainformant.networks.analysis.graph import BiologicalNetwork
from metainformant.networks.analysis.pathway import (
    PathwayNetwork,
    load_pathway_database,
    network_enrichment_analysis,
    pathway_enrichment,
)


class TestPathwayNetwork:
    """Test PathwayNetwork functionality."""

    def setup_method(self):
        """Set up test pathway network."""
        self.pathway_net = PathwayNetwork("test_pathway")

        # Add pathways
        self.pathway_net.add_pathway(
            "pathway_001",
            ["GENE1", "GENE2", "GENE3", "GENE4"],
            metadata={"name": "Test Pathway 1", "description": "First test pathway"},
        )

        self.pathway_net.add_pathway(
            "pathway_002",
            ["GENE3", "GENE4", "GENE5", "GENE6"],
            metadata={"name": "Test Pathway 2", "description": "Second test pathway"},
        )

        self.pathway_net.add_pathway(
            "pathway_003", ["GENE7", "GENE8"], metadata={"name": "Test Pathway 3", "description": "Third test pathway"}
        )

    def test_pathway_network_initialization(self):
        """Test PathwayNetwork initialization."""
        pathway_net = PathwayNetwork("test_network")

        assert pathway_net.name == "test_network"
        assert len(pathway_net.pathways) == 0
        assert len(pathway_net.gene_pathways) == 0

    def test_add_pathway(self):
        """Test adding pathways."""
        pathway_net = PathwayNetwork()

        # Add pathway with metadata
        pathway_net.add_pathway(
            "test_pathway", ["BRCA1", "BRCA2", "TP53"], metadata={"name": "DNA Repair", "category": "cancer"}
        )

        assert "test_pathway" in pathway_net.pathways
        assert pathway_net.pathways["test_pathway"] == {"BRCA1", "BRCA2", "TP53"}

        # Check gene-pathway mappings
        assert "test_pathway" in pathway_net.gene_pathways["BRCA1"]
        assert "test_pathway" in pathway_net.gene_pathways["BRCA2"]
        assert "test_pathway" in pathway_net.gene_pathways["TP53"]

        # Check metadata
        assert pathway_net.pathway_metadata["test_pathway"]["name"] == "DNA Repair"
        assert pathway_net.pathway_metadata["test_pathway"]["category"] == "cancer"

    def test_get_pathway_genes(self):
        """Test getting genes in a pathway."""
        genes = self.pathway_net.get_pathway_genes("pathway_001")
        expected_genes = {"GENE1", "GENE2", "GENE3", "GENE4"}
        assert genes == expected_genes

        # Test non-existent pathway
        empty_genes = self.pathway_net.get_pathway_genes("nonexistent")
        assert empty_genes == set()

    def test_get_gene_pathways(self):
        """Test getting pathways for a gene."""
        # GENE3 and GENE4 are in both pathway_001 and pathway_002
        pathways_gene3 = self.pathway_net.get_gene_pathways("GENE3")
        assert pathways_gene3 == {"pathway_001", "pathway_002"}

        pathways_gene4 = self.pathway_net.get_gene_pathways("GENE4")
        assert pathways_gene4 == {"pathway_001", "pathway_002"}

        # GENE1 is only in pathway_001
        pathways_gene1 = self.pathway_net.get_gene_pathways("GENE1")
        assert pathways_gene1 == {"pathway_001"}

        # Non-existent gene
        pathways_nonexistent = self.pathway_net.get_gene_pathways("NONEXISTENT")
        assert pathways_nonexistent == set()

    def test_pathway_similarity(self):
        """Test pathway similarity calculation."""
        # pathway_001 and pathway_002 share GENE3 and GENE4
        similarity = self.pathway_net.pathway_similarity("pathway_001", "pathway_002")

        # Jaccard similarity: |intersection| / |union|
        # intersection = {GENE3, GENE4} = 2 genes
        # union = {GENE1, GENE2, GENE3, GENE4, GENE5, GENE6} = 6 genes
        # similarity = 2/6 = 0.333...
        expected_similarity = 2.0 / 6.0
        assert abs(similarity - expected_similarity) < 0.01

        # pathway_001 and pathway_003 share no genes
        similarity_no_overlap = self.pathway_net.pathway_similarity("pathway_001", "pathway_003")
        assert similarity_no_overlap == 0.0

        # Self-similarity should be 1.0
        self_similarity = self.pathway_net.pathway_similarity("pathway_001", "pathway_001")
        assert self_similarity == 1.0

    def test_get_pathway_statistics(self):
        """Test pathway statistics calculation."""
        stats = self.pathway_net.get_statistics()

        # Should return reasonable statistics
        assert stats["num_pathways"] == 3
        assert stats["num_genes"] == 8  # GENE1-GENE8 (all unique)
        assert stats["avg_pathway_size"] == (4 + 4 + 2) / 3  # Average size
        assert stats["max_pathway_size"] == 4
        assert stats["min_pathway_size"] == 2

    def test_find_overlapping_pathways(self):
        """Test finding pathways that share genes."""
        overlapping = self.pathway_net.find_overlapping_pathways("pathway_001")

        # pathway_001 overlaps with pathway_002 (shares GENE3, GENE4)
        assert "pathway_002" in overlapping
        assert "pathway_003" not in overlapping  # No shared genes

        # Check overlap details
        overlap_info = overlapping["pathway_002"]
        assert overlap_info["shared_genes"] == {"GENE3", "GENE4"}
        assert overlap_info["overlap_size"] == 2


class TestLoadPathwayDatabase:
    """Test pathway database loading functionality."""

    def test_load_from_dict(self):
        """Test loading pathways from dictionary format."""
        pathway_data = {
            "KEGG:hsa00010": {"name": "Glycolysis", "genes": ["HK1", "GPI", "PFKL", "ALDOA"], "category": "Metabolism"},
            "KEGG:hsa00020": {"name": "TCA Cycle", "genes": ["CS", "ACO2", "IDH1", "OGDH"], "category": "Metabolism"},
        }

        pathway_network = load_pathway_database(pathway_data)

        assert isinstance(pathway_network, PathwayNetwork)
        assert len(pathway_network.pathways) == 2
        assert "KEGG:hsa00010" in pathway_network.pathways
        assert "KEGG:hsa00020" in pathway_network.pathways

        # Check pathway content
        glycolysis_genes = pathway_network.get_pathway_genes("KEGG:hsa00010")
        assert glycolysis_genes == {"HK1", "GPI", "PFKL", "ALDOA"}

        # Check metadata
        assert pathway_network.pathway_metadata["KEGG:hsa00010"]["name"] == "Glycolysis"
        assert pathway_network.pathway_metadata["KEGG:hsa00010"]["category"] == "Metabolism"

    def test_load_empty_database(self):
        """Test loading empty pathway database."""
        empty_data = {}
        pathway_network = load_pathway_database(empty_data)

        assert isinstance(pathway_network, PathwayNetwork)
        assert len(pathway_network.pathways) == 0

    def test_load_minimal_pathways(self):
        """Test loading pathways with minimal information."""
        minimal_data = {
            "minimal_pathway": {
                "genes": ["GENE_A", "GENE_B"]
                # Missing name and other metadata
            }
        }

        pathway_network = load_pathway_database(minimal_data)

        assert len(pathway_network.pathways) == 1
        assert pathway_network.get_pathway_genes("minimal_pathway") == {"GENE_A", "GENE_B"}

    def test_load_invalid_format(self):
        """Test error handling for invalid pathway format."""
        invalid_data = {
            "invalid_pathway": {
                "genes": "not_a_list",  # Should be list
            }
        }

        # Should handle gracefully or raise informative error
        try:
            pathway_network = load_pathway_database(invalid_data)
            assert isinstance(pathway_network, PathwayNetwork)
        except (ValueError, TypeError):
            # Acceptable to raise error for invalid format
            pass


class TestPathwayEnrichment:
    """Test pathway enrichment analysis."""

    def setup_method(self):
        """Set up test data for enrichment analysis."""
        # Create test pathway network
        self.pathway_network = PathwayNetwork()

        # Add pathways with different sizes
        self.pathway_network.add_pathway(
            "cell_cycle",
            [f"CELL_CYCLE_GENE_{i}" for i in range(10)],
            metadata={"name": "Cell Cycle", "category": "cell_division"},
        )

        self.pathway_network.add_pathway(
            "apoptosis",
            [f"APOPTOSIS_GENE_{i}" for i in range(8)],
            metadata={"name": "Apoptosis", "category": "cell_death"},
        )

        self.pathway_network.add_pathway(
            "metabolism",
            [f"METABOLISM_GENE_{i}" for i in range(15)],
            metadata={"name": "Metabolism", "category": "energy"},
        )

        # Test gene set (differentially expressed genes)
        # Include 6 cell cycle genes, 2 apoptosis genes, 1 metabolism gene
        self.test_gene_set = (
            [f"CELL_CYCLE_GENE_{i}" for i in range(6)]
            + [f"APOPTOSIS_GENE_{i}" for i in range(2)]
            + [f"METABOLISM_GENE_0"]
        )

    def test_pathway_enrichment_basic(self):
        """Test basic pathway enrichment analysis."""
        results = pathway_enrichment(gene_list=self.test_gene_set, pathway_network=self.pathway_network)

        # Should return results for all pathways
        assert isinstance(results, dict)
        assert len(results) == 3
        assert all(pathway_id in results for pathway_id in ["cell_cycle", "apoptosis", "metabolism"])

        # Check result structure
        for pathway_id, result in results.items():
            assert "overlap_size" in result
            assert "pathway_size" in result
            assert "p_value" in result
            assert "enrichment_ratio" in result
            assert "overlapping_genes" in result

        # Cell cycle should be most enriched (6/10 genes in test set)
        cell_cycle_result = results["cell_cycle"]
        assert cell_cycle_result["overlap_size"] == 6
        assert cell_cycle_result["pathway_size"] == 10
        assert cell_cycle_result["enrichment_ratio"] > 1.0  # Should be enriched

    def test_pathway_enrichment_significance(self):
        """Test statistical significance of enrichment."""
        results = pathway_enrichment(gene_list=self.test_gene_set, pathway_network=self.pathway_network)

        # Check p-values
        for pathway_id, result in results.items():
            assert 0.0 <= result["p_value"] <= 1.0

        # Cell cycle should have lowest p-value (most significant)
        p_values = {pid: res["p_value"] for pid, res in results.items()}
        most_significant = min(p_values, key=p_values.get)
        assert most_significant == "cell_cycle"

    def test_pathway_enrichment_no_overlap(self):
        """Test enrichment when gene set has no overlap with pathways."""
        # Gene set with no overlap
        no_overlap_genes = [f"UNRELATED_GENE_{i}" for i in range(5)]

        results = pathway_enrichment(gene_list=no_overlap_genes, pathway_network=self.pathway_network)

        # All pathways should have zero overlap
        for result in results.values():
            assert result["overlap_size"] == 0
            assert result["enrichment_ratio"] == 0.0
            # p-value should be high (not significant)
            assert result["p_value"] >= 0.5

    def test_pathway_enrichment_empty_gene_set(self):
        """Test enrichment with empty gene set."""
        results = pathway_enrichment(gene_list=[], pathway_network=self.pathway_network)

        # Should handle empty gene set gracefully
        assert isinstance(results, dict)
        assert len(results) == 3
        for result in results.values():
            assert result["overlap_size"] == 0


class TestNetworkEnrichmentAnalysis:
    """Test network-based enrichment analysis."""

    def setup_method(self):
        """Set up network and pathway data."""
        # Create protein interaction network
        self.ppi_network = BiologicalNetwork()

        # Add proteins and interactions
        proteins = ["PROTEIN_A", "PROTEIN_B", "PROTEIN_C", "PROTEIN_D", "PROTEIN_E"]
        for protein in proteins:
            self.ppi_network.add_node(protein)

        # Add interactions
        interactions = [
            ("PROTEIN_A", "PROTEIN_B", 0.9),
            ("PROTEIN_B", "PROTEIN_C", 0.8),
            ("PROTEIN_C", "PROTEIN_D", 0.7),
            ("PROTEIN_A", "PROTEIN_E", 0.6),
        ]

        for p1, p2, conf in interactions:
            self.ppi_network.add_edge(p1, p2, weight=conf)

        # Create pathway network
        self.pathway_network = PathwayNetwork()
        self.pathway_network.add_pathway(
            "test_pathway", ["PROTEIN_A", "PROTEIN_B", "PROTEIN_C"], metadata={"name": "Test Pathway"}
        )

        # Test gene set
        self.test_genes = ["PROTEIN_A", "PROTEIN_B"]

    def test_network_enrichment_basic(self):
        """Test basic network enrichment analysis."""
        results = network_enrichment_analysis(
            gene_list=self.test_genes, ppi_network=self.ppi_network, pathway_network=self.pathway_network
        )

        # Should return enrichment results
        assert isinstance(results, dict)
        assert "pathway_enrichment" in results
        assert "network_connectivity" in results

        # Pathway enrichment should show overlap
        pathway_results = results["pathway_enrichment"]
        assert "test_pathway" in pathway_results
        assert pathway_results["test_pathway"]["overlap_size"] == 2

        # Network connectivity should be positive
        network_results = results["network_connectivity"]
        assert network_results["observed_edges"] >= 0
        assert network_results["connectivity_score"] >= 0.0


class TestPathwayNetworkIntegration:
    """Integration tests for pathway analysis."""

    def test_complete_pathway_analysis_workflow(self):
        """Test complete pathway analysis workflow."""
        # 1. Load pathway database
        pathway_data = {
            "glycolysis": {
                "name": "Glycolysis",
                "genes": ["HK1", "GPI", "PFKL", "ALDOA", "TPI1"],
                "category": "metabolism",
            },
            "gluconeogenesis": {
                "name": "Gluconeogenesis",
                "genes": ["G6PC", "FBP1", "PCK1", "HK1"],  # HK1 shared with glycolysis
                "category": "metabolism",
            },
        }

        pathway_network = load_pathway_database(pathway_data)

        # 2. Perform enrichment analysis
        test_genes = ["HK1", "GPI", "PFKL"]  # Enriched in glycolysis

        enrichment_results = pathway_enrichment(gene_list=test_genes, pathway_network=pathway_network)

        # 3. Analyze pathway similarity
        similarity = pathway_network.pathway_similarity("glycolysis", "gluconeogenesis")

        # 4. Get network statistics
        stats = pathway_network.get_statistics()

        # 5. Verify integrated results
        assert len(enrichment_results) == 2
        assert len(pathway_network.pathways) == 2

        # Glycolysis should be enriched
        glycolysis_enrichment = enrichment_results["glycolysis"]
        assert glycolysis_enrichment["overlap_size"] == 3
        assert glycolysis_enrichment["enrichment_ratio"] > 1.0

        # Should detect similarity (shared HK1 gene)
        assert similarity > 0.0

        # Statistics should be reasonable
        assert stats["num_pathways"] == 2
        assert stats["num_genes"] >= 7  # All unique genes


class TestPathwayEdgeCases:
    """Test edge cases and error conditions."""

    def test_empty_pathway_network(self):
        """Test operations on empty pathway network."""
        empty_network = PathwayNetwork()

        stats = empty_network.get_statistics()
        assert stats["num_pathways"] == 0
        assert stats["num_genes"] == 0

        # Enrichment should handle empty network
        results = pathway_enrichment(gene_list=["GENE1", "GENE2"], pathway_network=empty_network)
        assert isinstance(results, dict)
        assert len(results) == 0

    def test_duplicate_pathways(self):
        """Test handling of duplicate pathway addition."""
        pathway_network = PathwayNetwork()

        # Add same pathway twice
        pathway_network.add_pathway("dup_pathway", ["GENE1", "GENE2"])
        pathway_network.add_pathway("dup_pathway", ["GENE3", "GENE4"])  # Should overwrite

        assert len(pathway_network.pathways) == 1
        # Should have the latest version
        assert pathway_network.get_pathway_genes("dup_pathway") == {"GENE3", "GENE4"}

    def test_very_large_pathway(self):
        """Test handling of very large pathways."""
        pathway_network = PathwayNetwork()

        # Add pathway with many genes
        large_genes = [f"GENE_{i}" for i in range(1000)]
        pathway_network.add_pathway("large_pathway", large_genes)

        assert len(pathway_network.get_pathway_genes("large_pathway")) == 1000

        # Statistics should handle large pathway
        stats = pathway_network.get_statistics()
        assert stats["max_pathway_size"] == 1000

    def test_special_characters_in_names(self):
        """Test handling of special characters in pathway/gene names."""
        pathway_network = PathwayNetwork()

        # Pathway and gene names with special characters
        pathway_network.add_pathway(
            "pathway:with-special_chars.123",
            ["GENE-1", "GENE_2", "GENE.3", "GENE@4"],
            metadata={"name": "Pathway with special chars!"},
        )

        assert "pathway:with-special_chars.123" in pathway_network.pathways
        genes = pathway_network.get_pathway_genes("pathway:with-special_chars.123")
        assert "GENE-1" in genes
        assert "GENE_2" in genes
        assert "GENE.3" in genes
        assert "GENE@4" in genes


class TestNewPathwayFunctions:
    """Test new pathway functions."""

    def test_pathway_similarity(self):
        """Test pathway similarity calculation."""
        from metainformant.networks.analysis.pathway import pathway_similarity

        path1 = {"GENE1", "GENE2", "GENE3"}
        path2 = {"GENE2", "GENE3", "GENE4"}

        similarity = pathway_similarity(path1, path2, method="jaccard")
        assert similarity > 0.0
        assert similarity <= 1.0

        overlap = pathway_similarity(path1, path2, method="overlap")
        assert overlap > 0.0

        dice = pathway_similarity(path1, path2, method="dice")
        assert dice > 0.0

    def test_pathway_activity_score(self):
        """Test pathway activity scoring."""
        from metainformant.networks.analysis.pathway import pathway_activity_score

        pn = PathwayNetwork()
        pn.add_pathway("path1", ["GENE1", "GENE2", "GENE3"])

        expression = {"GENE1": 10.5, "GENE2": 8.2, "GENE3": 12.1}

        score = pathway_activity_score(pn, "path1", expression, method="mean")
        assert score > 0.0

        score_max = pathway_activity_score(pn, "path1", expression, method="max")
        assert score_max == 12.1

        score_sum = pathway_activity_score(pn, "path1", expression, method="sum")
        assert score_sum > score
