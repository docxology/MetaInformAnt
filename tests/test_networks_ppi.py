"""Tests for protein-protein interaction network functionality.

Real implementation testing for PPI network analysis methods.
No mocking used - all tests use real computational methods and data.
"""

from __future__ import annotations

from typing import Dict, List, Set

import numpy as np
import pandas as pd
import pytest

from metainformant.networks.graph import BiologicalNetwork
from metainformant.networks.ppi import (
    ProteinNetwork,
    functional_enrichment_ppi,
    load_string_interactions,
    predict_interactions,
)


class TestProteinNetwork:
    """Test ProteinNetwork functionality."""

    def setup_method(self):
        """Set up test protein network."""
        self.ppi_network = ProteinNetwork("test_ppi")

        # Add protein interactions
        interactions = [
            ("P1", "P2", 0.9, ["experimental"]),
            ("P2", "P3", 0.8, ["computational"]),
            ("P3", "P4", 0.7, ["experimental"]),
            ("P1", "P5", 0.6, ["database"]),
            ("P4", "P5", 0.5, ["textmining"]),
        ]

        for p1, p2, confidence, evidence_types in interactions:
            self.ppi_network.add_interaction(p1, p2, confidence, evidence_types)

        # Add protein metadata
        proteins_metadata = [
            ("P1", {"name": "Protein1", "length": 300, "function": "kinase"}),
            ("P2", {"name": "Protein2", "length": 250, "function": "phosphatase"}),
            ("P3", {"name": "Protein3", "length": 400, "function": "transcription_factor"}),
            ("P4", {"name": "Protein4", "length": 180, "function": "transporter"}),
            ("P5", {"name": "Protein5", "length": 500, "function": "enzyme"}),
        ]

        for protein_id, metadata in proteins_metadata:
            self.ppi_network.add_protein_metadata(protein_id, **metadata)

    def test_protein_network_initialization(self):
        """Test ProteinNetwork initialization."""
        network = ProteinNetwork("test_network")

        assert network.name == "test_network"
        assert len(network.interactions) == 0
        assert len(network.proteins) == 0

    def test_add_interaction(self):
        """Test adding protein interactions."""
        network = ProteinNetwork()

        # Add interaction with metadata
        network.add_interaction(
            "BRCA1",
            "BRCA2",
            confidence=0.95,
            evidence_types=["yeast_two_hybrid", "coimmunoprecipitation"],
            pubmed_id="12345678",
        )

        assert len(network.interactions) == 1
        assert "BRCA1" in network.proteins
        assert "BRCA2" in network.proteins

        # Check interaction details
        interaction = network.interactions[0]
        assert interaction[0] == "BRCA1"
        assert interaction[1] == "BRCA2"
        assert interaction[2]["confidence"] == 0.95
        assert "yeast_two_hybrid" in interaction[2]["evidence_types"]
        assert interaction[2]["pubmed_id"] == "12345678"

    def test_add_protein_metadata(self):
        """Test adding protein metadata."""
        network = ProteinNetwork()

        network.add_protein_metadata(
            "TP53",
            name="Tumor protein p53",
            uniprot_id="P04637",
            length=393,
            molecular_weight=43653,
            function="transcription factor",
        )

        assert "TP53" in network.protein_metadata
        metadata = network.protein_metadata["TP53"]
        assert metadata["name"] == "Tumor protein p53"
        assert metadata["length"] == 393
        assert metadata["function"] == "transcription factor"

    def test_get_protein_partners(self):
        """Test getting interaction partners of a protein."""
        partners_p1 = self.ppi_network.get_protein_partners("P1")
        expected_partners_p1 = {"P2", "P5"}
        assert set(partners_p1) == expected_partners_p1

        partners_p3 = self.ppi_network.get_protein_partners("P3")
        expected_partners_p3 = {"P2", "P4"}
        assert set(partners_p3) == expected_partners_p3

        # Test protein with no partners
        partners_nonexistent = self.ppi_network.get_protein_partners("NONEXISTENT")
        assert len(partners_nonexistent) == 0

    def test_filter_by_confidence(self):
        """Test filtering interactions by confidence threshold."""
        # High confidence threshold
        high_conf_network = self.ppi_network.filter_by_confidence(threshold=0.8)

        # Should include P1-P2 (0.9) and P2-P3 (0.8)
        high_conf_interactions = high_conf_network.interactions
        assert len(high_conf_interactions) == 2

        # Check specific interactions
        interaction_pairs = [(i[0], i[1]) for i in high_conf_interactions]
        assert ("P1", "P2") in interaction_pairs or ("P2", "P1") in interaction_pairs
        assert ("P2", "P3") in interaction_pairs or ("P3", "P2") in interaction_pairs

    def test_filter_by_evidence(self):
        """Test filtering interactions by evidence type."""
        # Get only experimental interactions
        exp_network = self.ppi_network.filter_by_evidence("experimental")

        # Should include P1-P2 and P3-P4 (both have experimental evidence)
        exp_interactions = exp_network.interactions
        assert len(exp_interactions) == 2

        # Check evidence types
        for interaction in exp_interactions:
            assert "experimental" in interaction[2]["evidence_types"]

    def test_get_network_statistics(self):
        """Test network statistics calculation."""
        stats = self.ppi_network.get_network_statistics()

        # Should include various metrics
        expected_keys = ["num_proteins", "num_interactions", "avg_confidence", "density"]
        for key in expected_keys:
            assert key in stats

        assert stats["num_proteins"] == 5
        assert stats["num_interactions"] == 5
        assert 0.0 <= stats["avg_confidence"] <= 1.0
        assert 0.0 <= stats["density"] <= 1.0

    def test_to_biological_network(self):
        """Test conversion to BiologicalNetwork."""
        bio_network = self.ppi_network.to_biological_network()

        assert isinstance(bio_network, BiologicalNetwork)
        assert bio_network.num_nodes() == 5
        assert bio_network.num_edges() == 5

        # Check that interactions are preserved
        assert bio_network.has_edge("P1", "P2")
        assert bio_network.has_edge("P2", "P3")
        assert bio_network.has_edge("P3", "P4")


class TestLoadStringInteractions:
    """Test STRING database interaction loading."""

    def test_load_string_from_dataframe(self):
        """Test loading STRING interactions from DataFrame."""
        # Create sample STRING data
        interactions_data = pd.DataFrame(
            {
                "protein1": ["ENSP1", "ENSP2", "ENSP3"],
                "protein2": ["ENSP2", "ENSP3", "ENSP4"],
                "combined_score": [900, 750, 600],
                "experimental": [200, 100, 0],
                "database": [300, 400, 200],
                "textmining": [400, 250, 400],
            }
        )

        proteins_data = pd.DataFrame(
            {
                "protein_id": ["ENSP1", "ENSP2", "ENSP3", "ENSP4"],
                "gene_name": ["GENE1", "GENE2", "GENE3", "GENE4"],
                "protein_name": ["Protein 1", "Protein 2", "Protein 3", "Protein 4"],
            }
        )

        ppi_network = load_string_interactions(
            interactions_df=interactions_data, proteins_df=proteins_data, confidence_threshold=700
        )

        # Should load interactions above threshold
        assert isinstance(ppi_network, ProteinNetwork)
        # 900 and 750 are above 700, 600 is below
        assert len(ppi_network.interactions) == 2
        assert len(ppi_network.proteins) >= 3

        # Check protein metadata loaded
        assert "ENSP1" in ppi_network.protein_metadata
        metadata = ppi_network.protein_metadata["ENSP1"]
        assert metadata["gene_name"] == "GENE1"
        assert metadata["protein_name"] == "Protein 1"

    def test_load_string_confidence_filtering(self):
        """Test confidence threshold filtering during loading."""
        interactions_data = pd.DataFrame(
            {
                "protein1": ["A", "B", "C"],
                "protein2": ["B", "C", "D"],
                "combined_score": [950, 600, 200],  # High, medium, low
            }
        )

        # Load with high threshold
        high_conf_network = load_string_interactions(interactions_df=interactions_data, confidence_threshold=800)
        assert len(high_conf_network.interactions) == 1  # Only A-B (950)

        # Load with medium threshold
        med_conf_network = load_string_interactions(interactions_df=interactions_data, confidence_threshold=500)
        assert len(med_conf_network.interactions) == 2  # A-B and B-C

        # Load with low threshold
        low_conf_network = load_string_interactions(interactions_df=interactions_data, confidence_threshold=100)
        assert len(low_conf_network.interactions) == 3  # All interactions

    def test_load_string_empty_data(self):
        """Test loading empty STRING data."""
        empty_interactions = pd.DataFrame(columns=["protein1", "protein2", "combined_score"])

        ppi_network = load_string_interactions(interactions_df=empty_interactions)

        assert isinstance(ppi_network, ProteinNetwork)
        assert len(ppi_network.interactions) == 0
        assert len(ppi_network.proteins) == 0

    def test_load_string_invalid_data(self):
        """Test error handling for invalid STRING data."""
        invalid_data = pd.DataFrame(
            {
                "protein1": ["A"],
                # Missing protein2 and combined_score columns
            }
        )

        # Should handle gracefully or raise informative error
        try:
            ppi_network = load_string_interactions(interactions_df=invalid_data)
            assert isinstance(ppi_network, ProteinNetwork)
        except (ValueError, KeyError):
            # Acceptable to raise error for invalid format
            pass


class TestPredictInteractions:
    """Test protein interaction prediction."""

    def setup_method(self):
        """Set up test data for interaction prediction."""
        # Create known PPI network
        self.known_network = ProteinNetwork()

        # Add known interactions
        known_interactions = [
            ("P1", "P2", 0.9, ["experimental"]),
            ("P2", "P3", 0.8, ["experimental"]),
            ("P3", "P4", 0.85, ["experimental"]),
            ("P4", "P5", 0.7, ["database"]),
            ("P1", "P6", 0.6, ["textmining"]),
        ]

        for p1, p2, conf, evidence in known_interactions:
            self.known_network.add_interaction(p1, p2, conf, evidence)

        # Target proteins for prediction
        self.target_proteins = ["P7", "P8"]

    def test_predict_interactions_basic(self):
        """Test basic interaction prediction."""
        predictions = predict_interactions(
            target_proteins=self.target_proteins, known_network=self.known_network, method="similarity"
        )

        # Should return prediction results
        assert isinstance(predictions, dict)

        # Should have predictions for target proteins
        for protein in self.target_proteins:
            if protein in predictions:
                protein_predictions = predictions[protein]
                assert isinstance(protein_predictions, list)

                # Check prediction format
                for prediction in protein_predictions:
                    assert "partner" in prediction
                    assert "confidence" in prediction
                    assert 0.0 <= prediction["confidence"] <= 1.0

    def test_predict_interactions_different_methods(self):
        """Test interaction prediction with different methods."""
        methods = ["similarity", "coexpression", "domain"]

        for method in methods:
            try:
                predictions = predict_interactions(
                    target_proteins=["P7"], known_network=self.known_network, method=method
                )

                # All methods should return valid predictions
                assert isinstance(predictions, dict)

            except ValueError as e:
                if "Unknown method" in str(e):
                    # Some methods might not be implemented
                    continue
                else:
                    raise

    def test_predict_interactions_confidence_threshold(self):
        """Test prediction with confidence threshold."""
        # High threshold (strict predictions)
        strict_predictions = predict_interactions(
            target_proteins=["P7"], known_network=self.known_network, confidence_threshold=0.8
        )

        # Low threshold (permissive predictions)
        permissive_predictions = predict_interactions(
            target_proteins=["P7"], known_network=self.known_network, confidence_threshold=0.3
        )

        # Should return predictions (might be empty with high threshold)
        assert isinstance(strict_predictions, dict)
        assert isinstance(permissive_predictions, dict)

    def test_predict_interactions_empty_network(self):
        """Test prediction with empty known network."""
        empty_network = ProteinNetwork()

        predictions = predict_interactions(target_proteins=["P1", "P2"], known_network=empty_network)

        # Should handle empty network gracefully
        assert isinstance(predictions, dict)


class TestFunctionalEnrichmentPPI:
    """Test PPI-based functional enrichment analysis."""

    def setup_method(self):
        """Set up PPI network with functional annotations."""
        self.ppi_network = ProteinNetwork()

        # Add proteins with functional annotations
        proteins_functions = [
            ("A1", {"function": "DNA repair", "complex": "BRCA1_complex"}),
            ("A2", {"function": "DNA repair", "complex": "BRCA1_complex"}),
            ("A3", {"function": "DNA repair", "complex": "BRCA2_complex"}),
            ("B1", {"function": "cell cycle", "complex": "cyclin_complex"}),
            ("B2", {"function": "cell cycle", "complex": "cyclin_complex"}),
            ("C1", {"function": "metabolism", "complex": "enzyme_complex"}),
        ]

        for protein_id, metadata in proteins_functions:
            self.ppi_network.add_protein_metadata(protein_id, **metadata)

        # Add interactions (proteins with similar functions interact more)
        interactions = [
            # DNA repair proteins interact
            ("A1", "A2", 0.9, ["experimental"]),
            ("A1", "A3", 0.8, ["experimental"]),
            ("A2", "A3", 0.85, ["experimental"]),
            # Cell cycle proteins interact
            ("B1", "B2", 0.9, ["experimental"]),
            # Cross-functional interactions (weaker)
            ("A1", "B1", 0.4, ["database"]),
            ("B2", "C1", 0.3, ["textmining"]),
        ]

        for p1, p2, conf, evidence in interactions:
            self.ppi_network.add_interaction(p1, p2, conf, evidence)

    def test_functional_enrichment_basic(self):
        """Test basic functional enrichment analysis."""
        # Test set enriched for DNA repair proteins
        test_proteins = ["A1", "A2", "A3"]

        enrichment_results = functional_enrichment_ppi(
            protein_list=test_proteins, ppi_network=self.ppi_network, function_key="function"
        )

        # Should detect DNA repair enrichment
        assert isinstance(enrichment_results, dict)
        assert "DNA repair" in enrichment_results

        # Check enrichment details
        dna_repair_result = enrichment_results["DNA repair"]
        assert dna_repair_result["count"] == 3
        assert dna_repair_result["enrichment_ratio"] > 1.0
        assert dna_repair_result["p_value"] <= 0.05  # Should be significant

    def test_functional_enrichment_complex_analysis(self):
        """Test enrichment analysis for protein complexes."""
        test_proteins = ["A1", "A2", "B1", "B2"]  # Two complexes represented

        enrichment_results = functional_enrichment_ppi(
            protein_list=test_proteins, ppi_network=self.ppi_network, function_key="complex"
        )

        # Should detect enrichment for both complexes
        assert isinstance(enrichment_results, dict)
        assert "BRCA1_complex" in enrichment_results
        assert "cyclin_complex" in enrichment_results

        # Both should show enrichment
        for complex_name, result in enrichment_results.items():
            if complex_name in ["BRCA1_complex", "cyclin_complex"]:
                assert result["count"] >= 2
                assert result["enrichment_ratio"] >= 1.0

    def test_functional_enrichment_no_annotation(self):
        """Test enrichment when proteins lack functional annotations."""
        # Test set with proteins that have no function annotation
        test_proteins = ["UNANNOTATED1", "UNANNOTATED2"]

        # Add these proteins to network without functional metadata
        for protein in test_proteins:
            self.ppi_network.add_protein_metadata(protein)

        enrichment_results = functional_enrichment_ppi(
            protein_list=test_proteins, ppi_network=self.ppi_network, function_key="function"
        )

        # Should handle gracefully
        assert isinstance(enrichment_results, dict)
        # Might have empty results or "unknown" category

    def test_functional_enrichment_statistical_significance(self):
        """Test statistical significance calculation."""
        # Highly enriched set (all DNA repair)
        highly_enriched = ["A1", "A2", "A3"]

        # Moderately enriched set (mixed functions)
        moderately_enriched = ["A1", "B1", "C1"]

        high_results = functional_enrichment_ppi(
            protein_list=highly_enriched, ppi_network=self.ppi_network, function_key="function"
        )

        moderate_results = functional_enrichment_ppi(
            protein_list=moderately_enriched, ppi_network=self.ppi_network, function_key="function"
        )

        # Highly enriched should have lower p-values
        if "DNA repair" in high_results and "DNA repair" in moderate_results:
            assert high_results["DNA repair"]["p_value"] <= moderate_results["DNA repair"]["p_value"]


class TestPPINetworkIntegration:
    """Integration tests for PPI network analysis."""

    def test_complete_ppi_analysis_workflow(self):
        """Test complete PPI analysis workflow."""
        # 1. Create comprehensive PPI network
        ppi_network = ProteinNetwork("comprehensive_ppi")

        # Add proteins with functional annotations
        proteins = [
            ("BRCA1", {"function": "DNA repair", "length": 1863}),
            ("BRCA2", {"function": "DNA repair", "length": 3418}),
            ("TP53", {"function": "tumor suppressor", "length": 393}),
            ("MDM2", {"function": "ubiquitin ligase", "length": 491}),
            ("ATM", {"function": "kinase", "length": 3056}),
        ]

        for protein_id, metadata in proteins:
            ppi_network.add_protein_metadata(protein_id, **metadata)

        # Add biologically relevant interactions
        interactions = [
            ("BRCA1", "BRCA2", 0.8, ["experimental"]),
            ("TP53", "MDM2", 0.95, ["experimental"]),  # Well-known interaction
            ("ATM", "BRCA1", 0.7, ["experimental"]),
            ("ATM", "TP53", 0.85, ["experimental"]),
            ("BRCA1", "TP53", 0.6, ["database"]),
        ]

        for p1, p2, conf, evidence in interactions:
            ppi_network.add_interaction(p1, p2, conf, evidence)

        # 2. Network analysis
        stats = ppi_network.get_network_statistics()

        # 3. Functional enrichment
        dna_proteins = ["BRCA1", "BRCA2", "ATM"]
        enrichment_results = functional_enrichment_ppi(
            protein_list=dna_proteins, ppi_network=ppi_network, function_key="function"
        )

        # 4. Interaction prediction
        predictions = predict_interactions(target_proteins=["NEW_PROTEIN"], known_network=ppi_network)

        # 5. Convert to biological network for graph analysis
        bio_network = ppi_network.to_biological_network()

        # 6. Verify integrated results
        assert stats["num_proteins"] == 5
        assert stats["num_interactions"] == 5
        assert 0.0 < stats["density"] < 1.0

        # Should detect DNA repair enrichment
        if "DNA repair" in enrichment_results:
            assert enrichment_results["DNA repair"]["count"] >= 2

        # Predictions should be valid
        assert isinstance(predictions, dict)

        # Biological network should preserve interactions
        assert isinstance(bio_network, BiologicalNetwork)
        assert bio_network.num_nodes() == 5


class TestPPIEdgeCases:
    """Test edge cases and error conditions."""

    def test_empty_ppi_network(self):
        """Test operations on empty PPI network."""
        empty_network = ProteinNetwork()

        stats = empty_network.get_network_statistics()
        assert stats["num_proteins"] == 0
        assert stats["num_interactions"] == 0

        # Should handle conversion to biological network
        bio_network = empty_network.to_biological_network()
        assert bio_network.num_nodes() == 0

    def test_duplicate_interactions(self):
        """Test handling of duplicate interactions."""
        network = ProteinNetwork()

        # Add same interaction twice with different confidence
        network.add_interaction("A", "B", 0.7, ["experimental"])
        network.add_interaction("A", "B", 0.9, ["computational"])  # Should add as separate

        # Both interactions should be stored
        assert len(network.interactions) == 2

    def test_self_interactions(self):
        """Test handling of protein self-interactions."""
        network = ProteinNetwork()

        # Add self-interaction (homodimer)
        network.add_interaction("SELF_PROTEIN", "SELF_PROTEIN", 0.8, ["experimental"])

        assert len(network.interactions) == 1
        assert "SELF_PROTEIN" in network.proteins

        # Should appear as partner
        partners = network.get_protein_partners("SELF_PROTEIN")
        assert "SELF_PROTEIN" in partners

    def test_invalid_confidence_values(self):
        """Test handling of invalid confidence values."""
        network = ProteinNetwork()

        # Negative confidence
        network.add_interaction("A", "B", -0.5, ["experimental"])

        # Very high confidence (>1.0)
        network.add_interaction("C", "D", 2.0, ["experimental"])

        # Should store as-is (validation might be in application layer)
        assert len(network.interactions) == 2

    @pytest.mark.slow
    def test_large_ppi_network(self):
        """Test performance with large PPI network.
        
        This test creates a moderately-sized network to verify performance.
        Marked as slow due to computational complexity of network statistics.
        """
        large_network = ProteinNetwork()

        # Add many proteins and interactions (reduced from 100/200 to 50/100 for faster execution)
        proteins = [f"PROTEIN_{i}" for i in range(50)]

        # Add random interactions
        np.random.seed(42)
        for _ in range(100):
            p1 = np.random.choice(proteins)
            p2 = np.random.choice(proteins)
            if p1 != p2:
                conf = np.random.random()
                large_network.add_interaction(p1, p2, conf, ["computational"])

        # Should handle large network - use get_network_statistics for basic stats
        stats = large_network.get_network_statistics()
        assert stats["num_proteins"] <= 50
        assert stats["num_interactions"] <= 100
        assert stats["num_proteins"] > 0
        assert stats["num_interactions"] > 0


class TestNewPPIMethods:
    """Test new PPI methods and functions."""

    def test_get_protein_partners(self):
        """Test get_protein_partners method."""
        ppi = ProteinNetwork()
        ppi.add_interaction("P1", "P2", confidence=0.9)
        ppi.add_interaction("P1", "P3", confidence=0.5)
        ppi.add_interaction("P1", "P4", confidence=0.8)

        partners = ppi.get_protein_partners("P1", min_confidence=0.7)
        assert len(partners) == 2  # P2 and P4 above threshold
        assert "P2" in partners
        assert "P4" in partners

    def test_protein_similarity(self):
        """Test protein similarity calculation."""
        from metainformant.networks.ppi import protein_similarity

        ppi = ProteinNetwork()
        ppi.add_interaction("P1", "P2")
        ppi.add_interaction("P1", "P3")
        ppi.add_interaction("P2", "P3")  # Common neighbor

        sim = protein_similarity(ppi, "P1", "P2")
        assert "jaccard_similarity" in sim
        assert "common_neighbors" in sim
        assert sim["common_neighbors"] >= 0

    def test_detect_complexes(self):
        """Test protein complex detection."""
        from metainformant.networks.ppi import detect_complexes

        ppi = ProteinNetwork()
        # Create dense subgraph (complex)
        for i in range(3):
            for j in range(i + 1, 3):
                ppi.add_interaction(f"P{i}", f"P{j}", confidence=0.8)

        complexes = detect_complexes(ppi, min_confidence=0.7, min_size=3)
        assert len(complexes) > 0

    def test_export_to_string_format(self, tmp_path):
        """Test STRING format export."""
        from metainformant.networks.ppi import export_to_string_format

        ppi = ProteinNetwork()
        ppi.add_interaction("P1", "P2", confidence=0.8)
        ppi.add_interaction("P3", "P4", confidence=0.3)

        output_file = tmp_path / "string_output.tsv"
        export_to_string_format(ppi, str(output_file), score_threshold=400)

        # Check file was created
        assert output_file.exists()
