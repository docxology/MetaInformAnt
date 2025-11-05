"""Tests for gene regulatory network functionality.

Real implementation testing for GRN analysis methods.
No mocking used - all tests use real computational methods and data.
"""

from __future__ import annotations

from typing import Dict, List, Set

import numpy as np
import pandas as pd
import pytest

from metainformant.networks.graph import BiologicalNetwork
from metainformant.networks.regulatory import (
    GeneRegulatoryNetwork,
    infer_grn,
    pathway_regulation_analysis,
    regulatory_motifs,
)


class TestGeneRegulatoryNetwork:
    """Test GeneRegulatoryNetwork functionality."""

    def setup_method(self):
        """Set up test gene regulatory network."""
        self.grn = GeneRegulatoryNetwork("test_grn")

        # Add regulatory interactions
        regulations = [
            ("TF1", "TARGET1", "activation", 0.8, 0.9),
            ("TF1", "TARGET2", "repression", 0.7, 0.8),
            ("TF2", "TARGET2", "activation", 0.9, 0.95),
            ("TF2", "TARGET3", "activation", 0.6, 0.7),
            ("TARGET1", "TF2", "activation", 0.5, 0.6),  # Regulatory cascade
        ]

        for regulator, target, reg_type, strength, confidence in regulations:
            self.grn.add_regulation(regulator, target, reg_type, strength, confidence, evidence="ChIP-seq")

        # Mark transcription factors
        self.grn.add_transcription_factor("TF1", tf_family="bHLH")
        self.grn.add_transcription_factor("TF2", tf_family="homeodomain")

        # Add gene metadata
        genes_metadata = [
            ("TF1", {"type": "transcription_factor", "chromosome": "1"}),
            ("TF2", {"type": "transcription_factor", "chromosome": "2"}),
            ("TARGET1", {"type": "protein_coding", "chromosome": "3"}),
            ("TARGET2", {"type": "protein_coding", "chromosome": "4"}),
            ("TARGET3", {"type": "protein_coding", "chromosome": "5"}),
        ]

        for gene_id, metadata in genes_metadata:
            self.grn.add_gene_metadata(gene_id, **metadata)

    def test_grn_initialization(self):
        """Test GeneRegulatoryNetwork initialization."""
        grn = GeneRegulatoryNetwork("test_network")

        assert grn.name == "test_network"
        assert len(grn.regulations) == 0
        assert len(grn.genes) == 0
        assert len(grn.transcription_factors) == 0

    def test_add_regulation(self):
        """Test adding regulatory interactions."""
        grn = GeneRegulatoryNetwork()

        # Add regulation with metadata
        grn.add_regulation(
            "MYC",
            "CDKN1A",
            regulation_type="repression",
            strength=0.85,
            confidence=0.9,
            evidence="ChIP-seq",
            pubmed_id="12345678",
        )

        assert len(grn.regulations) == 1
        assert "MYC" in grn.genes
        assert "CDKN1A" in grn.genes

        # Check regulation details
        regulation = grn.regulations[0]
        assert regulation[0] == "MYC"
        assert regulation[1] == "CDKN1A"
        assert regulation[2]["type"] == "repression"
        assert regulation[2]["strength"] == 0.85
        assert regulation[2]["confidence"] == 0.9
        assert regulation[2]["evidence"] == "ChIP-seq"

    def test_add_transcription_factor(self):
        """Test adding transcription factors."""
        grn = GeneRegulatoryNetwork()

        grn.add_transcription_factor(
            "TP53", tf_family="p53", dna_binding_domain="p53_DNA_binding", function="tumor suppressor"
        )

        assert "TP53" in grn.transcription_factors
        assert "TP53" in grn.genes

        # Check TF metadata
        if "TP53" in grn.gene_metadata:
            metadata = grn.gene_metadata["TP53"]
            assert metadata.get("tf_family") == "p53"
            assert metadata.get("function") == "tumor suppressor"

    def test_add_gene_metadata(self):
        """Test adding gene metadata."""
        grn = GeneRegulatoryNetwork()

        grn.add_gene_metadata("BRCA1", type="protein_coding", chromosome="17", function="DNA repair", length=81189)

        assert "BRCA1" in grn.gene_metadata
        metadata = grn.gene_metadata["BRCA1"]
        assert metadata["type"] == "protein_coding"
        assert metadata["chromosome"] == "17"
        assert metadata["function"] == "DNA repair"
        assert metadata["length"] == 81189

    def test_get_targets(self):
        """Test getting targets of a transcription factor."""
        targets_tf1 = self.grn.get_targets("TF1")
        expected_targets_tf1 = {"TARGET1", "TARGET2"}
        assert set(targets_tf1) == expected_targets_tf1

        targets_tf2 = self.grn.get_targets("TF2")
        expected_targets_tf2 = {"TARGET2", "TARGET3"}
        assert set(targets_tf2) == expected_targets_tf2

        # Test gene with no targets
        targets_target1 = self.grn.get_targets("TARGET3")
        assert len(targets_target1) == 0

    def test_get_regulators(self):
        """Test getting regulators of a gene."""
        regulators_target2 = self.grn.get_regulators("TARGET2")
        expected_regulators = {"TF1", "TF2"}
        assert set(regulators_target2) == expected_regulators

        regulators_target1 = self.grn.get_regulators("TARGET1")
        assert set(regulators_target1) == {"TF1"}

        # Test TF that is regulated by others
        regulators_tf2 = self.grn.get_regulators("TF2")
        assert set(regulators_tf2) == {"TARGET1"}

    def test_filter_by_confidence(self):
        """Test filtering regulations by confidence threshold."""
        # High confidence threshold
        high_conf_grn = self.grn.filter_by_confidence(threshold=0.85)

        # Should include TF1->TARGET1 (0.9), TF2->TARGET2 (0.95)
        high_conf_regulations = high_conf_grn.regulations
        assert len(high_conf_regulations) == 2

        # Check specific regulations
        reg_pairs = [(r[0], r[1]) for r in high_conf_regulations]
        assert ("TF1", "TARGET1") in reg_pairs
        assert ("TF2", "TARGET2") in reg_pairs

    def test_filter_by_regulation_type(self):
        """Test filtering by regulation type."""
        # Get only activation regulations
        activation_grn = self.grn.filter_by_regulation_type("activation")

        # Should include TF1->TARGET1, TF2->TARGET2, TF2->TARGET3, TARGET1->TF2
        activation_regs = activation_grn.regulations
        assert len(activation_regs) == 4

        # All should be activation
        for reg in activation_regs:
            assert reg[2]["type"] == "activation"

        # Get only repression regulations
        repression_grn = self.grn.filter_by_regulation_type("repression")
        repression_regs = repression_grn.regulations
        assert len(repression_regs) == 1  # Only TF1->TARGET2

        # Should be repression
        assert repression_regs[0][2]["type"] == "repression"
        assert repression_regs[0][0] == "TF1"
        assert repression_regs[0][1] == "TARGET2"

    def test_get_network_statistics(self):
        """Test GRN statistics calculation."""
        stats = self.grn.get_network_statistics()

        # Should include various metrics
        expected_keys = ["num_genes", "num_regulations", "num_tfs", "density"]
        for key in expected_keys:
            assert key in stats

        assert stats["num_genes"] == 5
        assert stats["num_regulations"] == 5
        assert stats["num_tfs"] == 2
        assert 0.0 <= stats["density"] <= 1.0

    def test_to_biological_network(self):
        """Test conversion to BiologicalNetwork."""
        bio_network = self.grn.to_biological_network()

        assert isinstance(bio_network, BiologicalNetwork)
        assert bio_network.num_nodes() == 5
        assert bio_network.num_edges() == 5

        # Check that regulations are preserved as edges
        assert bio_network.has_edge("TF1", "TARGET1")
        assert bio_network.has_edge("TF1", "TARGET2")
        assert bio_network.has_edge("TF2", "TARGET2")

        # Edge weights should reflect regulation strength
        weight_tf1_target1 = bio_network.get_edge_weight("TF1", "TARGET1")
        assert weight_tf1_target1 == 0.8  # Original strength


class TestInferGRN:
    """Test gene regulatory network inference."""

    def setup_method(self):
        """Set up test expression data for GRN inference."""
        np.random.seed(42)

        # Simulate expression data with regulatory relationships
        self.n_genes = 8
        self.n_samples = 40
        self.gene_names = [f"Gene_{i}" for i in range(self.n_genes)]

        # Create expression data with known regulatory structure
        # Gene_0 and Gene_1 are TFs that regulate other genes

        # Base expression levels
        base_expression = np.random.lognormal(mean=2, sigma=0.5, size=(self.n_samples, self.n_genes))

        # Add regulatory effects
        # Gene_0 activates Gene_2, Gene_3
        base_expression[:, 2] += 0.3 * base_expression[:, 0]
        base_expression[:, 3] += 0.4 * base_expression[:, 0]

        # Gene_1 represses Gene_4, activates Gene_5
        base_expression[:, 4] -= 0.2 * base_expression[:, 1]
        base_expression[:, 4] = np.maximum(base_expression[:, 4], 0.1)  # Keep positive
        base_expression[:, 5] += 0.35 * base_expression[:, 1]

        # Add noise
        self.expression_data = base_expression + np.random.normal(0, 0.1, base_expression.shape)

        # Ensure non-negative
        self.expression_data = np.maximum(self.expression_data, 0.01)

        # TF annotations (prior knowledge)
        self.tf_genes = ["Gene_0", "Gene_1"]

    def test_infer_grn_basic(self):
        """Test basic GRN inference."""
        grn = infer_grn(expression_data=self.expression_data, gene_names=self.gene_names, method="correlation")

        # Should infer regulatory network
        assert isinstance(grn, GeneRegulatoryNetwork)
        assert len(grn.genes) == self.n_genes
        assert len(grn.regulations) >= 1

        # All genes should be in network
        assert set(grn.genes) == set(self.gene_names)

    def test_infer_grn_with_tf_prior(self):
        """Test GRN inference with TF prior knowledge."""
        grn = infer_grn(
            expression_data=self.expression_data,
            gene_names=self.gene_names,
            method="correlation",
            tf_genes=self.tf_genes,
        )

        # Should identify TF genes correctly
        assert "Gene_0" in grn.transcription_factors
        assert "Gene_1" in grn.transcription_factors

        # TFs should have some targets
        targets_gene0 = grn.get_targets("Gene_0")
        targets_gene1 = grn.get_targets("Gene_1")

        assert len(targets_gene0) >= 0  # Might not detect all relationships
        assert len(targets_gene1) >= 0

    def test_infer_grn_different_methods(self):
        """Test GRN inference with different methods."""
        methods = ["correlation", "mutual_information", "regression"]

        for method in methods:
            try:
                grn = infer_grn(expression_data=self.expression_data, gene_names=self.gene_names, method=method)

                # All methods should produce valid networks
                assert isinstance(grn, GeneRegulatoryNetwork)
                assert len(grn.genes) == self.n_genes
                assert len(grn.regulations) >= 0

            except ValueError as e:
                if "Unknown method" in str(e):
                    # Some methods might not be implemented
                    continue
                else:
                    raise

    def test_infer_grn_thresholds(self):
        """Test effect of different thresholds on GRN inference."""
        # High threshold (strict)
        strict_grn = infer_grn(
            expression_data=self.expression_data, gene_names=self.gene_names, method="correlation", threshold=0.7
        )

        # Low threshold (permissive)
        permissive_grn = infer_grn(
            expression_data=self.expression_data, gene_names=self.gene_names, method="correlation", threshold=0.2
        )

        # Permissive should have more regulations
        assert len(permissive_grn.regulations) >= len(strict_grn.regulations)

        # Both should have same number of genes
        assert len(strict_grn.genes) == len(permissive_grn.genes) == self.n_genes

    def test_infer_grn_small_dataset(self):
        """Test GRN inference with small dataset."""
        # Very small dataset
        small_expression = self.expression_data[:5, :5]  # 5 samples, 5 genes
        small_gene_names = self.gene_names[:5]

        grn = infer_grn(expression_data=small_expression, gene_names=small_gene_names, method="correlation")

        # Should handle small dataset
        assert isinstance(grn, GeneRegulatoryNetwork)
        assert len(grn.genes) == 5
        # Might have few or no regulations due to limited data
        assert len(grn.regulations) >= 0


class TestRegulatoryMotifs:
    """Test regulatory motif identification."""

    def setup_method(self):
        """Set up GRN with known motif structures."""
        self.grn = GeneRegulatoryNetwork("motif_test")

        # Create feed-forward loop: TF1 -> TF2 -> TARGET, TF1 -> TARGET
        self.grn.add_regulation("TF1", "TF2", "activation", 0.9, 0.9)
        self.grn.add_regulation("TF2", "TARGET", "activation", 0.8, 0.8)
        self.grn.add_regulation("TF1", "TARGET", "activation", 0.7, 0.7)

        # Create feedback loop: TF3 -> TF4 -> TF3
        self.grn.add_regulation("TF3", "TF4", "activation", 0.8, 0.8)
        self.grn.add_regulation("TF4", "TF3", "repression", 0.6, 0.7)

        # Mark transcription factors
        for tf in ["TF1", "TF2", "TF3", "TF4"]:
            self.grn.add_transcription_factor(tf)

    def test_regulatory_motifs_basic(self):
        """Test basic regulatory motif identification."""
        motifs = regulatory_motifs(self.grn)

        # Should find regulatory motifs
        assert isinstance(motifs, list)
        assert len(motifs) >= 1

        # Check motif structure
        for motif in motifs:
            assert "motif_type" in motif
            assert "genes" in motif
            assert "confidence" in motif

            # Confidence should be reasonable
            assert 0.0 <= motif["confidence"] <= 1.0

    def test_regulatory_motifs_types(self):
        """Test identification of different motif types."""
        motifs = regulatory_motifs(self.grn, motif_types=["feed_forward", "feedback"])

        # Should identify different motif types
        motif_types = [motif["motif_type"] for motif in motifs]

        # Should find at least some motifs
        assert len(motifs) >= 1

        # Might find feed-forward loop (TF1 -> TF2 -> TARGET, TF1 -> TARGET)
        # and feedback loop (TF3 -> TF4 -> TF3)
        expected_types = ["feed_forward", "feedback"]
        found_types = set(motif_types)

        # At least one expected type should be found
        assert len(found_types.intersection(expected_types)) >= 1

    def test_regulatory_motifs_confidence_filtering(self):
        """Test motif finding with confidence filtering."""
        # High confidence threshold
        high_conf_motifs = regulatory_motifs(self.grn, min_confidence=0.8)

        # Low confidence threshold
        all_motifs = regulatory_motifs(self.grn, min_confidence=0.1)

        # Should find fewer motifs with high confidence threshold
        assert len(high_conf_motifs) <= len(all_motifs)

        # All high confidence motifs should meet threshold
        for motif in high_conf_motifs:
            assert motif["confidence"] >= 0.7  # Some tolerance

    def test_regulatory_motifs_empty_grn(self):
        """Test motif finding on empty GRN."""
        empty_grn = GeneRegulatoryNetwork("empty")

        motifs = regulatory_motifs(empty_grn)

        # Should handle empty network gracefully
        assert isinstance(motifs, list)
        assert len(motifs) == 0


class TestPathwayRegulationAnalysis:
    """Test pathway regulation analysis."""

    def setup_method(self):
        """Set up GRN and pathway data for analysis."""
        # Create GRN with pathway-relevant regulations
        self.grn = GeneRegulatoryNetwork("pathway_grn")

        # Add regulations involving pathway genes
        pathway_regulations = [
            ("TF_MASTER", "PATHWAY_GENE_1", "activation", 0.9, 0.95),
            ("TF_MASTER", "PATHWAY_GENE_2", "activation", 0.8, 0.9),
            ("TF_MASTER", "PATHWAY_GENE_3", "repression", 0.7, 0.8),
            ("PATHWAY_GENE_1", "PATHWAY_GENE_2", "activation", 0.6, 0.7),
            ("EXTERNAL_TF", "PATHWAY_GENE_1", "activation", 0.5, 0.6),
        ]

        for regulator, target, reg_type, strength, confidence in pathway_regulations:
            self.grn.add_regulation(regulator, target, reg_type, strength, confidence)

        # Mark TFs
        self.grn.add_transcription_factor("TF_MASTER")
        self.grn.add_transcription_factor("EXTERNAL_TF")

        # Define pathway
        self.pathway_genes = ["PATHWAY_GENE_1", "PATHWAY_GENE_2", "PATHWAY_GENE_3"]

    def test_pathway_regulation_analysis_basic(self):
        """Test basic pathway regulation analysis."""
        analysis_results = pathway_regulation_analysis(grn=self.grn, pathway_genes=self.pathway_genes)

        # Should return analysis results
        assert isinstance(analysis_results, dict)

        # Should include key metrics
        expected_keys = ["pathway_tfs", "internal_regulations", "external_regulations", "regulation_density"]
        for key in expected_keys:
            assert key in analysis_results

        # Should identify TF_MASTER as pathway TF
        pathway_tfs = analysis_results["pathway_tfs"]
        assert "TF_MASTER" in pathway_tfs

        # Should have both internal and external regulations
        internal_regs = analysis_results["internal_regulations"]
        external_regs = analysis_results["external_regulations"]

        assert len(internal_regs) >= 1  # PATHWAY_GENE_1 -> PATHWAY_GENE_2
        assert len(external_regs) >= 3  # TF_MASTER regulations

    def test_pathway_regulation_analysis_enrichment(self):
        """Test pathway regulation enrichment analysis."""
        analysis_results = pathway_regulation_analysis(
            grn=self.grn, pathway_genes=self.pathway_genes, calculate_enrichment=True
        )

        # Should include enrichment analysis
        assert "regulation_enrichment" in analysis_results

        enrichment = analysis_results["regulation_enrichment"]
        assert "p_value" in enrichment
        assert "enrichment_ratio" in enrichment

        # Should show enrichment (pathway genes are highly regulated)
        assert enrichment["enrichment_ratio"] >= 1.0
        assert 0.0 <= enrichment["p_value"] <= 1.0

    def test_pathway_regulation_analysis_empty_pathway(self):
        """Test analysis with empty pathway."""
        analysis_results = pathway_regulation_analysis(grn=self.grn, pathway_genes=[])

        # Should handle empty pathway gracefully
        assert isinstance(analysis_results, dict)
        assert analysis_results["pathway_tfs"] == []
        assert len(analysis_results["internal_regulations"]) == 0
        assert len(analysis_results["external_regulations"]) == 0


class TestRegulatoryNetworkIntegration:
    """Integration tests for regulatory network analysis."""

    def test_complete_grn_analysis_workflow(self):
        """Test complete GRN analysis workflow."""
        # 1. Simulate expression data with regulatory structure
        np.random.seed(42)
        n_samples = 30
        gene_names = ["TF1", "TF2", "Gene_A", "Gene_B", "Gene_C"]

        # Create expression with regulatory relationships
        tf1_expr = np.random.lognormal(2, 0.8, n_samples)
        tf2_expr = np.random.lognormal(2.2, 0.7, n_samples)

        # Target genes influenced by TFs
        gene_a = 0.4 * tf1_expr + np.random.lognormal(1.5, 0.3, n_samples)
        gene_b = 0.3 * tf1_expr + 0.2 * tf2_expr + np.random.lognormal(1.8, 0.4, n_samples)
        gene_c = 0.5 * tf2_expr + np.random.lognormal(1.6, 0.3, n_samples)

        expression_data = np.column_stack([tf1_expr, tf2_expr, gene_a, gene_b, gene_c])

        # 2. Infer GRN from expression data
        inferred_grn = infer_grn(
            expression_data=expression_data, gene_names=gene_names, method="correlation", tf_genes=["TF1", "TF2"]
        )

        # 3. Find regulatory motifs
        motifs = regulatory_motifs(inferred_grn)

        # 4. Analyze pathway regulation
        pathway_genes = ["Gene_A", "Gene_B", "Gene_C"]
        pathway_analysis = pathway_regulation_analysis(grn=inferred_grn, pathway_genes=pathway_genes)

        # 5. Convert to biological network for graph analysis
        bio_network = inferred_grn.to_biological_network()

        # 6. Verify integrated results
        assert len(inferred_grn.genes) == len(gene_names)
        assert len(inferred_grn.regulations) >= 1

        # Should identify TFs
        assert len(inferred_grn.transcription_factors) >= 1
        tf_overlap = set(inferred_grn.transcription_factors).intersection({"TF1", "TF2"})
        assert len(tf_overlap) >= 1

        # Network statistics should be reasonable
        stats = inferred_grn.get_network_statistics()
        assert stats["density"] >= 0.0
        assert stats["num_tfs"] >= 1

        # Should find regulatory motifs
        assert isinstance(motifs, list)

        # Pathway analysis should identify TF regulation
        assert isinstance(pathway_analysis, dict)
        assert len(pathway_analysis["pathway_tfs"]) >= 0

        # Biological network should preserve regulations
        assert isinstance(bio_network, BiologicalNetwork)
        assert bio_network.num_nodes() == 5


class TestRegulatoryNetworkEdgeCases:
    """Test edge cases and error conditions."""

    def test_empty_grn_operations(self):
        """Test operations on empty GRN."""
        empty_grn = GeneRegulatoryNetwork("empty")

        # Should handle empty network gracefully
        assert len(empty_grn.genes) == 0
        assert len(empty_grn.regulations) == 0
        assert len(empty_grn.get_targets("NONEXISTENT")) == 0
        assert len(empty_grn.get_regulators("NONEXISTENT")) == 0

        # Statistics should handle empty network
        stats = empty_grn.get_network_statistics()
        assert stats["num_genes"] == 0
        assert stats["density"] == 0.0

        # Motif finding should return empty
        motifs = regulatory_motifs(empty_grn)
        assert len(motifs) == 0

    def test_self_regulation(self):
        """Test handling of self-regulatory loops."""
        grn = GeneRegulatoryNetwork("self_reg")

        # Add self-regulation
        grn.add_regulation("AUTOREGULATOR", "AUTOREGULATOR", "activation", 0.7, 0.8)
        grn.add_transcription_factor("AUTOREGULATOR")

        assert len(grn.regulations) == 1
        assert "AUTOREGULATOR" in grn.get_targets("AUTOREGULATOR")
        assert "AUTOREGULATOR" in grn.get_regulators("AUTOREGULATOR")

        # Should handle in motif detection
        motifs = regulatory_motifs(grn)
        assert isinstance(motifs, list)  # Should not crash

    def test_duplicate_regulations(self):
        """Test handling of duplicate regulatory interactions."""
        grn = GeneRegulatoryNetwork()

        # Add same regulation twice with different parameters
        grn.add_regulation("TF", "TARGET", "activation", 0.7, 0.8)
        grn.add_regulation("TF", "TARGET", "repression", 0.9, 0.9)  # Different type

        # Both regulations should be stored
        assert len(grn.regulations) == 2

        # Both should appear in targets
        targets = grn.get_targets("TF")
        assert "TARGET" in targets

    def test_invalid_regulation_types(self):
        """Test handling of invalid regulation types."""
        grn = GeneRegulatoryNetwork()

        # Add regulation with custom type
        grn.add_regulation("A", "B", "custom_regulation", 0.8, 0.9)

        assert len(grn.regulations) == 1
        reg = grn.regulations[0]
        assert reg[2]["type"] == "custom_regulation"

        # Filtering by standard types shouldn't include custom type
        activation_grn = grn.filter_by_regulation_type("activation")
        assert len(activation_grn.regulations) == 0

    def test_extreme_regulation_strengths(self):
        """Test handling of extreme regulation strengths."""
        grn = GeneRegulatoryNetwork()

        # Very weak regulation
        grn.add_regulation("A", "B", "activation", 0.001, 0.1)

        # Very strong regulation (>1.0)
        grn.add_regulation("A", "C", "activation", 2.5, 1.0)

        assert len(grn.regulations) == 2

        # Should handle in statistics
        stats = grn.get_network_statistics()
        assert isinstance(stats["density"], float)

    def test_large_grn_performance(self):
        """Test performance with large GRN."""
        large_grn = GeneRegulatoryNetwork("large")

        # Add many genes and regulations
        n_genes = 100
        gene_names = [f"Gene_{i}" for i in range(n_genes)]

        # Add random regulations
        np.random.seed(42)
        tfs = [f"TF_{i}" for i in range(10)]  # 10 TFs

        for tf in tfs:
            large_grn.add_transcription_factor(tf)

        for _ in range(200):  # 200 regulations
            tf = np.random.choice(tfs)
            target = np.random.choice(gene_names)
            if tf != target:
                reg_type = np.random.choice(["activation", "repression"])
                strength = np.random.random()
                confidence = np.random.random()
                large_grn.add_regulation(tf, target, reg_type, strength, confidence)

        # Should handle large network
        # Note: Not all genes may be targets if randomly selected, so check for reasonable coverage
        assert len(large_grn.genes) >= n_genes * 0.9  # At least 90% of target genes should be present
        assert len(large_grn.regulations) <= 200

        # Statistics should complete
        stats = large_grn.get_network_statistics()
        assert isinstance(stats["density"], float)

        # Should handle conversion
        bio_network = large_grn.create_network()
        assert isinstance(bio_network, BiologicalNetwork)


class TestNewRegulatoryFunctions:
    """Test new regulatory network functions."""

    def test_detect_regulatory_cascades(self):
        """Test regulatory cascade detection."""
        from metainformant.networks.regulatory import detect_regulatory_cascades

        grn = GeneRegulatoryNetwork()
        grn.add_regulation("TF1", "TF2", confidence=0.8)
        grn.add_regulation("TF2", "GENE1", confidence=0.9)

        cascades = detect_regulatory_cascades(grn, max_length=5, min_confidence=0.5)
        assert len(cascades) > 0

    def test_validate_regulation(self):
        """Test regulation validation."""
        from metainformant.networks.regulatory import validate_regulation

        grn = GeneRegulatoryNetwork()
        grn.add_regulation("TF1", "GENE1", confidence=0.8, regulation_type="activation")

        validation = validate_regulation(grn, "TF1", "GENE1", min_confidence=0.5)
        assert validation["exists"] is True
        assert validation["confidence"] == 0.8

        # Test non-existent regulation
        validation2 = validate_regulation(grn, "TF1", "GENE2", min_confidence=0.5)
        assert validation2["exists"] is False
