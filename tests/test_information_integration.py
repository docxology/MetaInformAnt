"""Integration tests for information module with other modules.

Tests verify that information theory methods work correctly when
integrated with DNA, RNA, single-cell, multi-omics, and ML modules.
"""

from __future__ import annotations

import numpy as np
import pytest

from metainformant.information import (
    information_profile,
    mutual_information,
    shannon_entropy_from_counts,
)
from metainformant.information.integration import (
    dna_integration,
    ml_integration,
    multiomics_integration,
    rna_integration,
    singlecell_integration,
)


class TestDNAIntegration:
    """Integration tests with DNA module."""

    def test_dna_sequence_entropy(self):
        """Test entropy calculation with DNA sequences."""
        from metainformant.dna import sequences

        # Create test sequences
        test_seqs = {
            "seq1": "ATCGATCG",
            "seq2": "AAAAA",
            "seq3": "ATCGATCG",
        }

        # Use information module functions
        seq_list = list(test_seqs.values())
        profile = information_profile(seq_list, k=1)

        assert profile["entropy"] > 0.0
        assert "A" in profile["kmer_frequencies"]

    def test_dna_integration_function(self):
        """Test DNA integration wrapper function."""
        sequences = {"seq1": "ATCG", "seq2": "AAAA"}
        results = dna_integration(sequences, k=1, analysis_type="entropy")
        assert "entropy_analysis" in results
        assert len(results["entropy_analysis"]) == 2

    def test_dna_integration_profile(self):
        """Test DNA integration with profile analysis."""
        sequences = {"seq1": "ATCGATCG", "seq2": "AAAA"}
        results = dna_integration(sequences, k=2, analysis_type="profile")
        assert "profile" in results


class TestRNAIntegration:
    """Integration tests with RNA module."""

    def test_rna_expression_entropy(self):
        """Test entropy calculation with RNA expression data."""
        import numpy as np

        # Simulate expression data
        expression = np.random.randn(100, 50)  # 100 samples, 50 genes

        # Calculate entropy for each gene
        gene_entropies = []
        for i in range(expression.shape[1]):
            gene_expr = expression[:, i]
            from collections import Counter

            counts = Counter(gene_expr)
            entropy = shannon_entropy_from_counts(counts)
            gene_entropies.append(entropy)

        assert len(gene_entropies) == 50
        assert all(e >= 0.0 for e in gene_entropies)

    def test_rna_integration_function(self):
        """Test RNA integration wrapper function."""
        expression = np.random.randn(100, 50)
        results = rna_integration(expression, method="entropy")
        assert "gene_entropies" in results
        assert len(results["gene_entropies"]) == 50

    def test_rna_mutual_information(self):
        """Test MI calculation between genes."""
        expression = np.random.randn(100, 50)
        results = rna_integration(expression, method="mutual_information")
        assert "mi_matrix" in results
        assert len(results["mi_matrix"]) == 50


class TestSingleCellIntegration:
    """Integration tests with single-cell module."""

    def test_singlecell_integration(self):
        """Test single-cell integration function."""
        count_matrix = np.random.randint(0, 100, (100, 50))  # 100 cells, 50 genes
        cell_types = ["TypeA"] * 50 + ["TypeB"] * 50

        results = singlecell_integration(
            count_matrix, cell_types=cell_types, method="cell_type_entropy"
        )
        assert "cell_type_entropy" in results
        assert results["num_cell_types"] == 2

    def test_singlecell_gene_entropy(self):
        """Test gene entropy calculation in single-cell data."""
        count_matrix = np.random.randint(0, 100, (100, 50))
        results = singlecell_integration(count_matrix, method="gene_entropy")
        assert "gene_entropies" in results
        assert len(results["gene_entropies"]) == 50


class TestMultiOmicsIntegration:
    """Integration tests with multi-omics module."""

    def test_multiomics_integration(self):
        """Test multi-omics integration function."""
        genomics = np.random.randn(100, 50)
        transcriptomics = np.random.randn(100, 50)
        proteomics = np.random.randn(100, 50)

        results = multiomics_integration(
            genomics_data=genomics,
            transcriptomics_data=transcriptomics,
            proteomics_data=proteomics,
            method="platform_entropy",
        )

        assert "genomics_entropy" in results
        assert "transcriptomics_entropy" in results
        assert "proteomics_entropy" in results

    def test_multiomics_cross_platform_mi(self):
        """Test cross-platform mutual information."""
        genomics = np.random.randn(100, 50)
        transcriptomics = np.random.randn(100, 50)

        results = multiomics_integration(
            genomics_data=genomics,
            transcriptomics_data=transcriptomics,
            method="cross_platform_mi",
        )

        # Results always in matrix format
        assert "method" in results
        assert "genomics_transcriptomics_mi_matrix" in results
        assert "genomics_transcriptomics_mean_mi" in results
        assert "genomics_transcriptomics_max_mi" in results
        assert "genomics_transcriptomics_min_mi" in results
        assert len(results["genomics_transcriptomics_mi_matrix"]) == 1
        assert len(results["genomics_transcriptomics_mi_matrix"][0]) == 1

    def test_multiomics_feature_selection(self):
        """Test multi-omics integration with feature selection."""
        genomics = np.random.randn(100, 50)
        transcriptomics = np.random.randn(100, 50)

        # Test with specific feature indices
        results = multiomics_integration(
            genomics_data=genomics,
            transcriptomics_data=transcriptomics,
            method="cross_platform_mi",
            feature_indices={"genomics": [0, 1, 2], "transcriptomics": [0, 3, 5]},
        )

        assert "method" in results
        # Should have matrix for multiple features
        assert "genomics_transcriptomics_mi_matrix" in results
        assert len(results["genomics_transcriptomics_mi_matrix"]) == 3
        assert len(results["genomics_transcriptomics_mi_matrix"][0]) == 3
        assert "genomics_transcriptomics_mean_mi" in results
        assert "genomics_transcriptomics_max_mi" in results
        assert "genomics_transcriptomics_min_mi" in results
        assert "genomics_feature_indices" in results
        assert "transcriptomics_feature_indices" in results

    def test_multiomics_feature_selection_single(self):
        """Test multi-omics integration with single feature indices."""
        genomics = np.random.randn(100, 50)
        transcriptomics = np.random.randn(100, 50)

        # Test with single integer (should be converted to list)
        results = multiomics_integration(
            genomics_data=genomics,
            transcriptomics_data=transcriptomics,
            method="cross_platform_mi",
            feature_indices={"genomics": 0, "transcriptomics": 0},
        )

        # Always returns matrix format
        assert "method" in results
        assert "genomics_transcriptomics_mi_matrix" in results
        assert len(results["genomics_transcriptomics_mi_matrix"]) == 1
        assert len(results["genomics_transcriptomics_mi_matrix"][0]) == 1
        assert "genomics_transcriptomics_mean_mi" in results


class TestMLIntegration:
    """Integration tests with ML module."""

    def test_ml_feature_selection_mi(self):
        """Test feature selection using mutual information."""
        X = np.random.randn(100, 50)
        y = np.random.randint(0, 2, 100)

        results = ml_integration(X, y, method="feature_mi")
        assert "feature_mis" in results
        assert "top_features" in results
        assert len(results["top_features"]) <= 10

    def test_ml_feature_entropy(self):
        """Test feature entropy calculation."""
        X = np.random.randn(100, 50)
        y = np.random.randint(0, 2, 100)

        results = ml_integration(X, y, method="feature_entropy")
        assert "feature_entropies" in results
        assert len(results["feature_entropies"]) == 50


class TestVisualizationIntegration:
    """Integration tests with visualization module."""

    def test_entropy_distribution_plot(self):
        """Test entropy distribution plotting."""
        try:
            from metainformant.information.visualization import plot_entropy_distribution

            entropies = [2.5, 3.1, 2.8, 3.0, 2.9]
            result = plot_entropy_distribution(
                entropies, output_path="output/information/test_entropy.png"
            )
            assert "output_path" in result
            assert result["num_values"] == 5
        except ImportError:
            pytest.skip("Visualization module not available")

    def test_mi_matrix_plot(self):
        """Test MI matrix plotting."""
        try:
            from metainformant.information.visualization import plot_mutual_information_matrix

            mi_matrix = np.random.rand(10, 10)
            result = plot_mutual_information_matrix(
                mi_matrix, output_path="output/information/test_mi.png"
            )
            assert "output_path" in result
        except ImportError:
            pytest.skip("Visualization module not available")


class TestNetworkIntegration:
    """Integration tests with networks module."""

    def test_network_entropy(self):
        """Test network entropy calculation."""
        try:
            import networkx as nx
            from metainformant.information.networks import network_entropy

            G = nx.karate_club_graph()
            entropy = network_entropy(G)
            assert entropy >= 0.0

            # Test with attribute
            for node in G.nodes():
                G.nodes[node]["club"] = "A" if node < 17 else "B"
            entropy_attr = network_entropy(G, attribute="club")
            assert entropy_attr >= 0.0
        except ImportError:
            pytest.skip("NetworkX not available")

    def test_information_flow(self):
        """Test information flow calculation."""
        try:
            import networkx as nx
            from metainformant.information.networks import information_flow

            G = nx.karate_club_graph()
            flow = information_flow(G)
            assert "path_length_entropy" in flow
        except ImportError:
            pytest.skip("NetworkX not available")


class TestCrossModuleIntegration:
    """Tests for cross-module integration patterns."""

    def test_dna_to_network_information(self):
        """Test information flow from DNA analysis to network analysis."""
        from metainformant.information import information_profile
        from metainformant.information.networks import network_entropy

        try:
            import networkx as nx

            # Analyze DNA sequences
            sequences = ["ATCGATCG", "AAAA", "ATCGATCG"]
            profile = information_profile(sequences, k=1)

            # Create network with entropy information
            G = nx.Graph()
            G.add_node("seq1", entropy=profile["entropy"])
            G.add_node("seq2", entropy=profile["entropy"] * 0.5)

            # Calculate network entropy
            entropy = network_entropy(G, attribute="entropy")
            assert entropy >= 0.0
        except ImportError:
            pytest.skip("NetworkX not available")

    def test_information_to_visualization_pipeline(self):
        """Test pipeline from information analysis to visualization."""
        try:
            from metainformant.information import information_profile
            from metainformant.information.visualization import plot_information_profile

            sequences = ["ATCGATCG", "AAAA"]
            profile = information_profile(sequences, k=1)

            plot_result = plot_information_profile(
                profile, output_path="output/information/test_pipeline.png"
            )
            assert "output_path" in plot_result
        except ImportError:
            pytest.skip("Visualization module not available")

