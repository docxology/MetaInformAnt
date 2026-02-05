"""Integration tests for information module with other modules.

Tests verify that information theory methods work correctly when
integrated with DNA, RNA, single-cell, multi-omics, and ML modules.
"""

from __future__ import annotations

from collections import Counter

import numpy as np
import pytest

from metainformant.information import (
    information_profile,
    mutual_information,
    shannon_entropy,
    shannon_entropy_from_counts,
)
from metainformant.information.integration.integration import (
    dna_integration,
    ml_integration,
    multiomics_integration,
    rna_integration,
    singlecell_integration,
)

# ============================================================
# DNA Integration Tests
# ============================================================


class TestDNAIntegration:
    """Integration tests with DNA module."""

    def test_dna_integration_basic(self) -> None:
        """Test DNA integration with default settings."""
        sequences = ["ATCGATCG", "GCTAGCTA", "AAATTTT"]
        results = dna_integration(sequences)
        assert results["n_sequences"] == 3
        assert results["k_values"] == [1]  # default k=1 -> k_values=[1]
        assert results["analysis_type"] == "entropy"
        assert "entropy_analysis" in results
        assert "integrated_metrics" in results

    def test_dna_integration_custom_k(self) -> None:
        """Test DNA integration with custom k_values."""
        sequences = ["ATCGATCG", "GCTAGCTA"]
        results = dna_integration(sequences, k_values=[1, 2])
        assert results["k_values"] == [1, 2]
        assert results["n_sequences"] == 2

    def test_dna_integration_single_k(self) -> None:
        """Test DNA integration with single k parameter."""
        sequences = ["ATCGATCG"]
        results = dna_integration(sequences, k=2)
        assert results["k_values"] == [2]

    def test_dna_integration_entropy_metrics(self) -> None:
        """Test that integrated metrics include entropy statistics."""
        sequences = ["ATCGATCG", "AAAAAAAA"]
        results = dna_integration(sequences)
        metrics = results["integrated_metrics"]
        assert "mean_entropy" in metrics
        assert "std_entropy" in metrics
        assert metrics["mean_entropy"] >= 0.0

    def test_dna_integration_profile_analysis(self) -> None:
        """Test DNA integration with profile analysis type."""
        sequences = ["ATCGATCG", "GCTAGCTA", "ATCGATCG"]
        results = dna_integration(sequences, analysis_type="profile", k=1)
        assert "profile" in results

    def test_dna_integration_empty_raises(self) -> None:
        """Test that empty sequences raises ValueError."""
        with pytest.raises(ValueError, match="No sequences"):
            dna_integration([])

    def test_dna_integration_dict_input(self) -> None:
        """Test DNA integration accepts dict input."""
        sequences = {"seq1": "ATCGATCG", "seq2": "GCTAGCTA"}
        results = dna_integration(sequences)
        assert results["n_sequences"] == 2
        assert "seq1" in results["entropy_analysis"]

    def test_dna_information_profile_integration(self) -> None:
        """Test information profile on aligned DNA sequences."""
        sequences = ["ATCGATCG", "ATCGATCG", "ATCGATCG"]
        profile = information_profile(sequences, k=1)
        assert "entropy" in profile
        assert "positions" in profile
        assert profile["entropy"] >= 0.0


# ============================================================
# RNA Integration Tests
# ============================================================


class TestRNAIntegration:
    """Integration tests with RNA expression data."""

    def test_rna_integration_basic(self) -> None:
        """Test RNA integration with default settings.

        Convention: rows = samples, columns = genes.
        """
        expression = np.random.randn(50, 20)  # 50 samples, 20 genes
        results = rna_integration(expression)
        assert results["n_samples"] == 50
        assert results["n_genes"] == 20
        assert results["method"] == "entropy"
        assert "gene_entropies" in results
        assert "integrated_metrics" in results

    def test_rna_integration_normalized(self) -> None:
        """Test RNA integration with normalization enabled."""
        expression = np.abs(np.random.randn(30, 10))
        results = rna_integration(expression, normalize=True)
        assert "integrated_metrics" in results
        metrics = results["integrated_metrics"]
        assert "gene_entropy_mean" in metrics

    def test_rna_integration_no_normalize(self) -> None:
        """Test RNA integration without normalization."""
        expression = np.random.randn(20, 10)  # 20 samples, 10 genes
        results = rna_integration(expression, normalize=False)
        assert results["n_samples"] == 20
        assert results["n_genes"] == 10

    def test_rna_integration_mi_method(self) -> None:
        """Test RNA integration with mutual information method."""
        expression = np.random.randn(30, 5)  # Small for speed
        results = rna_integration(expression, method="mutual_information")
        assert results["method"] == "mutual_information"
        assert "mi_matrix" in results
        assert len(results["mi_matrix"]) == 5

    def test_rna_entropy_from_counts_integration(self) -> None:
        """Test that entropy from counts works with discretized expression."""
        expression = np.random.randint(0, 10, size=100)
        counts = Counter(expression.tolist())
        entropy = shannon_entropy_from_counts(counts)
        assert entropy >= 0.0


# ============================================================
# Single-Cell Integration Tests
# ============================================================


class TestSingleCellIntegration:
    """Integration tests with single-cell data."""

    def test_singlecell_integration_with_anndata(self) -> None:
        """Test single-cell integration with AnnData object."""
        try:
            import anndata

            adata = anndata.AnnData(X=np.random.rand(100, 50))
            results = singlecell_integration(adata)
            assert results["n_cells"] == 100
            assert results["n_genes"] == 50
        except ImportError:
            pytest.skip("anndata not available")

    def test_singlecell_integration_mock_adata(self) -> None:
        """Test single-cell integration with mock AnnData structure."""

        class MockAnnData:
            def __init__(self, X: np.ndarray) -> None:
                self.X = X
                self.n_obs = X.shape[0]
                self.n_vars = X.shape[1]

        mock = MockAnnData(np.random.rand(50, 20))
        results = singlecell_integration(mock)
        assert results["n_cells"] == 50
        assert results["n_genes"] == 20
        assert results["method"] == "gene_entropy"
        assert "gene_entropies" in results

    def test_singlecell_cell_type_entropy(self) -> None:
        """Test single-cell cell type entropy analysis."""
        data = np.random.randint(0, 100, (100, 20))
        cell_types = ["TypeA"] * 50 + ["TypeB"] * 50
        results = singlecell_integration(data, cell_types=cell_types, method="cell_type_entropy")
        assert "cell_type_entropy" in results
        assert results["num_cell_types"] == 2
        assert "per_type_entropy" in results


# ============================================================
# Multi-Omics Integration Tests
# ============================================================


class TestMultiOmicsIntegration:
    """Integration tests with multi-omics data.

    Note: multiomics_integration uses keyword-only arguments.
    """

    def test_multiomics_platform_entropy(self) -> None:
        """Test multi-omics with platform entropy method."""
        omics = {
            "genomics": np.random.randn(100, 20),
            "transcriptomics": np.random.randn(100, 15),
        }
        results = multiomics_integration(omics_data=omics)
        assert "genomics" in results["omics_types"]
        assert "transcriptomics" in results["omics_types"]
        assert "genomics_entropy" in results
        assert "transcriptomics_entropy" in results

    def test_multiomics_cross_platform_mi(self) -> None:
        """Test cross-platform mutual information."""
        omics = {
            "genomics": np.random.randn(100, 10),
            "transcriptomics": np.random.randn(100, 10),
        }
        results = multiomics_integration(omics_data=omics, method="cross_platform_mi")
        assert results["method"] == "cross_platform_mi"
        assert "genomics_transcriptomics_mi_matrix" in results

    def test_multiomics_kwargs_input(self) -> None:
        """Test multi-omics with named kwargs for platforms."""
        results = multiomics_integration(
            genomics_data=np.random.randn(100, 10),
            transcriptomics_data=np.random.randn(100, 10),
        )
        assert "genomics" in results["omics_types"]
        assert "transcriptomics" in results["omics_types"]

    def test_multiomics_feature_indices(self) -> None:
        """Test cross-platform MI with specific feature indices."""
        omics = {
            "genomics": np.random.randn(100, 20),
            "transcriptomics": np.random.randn(100, 20),
        }
        results = multiomics_integration(
            omics_data=omics,
            method="cross_platform_mi",
            feature_indices={"genomics": [0, 1, 2], "transcriptomics": [0, 3, 5]},
        )
        assert "genomics_transcriptomics_mi_matrix" in results
        assert len(results["genomics_transcriptomics_mi_matrix"]) == 3
        assert len(results["genomics_transcriptomics_mi_matrix"][0]) == 3

    def test_multiomics_empty(self) -> None:
        """Test multi-omics with no data."""
        results = multiomics_integration()
        assert results["omics_types"] == []


# ============================================================
# ML Integration Tests
# ============================================================


class TestMLIntegration:
    """Integration tests with ML feature analysis."""

    def test_ml_feature_mi(self) -> None:
        """Test ML integration with feature mutual information."""
        X = np.random.randn(100, 20)
        y = np.random.randint(0, 2, 100)
        results = ml_integration(X, y, method="feature_mi")
        assert results["n_features"] == 20
        assert results["n_samples"] == 100
        assert "feature_mis" in results
        assert "top_features" in results
        assert len(results["top_features"]) <= 10

    def test_ml_mutual_info_alias(self) -> None:
        """Test that 'mutual_info' works as alias for 'feature_mi'."""
        X = np.random.randn(100, 10)
        y = np.random.randint(0, 2, 100)
        results = ml_integration(X, y, method="mutual_info")
        assert "feature_mis" in results

    def test_ml_feature_entropy(self) -> None:
        """Test ML integration with feature entropy."""
        X = np.random.randn(100, 15)
        y = np.random.randint(0, 3, 100)
        results = ml_integration(X, y, method="feature_entropy")
        assert "feature_entropies" in results
        assert "feature_scores" in results
        assert "mean_feature_entropy" in results["feature_scores"]

    def test_ml_entropy_alias(self) -> None:
        """Test that 'entropy' works as alias for 'feature_entropy'."""
        X = np.random.randn(100, 10)
        y = np.random.randint(0, 3, 100)
        results = ml_integration(X, y, method="entropy")
        assert "feature_entropies" in results

    def test_ml_correlation(self) -> None:
        """Test ML integration with correlation method."""
        X = np.random.randn(100, 10)
        y = np.random.randn(100)
        results = ml_integration(X, y, method="correlation")
        assert results["method"] == "correlation"
        assert "feature_scores" in results


# ============================================================
# Network Integration Tests
# ============================================================


class TestNetworkIntegration:
    """Integration tests with network information functions."""

    def test_network_entropy(self) -> None:
        """Test network entropy calculation."""
        try:
            import networkx as nx

            from metainformant.information.integration.networks import network_entropy

            G = nx.karate_club_graph()
            entropy = network_entropy(G)
            assert entropy >= 0.0
        except ImportError:
            pytest.skip("NetworkX not available")

    def test_network_entropy_with_attribute(self) -> None:
        """Test network entropy with node attribute."""
        try:
            import networkx as nx

            from metainformant.information.integration.networks import network_entropy

            G = nx.karate_club_graph()
            for node in G.nodes():
                G.nodes[node]["club"] = "A" if node < 17 else "B"
            entropy_attr = network_entropy(G, attribute="club")
            assert entropy_attr >= 0.0
        except ImportError:
            pytest.skip("NetworkX not available")

    def test_network_entropy_empty_graph(self) -> None:
        """Test network entropy on empty graph returns 0."""
        try:
            import networkx as nx

            from metainformant.information.integration.networks import network_entropy

            G = nx.Graph()
            entropy = network_entropy(G)
            assert entropy == 0.0
        except ImportError:
            pytest.skip("NetworkX not available")

    def test_information_flow(self) -> None:
        """Test information flow calculation."""
        try:
            import networkx as nx

            from metainformant.information.integration.networks import information_flow

            G = nx.karate_club_graph()
            flow = information_flow(G)
            assert "flow_matrix" in flow
            assert "method" in flow
        except ImportError:
            pytest.skip("NetworkX not available")

    def test_network_information_centrality(self) -> None:
        """Test network information centrality returns per-node scores."""
        try:
            import networkx as nx

            from metainformant.information.integration.networks import (
                network_information_centrality,
            )

            G = nx.karate_club_graph()
            centrality = network_information_centrality(G)
            # Returns flat dict: {node_id: score}
            assert isinstance(centrality, dict)
            assert len(centrality) == len(G.nodes())
            assert all(isinstance(v, (int, float, np.floating)) for v in centrality.values())
        except ImportError:
            pytest.skip("NetworkX not available")


# ============================================================
# Cross-Module Pipeline Tests
# ============================================================


class TestCrossModulePipelines:
    """Tests for pipelines that combine multiple modules."""

    def test_dna_to_entropy_pipeline(self) -> None:
        """Test pipeline: DNA sequences -> entropy analysis."""
        sequences = ["ATCGATCG", "GCTAGCTA"]
        dna_results = dna_integration(sequences)
        assert dna_results["n_sequences"] == 2
        metrics = dna_results.get("integrated_metrics", {})
        if metrics and "mean_entropy" in metrics:
            assert metrics["mean_entropy"] >= 0.0

    def test_rna_to_ml_pipeline(self) -> None:
        """Test pipeline: RNA expression -> ML feature analysis."""
        expression = np.abs(np.random.randn(50, 20))  # 50 samples, 20 genes
        rna_results = rna_integration(expression)
        assert "gene_entropies" in rna_results

        # Use same data for ML
        y = np.random.randint(0, 2, 50)
        ml_results = ml_integration(expression, y, method="feature_entropy")
        assert "feature_entropies" in ml_results

    def test_multiomics_pipeline(self) -> None:
        """Test full multi-omics integration pipeline."""
        results = multiomics_integration(
            omics_data={
                "genomics": np.random.randn(100, 10),
                "transcriptomics": np.random.randn(100, 10),
            },
        )
        assert len(results["omics_types"]) == 2

    def test_mutual_information_integration(self) -> None:
        """Test mutual information between paired sequences."""
        x = list("ATCGATCG")
        y = list("ATCGATCG")
        mi = mutual_information(x, y)
        assert mi >= 0.0

        x2 = list("ATCGATCG")
        y2 = list("GCTAGCTA")
        mi2 = mutual_information(x2, y2)
        assert mi2 >= 0.0
        assert mi >= mi2
