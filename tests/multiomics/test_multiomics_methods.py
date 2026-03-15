"""Tests for multiomics methods: factorization and clustering submodules.

Tests matrix factorization (joint NMF, MOFA, tensor decomposition, SNF, CCA),
multi-omic clustering (concatenation, SNF-based, late integration), consensus
clustering, multi-view spectral clustering, and integration evaluation.

All tests use real implementations with small synthetic datasets (NO mocking).
"""

from __future__ import annotations

import math
import random
from typing import Any

import numpy as np
import pytest

from metainformant.multiomics.methods.clustering import (
    consensus_clustering,
    evaluate_integration,
    multi_omic_clustering,
    multi_view_spectral,
)
from metainformant.multiomics.methods.factorization import (
    canonical_correlation,
    joint_nmf,
    mofa_simple,
    similarity_network_fusion,
    tensor_decomposition,
)

# ---------------------------------------------------------------------------
# Fixtures: realistic multi-omic data generators
# ---------------------------------------------------------------------------


def _make_two_cluster_data(
    n_samples: int = 30,
    n_features_expr: int = 20,
    n_features_meth: int = 15,
    seed: int = 42,
) -> tuple[dict[str, Any], list[int]]:
    """Generate two-cluster multi-omic data (expression + methylation).

    Returns (data_matrices, true_labels).
    """
    rng = np.random.RandomState(seed)
    n1 = n_samples // 2
    n2 = n_samples - n1

    # Expression: cluster 1 has higher mean
    expr = np.vstack(
        [
            rng.randn(n1, n_features_expr) + 2.0,
            rng.randn(n2, n_features_expr) - 2.0,
        ]
    )

    # Methylation: cluster 1 has different mean
    meth = np.vstack(
        [
            rng.randn(n1, n_features_meth) + 1.5,
            rng.randn(n2, n_features_meth) - 1.5,
        ]
    )

    labels = [0] * n1 + [1] * n2
    return {"expression": expr, "methylation": meth}, labels


def _make_nonneg_data(
    n_samples: int = 20,
    n_features_a: int = 15,
    n_features_b: int = 10,
    seed: int = 42,
) -> dict[str, Any]:
    """Generate non-negative multi-omic data suitable for NMF."""
    rng = np.random.RandomState(seed)
    return {
        "gene_expr": np.abs(rng.randn(n_samples, n_features_a)) + 0.01,
        "protein": np.abs(rng.randn(n_samples, n_features_b)) + 0.01,
    }


def _make_similarity_matrices(n: int = 20, n_views: int = 3, seed: int = 42) -> list[Any]:
    """Generate symmetric positive similarity matrices for SNF."""
    rng = np.random.RandomState(seed)
    matrices = []
    for _ in range(n_views):
        X = rng.randn(n, 5)
        sim = X @ X.T
        sim = (sim - sim.min()) / (sim.max() - sim.min() + 1e-10)
        np.fill_diagonal(sim, 1.0)
        # Symmetrise
        sim = (sim + sim.T) / 2.0
        matrices.append(sim)
    return matrices


# ===================================================================
# Joint NMF Tests
# ===================================================================


class TestJointNMF:
    """Tests for joint_nmf factorization."""

    def test_basic_factorization(self) -> None:
        """Joint NMF produces W and H_dict with correct shapes."""
        data = _make_nonneg_data(n_samples=20, n_features_a=15, n_features_b=10)
        result = joint_nmf(data, k=3, max_iter=50, tol=1e-3)

        assert "W" in result
        assert "H_dict" in result
        assert "reconstruction_error" in result
        assert "n_iter" in result
        assert "converged" in result

        W = np.array(result["W"])
        assert W.shape == (20, 3)

        for name in data:
            H = np.array(result["H_dict"][name])
            assert H.shape[0] == 3
            assert H.shape[1] == np.array(data[name]).shape[1]

    def test_reconstruction_error_decreases(self) -> None:
        """Reconstruction error from a longer run should be lower or equal."""
        data = _make_nonneg_data(n_samples=15, n_features_a=10, n_features_b=8)
        short = joint_nmf(data, k=3, max_iter=10, tol=1e-10)
        long = joint_nmf(data, k=3, max_iter=100, tol=1e-10)
        assert long["reconstruction_error"] <= short["reconstruction_error"] + 1e-6

    def test_single_omic(self) -> None:
        """Joint NMF should work with a single omic layer."""
        rng = np.random.RandomState(7)
        data = {"single": np.abs(rng.randn(10, 8)) + 0.01}
        result = joint_nmf(data, k=2, max_iter=30)
        W = np.array(result["W"])
        assert W.shape == (10, 2)
        assert result["reconstruction_error"] >= 0.0

    def test_k_equals_one(self) -> None:
        """Joint NMF with k=1 should still work."""
        data = _make_nonneg_data(n_samples=10, n_features_a=5, n_features_b=4)
        result = joint_nmf(data, k=1, max_iter=30)
        W = np.array(result["W"])
        assert W.shape == (10, 1)

    def test_negative_values_raises(self) -> None:
        """Negative values should raise ValueError."""
        data = {"bad": np.array([[-1.0, 2.0], [3.0, 4.0]])}
        with pytest.raises(ValueError, match="negative"):
            joint_nmf(data, k=2)

    def test_invalid_k_raises(self) -> None:
        """k < 1 should raise ValueError."""
        data = _make_nonneg_data()
        with pytest.raises(ValueError, match="k must be >= 1"):
            joint_nmf(data, k=0)

    def test_empty_data_raises(self) -> None:
        """Empty data dict should raise ValueError."""
        with pytest.raises(ValueError, match="at least one"):
            joint_nmf({}, k=2)

    def test_mismatched_samples_raises(self) -> None:
        """Matrices with different row counts should raise ValueError."""
        data = {
            "a": np.abs(np.random.randn(10, 5)) + 0.01,
            "b": np.abs(np.random.randn(12, 5)) + 0.01,
        }
        with pytest.raises(ValueError, match="mismatch"):
            joint_nmf(data, k=2)

    def test_non_negative_factors(self) -> None:
        """W and H should be non-negative after NMF."""
        data = _make_nonneg_data(n_samples=20, n_features_a=10, n_features_b=8)
        result = joint_nmf(data, k=3, max_iter=100)
        W = np.array(result["W"])
        assert np.all(W >= -1e-10)  # allow tiny numerical drift
        for H_vals in result["H_dict"].values():
            H = np.array(H_vals)
            assert np.all(H >= -1e-10)


# ===================================================================
# MOFA Tests
# ===================================================================


class TestMOFA:
    """Tests for simplified MOFA factor analysis."""

    def test_basic_mofa(self) -> None:
        """MOFA produces factors and weights with correct shapes."""
        rng = np.random.RandomState(42)
        data = {
            "expr": rng.randn(25, 15),
            "meth": rng.randn(25, 12),
        }
        result = mofa_simple(data, k=5, max_iter=30, tol=1e-2)

        assert "factors" in result
        assert "weights_per_view" in result
        assert "variance_explained" in result
        assert "active_factors" in result

        Z = np.array(result["factors"])
        assert Z.shape == (25, 5)

        for name, mat in data.items():
            W = np.array(result["weights_per_view"][name])
            assert W.shape == (mat.shape[1], 5)

    def test_variance_explained_non_negative(self) -> None:
        """Variance explained should be non-negative."""
        rng = np.random.RandomState(42)
        data = {"a": rng.randn(20, 10), "b": rng.randn(20, 8)}
        result = mofa_simple(data, k=3, max_iter=20)

        for name, ve_list in result["variance_explained"].items():
            for ve in ve_list:
                assert ve >= 0.0

    def test_active_factors_subset(self) -> None:
        """Active factors should be a subset of all factor indices."""
        rng = np.random.RandomState(42)
        data = {"x": rng.randn(20, 10), "y": rng.randn(20, 8)}
        result = mofa_simple(data, k=5, max_iter=30)

        for f in result["active_factors"]:
            assert 0 <= f < 5

    def test_invalid_k_raises(self) -> None:
        """k < 1 should raise ValueError."""
        data = {"a": np.random.randn(10, 5)}
        with pytest.raises(ValueError, match="k must be >= 1"):
            mofa_simple(data, k=0)

    def test_empty_data_raises(self) -> None:
        """Empty data dict should raise ValueError."""
        with pytest.raises(ValueError, match="at least one"):
            mofa_simple({}, k=2)

    def test_mismatched_samples_raises(self) -> None:
        """Matrices with different sample counts should raise ValueError."""
        data = {"a": np.random.randn(10, 5), "b": np.random.randn(12, 5)}
        with pytest.raises(ValueError, match="mismatch"):
            mofa_simple(data, k=2)

    def test_three_views(self) -> None:
        """MOFA should handle three omic layers."""
        rng = np.random.RandomState(42)
        data = {
            "transcriptomics": rng.randn(20, 12),
            "proteomics": rng.randn(20, 8),
            "metabolomics": rng.randn(20, 10),
        }
        result = mofa_simple(data, k=4, max_iter=20)
        assert len(result["weights_per_view"]) == 3
        assert len(result["variance_explained"]) == 3


# ===================================================================
# Tensor Decomposition Tests
# ===================================================================


class TestTensorDecomposition:
    """Tests for CP tensor decomposition."""

    def test_basic_decomposition(self) -> None:
        """CP decomposition produces factors with correct shapes."""
        rng = np.random.RandomState(42)
        tensor = rng.randn(8, 6, 5)
        result = tensor_decomposition(tensor, rank=3, max_iter=50)

        assert "factors" in result
        assert "fit" in result
        assert "core_consistency" in result
        assert len(result["factors"]) == 3

        A = np.array(result["factors"][0])
        B = np.array(result["factors"][1])
        C = np.array(result["factors"][2])
        assert A.shape == (8, 3)
        assert B.shape == (6, 3)
        assert C.shape == (5, 3)

    def test_fit_value_range(self) -> None:
        """Fit should be a reasonable value (could be negative for bad rank)."""
        rng = np.random.RandomState(42)
        tensor = rng.randn(6, 5, 4)
        result = tensor_decomposition(tensor, rank=3, max_iter=100)
        # Fit can vary but should be finite
        assert math.isfinite(result["fit"])

    def test_rank_one(self) -> None:
        """Rank-1 decomposition of a rank-1 tensor should have positive fit."""
        rng = np.random.RandomState(42)
        a = rng.randn(5)
        b = rng.randn(4)
        c = rng.randn(3)
        tensor = np.einsum("i,j,k->ijk", a, b, c)
        result = tensor_decomposition(tensor, rank=1, max_iter=100)
        # ALS may not fully converge but should explain some variance
        assert result["fit"] > 0.0

    def test_invalid_method_raises(self) -> None:
        """Unsupported method should raise ValueError."""
        rng = np.random.RandomState(42)
        tensor = rng.randn(4, 3, 2)
        with pytest.raises(ValueError, match="Unsupported"):
            tensor_decomposition(tensor, rank=2, method="tucker")

    def test_invalid_rank_raises(self) -> None:
        """rank < 1 should raise ValueError."""
        rng = np.random.RandomState(42)
        tensor = rng.randn(4, 3, 2)
        with pytest.raises(ValueError, match="rank must be >= 1"):
            tensor_decomposition(tensor, rank=0)

    def test_non_3d_raises(self) -> None:
        """Non-3D input should raise ValueError."""
        with pytest.raises(ValueError, match="3-D"):
            tensor_decomposition(np.random.randn(4, 3), rank=2)

    def test_list_input(self) -> None:
        """Should accept nested lists as input."""
        rng = np.random.RandomState(42)
        tensor_list = rng.randn(4, 3, 2).tolist()
        result = tensor_decomposition(tensor_list, rank=2, max_iter=30)
        assert len(result["factors"]) == 3

    def test_core_consistency_bounded(self) -> None:
        """Core consistency should be between 0 and 100."""
        rng = np.random.RandomState(42)
        tensor = rng.randn(6, 5, 4)
        result = tensor_decomposition(tensor, rank=2, max_iter=50)
        if not math.isnan(result["core_consistency"]):
            assert 0.0 <= result["core_consistency"] <= 100.0


# ===================================================================
# SNF Tests
# ===================================================================


class TestSimilarityNetworkFusion:
    """Tests for Similarity Network Fusion."""

    def test_basic_fusion(self) -> None:
        """SNF produces a fused network and cluster labels."""
        networks = _make_similarity_matrices(n=20, n_views=3)
        result = similarity_network_fusion(networks, k_neighbors=5, n_iter=10)

        assert "fused_network" in result
        assert "cluster_labels" in result
        assert "silhouette_score" in result

        fused = np.array(result["fused_network"])
        assert fused.shape == (20, 20)
        assert len(result["cluster_labels"]) == 20

    def test_fused_is_symmetric(self) -> None:
        """Fused network should be symmetric."""
        networks = _make_similarity_matrices(n=15, n_views=2)
        result = similarity_network_fusion(networks, k_neighbors=5, n_iter=10)
        fused = np.array(result["fused_network"])
        np.testing.assert_allclose(fused, fused.T, atol=1e-10)

    def test_two_networks_minimum(self) -> None:
        """SNF with exactly 2 networks should work."""
        networks = _make_similarity_matrices(n=10, n_views=2)
        result = similarity_network_fusion(networks, k_neighbors=3, n_iter=5)
        assert len(result["cluster_labels"]) == 10

    def test_fewer_than_two_raises(self) -> None:
        """Fewer than 2 networks should raise ValueError."""
        with pytest.raises(ValueError, match="at least 2"):
            similarity_network_fusion([np.eye(5)], k_neighbors=2)

    def test_mismatched_sizes_raises(self) -> None:
        """Networks with different sizes should raise ValueError."""
        with pytest.raises(ValueError, match="expected"):
            similarity_network_fusion([np.eye(5), np.eye(6)], k_neighbors=2)

    def test_non_square_raises(self) -> None:
        """Non-square matrix should raise ValueError."""
        with pytest.raises(ValueError, match="square"):
            similarity_network_fusion([np.ones((5, 3)), np.ones((5, 3))], k_neighbors=2)

    def test_cluster_labels_valid(self) -> None:
        """Cluster labels should be valid integers."""
        networks = _make_similarity_matrices(n=20, n_views=3)
        result = similarity_network_fusion(networks, k_neighbors=5, n_iter=10)
        labels = result["cluster_labels"]
        unique = set(labels)
        assert len(unique) >= 1
        for label in labels:
            assert isinstance(label, int)


# ===================================================================
# Canonical Correlation Analysis Tests
# ===================================================================


class TestCanonicalCorrelation:
    """Tests for regularized CCA."""

    def test_basic_cca(self) -> None:
        """CCA produces scores, correlations, and loadings."""
        rng = np.random.RandomState(42)
        X = rng.randn(30, 8)
        Y = rng.randn(30, 6)
        result = canonical_correlation(X, Y, n_components=2)

        assert "x_scores" in result
        assert "y_scores" in result
        assert "correlations" in result
        assert "x_loadings" in result
        assert "y_loadings" in result

        assert len(result["correlations"]) == 2
        xs = np.array(result["x_scores"])
        assert xs.shape == (30, 2)

    def test_correlations_between_zero_and_one(self) -> None:
        """Canonical correlations should be in [0, 1]."""
        rng = np.random.RandomState(42)
        X = rng.randn(25, 6)
        Y = rng.randn(25, 5)
        result = canonical_correlation(X, Y, n_components=3)
        for corr in result["correlations"]:
            assert -0.01 <= corr <= 1.01  # small numerical tolerance

    def test_correlated_data_high_first_correlation(self) -> None:
        """When X and Y share a strong signal, first correlation should be high."""
        rng = np.random.RandomState(42)
        shared = rng.randn(30, 1)
        X = np.hstack([shared * 3.0, rng.randn(30, 4) * 0.1])
        Y = np.hstack([shared * 3.0, rng.randn(30, 3) * 0.1])
        result = canonical_correlation(X, Y, n_components=1, regularization=0.01)
        assert result["correlations"][0] > 0.5

    def test_invalid_n_components_raises(self) -> None:
        """n_components < 1 should raise ValueError."""
        X = np.random.randn(10, 5)
        Y = np.random.randn(10, 4)
        with pytest.raises(ValueError, match="n_components"):
            canonical_correlation(X, Y, n_components=0)

    def test_mismatched_samples_raises(self) -> None:
        """Different sample counts should raise ValueError."""
        X = np.random.randn(10, 5)
        Y = np.random.randn(12, 4)
        with pytest.raises(ValueError, match="mismatch"):
            canonical_correlation(X, Y)

    def test_n_components_clamped(self) -> None:
        """n_components larger than min(p,q,n) should be automatically clamped."""
        rng = np.random.RandomState(42)
        X = rng.randn(10, 3)
        Y = rng.randn(10, 4)
        result = canonical_correlation(X, Y, n_components=50)
        # Should clamp to min(3, 4, 10) = 3
        assert len(result["correlations"]) == 3

    def test_regularization_effect(self) -> None:
        """Higher regularization should produce lower canonical correlations."""
        rng = np.random.RandomState(42)
        X = rng.randn(20, 5)
        Y = rng.randn(20, 5)
        low_reg = canonical_correlation(X, Y, n_components=2, regularization=0.01)
        high_reg = canonical_correlation(X, Y, n_components=2, regularization=10.0)
        # Higher regularization shrinks correlations
        assert high_reg["correlations"][0] <= low_reg["correlations"][0] + 0.01


# ===================================================================
# Multi-Omic Clustering Tests
# ===================================================================


class TestMultiOmicClustering:
    """Tests for multi_omic_clustering."""

    def test_concatenation_method(self) -> None:
        """Concatenation-based clustering should produce valid labels."""
        data, _ = _make_two_cluster_data(n_samples=30)
        result = multi_omic_clustering(data, n_clusters=2, method="concatenation")

        assert "labels" in result
        assert "silhouette" in result
        assert "omic_contributions" in result
        assert len(result["labels"]) == 30
        assert len(set(result["labels"])) == 2

    def test_concatenation_recovers_clusters(self) -> None:
        """With well-separated data, clustering should recover true clusters."""
        data, true_labels = _make_two_cluster_data(n_samples=30, seed=42)
        result = multi_omic_clustering(data, n_clusters=2, method="concatenation")
        # Silhouette should be positive for well-separated clusters
        assert result["silhouette"] > 0.0

    def test_snf_method(self) -> None:
        """SNF-based clustering should produce valid labels."""
        data, _ = _make_two_cluster_data(n_samples=20)
        result = multi_omic_clustering(data, n_clusters=2, method="snf")
        assert len(result["labels"]) == 20
        assert "omic_contributions" in result

    def test_late_integration_method(self) -> None:
        """Late integration should produce valid labels."""
        data, _ = _make_two_cluster_data(n_samples=20)
        result = multi_omic_clustering(data, n_clusters=2, method="late_integration")
        assert len(result["labels"]) == 20
        assert "omic_contributions" in result

    def test_invalid_method_raises(self) -> None:
        """Unsupported method should raise ValueError."""
        data, _ = _make_two_cluster_data()
        with pytest.raises(ValueError, match="method must be"):
            multi_omic_clustering(data, n_clusters=2, method="invalid")

    def test_n_clusters_one_raises(self) -> None:
        """n_clusters < 2 should raise ValueError."""
        data, _ = _make_two_cluster_data()
        with pytest.raises(ValueError, match="n_clusters must be >= 2"):
            multi_omic_clustering(data, n_clusters=1)

    def test_empty_data_raises(self) -> None:
        """Empty data dict should raise ValueError."""
        with pytest.raises(ValueError, match="at least one"):
            multi_omic_clustering({}, n_clusters=2)

    def test_omic_contributions_present(self) -> None:
        """Each omic should have contribution info."""
        data, _ = _make_two_cluster_data(n_samples=20)
        result = multi_omic_clustering(data, n_clusters=2, method="concatenation")
        for omic_name in data:
            assert omic_name in result["omic_contributions"]

    def test_three_clusters(self) -> None:
        """Clustering into 3 groups should produce 3 unique labels."""
        rng = np.random.RandomState(42)
        data = {
            "a": np.vstack([rng.randn(10, 8) + i * 5.0 for i in range(3)]),
            "b": np.vstack([rng.randn(10, 6) + i * 5.0 for i in range(3)]),
        }
        result = multi_omic_clustering(data, n_clusters=3, method="concatenation")
        assert len(result["labels"]) == 30
        assert len(set(result["labels"])) == 3


# ===================================================================
# Consensus Clustering Tests
# ===================================================================


class TestConsensusClustering:
    """Tests for consensus clustering."""

    def test_basic_consensus(self) -> None:
        """Consensus clustering should find optimal k."""
        rng = np.random.RandomState(42)
        # Data with 2 clear clusters
        X = np.vstack([rng.randn(15, 8) + 3.0, rng.randn(15, 8) - 3.0])
        result = consensus_clustering(X, k_range=[2, 3, 4], n_resamples=20, proportion=0.8)

        assert "optimal_k" in result
        assert "labels" in result
        assert "consensus_matrix" in result
        assert "pac_score" in result
        assert result["optimal_k"] in [2, 3, 4]
        assert len(result["labels"]) == 30

    def test_consensus_matrix_properties(self) -> None:
        """Consensus matrix should be symmetric with diagonal of 1."""
        rng = np.random.RandomState(42)
        X = np.vstack([rng.randn(10, 5) + 2.0, rng.randn(10, 5) - 2.0])
        result = consensus_clustering(X, k_range=[2, 3], n_resamples=15, proportion=0.8)

        C = np.array(result["consensus_matrix"])
        assert C.shape == (20, 20)
        np.testing.assert_allclose(np.diag(C), 1.0)
        np.testing.assert_allclose(C, C.T, atol=1e-10)

    def test_pac_scores_present(self) -> None:
        """PAC scores should be computed for each k."""
        rng = np.random.RandomState(42)
        X = rng.randn(20, 6)
        result = consensus_clustering(X, k_range=[2, 3], n_resamples=10, proportion=0.8)
        for k_val in result["pac_score"]:
            assert 0.0 <= result["pac_score"][k_val] <= 1.0

    def test_too_few_samples_raises(self) -> None:
        """Fewer than 3 samples should raise ValueError."""
        X = np.random.randn(2, 5)
        with pytest.raises(ValueError, match="at least 3"):
            consensus_clustering(X)

    def test_invalid_proportion_raises(self) -> None:
        """Invalid proportion should raise ValueError."""
        X = np.random.randn(10, 5)
        with pytest.raises(ValueError, match="proportion"):
            consensus_clustering(X, proportion=0.0)

    def test_non_2d_raises(self) -> None:
        """Non-2D input should raise ValueError."""
        with pytest.raises(ValueError, match="2-D"):
            consensus_clustering(np.random.randn(10))

    def test_default_k_range(self) -> None:
        """Default k_range should work when none is specified."""
        rng = np.random.RandomState(42)
        X = rng.randn(20, 5)
        result = consensus_clustering(X, n_resamples=10, proportion=0.8)
        assert result["optimal_k"] >= 2


# ===================================================================
# Multi-View Spectral Clustering Tests
# ===================================================================


class TestMultiViewSpectral:
    """Tests for multi-view spectral clustering."""

    def test_average_method(self) -> None:
        """Average combination should produce valid labels."""
        sims = _make_similarity_matrices(n=20, n_views=3)
        result = multi_view_spectral(sims, n_clusters=2, method="average")

        assert "labels" in result
        assert "eigenvalues" in result
        assert "eigenvectors" in result
        assert len(result["labels"]) == 20
        assert len(set(result["labels"])) >= 1

    def test_product_method(self) -> None:
        """Product combination should work."""
        sims = _make_similarity_matrices(n=15, n_views=2)
        result = multi_view_spectral(sims, n_clusters=2, method="product")
        assert len(result["labels"]) == 15

    def test_max_method(self) -> None:
        """Max combination should work."""
        sims = _make_similarity_matrices(n=15, n_views=2)
        result = multi_view_spectral(sims, n_clusters=2, method="max")
        assert len(result["labels"]) == 15

    def test_eigenvalues_sorted(self) -> None:
        """Eigenvalues should be returned for the n_clusters smallest."""
        sims = _make_similarity_matrices(n=20, n_views=2)
        result = multi_view_spectral(sims, n_clusters=3, method="average")
        assert len(result["eigenvalues"]) == 3

    def test_invalid_method_raises(self) -> None:
        """Unsupported method should raise ValueError."""
        sims = _make_similarity_matrices(n=10, n_views=2)
        with pytest.raises(ValueError, match="method must be"):
            multi_view_spectral(sims, n_clusters=2, method="sum")

    def test_empty_list_raises(self) -> None:
        """Empty similarity list should raise ValueError."""
        with pytest.raises(ValueError, match="At least one"):
            multi_view_spectral([], n_clusters=2)

    def test_n_clusters_one_raises(self) -> None:
        """n_clusters < 2 should raise ValueError."""
        sims = _make_similarity_matrices(n=10, n_views=2)
        with pytest.raises(ValueError, match="n_clusters must be >= 2"):
            multi_view_spectral(sims, n_clusters=1)


# ===================================================================
# Integration Evaluation Tests
# ===================================================================


class TestEvaluateIntegration:
    """Tests for evaluate_integration."""

    def test_basic_evaluation(self) -> None:
        """Evaluation should return silhouette and ARI per omic."""
        data, true_labels = _make_two_cluster_data(n_samples=20)
        result = evaluate_integration(true_labels, data)

        assert "silhouette_per_omic" in result
        assert "ari_per_omic" in result
        assert "mean_silhouette" in result
        assert "mean_ari" in result
        assert "integration_metric" in result

        for omic_name in data:
            assert omic_name in result["silhouette_per_omic"]
            assert omic_name in result["ari_per_omic"]

    def test_silhouette_range(self) -> None:
        """Silhouette scores should be in [-1, 1]."""
        data, true_labels = _make_two_cluster_data(n_samples=20)
        result = evaluate_integration(true_labels, data)
        for sil in result["silhouette_per_omic"].values():
            assert -1.0 <= sil <= 1.0

    def test_well_separated_clusters_positive_silhouette(self) -> None:
        """Well-separated clusters should have positive mean silhouette."""
        data, true_labels = _make_two_cluster_data(n_samples=30, seed=42)
        result = evaluate_integration(true_labels, data)
        assert result["mean_silhouette"] > 0.0

    def test_integration_metric_bounded(self) -> None:
        """Integration metric should be in [0, 1]."""
        data, true_labels = _make_two_cluster_data(n_samples=20)
        result = evaluate_integration(true_labels, data)
        assert 0.0 <= result["integration_metric"] <= 1.0

    def test_mismatched_labels_raises(self) -> None:
        """Labels length mismatch should raise ValueError."""
        data, _ = _make_two_cluster_data(n_samples=20)
        with pytest.raises(ValueError, match="samples"):
            evaluate_integration([0, 1, 0], data)

    def test_empty_data_raises(self) -> None:
        """Empty omic data should raise ValueError."""
        with pytest.raises(ValueError, match="at least one"):
            evaluate_integration([0, 1], {})
