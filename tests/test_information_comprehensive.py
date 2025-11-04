"""Comprehensive tests for information module.

Tests cover syntactic and semantic information theory methods,
analysis functions, and integration with other modules.
"""

from __future__ import annotations

import math
from collections import Counter

import numpy as np
import pytest

from metainformant.information import (
    analyze_sequence_information,
    compare_sequences_information,
    conditional_entropy,
    conditional_mutual_information,
    cross_entropy,
    information_content,
    information_profile,
    information_signature,
    joint_entropy,
    kl_divergence,
    mutual_information,
    semantic_entropy,
    semantic_similarity,
    semantic_similarity_matrix,
    shannon_entropy,
    total_correlation,
    transfer_entropy,
)
from metainformant.information.syntactic import shannon_entropy_from_counts


class TestSyntacticInformation:
    """Tests for syntactic information theory methods."""

    def test_shannon_entropy_uniform(self):
        """Test Shannon entropy with uniform distribution."""
        probs = [0.25, 0.25, 0.25, 0.25]
        entropy = shannon_entropy(probs)
        assert abs(entropy - 2.0) < 1e-10  # Maximum entropy for 4 outcomes

    def test_shannon_entropy_certainty(self):
        """Test Shannon entropy with certain outcome."""
        probs = [1.0, 0.0, 0.0]
        entropy = shannon_entropy(probs)
        assert abs(entropy) < 1e-10  # Zero entropy for certainty

    def test_shannon_entropy_from_counts(self):
        """Test entropy calculation from counts."""
        counts = {"A": 50, "T": 30, "G": 20}
        entropy = shannon_entropy_from_counts(counts)
        assert entropy > 0.0
        assert entropy < 2.0  # Less than maximum for 3 symbols

    def test_mutual_information_perfect_correlation(self):
        """Test mutual information with perfect correlation."""
        x = [0, 1, 0, 1]
        y = [0, 1, 0, 1]  # Perfect correlation
        mi = mutual_information(x, y)
        assert abs(mi - 1.0) < 1e-10  # Should equal H(X) = 1.0

    def test_mutual_information_independent(self):
        """Test mutual information with independent variables."""
        x = [0, 0, 1, 1]
        y = [0, 1, 0, 1]  # Independent
        mi = mutual_information(x, y)
        assert abs(mi) < 1e-10  # Should be approximately zero

    def test_conditional_entropy(self):
        """Test conditional entropy calculation."""
        x = [0, 1, 0, 1]
        y = [0, 0, 1, 1]  # X = Y
        h_x_given_y = conditional_entropy(x, y)
        assert abs(h_x_given_y) < 1e-10  # No uncertainty given Y

    def test_kl_divergence_identical(self):
        """Test KL divergence with identical distributions."""
        p = [0.5, 0.3, 0.2]
        q = [0.5, 0.3, 0.2]
        kl = kl_divergence(p, q)
        assert abs(kl) < 1e-10  # Should be zero

    def test_kl_divergence_different(self):
        """Test KL divergence with different distributions."""
        p = [0.5, 0.5]
        q = [0.9, 0.1]
        kl = kl_divergence(p, q)
        assert kl > 0.0  # Should be positive

    def test_cross_entropy(self):
        """Test cross-entropy calculation."""
        p = [0.5, 0.5]
        q = [0.5, 0.5]
        ce = cross_entropy(p, q)
        assert ce > 0.0
        assert abs(ce - 1.0) < 1e-10  # Should equal entropy of p

    def test_total_correlation(self):
        """Test total correlation (multivariate mutual information)."""
        x = [0, 1, 0, 1]
        y = [0, 1, 0, 1]  # Perfect correlation
        tc = total_correlation([x, y])
        assert tc > 0.0  # Should have positive total correlation

    def test_transfer_entropy(self):
        """Test transfer entropy calculation."""
        x = [0, 1, 0, 1, 0, 1]
        y = [0, 0, 1, 1, 0, 0]
        te = transfer_entropy(x, y, lag=1)
        assert te >= 0.0  # Should be non-negative


class TestSemanticInformation:
    """Tests for semantic information theory methods."""

    def test_information_content(self):
        """Test information content calculation."""
        term_freqs = {"common": 100, "rare": 1}
        ic_common = information_content(term_freqs, "common")
        ic_rare = information_content(term_freqs, "rare")
        assert ic_rare > ic_common  # Rare terms have higher IC

    def test_information_content_not_found(self):
        """Test information content for non-existent term."""
        term_freqs = {"A": 10, "B": 20}
        ic = information_content(term_freqs, "C")
        assert ic == 0.0

    def test_semantic_entropy(self):
        """Test semantic entropy calculation."""
        annotations = {
            "entity1": {"GO:0008150", "GO:0003674"},
            "entity2": {"GO:0008150"},
            "entity3": {"GO:0003674"},
        }
        entropy = semantic_entropy(annotations)
        assert entropy > 0.0

    def test_semantic_similarity(self):
        """Test semantic similarity calculation."""
        term_ic = {"A": 2.0, "B": 2.0, "C": 1.0}
        hierarchy = {"A": {"C"}, "B": {"C"}}
        sim = semantic_similarity("A", "B", term_ic, hierarchy)
        assert 0.0 <= sim <= 1.0

    def test_semantic_similarity_identical(self):
        """Test semantic similarity for identical terms."""
        term_ic = {"A": 2.0}
        sim = semantic_similarity("A", "A", term_ic)
        assert abs(sim - 1.0) < 1e-10

    def test_semantic_similarity_matrix(self):
        """Test semantic similarity matrix."""
        terms = ["A", "B", "C"]
        term_ic = {"A": 2.0, "B": 2.0, "C": 1.0}
        matrix = semantic_similarity_matrix(terms, term_ic)
        assert matrix.shape == (3, 3)
        assert np.allclose(matrix, matrix.T)  # Should be symmetric


class TestAnalysisFunctions:
    """Tests for high-level analysis functions."""

    def test_information_profile(self):
        """Test information profile calculation."""
        sequences = ["ATCG", "ATCG", "AAAA"]
        profile = information_profile(sequences, k=1)
        assert "entropy" in profile
        assert "kmer_frequencies" in profile
        assert "unique_kmers" in profile
        assert profile["entropy"] > 0.0

    def test_information_profile_empty(self):
        """Test information profile with empty sequences."""
        profile = information_profile([], k=1)
        assert profile["entropy"] == 0.0

    def test_information_signature_entropy(self):
        """Test information signature with entropy method."""
        data = np.random.randn(100, 10)
        signature = information_signature(data, method="entropy")
        assert "signature" in signature
        assert "method" in signature
        assert signature["method"] == "entropy"

    def test_analyze_sequence_information(self):
        """Test sequence information analysis."""
        sequence = "ATCGATCG"
        analysis = analyze_sequence_information(sequence, k_values=[1, 2])
        assert "sequence_length" in analysis
        assert "kmer_analyses" in analysis
        assert 1 in analysis["kmer_analyses"]

    def test_compare_sequences_information(self):
        """Test sequence comparison."""
        seq1 = "ATCGATCG"
        seq2 = "ATCGATCG"
        comparison = compare_sequences_information(seq1, seq2, k=1)
        assert "entropy_1" in comparison
        assert "entropy_2" in comparison
        assert "kl_divergence" in comparison
        assert abs(comparison["kl_divergence"]) < 1e-10  # Identical sequences

    def test_compare_sequences_different(self):
        """Test comparison of different sequences."""
        seq1 = "ATCGATCG"
        seq2 = "AAAAAAAA"
        comparison = compare_sequences_information(seq1, seq2, k=1)
        assert comparison["kl_divergence"] > 0.0  # Different sequences


class TestEdgeCases:
    """Tests for edge cases and error handling."""

    def test_shannon_entropy_empty(self):
        """Test entropy with empty input."""
        entropy = shannon_entropy([])
        assert entropy == 0.0

    def test_mutual_information_length_mismatch(self):
        """Test mutual information with mismatched lengths."""
        x = [0, 1, 0]
        y = [0, 1]
        with pytest.raises(ValueError):
            mutual_information(x, y)

    def test_kl_divergence_length_mismatch(self):
        """Test KL divergence with mismatched lengths."""
        p = [0.5, 0.5]
        q = [0.33, 0.33, 0.34]
        with pytest.raises(ValueError):
            kl_divergence(p, q)

    def test_conditional_mutual_information_length_mismatch(self):
        """Test conditional MI with mismatched lengths."""
        x = [0, 1, 0]
        y = [0, 1]
        z = [0, 1, 0]
        with pytest.raises(ValueError):
            conditional_mutual_information(x, y, z)

    def test_transfer_entropy_short_sequence(self):
        """Test transfer entropy with short sequence."""
        x = [0, 1]
        y = [0, 1]
        te = transfer_entropy(x, y, lag=1)
        assert te == 0.0  # Should handle gracefully

    def test_total_correlation_single_variable(self):
        """Test total correlation with single variable."""
        x = [0, 1, 0, 1]
        tc = total_correlation([x])
        assert tc == 0.0  # Single variable has no correlation


class TestIntegration:
    """Tests for integration with other modules."""

    def test_information_profile_with_dna_sequences(self):
        """Test information profile with DNA-like sequences."""
        sequences = ["ATCGATCG", "GCTAGCTA", "ATCGATCG"]
        profile = information_profile(sequences, k=2)
        assert profile["entropy"] > 0.0
        assert "AT" in profile["kmer_frequencies"] or "TC" in profile["kmer_frequencies"]

    def test_conditional_entropy_properties(self):
        """Test properties of conditional entropy."""
        x = [0, 1, 0, 1, 0, 1]
        y = [0, 1, 0, 1, 0, 1]
        
        h_x = shannon_entropy_from_counts(Counter(x))
        h_x_given_y = conditional_entropy(x, y)
        
        # H(X|Y) <= H(X)
        assert h_x_given_y <= h_x + 1e-10

    def test_mutual_information_properties(self):
        """Test properties of mutual information."""
        x = [0, 1, 0, 1, 0, 1]
        y = [0, 1, 0, 1, 0, 1]
        
        h_x = shannon_entropy_from_counts(Counter(x))
        h_y = shannon_entropy_from_counts(Counter(y))
        mi = mutual_information(x, y)
        
        # I(X; Y) <= min(H(X), H(Y))
        assert mi <= min(h_x, h_y) + 1e-10
        # I(X; Y) >= 0
        assert mi >= 0.0


class TestAdditionalMethods:
    """Tests for additional information theory methods."""

    def test_jensen_shannon_divergence(self):
        """Test Jensen-Shannon divergence."""
        from metainformant.information.syntactic import jensen_shannon_divergence

        p = [0.5, 0.5]
        q = [0.5, 0.5]
        js = jensen_shannon_divergence(p, q)
        assert abs(js) < 1e-10  # Identical distributions

    def test_renyi_entropy(self):
        """Test Rényi entropy calculation."""
        from metainformant.information.syntactic import renyi_entropy, shannon_entropy

        probs = [0.25, 0.25, 0.25, 0.25]
        # α=1 should equal Shannon entropy
        h_renyi_1 = renyi_entropy(probs, alpha=1.0)
        h_shannon = shannon_entropy(probs)
        assert abs(h_renyi_1 - h_shannon) < 1e-10

        # α=2 (collision entropy)
        h_renyi_2 = renyi_entropy(probs, alpha=2.0)
        assert h_renyi_2 > 0.0

    def test_tsallis_entropy(self):
        """Test Tsallis entropy calculation."""
        from metainformant.information.syntactic import shannon_entropy, tsallis_entropy

        probs = [0.5, 0.5]
        # q=1 should equal Shannon entropy
        h_tsallis_1 = tsallis_entropy(probs, q=1.0)
        h_shannon = shannon_entropy(probs)
        assert abs(h_tsallis_1 - h_shannon) < 1e-10

        # q=2
        h_tsallis_2 = tsallis_entropy(probs, q=2.0)
        assert h_tsallis_2 > 0.0

    def test_normalized_mutual_information(self):
        """Test normalized mutual information."""
        from metainformant.information.syntactic import normalized_mutual_information

        x = [0, 1, 0, 1]
        y = [0, 1, 0, 1]  # Perfect correlation
        nmi = normalized_mutual_information(x, y)
        assert abs(nmi - 1.0) < 1e-10  # Should be 1.0

    def test_information_coefficient(self):
        """Test information coefficient."""
        from metainformant.information.syntactic import information_coefficient

        x = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
        y = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]  # Perfect correlation
        ic = information_coefficient(x, y)
        assert 0.0 <= ic <= 1.0


class TestContinuousMethods:
    """Tests for continuous information theory methods."""

    def test_differential_entropy(self):
        """Test differential entropy calculation."""
        from metainformant.information.continuous import differential_entropy

        import numpy as np

        samples = np.random.normal(0, 1, 1000)
        h = differential_entropy(samples, method="histogram")
        assert h > 0.0

    def test_mutual_information_continuous(self):
        """Test continuous mutual information."""
        from metainformant.information.continuous import mutual_information_continuous

        import numpy as np

        x = np.random.randn(1000)
        y = x + np.random.randn(1000) * 0.1  # Strong correlation
        mi = mutual_information_continuous(x, y)
        assert mi >= 0.0

    def test_entropy_estimation(self):
        """Test entropy estimation methods."""
        from metainformant.information.continuous import entropy_estimation

        import numpy as np

        samples = np.random.normal(0, 1, 1000)
        h_plugin = entropy_estimation(samples, method="plugin")
        assert h_plugin > 0.0

        h_mm = entropy_estimation(samples, method="miller_madow")
        assert h_mm > 0.0


class TestEstimationMethods:
    """Tests for estimation and bias correction methods."""

    def test_entropy_estimator(self):
        """Test entropy estimation with various methods."""
        from metainformant.information.estimation import entropy_estimator

        counts = {"A": 50, "T": 30, "G": 20}
        h_plugin = entropy_estimator(counts, method="plugin")
        h_mm = entropy_estimator(counts, method="miller_madow")
        assert h_plugin > 0.0
        assert h_mm > 0.0

    def test_mutual_information_estimator(self):
        """Test MI estimation with bias correction."""
        from metainformant.information.estimation import mutual_information_estimator

        x = [0, 1, 0, 1, 0, 1]
        y = [0, 1, 0, 1, 0, 1]
        mi = mutual_information_estimator(x, y, method="plugin", bias_correction=True)
        assert mi >= 0.0

    def test_bias_correction(self):
        """Test bias correction function."""
        from metainformant.information.estimation import bias_correction

        corrected = bias_correction(1.5, n_samples=100, n_bins=10, measure="entropy")
        assert corrected > 1.5  # Should increase for entropy

        corrected_mi = bias_correction(
            0.5, n_samples=100, n_bins=10, measure="mutual_information"
        )
        assert corrected_mi <= 0.5  # Should decrease for MI


class TestWorkflowFunctions:
    """Tests for workflow functions."""

    def test_batch_entropy_analysis(self):
        """Test batch entropy analysis."""
        from metainformant.information.workflows import batch_entropy_analysis

        sequences = ["ATCG", "ATCG", "AAAA"]
        results = batch_entropy_analysis(sequences, k=1)
        assert "sequences" in results
        assert "summary" in results
        assert len(results["sequences"]) == 3

    def test_information_workflow(self):
        """Test complete information workflow."""
        from metainformant.information.workflows import information_workflow

        sequences = ["ATCGATCG", "AAAA"]
        results = information_workflow(sequences, k_values=[1, 2])
        assert "profiles" in results
        assert 1 in results["profiles"]

    def test_compare_datasets(self):
        """Test dataset comparison."""
        from metainformant.information.workflows import compare_datasets

        dataset1 = ["ATCG", "ATCG"]
        dataset2 = ["AAAA", "AAAA"]
        comparison = compare_datasets(dataset1, dataset2, k=1, method="entropy")
        assert "method" in comparison
        assert comparison["method"] == "entropy"


class TestIntegrationFunctions:
    """Tests for integration functions."""

    def test_dna_integration(self):
        """Test DNA integration."""
        from metainformant.information.integration import dna_integration

        sequences = {"seq1": "ATCG", "seq2": "AAAA"}
        results = dna_integration(sequences, k=1, analysis_type="entropy")
        assert "entropy_analysis" in results

    def test_rna_integration(self):
        """Test RNA integration."""
        from metainformant.information.integration import rna_integration

        import numpy as np

        expression = np.random.randn(100, 50)
        results = rna_integration(expression, method="entropy")
        assert "gene_entropies" in results

    def test_ml_integration(self):
        """Test ML integration."""
        from metainformant.information.integration import ml_integration

        import numpy as np

        X = np.random.randn(100, 50)
        y = np.random.randint(0, 2, 100)
        results = ml_integration(X, y, method="feature_mi")
        assert "feature_mis" in results

