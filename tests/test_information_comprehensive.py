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

