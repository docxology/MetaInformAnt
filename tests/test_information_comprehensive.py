"""Comprehensive tests for the information theory module.

Tests cover all submodules: syntactic, semantic, continuous, estimation,
analysis, advanced measures, network analysis, integration, and workflows.
All tests use real implementations per NO_MOCKING_POLICY.
"""

from __future__ import annotations

import math
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

# ============================================================
# Syntactic information tests
# ============================================================


class TestShannonEntropy:
    """Tests for Shannon entropy and related discrete measures."""

    def test_uniform_distribution_entropy(self) -> None:
        from metainformant.information.metrics.core.syntactic import shannon_entropy

        probs = [0.25, 0.25, 0.25, 0.25]
        h = shannon_entropy(probs)
        assert abs(h - 2.0) < 1e-10

    def test_certain_outcome_zero_entropy(self) -> None:
        from metainformant.information.metrics.core.syntactic import shannon_entropy

        probs = [1.0, 0.0, 0.0]
        h = shannon_entropy(probs)
        assert abs(h) < 1e-10

    def test_binary_distribution(self) -> None:
        from metainformant.information.metrics.core.syntactic import shannon_entropy

        probs = [0.5, 0.5]
        h = shannon_entropy(probs)
        assert abs(h - 1.0) < 1e-10

    def test_nats_base(self) -> None:
        from metainformant.information.metrics.core.syntactic import shannon_entropy

        probs = [0.5, 0.5]
        h_nats = shannon_entropy(probs, base=math.e)
        assert abs(h_nats - math.log(2)) < 1e-10

    def test_negative_probs_raises(self) -> None:
        from metainformant.information.metrics.core.syntactic import shannon_entropy

        with pytest.raises(ValueError, match="negative"):
            shannon_entropy([-0.1, 0.5, 0.6])

    def test_probs_not_summing_to_one_raises(self) -> None:
        from metainformant.information.metrics.core.syntactic import shannon_entropy

        with pytest.raises(ValueError, match="sum to 1"):
            shannon_entropy([0.3, 0.3])

    def test_entropy_from_counts_dict(self) -> None:
        from metainformant.information.metrics.core.syntactic import shannon_entropy_from_counts

        counts = {"A": 50, "T": 30, "G": 20}
        h = shannon_entropy_from_counts(counts)
        assert 0.0 < h < math.log2(3) + 0.01

    def test_entropy_from_counts_list(self) -> None:
        from metainformant.information.metrics.core.syntactic import shannon_entropy_from_counts

        h = shannon_entropy_from_counts([10, 10, 10, 10])
        assert abs(h - 2.0) < 1e-10

    def test_entropy_from_counts_all_zero(self) -> None:
        from metainformant.information.metrics.core.syntactic import shannon_entropy_from_counts

        h = shannon_entropy_from_counts([0, 0, 0])
        assert h == 0.0


class TestJointAndConditionalEntropy:
    """Tests for joint entropy, conditional entropy, and mutual information."""

    def test_joint_entropy_independent(self) -> None:
        from metainformant.information.metrics.core.syntactic import joint_entropy

        x = [0, 0, 1, 1]
        y = [0, 1, 0, 1]
        h_xy = joint_entropy(x, y)
        # For independent uniform binary: H(X,Y) = 2.0
        assert abs(h_xy - 2.0) < 1e-10

    def test_joint_entropy_identical(self) -> None:
        from metainformant.information.metrics.core.syntactic import joint_entropy

        x = [0, 1, 0, 1]
        y = [0, 1, 0, 1]
        h_xy = joint_entropy(x, y)
        # H(X,Y) = H(X) = 1.0 when X == Y
        assert abs(h_xy - 1.0) < 1e-10

    def test_conditional_entropy_identical(self) -> None:
        from metainformant.information.metrics.core.syntactic import conditional_entropy

        x = [0, 1, 0, 1]
        y = [0, 1, 0, 1]
        h_x_given_y = conditional_entropy(x, y)
        assert abs(h_x_given_y) < 1e-10

    def test_conditional_entropy_independent(self) -> None:
        from metainformant.information.metrics.core.syntactic import conditional_entropy, shannon_entropy_from_counts

        x = [0, 0, 1, 1]
        y = [0, 1, 0, 1]
        h_x_given_y = conditional_entropy(x, y)
        h_x = shannon_entropy_from_counts(Counter(x))
        # H(X|Y) = H(X) when independent
        assert abs(h_x_given_y - h_x) < 0.1

    def test_mutual_information_perfect(self) -> None:
        from metainformant.information.metrics.core.syntactic import mutual_information

        x = [0, 1, 0, 1]
        y = [0, 1, 0, 1]
        mi = mutual_information(x, y)
        assert abs(mi - 1.0) < 1e-10

    def test_mutual_information_independent(self) -> None:
        from metainformant.information.metrics.core.syntactic import mutual_information

        x = [0, 0, 1, 1]
        y = [0, 1, 0, 1]
        mi = mutual_information(x, y)
        assert abs(mi) < 1e-10

    def test_mutual_information_non_negative(self) -> None:
        from metainformant.information.metrics.core.syntactic import mutual_information

        rng = np.random.default_rng(42)
        x = rng.integers(0, 5, size=100).tolist()
        y = rng.integers(0, 5, size=100).tolist()
        mi = mutual_information(x, y)
        assert mi >= -1e-10

    def test_mutual_information_length_mismatch(self) -> None:
        from metainformant.information.metrics.core.syntactic import mutual_information

        with pytest.raises(ValueError):
            mutual_information([0, 1], [0, 1, 0])

    def test_conditional_mutual_information(self) -> None:
        from metainformant.information.metrics.core.syntactic import conditional_mutual_information

        x = [0, 1, 0, 1, 0, 1]
        y = [0, 1, 0, 1, 0, 1]
        z = [0, 0, 0, 1, 1, 1]
        cmi = conditional_mutual_information(x, y, z)
        assert cmi >= -1e-10

    def test_conditional_mutual_information_length_mismatch(self) -> None:
        from metainformant.information.metrics.core.syntactic import conditional_mutual_information

        with pytest.raises(ValueError):
            conditional_mutual_information([0, 1, 0], [0, 1], [0, 1, 0])


class TestDivergences:
    """Tests for KL divergence, cross-entropy, JSD."""

    def test_kl_divergence_identical(self) -> None:
        from metainformant.information.metrics.core.syntactic import kl_divergence

        p = [0.5, 0.3, 0.2]
        kl = kl_divergence(p, p)
        assert abs(kl) < 1e-10

    def test_kl_divergence_positive(self) -> None:
        from metainformant.information.metrics.core.syntactic import kl_divergence

        p = [0.5, 0.5]
        q = [0.9, 0.1]
        kl = kl_divergence(p, q)
        assert kl > 0.0

    def test_kl_divergence_infinite(self) -> None:
        from metainformant.information.metrics.core.syntactic import kl_divergence

        p = [0.5, 0.5]
        q = [1.0, 0.0]
        kl = kl_divergence(p, q)
        assert kl == float("inf")

    def test_kl_divergence_length_mismatch(self) -> None:
        from metainformant.information.metrics.core.syntactic import kl_divergence

        with pytest.raises(ValueError):
            kl_divergence([0.5, 0.5], [0.33, 0.33, 0.34])

    def test_cross_entropy_self(self) -> None:
        from metainformant.information.metrics.core.syntactic import cross_entropy, shannon_entropy

        p = [0.5, 0.5]
        ce = cross_entropy(p, p)
        h = shannon_entropy(p)
        assert abs(ce - h) < 1e-10

    def test_cross_entropy_greater_than_entropy(self) -> None:
        from metainformant.information.metrics.core.syntactic import cross_entropy, shannon_entropy

        p = [0.5, 0.5]
        q = [0.9, 0.1]
        ce = cross_entropy(p, q)
        h = shannon_entropy(p)
        # H(p,q) >= H(p) always (Gibbs inequality)
        assert ce >= h - 1e-10

    def test_jensen_shannon_divergence_identical(self) -> None:
        from metainformant.information.metrics.core.syntactic import jensen_shannon_divergence

        p = [0.5, 0.3, 0.2]
        jsd = jensen_shannon_divergence(p, p)
        assert abs(jsd) < 1e-10

    def test_jensen_shannon_divergence_bounded(self) -> None:
        from metainformant.information.metrics.core.syntactic import jensen_shannon_divergence

        p = [1.0, 0.0]
        q = [0.0, 1.0]
        jsd = jensen_shannon_divergence(p, q)
        # JSD is bounded by log(2) = 1 bit
        assert 0.0 <= jsd <= 1.0 + 1e-10


class TestGeneralizedEntropies:
    """Tests for Renyi, Tsallis, NMI, IC."""

    def test_renyi_entropy_alpha2(self) -> None:
        from metainformant.information.metrics.core.syntactic import renyi_entropy

        probs = [0.25, 0.25, 0.25, 0.25]
        h = renyi_entropy(probs, alpha=2.0)
        # For uniform: Renyi = Shannon = 2.0
        assert abs(h - 2.0) < 1e-10

    def test_renyi_entropy_alpha0_hartley(self) -> None:
        from metainformant.information.metrics.core.syntactic import renyi_entropy

        probs = [0.3, 0.3, 0.2, 0.2]
        h = renyi_entropy(probs, alpha=0.0)
        # Hartley entropy = log2(4) = 2.0
        assert abs(h - 2.0) < 1e-10

    def test_renyi_entropy_alpha1_raises(self) -> None:
        from metainformant.information.metrics.core.syntactic import renyi_entropy

        with pytest.raises(ValueError):
            renyi_entropy([0.5, 0.5], alpha=1.0)

    def test_tsallis_entropy_positive(self) -> None:
        from metainformant.information.metrics.core.syntactic import tsallis_entropy

        # Tsallis entropy: S_q = (1 - sum(p^q)) / (q - 1)
        # For q=2.0 with uniform [0.5, 0.5]:
        # S_2 = (1 - 0.5) / (2 - 1) = 0.5
        probs = [0.5, 0.5]
        h = tsallis_entropy(probs, q=2.0)
        assert np.isclose(h, 0.5)
        assert h > 0.0

        # For q=0.5: S_0.5 = (1 - 2*sqrt(0.5)) / (0.5 - 1) = (1 - 1.4142) / (-0.5) ≈ 0.8284
        h_low_q = tsallis_entropy(probs, q=0.5)
        assert h_low_q > 0.0

    def test_tsallis_entropy_q1_raises(self) -> None:
        from metainformant.information.metrics.core.syntactic import tsallis_entropy

        with pytest.raises(ValueError):
            tsallis_entropy([0.5, 0.5], q=1.0)

    def test_normalized_mutual_information_perfect(self) -> None:
        from metainformant.information.metrics.core.syntactic import normalized_mutual_information

        x = [0, 1, 0, 1]
        y = [0, 1, 0, 1]
        nmi = normalized_mutual_information(x, y, method="arithmetic")
        assert abs(nmi - 1.0) < 1e-10

    def test_normalized_mutual_information_geometric(self) -> None:
        from metainformant.information.metrics.core.syntactic import normalized_mutual_information

        x = [0, 1, 0, 1]
        y = [0, 1, 0, 1]
        nmi = normalized_mutual_information(x, y, method="geometric")
        assert abs(nmi - 1.0) < 1e-10

    def test_normalized_mutual_information_max(self) -> None:
        from metainformant.information.metrics.core.syntactic import normalized_mutual_information

        x = [0, 1, 0, 1]
        y = [0, 1, 0, 1]
        nmi = normalized_mutual_information(x, y, method="max")
        assert abs(nmi - 1.0) < 1e-10

    def test_normalized_mutual_information_invalid_method(self) -> None:
        from metainformant.information.metrics.core.syntactic import normalized_mutual_information

        with pytest.raises(ValueError):
            normalized_mutual_information([0, 1], [0, 1], method="invalid")

    def test_information_coefficient_range(self) -> None:
        from metainformant.information.metrics.core.syntactic import information_coefficient

        x = [0, 1, 0, 1]
        y = [0, 1, 0, 1]
        ic = information_coefficient(x, y)
        assert 0.0 <= ic <= 1.0 + 1e-10

    def test_total_correlation_perfect(self) -> None:
        from metainformant.information.metrics.core.syntactic import total_correlation

        x = [0, 1, 0, 1]
        y = [0, 1, 0, 1]
        tc = total_correlation([x, y])
        assert tc > 0.0

    def test_total_correlation_single(self) -> None:
        from metainformant.information.metrics.core.syntactic import total_correlation

        tc = total_correlation([[0, 1, 0, 1]])
        assert tc == 0.0

    def test_total_correlation_empty(self) -> None:
        from metainformant.information.metrics.core.syntactic import total_correlation

        tc = total_correlation([])
        assert tc == 0.0

    def test_transfer_entropy_basic(self) -> None:
        from metainformant.information.metrics.core.syntactic import transfer_entropy

        x = [0, 1, 0, 1, 0, 1]
        y = [0, 0, 1, 1, 0, 0]
        te = transfer_entropy(x, y, lag=1)
        assert te >= 0.0

    def test_transfer_entropy_short_raises(self) -> None:
        from metainformant.information.metrics.core.syntactic import transfer_entropy

        with pytest.raises(ValueError, match="too short"):
            transfer_entropy([0, 1], [0, 1], lag=1)

    def test_transfer_entropy_invalid_lag(self) -> None:
        from metainformant.information.metrics.core.syntactic import transfer_entropy

        with pytest.raises(ValueError, match="Lag"):
            transfer_entropy([0, 1, 0, 1], [0, 1, 0, 1], lag=0)


# ============================================================
# Semantic information tests
# ============================================================


class TestSemanticInformation:
    """Tests for semantic information measures."""

    def test_information_content_rare_higher(self) -> None:
        from metainformant.information.metrics.advanced.semantic import information_content

        freqs = {"common": 100, "rare": 1}
        ic_common = information_content(freqs, "common")
        ic_rare = information_content(freqs, "rare")
        assert ic_rare > ic_common

    def test_information_content_not_found(self) -> None:
        from metainformant.information.metrics.advanced.semantic import information_content

        with pytest.raises(ValueError, match="not found"):
            information_content({"A": 10}, "B")

    def test_information_content_from_annotations(self) -> None:
        from metainformant.information.metrics.advanced.semantic import information_content_from_annotations

        annotations = {
            "gene1": {"GO:001", "GO:002"},
            "gene2": {"GO:001"},
            "gene3": {"GO:003"},
        }
        ic = information_content_from_annotations(annotations, "GO:001")
        assert ic > 0.0

    def test_semantic_entropy_positive(self) -> None:
        from metainformant.information.metrics.advanced.semantic import semantic_entropy

        annotations = {
            "e1": {"GO:001", "GO:002"},
            "e2": {"GO:001"},
            "e3": {"GO:002"},
        }
        h = semantic_entropy(annotations)
        assert h > 0.0

    def test_semantic_entropy_empty(self) -> None:
        from metainformant.information.metrics.advanced.semantic import semantic_entropy

        assert semantic_entropy({}) == 0.0

    def test_semantic_similarity_resnik(self) -> None:
        from metainformant.information.metrics.advanced.semantic import semantic_similarity

        term_ic = {"A": 2.0, "B": 2.0, "root": 0.5}
        hierarchy = {"A": {"root"}, "B": {"root"}}
        sim = semantic_similarity("A", "B", term_ic, hierarchy, method="resnik")
        assert 0.0 <= sim <= 1.0

    def test_semantic_similarity_lin(self) -> None:
        from metainformant.information.metrics.advanced.semantic import semantic_similarity

        term_ic = {"A": 2.0, "B": 2.0, "root": 1.0}
        hierarchy = {"A": {"root"}, "B": {"root"}}
        sim = semantic_similarity("A", "B", term_ic, hierarchy, method="lin")
        assert 0.0 <= sim <= 1.0

    def test_semantic_similarity_jiang_conrath(self) -> None:
        from metainformant.information.metrics.advanced.semantic import semantic_similarity

        term_ic = {"A": 2.0, "B": 2.0, "root": 1.0}
        hierarchy = {"A": {"root"}, "B": {"root"}}
        sim = semantic_similarity("A", "B", term_ic, hierarchy, method="jiang_conrath")
        assert 0.0 <= sim <= 1.0

    def test_semantic_similarity_matrix_symmetric(self) -> None:
        from metainformant.information.metrics.advanced.semantic import semantic_similarity_matrix

        terms = ["A", "B", "C"]
        term_ic = {"A": 2.0, "B": 1.5, "C": 1.0, "root": 0.5}
        hierarchy = {"A": {"root"}, "B": {"root"}, "C": {"root"}}
        matrix = semantic_similarity_matrix(terms, term_ic, hierarchy)
        assert len(matrix) == 3
        assert len(matrix[0]) == 3
        # Symmetric
        for i in range(3):
            for j in range(3):
                assert abs(matrix[i][j] - matrix[j][i]) < 1e-10
            # Diagonal = 1
            assert abs(matrix[i][i] - 1.0) < 1e-10

    def test_semantic_distance_jiang_conrath(self) -> None:
        from metainformant.information.metrics.advanced.semantic import semantic_distance

        term_ic = {"A": 2.0, "B": 2.0, "root": 0.5}
        hierarchy = {"A": {"root"}, "B": {"root"}}
        d = semantic_distance("A", "B", term_ic, hierarchy, method="jiang_conrath")
        assert d >= 0.0

    def test_term_specificity(self) -> None:
        from metainformant.information.metrics.advanced.semantic import term_specificity

        term_ic = {"leaf": 5.0, "root": 0.5}
        hierarchy = {"leaf": {"root"}}
        spec = term_specificity("leaf", term_ic, hierarchy)
        assert 0.0 < spec <= 1.0

    def test_ontology_complexity(self) -> None:
        from metainformant.information.metrics.advanced.semantic import ontology_complexity

        hierarchy = {"A": {"root"}, "B": {"root"}, "C": {"A"}}
        result = ontology_complexity(hierarchy)
        assert "n_terms" in result
        assert result["n_terms"] == 3
        assert "avg_depth" in result

    def test_ontology_complexity_with_ic(self) -> None:
        from metainformant.information.metrics.advanced.semantic import ontology_complexity

        hierarchy = {"A": {"root"}, "B": {"root"}}
        term_ic = {"A": 2.0, "B": 1.5}
        result = ontology_complexity(hierarchy, term_ic=term_ic)
        assert "ic_mean" in result
        assert abs(result["ic_mean"] - 1.75) < 1e-10

    def test_term_redundancy(self) -> None:
        from metainformant.information.metrics.advanced.semantic import term_redundancy

        annotations = {
            "g1": {"A", "B"},
            "g2": {"A"},
            "g3": {"A", "C"},
        }
        r = term_redundancy(annotations, "A")
        # A appears 3 times, alone 1 time => redundancy = 1 - 1/3 = 0.667
        assert abs(r - (1.0 - 1 / 3)) < 1e-10

    def test_annotation_specificity(self) -> None:
        from metainformant.information.metrics.advanced.semantic import annotation_specificity

        annotations = {"g1": {"A", "B"}, "g2": {"C"}}
        term_ic = {"A": 1.0, "B": 3.0, "C": 2.0}
        result = annotation_specificity(annotations, term_ic)
        assert "g1" in result
        assert abs(result["g1"] - 2.0) < 1e-10  # avg(1, 3)
        assert abs(result["g2"] - 2.0) < 1e-10


# ============================================================
# Continuous information tests
# ============================================================


class TestContinuousInformation:
    """Tests for continuous (differential) entropy and MI."""

    def test_differential_entropy_histogram(self) -> None:
        from metainformant.information.metrics.core.continuous import differential_entropy

        rng = np.random.default_rng(42)
        samples = rng.normal(0, 1, 1000)
        h = differential_entropy(samples, method="histogram")
        # Theoretical: 0.5 * log2(2*pi*e) ≈ 2.047
        assert h > 0.0

    def test_differential_entropy_knn(self) -> None:
        from metainformant.information.metrics.core.continuous import differential_entropy

        rng = np.random.default_rng(42)
        samples = rng.normal(0, 1, 500)
        h = differential_entropy(samples, method="knn")
        assert isinstance(h, float)

    def test_differential_entropy_too_few_samples(self) -> None:
        from metainformant.information.metrics.core.continuous import differential_entropy

        with pytest.raises(ValueError, match="at least 10"):
            differential_entropy(np.array([1.0, 2.0, 3.0]), method="histogram")

    def test_mutual_information_continuous_correlated(self) -> None:
        from metainformant.information.metrics.core.continuous import mutual_information_continuous

        rng = np.random.default_rng(42)
        x = rng.normal(0, 1, 500)
        y = x + rng.normal(0, 0.1, 500)
        mi = mutual_information_continuous(x, y)
        assert mi > 0.0

    def test_mutual_information_continuous_length_mismatch(self) -> None:
        from metainformant.information.metrics.core.continuous import mutual_information_continuous

        with pytest.raises(ValueError):
            mutual_information_continuous(np.array([1.0, 2.0]), np.array([1.0]))

    def test_kl_divergence_continuous_similar(self) -> None:
        from metainformant.information.metrics.core.continuous import kl_divergence_continuous

        rng = np.random.default_rng(42)
        p = rng.normal(0, 1, 500)
        q = rng.normal(0, 1, 500)
        kl = kl_divergence_continuous(p, q, method="histogram")
        assert kl >= 0.0

    def test_entropy_estimation_alias(self) -> None:
        from metainformant.information.metrics.core.continuous import entropy_estimation

        rng = np.random.default_rng(42)
        samples = rng.normal(0, 1, 100)
        h = entropy_estimation(samples)
        assert isinstance(h, float)

    def test_copula_entropy(self) -> None:
        from metainformant.information.metrics.core.continuous import copula_entropy

        rng = np.random.default_rng(42)
        data = rng.normal(0, 1, (100, 3))
        ce = copula_entropy(data)
        assert isinstance(ce, float)

    def test_copula_entropy_1d_raises(self) -> None:
        from metainformant.information.metrics.core.continuous import copula_entropy

        with pytest.raises(ValueError):
            copula_entropy(np.array([[1.0], [2.0]]))

    def test_transfer_entropy_continuous(self) -> None:
        from metainformant.information.metrics.core.continuous import transfer_entropy_continuous

        rng = np.random.default_rng(42)
        x = rng.normal(0, 1, 200)
        y = np.roll(x, 1) + rng.normal(0, 0.5, 200)
        te = transfer_entropy_continuous(x, y, lag=1)
        assert te >= 0.0

    def test_conditional_entropy_continuous(self) -> None:
        from metainformant.information.metrics.core.continuous import conditional_entropy_continuous

        rng = np.random.default_rng(42)
        x = rng.normal(0, 1, 200)
        y = rng.normal(0, 1, 200)
        h = conditional_entropy_continuous(x, y)
        assert h >= 0.0

    def test_information_flow_network(self) -> None:
        from metainformant.information.metrics.core.continuous import information_flow_network

        rng = np.random.default_rng(42)
        data = rng.normal(0, 1, (3, 50))
        flow = information_flow_network(data, lag=1)
        assert flow.shape == (3, 3)
        # Diagonal should be zero (no self-flow)
        for i in range(3):
            assert flow[i, i] == 0.0


# ============================================================
# Estimation tests
# ============================================================


class TestEstimation:
    """Tests for entropy estimation and bias correction."""

    def test_plugin_estimator(self) -> None:
        from metainformant.information.metrics.core.estimation import entropy_estimator

        counts = {"A": 50, "T": 30, "G": 20}
        h = entropy_estimator(counts, method="plugin")
        assert h > 0.0

    def test_miller_madow_estimator(self) -> None:
        from metainformant.information.metrics.core.estimation import entropy_estimator

        counts = {"A": 50, "T": 30, "G": 20}
        h = entropy_estimator(counts, method="miller_madow")
        assert h > 0.0

    def test_chao_shen_estimator(self) -> None:
        from metainformant.information.metrics.core.estimation import entropy_estimator

        counts = {"A": 50, "T": 2, "G": 1, "C": 1}
        h = entropy_estimator(counts, method="chao_shen")
        assert h >= 0.0

    def test_jackknife_estimator(self) -> None:
        from metainformant.information.metrics.core.estimation import entropy_estimator

        counts = {"A": 50, "T": 30, "G": 20}
        h = entropy_estimator(counts, method="jackknife")
        assert h >= 0.0

    def test_invalid_method_raises(self) -> None:
        from metainformant.information.metrics.core.estimation import entropy_estimator

        with pytest.raises(ValueError, match="Unknown"):
            entropy_estimator({"A": 10}, method="bogus")

    def test_negative_counts_raises(self) -> None:
        from metainformant.information.metrics.core.estimation import entropy_estimator

        with pytest.raises(ValueError, match="negative"):
            entropy_estimator({"A": -1, "B": 5})

    def test_zero_counts(self) -> None:
        from metainformant.information.metrics.core.estimation import entropy_estimator

        h = entropy_estimator({"A": 0, "B": 0})
        assert h == 0.0

    def test_mutual_information_estimator(self) -> None:
        from metainformant.information.metrics.core.estimation import mutual_information_estimator

        x = [0, 1, 0, 1, 0, 1]
        y = [0, 1, 0, 1, 0, 1]
        mi = mutual_information_estimator(x, y, method="plugin")
        assert mi >= 0.0

    def test_kl_divergence_estimator(self) -> None:
        from metainformant.information.metrics.core.estimation import kl_divergence_estimator

        p = [0, 0, 1, 1, 0, 0]
        q = [0, 0, 1, 1, 0, 0]
        kl = kl_divergence_estimator(p, q)
        assert abs(kl) < 0.1

    def test_bias_correction(self) -> None:
        from metainformant.information.metrics.core.estimation import bias_correction

        corrected = bias_correction(entropy=2.0, sample_size=100, alphabet_size=4)
        # Correction = (4-1)/(2*100) = 0.015 => corrected = 2.0 - 0.015
        assert abs(corrected - 1.985) < 1e-10

    def test_bias_correction_invalid_params(self) -> None:
        from metainformant.information.metrics.core.estimation import bias_correction

        with pytest.raises(ValueError):
            bias_correction(2.0, sample_size=0, alphabet_size=4)

    def test_entropy_bootstrap_confidence(self) -> None:
        from metainformant.information.metrics.core.estimation import entropy_bootstrap_confidence

        counts = {"A": 50, "T": 30, "G": 20}
        result = entropy_bootstrap_confidence(counts, n_bootstraps=100, random_state=42)
        assert "entropy" in result
        assert "ci_lower" in result
        assert "ci_upper" in result
        assert result["ci_lower"] <= result["entropy"] <= result["ci_upper"]

    def test_panzeri_treves_correction(self) -> None:
        from metainformant.information.metrics.core.estimation import panzeri_treves_bias_correction

        corrected = panzeri_treves_bias_correction(entropy=2.0, sample_size=100, alphabet_size=4)
        assert corrected >= 0.0

    def test_entropy_rate_estimator(self) -> None:
        from metainformant.information.metrics.core.estimation import entropy_rate_estimator

        sequence = [0, 1, 0, 1, 0, 1, 0, 1, 0, 1]
        rate = entropy_rate_estimator(sequence, order=1)
        assert rate >= 0.0

    def test_entropy_rate_short_raises(self) -> None:
        from metainformant.information.metrics.core.estimation import entropy_rate_estimator

        with pytest.raises(ValueError, match="too short"):
            entropy_rate_estimator([0, 1], order=1)


# ============================================================
# Analysis tests
# ============================================================


class TestAnalysis:
    """Tests for high-level analysis functions."""

    def test_information_profile(self) -> None:
        from metainformant.information.metrics.analysis import information_profile

        sequences = ["ATCGATCG", "GCTAGCTA", "ATCGATCG"]
        result = information_profile(sequences, k=1)
        assert "profile" in result
        assert "statistics" in result
        assert "parameters" in result
        assert len(result["profile"]) == 8  # 8 positions for k=1

    def test_information_profile_k2(self) -> None:
        from metainformant.information.metrics.analysis import information_profile

        sequences = ["ATCGATCG", "GCTAGCTA", "ATCGATCG"]
        result = information_profile(sequences, k=2)
        assert len(result["profile"]) == 7  # 7 positions for k=2

    def test_information_profile_empty_raises(self) -> None:
        from metainformant.information.metrics.analysis import information_profile

        with pytest.raises(ValueError, match="empty"):
            information_profile([], k=1)

    def test_information_profile_unequal_raises(self) -> None:
        from metainformant.information.metrics.analysis import information_profile

        with pytest.raises(ValueError, match="same length"):
            information_profile(["ATCG", "AT"], k=1)

    def test_information_signature_entropy(self) -> None:
        from metainformant.information.metrics.analysis import information_signature

        rng = np.random.default_rng(42)
        data = rng.normal(0, 1, (100, 5))
        sig = information_signature(data, method="entropy")
        assert "method" in sig
        assert sig["method"] == "entropy"
        assert "features" in sig
        assert "summary" in sig
        assert len(sig["features"]) == 5

    def test_information_signature_mutual_info(self) -> None:
        from metainformant.information.metrics.analysis import information_signature

        rng = np.random.default_rng(42)
        data = rng.normal(0, 1, (50, 3))
        sig = information_signature(data, method="mutual_info")
        assert "mutual_information_matrix" in sig

    def test_information_signature_invalid_method(self) -> None:
        from metainformant.information.metrics.analysis import information_signature

        with pytest.raises(ValueError, match="Unknown method"):
            information_signature(np.random.randn(20, 3), method="invalid")

    def test_analyze_sequence_dna(self) -> None:
        from metainformant.information.metrics.analysis import analyze_sequence_information

        # Pure ATCG sequence — all chars are valid nucleotides → DNA
        result = analyze_sequence_information("ATCGATCGATCG", k_values=[1, 2])
        assert result["sequence_type"] == "DNA"
        assert result["sequence_length"] == 12
        assert "k1" in result["k_mer_analysis"]
        assert "k2" in result["k_mer_analysis"]

    def test_analyze_sequence_protein(self) -> None:
        from metainformant.information.metrics.analysis import analyze_sequence_information

        result = analyze_sequence_information("MKFLILLFNILCL", k_values=[1])
        assert result["sequence_type"] == "protein"

    def test_analyze_sequence_short_raises(self) -> None:
        from metainformant.information.metrics.analysis import analyze_sequence_information

        with pytest.raises(ValueError, match="too short"):
            analyze_sequence_information("ATCG")

    def test_analyze_sequence_complexity(self) -> None:
        from metainformant.information.metrics.analysis import analyze_sequence_information

        result = analyze_sequence_information("ATCGATCGATCG", methods=["complexity"])
        assert "complexity" in result["methods"]
        assert "linguistic_complexity" in result["methods"]["complexity"]
        assert "compression_ratio" in result["methods"]["complexity"]

    def test_compare_sequences_mutual_info(self) -> None:
        from metainformant.information.metrics.analysis import compare_sequences_information

        result = compare_sequences_information("ATCGATCG", "GCTAGCTA", k=1, method="mutual_info")
        assert "positional_mutual_info" in result
        assert "mean_mi" in result

    def test_compare_sequences_length_mismatch(self) -> None:
        from metainformant.information.metrics.analysis import compare_sequences_information

        with pytest.raises(ValueError, match="same length"):
            compare_sequences_information("ATCG", "AT", k=1)


# ============================================================
# Advanced measures tests
# ============================================================


class TestAdvancedMeasures:
    """Tests for Fisher information, variation of information, etc."""

    def test_fisher_information_normal(self) -> None:
        from metainformant.information.metrics.analysis import fisher_information

        rng = np.random.default_rng(42)
        samples = rng.normal(0, 2.0, 1000)
        fi = fisher_information(samples, method="parametric_normal")
        # I(mu) = 1/sigma^2 ≈ 0.25
        assert abs(fi - 0.25) < 0.1

    def test_fisher_information_empirical(self) -> None:
        from metainformant.information.metrics.analysis import fisher_information

        rng = np.random.default_rng(42)
        samples = rng.normal(0, 1, 500)
        fi = fisher_information(samples, method="empirical")
        assert fi > 0.0

    def test_fisher_information_too_few(self) -> None:
        from metainformant.information.metrics.analysis import fisher_information

        with pytest.raises(ValueError, match="at least 5"):
            fisher_information(np.array([1.0, 2.0]))

    def test_fisher_information_matrix_normal(self) -> None:
        from metainformant.information.metrics.analysis import fisher_information_matrix

        rng = np.random.default_rng(42)
        cov = np.array([[1.0, 0.5], [0.5, 2.0]])
        samples = rng.multivariate_normal([0, 0], cov, 500)
        fim = fisher_information_matrix(samples, method="parametric_normal")
        assert fim.shape == (2, 2)
        # FIM ≈ inv(cov)
        inv_cov = np.linalg.inv(cov)
        assert np.allclose(fim, inv_cov, atol=0.3)

    def test_variation_of_information_identical(self) -> None:
        from metainformant.information.metrics.analysis import variation_of_information

        x = [0, 1, 0, 1, 0]
        vi = variation_of_information(x, x)
        assert abs(vi) < 1e-10

    def test_variation_of_information_different(self) -> None:
        from metainformant.information.metrics.analysis import variation_of_information

        x = [0, 0, 1, 1]
        y = [0, 1, 0, 1]
        vi = variation_of_information(x, y)
        assert vi > 0.0

    def test_relative_information_gain(self) -> None:
        from metainformant.information.metrics.analysis import relative_information_gain

        prior = [0.5, 0.5]
        posterior = [0.9, 0.1]
        gain = relative_information_gain(prior, posterior)
        assert gain > 0.0

    def test_interaction_information(self) -> None:
        from metainformant.information.metrics.analysis import interaction_information

        x = [0, 1, 0, 1, 0, 1]
        y = [0, 1, 0, 1, 0, 1]
        z = [0, 0, 1, 1, 0, 0]
        ii = interaction_information(x, y, z)
        assert isinstance(ii, float)

    def test_binding_information(self) -> None:
        from metainformant.information.metrics.analysis import binding_information

        x = [0, 1, 0, 1]
        y = [0, 1, 0, 1]
        bi = binding_information(x, y)
        # Same as MI for pairs
        assert abs(bi - 1.0) < 1e-10

    def test_lautum_information_identical(self) -> None:
        from metainformant.information.metrics.analysis import lautum_information

        x = [0, 0, 1, 1]
        y = [0, 1, 0, 1]
        li = lautum_information(x, y)
        assert li >= 0.0

    def test_lautum_information_length_mismatch(self) -> None:
        from metainformant.information.metrics.analysis import lautum_information

        with pytest.raises(ValueError):
            lautum_information([0, 1], [0, 1, 0])


# ============================================================
# Workflow tests
# ============================================================


class TestWorkflows:
    """Tests for workflow functions."""

    def test_batch_entropy_analysis(self) -> None:
        from metainformant.information.workflow.workflows import batch_entropy_analysis

        sequences = ["ATCGATCG", "GCTAGCTA", "AAAAAAAA"]
        result = batch_entropy_analysis(sequences, k=1)
        assert "batch_info" in result
        assert "sequence_results" in result
        assert "summary" in result
        assert len(result["sequence_results"]) == 3

    def test_batch_entropy_analysis_output(self, tmp_path: Path) -> None:
        from metainformant.information.workflow.workflows import batch_entropy_analysis

        sequences = ["ATCGATCG", "GCTAGCTA"]
        result = batch_entropy_analysis(sequences, k=1, output_dir=tmp_path)
        assert (tmp_path / "batch_entropy_analysis.json").exists()

    def test_information_workflow(self) -> None:
        from metainformant.information.workflow.workflows import information_workflow

        # Sequences must be >= 10 chars for analyze_sequence_information
        sequences = ["ATCGATCGATCG", "GCTAGCTAGCTA", "ATCGATCGATCG"]
        result = information_workflow(sequences, k_values=[1, 2])
        assert result["workflow_status"] == "completed"
        assert len(result["sequence_analyses"]) == 3

    def test_compare_datasets(self) -> None:
        from metainformant.information.workflow.workflows import compare_datasets

        ds1 = ["ATCGATCG", "GCTAGCTA"]
        ds2 = ["AAAAAAAA", "TTTTTTTT"]
        result = compare_datasets(ds1, ds2, k=1)
        assert "comparisons" in result
        assert "dataset_info" in result

    def test_information_report_markdown(self, tmp_path: Path) -> None:
        from metainformant.information.workflow.workflows import information_report, information_workflow

        # Sequences must be >= 10 chars for analyze_sequence_information
        sequences = ["ATCGATCGATCG", "GCTAGCTAGCTA"]
        results = information_workflow(sequences, k_values=[1])
        output_path = tmp_path / "report.md"
        information_report(results, output_path=output_path, format="markdown")
        assert output_path.exists()
        content = output_path.read_text()
        assert "Information Analysis Report" in content

    def test_information_report_json(self, tmp_path: Path) -> None:
        from metainformant.information.workflow.workflows import information_report, information_workflow

        # Sequences must be >= 10 chars for analyze_sequence_information
        sequences = ["ATCGATCGATCG", "GCTAGCTAGCTA"]
        results = information_workflow(sequences, k_values=[1])
        output_path = tmp_path / "report.json"
        information_report(results, output_path=output_path, format="json")
        assert output_path.exists()


# ============================================================
# Integration tests
# ============================================================


class TestIntegrationFunctions:
    """Tests for cross-module integration functions."""

    def test_dna_integration(self) -> None:
        from metainformant.information.integration.integration import dna_integration

        sequences = ["ATCGATCG", "GCTAGCTA"]
        result = dna_integration(sequences)
        assert result["n_sequences"] == 2
        assert "integrated_metrics" in result

    def test_rna_integration(self) -> None:
        from metainformant.information.integration.integration import rna_integration

        rng = np.random.default_rng(42)
        expr = rng.random((10, 5))
        result = rna_integration(expr)
        assert "n_genes" in result
        assert "integrated_metrics" in result

    def test_multiomics_integration(self) -> None:
        from metainformant.information.integration.integration import multiomics_integration

        rng = np.random.default_rng(42)
        omics = {
            "genomics": rng.random((10, 5)),
            "transcriptomics": rng.random((10, 5)),
        }
        result = multiomics_integration(omics_data=omics)
        assert "omics_types" in result
        assert "genomics" in result["omics_types"]

    def test_ml_integration(self) -> None:
        from metainformant.information.integration.integration import ml_integration

        rng = np.random.default_rng(42)
        X = rng.random((50, 10))
        y = rng.integers(0, 2, 50)
        result = ml_integration(X, y, method="mutual_info")
        assert "feature_scores" in result


# ============================================================
# Network tests
# ============================================================


class TestNetworkInformation:
    """Tests for network information analysis."""

    @pytest.fixture
    def small_graph(self):
        try:
            import networkx as nx
        except ImportError:
            pytest.skip("networkx not available")
        G = nx.Graph()
        G.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)])
        return G

    def test_network_entropy_shannon(self, small_graph) -> None:
        from metainformant.information.integration.networks import network_entropy

        h = network_entropy(small_graph, method="shannon")
        assert h >= 0.0

    def test_network_entropy_von_neumann(self, small_graph) -> None:
        from metainformant.information.integration.networks import network_entropy

        h = network_entropy(small_graph, method="von_neumann")
        assert isinstance(h, (float, complex))

    def test_network_entropy_from_adjacency(self) -> None:
        from metainformant.information.integration.networks import network_entropy

        try:
            import networkx  # noqa: F401
        except ImportError:
            pytest.skip("networkx not available")
        adj = np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]])
        h = network_entropy(adj)
        assert h >= 0.0

    def test_information_flow_random_walk(self, small_graph) -> None:
        from metainformant.information.integration.networks import information_flow

        result = information_flow(small_graph, method="random_walk", steps=10)
        assert "flow_matrix" in result
        assert result["method"] == "random_walk"

    def test_information_flow_diffusion(self, small_graph) -> None:
        from metainformant.information.integration.networks import information_flow

        result = information_flow(small_graph, method="diffusion", steps=10)
        assert "node_flows" in result
        assert result["method"] == "diffusion"

    def test_network_information_centrality(self, small_graph) -> None:
        from metainformant.information.integration.networks import network_information_centrality

        centrality = network_information_centrality(small_graph, method="entropy")
        assert len(centrality) == 4
        for v in centrality.values():
            assert 0.0 <= v <= 1.0

    def test_network_motif_information(self, small_graph) -> None:
        from metainformant.information.integration.networks import network_motif_information

        result = network_motif_information(small_graph, motif_size=3, n_random=5)
        assert result["status"] == "completed"
        assert "motifs_found" in result
        assert "z_scores" in result
        assert "motif_entropy" in result

    def test_information_graph_distance_entropy(self) -> None:
        from metainformant.information.integration.networks import information_graph_distance

        try:
            import networkx as nx
        except ImportError:
            pytest.skip("networkx not available")
        G1 = nx.cycle_graph(5)
        G2 = nx.complete_graph(5)
        d = information_graph_distance(G1, G2, method="entropy")
        assert d >= 0.0

    def test_information_graph_distance_jsd(self) -> None:
        from metainformant.information.integration.networks import information_graph_distance

        try:
            import networkx as nx
        except ImportError:
            pytest.skip("networkx not available")
        G1 = nx.cycle_graph(5)
        G2 = nx.complete_graph(5)
        d = information_graph_distance(G1, G2, method="jensen_shannon")
        assert d >= 0.0


# ============================================================
# Module import and export tests
# ============================================================


class TestModuleExports:
    """Tests that all exported names are accessible."""

    def test_all_syntactic_exports(self) -> None:
        from metainformant.information.metrics.core.syntactic import conditional_entropy, conditional_mutual_information, cross_entropy, information_coefficient, jensen_shannon_divergence, joint_entropy, kl_divergence, mutual_information, normalized_mutual_information, renyi_entropy, shannon_entropy, shannon_entropy_from_counts, total_correlation, transfer_entropy, tsallis_entropy

        assert callable(shannon_entropy)
        assert callable(mutual_information)

    def test_all_semantic_exports(self) -> None:
        from metainformant.information.metrics.advanced.semantic import annotation_specificity, information_content, information_content_from_annotations, ontology_complexity, semantic_distance, semantic_entropy, semantic_similarity, semantic_similarity_matrix, term_redundancy, term_specificity

        assert callable(information_content)

    def test_all_continuous_exports(self) -> None:
        from metainformant.information.metrics.core.continuous import conditional_entropy_continuous, copula_entropy, differential_entropy, entropy_estimation, information_flow_network, kl_divergence_continuous, mutual_information_continuous, transfer_entropy_continuous

        assert callable(differential_entropy)

    def test_all_estimation_exports(self) -> None:
        from metainformant.information.metrics.core.estimation import bias_correction, entropy_bootstrap_confidence, entropy_estimator, entropy_rate_estimator, kl_divergence_estimator, mutual_information_estimator, panzeri_treves_bias_correction

        assert callable(entropy_estimator)

    def test_all_advanced_exports(self) -> None:
        from metainformant.information.metrics.analysis import binding_information, fisher_information, fisher_information_matrix, interaction_information, lautum_information, relative_information_gain, variation_of_information

        assert callable(fisher_information)

    def test_all_network_exports(self) -> None:
        from metainformant.information.integration.networks import information_community_detection, information_flow, information_graph_distance, network_entropy, network_information_centrality, network_motif_information

        assert callable(network_entropy)

    def test_all_workflow_exports(self) -> None:
        from metainformant.information.workflow.workflows import batch_entropy_analysis, compare_datasets, information_report, information_workflow

        assert callable(batch_entropy_analysis)
