"""Tests for information geometry and partial information decomposition modules.

Tests cover all functions in:
- metainformant.information.metrics.geometry
- metainformant.information.metrics.decomposition
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from metainformant.information.metrics.geometry import (
    channel_capacity,
    entropy_power_inequality,
    exponential_family_entropy,
    fisher_rao_distance,
    hellinger_distance,
    information_bottleneck,
    information_dimension,
    information_projection,
    natural_gradient,
    rate_distortion_function,
    statistical_divergence,
)
from metainformant.information.metrics.decomposition import (
    co_information,
    dual_total_correlation,
    o_information,
    partial_information_decomposition,
    redundant_information,
    synergistic_information,
    unique_information,
)


# ============================================================
# Fisher-Rao Distance Tests
# ============================================================


class TestFisherRaoDistance:
    """Tests for Fisher-Rao geodesic distance."""

    def test_identical_distributions(self) -> None:
        """Distance between identical distributions is 0."""
        p = [0.25, 0.25, 0.25, 0.25]
        assert fisher_rao_distance(p, p) == pytest.approx(0.0, abs=1e-10)

    def test_different_distributions(self) -> None:
        """Distance between different distributions is positive."""
        p = [0.5, 0.5]
        q = [0.9, 0.1]
        d = fisher_rao_distance(p, q)
        assert d > 0.0

    def test_symmetric(self) -> None:
        """Fisher-Rao distance is symmetric."""
        p = [0.3, 0.7]
        q = [0.6, 0.4]
        assert fisher_rao_distance(p, q) == pytest.approx(fisher_rao_distance(q, p), abs=1e-10)

    def test_triangle_inequality(self) -> None:
        """Fisher-Rao satisfies the triangle inequality."""
        p = [0.2, 0.8]
        q = [0.5, 0.5]
        r = [0.9, 0.1]
        d_pq = fisher_rao_distance(p, q)
        d_qr = fisher_rao_distance(q, r)
        d_pr = fisher_rao_distance(p, r)
        assert d_pr <= d_pq + d_qr + 1e-10

    def test_mismatched_lengths_raises(self) -> None:
        """Mismatched distribution lengths raise ValueError."""
        with pytest.raises(ValueError, match="same length"):
            fisher_rao_distance([0.5, 0.5], [0.3, 0.3, 0.4])


# ============================================================
# Natural Gradient Tests
# ============================================================


class TestNaturalGradient:
    """Tests for natural gradient computation."""

    def test_identity_fim(self) -> None:
        """With identity FIM, natural gradient equals Euclidean gradient."""
        grad = np.array([1.0, 2.0, 3.0])
        fim = np.eye(3)
        nat_grad = natural_gradient(grad, fim)
        np.testing.assert_allclose(nat_grad, grad, atol=1e-10)

    def test_scaled_fim(self) -> None:
        """With scaled identity FIM, natural gradient is inverse-scaled."""
        grad = np.array([2.0, 4.0])
        fim = 2.0 * np.eye(2)
        nat_grad = natural_gradient(grad, fim)
        np.testing.assert_allclose(nat_grad, [1.0, 2.0], atol=1e-10)

    def test_dimension_mismatch_raises(self) -> None:
        """Dimension mismatch raises ValueError."""
        with pytest.raises(ValueError, match="Dimension mismatch"):
            natural_gradient(np.array([1.0, 2.0]), np.eye(3))


# ============================================================
# Information Projection Tests
# ============================================================


class TestInformationProjection:
    """Tests for information projection (m-projection)."""

    def test_basic_projection(self) -> None:
        """Projection produces a valid distribution."""
        p = [0.25, 0.25, 0.25, 0.25]
        constraints = [([0, 1], [0.5, 0.5])]
        result = information_projection(p, constraints)
        assert "projected" in result
        assert "kl_divergence" in result
        assert "n_iterations" in result
        proj = result["projected"]
        assert pytest.approx(sum(proj), abs=1e-6) == 1.0
        assert all(x >= 0 for x in proj)

    def test_kl_divergence_non_negative(self) -> None:
        """KL divergence of projection is non-negative."""
        p = [0.1, 0.2, 0.3, 0.4]
        constraints = [([0, 1], [0.4, 0.4])]
        result = information_projection(p, constraints)
        assert result["kl_divergence"] >= 0.0

    def test_empty_constraints_raises(self) -> None:
        """Empty constraint set raises ValueError."""
        with pytest.raises(ValueError, match="must not be empty"):
            information_projection([0.5, 0.5], [])


# ============================================================
# Statistical (Alpha) Divergence Tests
# ============================================================


class TestStatisticalDivergence:
    """Tests for alpha-divergence."""

    def test_identical_distributions(self) -> None:
        """Alpha-divergence between identical distributions is 0."""
        p = [0.3, 0.7]
        assert statistical_divergence(p, p, alpha=0.5) == pytest.approx(0.0, abs=1e-10)

    def test_positive_for_different(self) -> None:
        """Alpha-divergence is positive for different distributions."""
        p = [0.5, 0.5]
        q = [0.9, 0.1]
        assert statistical_divergence(p, q, alpha=0.5) > 0.0

    def test_different_alphas(self) -> None:
        """Alpha-divergence varies with alpha parameter."""
        p = [0.3, 0.7]
        q = [0.6, 0.4]
        d1 = statistical_divergence(p, q, alpha=0.3)
        d2 = statistical_divergence(p, q, alpha=0.7)
        # Both should be positive but may differ
        assert d1 > 0.0
        assert d2 > 0.0


# ============================================================
# Exponential Family Entropy Tests
# ============================================================


class TestExponentialFamilyEntropy:
    """Tests for exponential family entropy."""

    def test_gaussian_entropy(self) -> None:
        """Gaussian N(0,1): H = 0.5 * ln(2*pi*e) using 2-parameter exponential form."""
        # 2-parameter Gaussian: eta = (mu/sigma^2, -1/(2*sigma^2)), T(x) = (x, x^2)
        # For N(0,1): eta = (0, -0.5), E[T(x)] = (0, 1)
        # A(eta) = -eta1^2/(4*eta2) - 0.5*ln(-2*eta2) + 0.5*ln(2*pi)
        # A = 0 - 0.5*ln(1) + 0.5*ln(2*pi) = 0.5*ln(2*pi)
        eta = [0.0, -0.5]
        e_t = [0.0, 1.0]
        a_eta = 0.5 * math.log(2 * math.pi)
        h = exponential_family_entropy(eta, e_t, a_eta)
        expected = 0.5 * math.log(2 * math.pi * math.e)
        assert h == pytest.approx(expected, abs=1e-6)

    def test_mismatched_lengths_raises(self) -> None:
        """Mismatched parameter lengths raise ValueError."""
        with pytest.raises(ValueError, match="same length"):
            exponential_family_entropy([1.0, 2.0], [1.0], 0.5)


# ============================================================
# Hellinger Distance Tests
# ============================================================


class TestHellingerDistance:
    """Tests for Hellinger distance."""

    def test_identical_distributions(self) -> None:
        """Hellinger distance between identical distributions is 0."""
        p = [0.25, 0.25, 0.25, 0.25]
        assert hellinger_distance(p, p) == pytest.approx(0.0, abs=1e-10)

    def test_range_zero_to_one(self) -> None:
        """Hellinger distance is in [0, 1]."""
        p = [0.1, 0.9]
        q = [0.8, 0.2]
        h = hellinger_distance(p, q)
        assert 0.0 <= h <= 1.0

    def test_symmetric(self) -> None:
        """Hellinger distance is symmetric."""
        p = [0.3, 0.7]
        q = [0.6, 0.4]
        assert hellinger_distance(p, q) == pytest.approx(hellinger_distance(q, p), abs=1e-10)

    def test_disjoint_support(self) -> None:
        """Hellinger distance for disjoint support is 1."""
        p = [1.0, 0.0]
        q = [0.0, 1.0]
        assert hellinger_distance(p, q) == pytest.approx(1.0, abs=1e-10)


# ============================================================
# Channel Capacity (Blahut-Arimoto) Tests
# ============================================================


class TestChannelCapacity:
    """Tests for channel capacity via Blahut-Arimoto."""

    def test_binary_symmetric_channel(self) -> None:
        """BSC(p): C = 1 - H(p) where H is binary entropy."""
        p_error = 0.1
        transition = np.array([
            [1 - p_error, p_error],
            [p_error, 1 - p_error],
        ])
        result = channel_capacity(transition)
        assert "capacity" in result
        assert "optimal_input" in result
        assert "n_iterations" in result
        # Capacity should be positive
        assert result["capacity"] > 0.0

    def test_noiseless_channel(self) -> None:
        """Noiseless channel has capacity = log(|X|) nats."""
        transition = np.eye(3)  # Perfect 3-symbol channel
        result = channel_capacity(transition)
        expected = math.log(3)
        assert result["capacity"] == pytest.approx(expected, rel=0.01)

    def test_optimal_input_is_distribution(self) -> None:
        """Optimal input sums to 1 and is non-negative."""
        transition = np.array([[0.8, 0.2], [0.3, 0.7]])
        result = channel_capacity(transition)
        p = result["optimal_input"]
        assert np.all(p >= 0)
        assert np.sum(p) == pytest.approx(1.0, abs=1e-6)


# ============================================================
# Rate-Distortion Tests
# ============================================================


class TestRateDistortion:
    """Tests for rate-distortion function."""

    def test_binary_source_hamming(self) -> None:
        """Binary source with Hamming distortion produces valid R(D) curve."""
        source_probs = [0.5, 0.5]
        distortion = np.array([[0, 1], [1, 0]], dtype=float)  # Hamming
        result = rate_distortion_function(source_probs, distortion, n_points=20)
        assert "rates" in result
        assert "distortions" in result
        assert "critical_distortion" in result
        assert len(result["rates"]) == 20
        # All rates should be non-negative
        assert all(r >= 0 for r in result["rates"])

    def test_dimension_mismatch_raises(self) -> None:
        """Mismatched source_probs and distortion_matrix raises ValueError."""
        with pytest.raises(ValueError, match="rows"):
            rate_distortion_function([0.5, 0.5], np.zeros((3, 2)))


# ============================================================
# Information Bottleneck Tests
# ============================================================


class TestInformationBottleneck:
    """Tests for Tishby's information bottleneck."""

    def test_basic_joint_distribution(self) -> None:
        """Information bottleneck on a simple joint distribution."""
        # Create a block-diagonal joint distribution
        joint = np.array([
            [0.2, 0.05, 0.0, 0.0],
            [0.05, 0.2, 0.0, 0.0],
            [0.0, 0.0, 0.2, 0.05],
            [0.0, 0.0, 0.05, 0.2],
        ])
        result = information_bottleneck(joint, beta=10.0, n_clusters=2)
        assert "compression_rate" in result
        assert "relevance" in result
        assert "cluster_assignments" in result
        assert "n_iterations" in result
        assert result["compression_rate"] >= 0.0
        assert result["relevance"] >= 0.0
        # Should find 2 clusters
        assert len(set(result["cluster_assignments"])) <= 2

    def test_invalid_beta_raises(self) -> None:
        """Negative beta raises ValueError."""
        joint = np.ones((2, 2)) / 4
        with pytest.raises(ValueError, match="beta must be positive"):
            information_bottleneck(joint, beta=-1.0)


# ============================================================
# Entropy Power Inequality Tests
# ============================================================


class TestEntropyPowerInequality:
    """Tests for entropy power inequality."""

    def test_gaussian_variances(self) -> None:
        """EPI with Gaussian variances: N(X) = sigma^2."""
        result = entropy_power_inequality([1.0, 2.0, 3.0])
        assert result["entropy_powers"] == [1.0, 2.0, 3.0]
        assert result["sum_entropy_power"] == pytest.approx(6.0)
        assert result["epi_bound"] == pytest.approx(6.0)
        assert result["epi_satisfied"] is True

    def test_single_variable(self) -> None:
        """Single variable still works."""
        result = entropy_power_inequality([5.0])
        assert result["entropy_powers"] == [5.0]
        assert result["epi_satisfied"] is True

    def test_non_positive_variance_raises(self) -> None:
        """Non-positive variance raises ValueError."""
        with pytest.raises(ValueError, match="positive"):
            entropy_power_inequality([1.0, -0.5])


# ============================================================
# Information Dimension Tests
# ============================================================


class TestInformationDimension:
    """Tests for Grassberger-Procaccia correlation dimension."""

    def test_2d_gaussian_cloud(self) -> None:
        """2D Gaussian cloud should have dimension close to 2."""
        rng = np.random.default_rng(42)
        samples = rng.normal(size=(500, 2))
        result = information_dimension(samples)
        assert "dimension" in result
        assert "r_squared" in result
        assert "r_values" in result
        assert "correlation_integral" in result
        # Dimension should be approximately 2
        assert 1.0 < result["dimension"] < 3.0

    def test_1d_line(self) -> None:
        """1D data should have dimension close to 1."""
        rng = np.random.default_rng(42)
        samples = rng.normal(size=(200, 1))
        result = information_dimension(samples)
        assert 0.5 < result["dimension"] < 1.5

    def test_too_few_samples_raises(self) -> None:
        """Fewer than 10 samples raises ValueError."""
        with pytest.raises(ValueError, match="at least 10"):
            information_dimension(np.array([[1, 2], [3, 4]]))


# ============================================================
# Co-Information Tests
# ============================================================


class TestCoInformation:
    """Tests for co-information (interaction information)."""

    def test_two_variables_equals_mi(self) -> None:
        """For two variables, co-information equals mutual information."""
        rng = np.random.default_rng(42)
        x = list(rng.choice(["A", "B"], size=200))
        y = x.copy()  # Perfect correlation
        ci = co_information([x, y])
        from metainformant.information import mutual_information

        mi = mutual_information(x, y)
        assert ci == pytest.approx(mi, abs=0.01)

    def test_independent_variables(self) -> None:
        """Independent variables have co-information near 0."""
        rng = np.random.default_rng(42)
        x = list(rng.choice([0, 1], size=500))
        y = list(rng.choice([0, 1], size=500))
        ci = co_information([x, y])
        assert abs(ci) < 0.1

    def test_three_redundant_variables(self) -> None:
        """Three redundant copies have positive co-information."""
        rng = np.random.default_rng(42)
        x = list(rng.choice([0, 1], size=500))
        ci = co_information([x, x, x])
        assert ci > 0.0

    def test_minimum_variables_raises(self) -> None:
        """Fewer than 2 variables raises ValueError."""
        with pytest.raises(ValueError, match="At least 2"):
            co_information([[1, 2, 3]])


# ============================================================
# Partial Information Decomposition Tests
# ============================================================


class TestPartialInformationDecomposition:
    """Tests for PID (MMI approach)."""

    def test_redundant_sources(self) -> None:
        """When both sources are identical, all info is redundant."""
        rng = np.random.default_rng(42)
        target = list(rng.choice([0, 1], size=500))
        # Both sources are copies of target
        result = partial_information_decomposition(target, target, target)
        assert result["redundancy"] > 0.0
        assert result["unique_x"] == pytest.approx(0.0, abs=0.01)
        assert result["unique_y"] == pytest.approx(0.0, abs=0.01)

    def test_xor_synergy(self) -> None:
        """XOR pattern: jointly informative but individually uninformative."""
        rng = np.random.default_rng(42)
        x = list(rng.choice([0, 1], size=500))
        y = list(rng.choice([0, 1], size=500))
        target = [a ^ b for a, b in zip(x, y)]
        result = partial_information_decomposition(x, y, target)
        # Synergy should be present
        assert result["total"] > 0.0
        assert "synergy" in result

    def test_pid_components_sum_to_total(self) -> None:
        """Redundancy + unique_x + unique_y + synergy ~= total."""
        rng = np.random.default_rng(42)
        x = list(rng.choice([0, 1, 2], size=300))
        y = list(rng.choice([0, 1, 2], size=300))
        target = list(rng.choice([0, 1], size=300))
        result = partial_information_decomposition(x, y, target)
        component_sum = result["redundancy"] + result["unique_x"] + result["unique_y"] + result["synergy"]
        assert component_sum == pytest.approx(result["total"], abs=0.05)

    def test_empty_raises(self) -> None:
        """Empty sequences raise ValueError."""
        with pytest.raises(ValueError, match="empty"):
            partial_information_decomposition([], [], [])


# ============================================================
# Unique Information Tests
# ============================================================


class TestUniqueInformation:
    """Tests for unique information."""

    def test_source_carries_unique_info(self) -> None:
        """Source has unique information when it predicts target better."""
        rng = np.random.default_rng(42)
        target = list(rng.choice([0, 1], size=500))
        source = target.copy()  # Perfect predictor
        other = list(rng.choice([0, 1], size=500))  # Random
        ui = unique_information(source, target, other)
        assert ui > 0.0

    def test_non_negative(self) -> None:
        """Unique information is non-negative."""
        rng = np.random.default_rng(42)
        x = list(rng.choice([0, 1], size=200))
        y = list(rng.choice([0, 1], size=200))
        t = list(rng.choice([0, 1], size=200))
        ui = unique_information(x, t, y)
        assert ui >= 0.0


# ============================================================
# Redundant Information Tests
# ============================================================


class TestRedundantInformation:
    """Tests for redundant information."""

    def test_identical_sources(self) -> None:
        """Identical sources share all information redundantly."""
        rng = np.random.default_rng(42)
        target = list(rng.choice([0, 1], size=500))
        red = redundant_information([target, target], target)
        # Redundancy should equal MI between target and itself
        from metainformant.information import mutual_information

        mi = mutual_information(target, target)
        assert red == pytest.approx(mi, abs=0.01)

    def test_independent_sources(self) -> None:
        """Independent random sources have low redundancy about random target."""
        rng = np.random.default_rng(42)
        x = list(rng.choice([0, 1], size=500))
        y = list(rng.choice([0, 1], size=500))
        target = list(rng.choice([0, 1], size=500))
        red = redundant_information([x, y], target)
        assert red >= 0.0

    def test_minimum_sources_raises(self) -> None:
        """Fewer than 2 sources raises ValueError."""
        with pytest.raises(ValueError, match="At least 2"):
            redundant_information([[1, 2]], [1, 2])


# ============================================================
# Synergistic Information Tests
# ============================================================


class TestSynergisticInformation:
    """Tests for synergistic information."""

    def test_xor_has_synergy(self) -> None:
        """XOR pattern should produce synergistic information."""
        rng = np.random.default_rng(42)
        x = list(rng.choice([0, 1], size=500))
        y = list(rng.choice([0, 1], size=500))
        target = [a ^ b for a, b in zip(x, y)]
        syn = synergistic_information([x, y], target)
        assert syn >= 0.0

    def test_non_negative(self) -> None:
        """Synergistic information is non-negative."""
        rng = np.random.default_rng(42)
        x = list(rng.choice([0, 1], size=200))
        y = list(rng.choice([0, 1], size=200))
        t = list(rng.choice([0, 1], size=200))
        syn = synergistic_information([x, y], t)
        assert syn >= 0.0


# ============================================================
# Dual Total Correlation Tests
# ============================================================


class TestDualTotalCorrelation:
    """Tests for dual total correlation."""

    def test_independent_variables(self) -> None:
        """Independent variables have DTC near 0."""
        rng = np.random.default_rng(42)
        x = list(rng.choice([0, 1], size=500))
        y = list(rng.choice([0, 1], size=500))
        dtc = dual_total_correlation([x, y])
        assert abs(dtc) < 0.1

    def test_identical_variables(self) -> None:
        """Identical variables have positive DTC."""
        rng = np.random.default_rng(42)
        x = list(rng.choice([0, 1], size=500))
        dtc = dual_total_correlation([x, x])
        assert dtc > 0.0

    def test_non_negative(self) -> None:
        """DTC is non-negative."""
        rng = np.random.default_rng(42)
        x = list(rng.choice([0, 1, 2], size=300))
        y = list(rng.choice([0, 1, 2], size=300))
        z = list(rng.choice([0, 1, 2], size=300))
        dtc = dual_total_correlation([x, y, z])
        assert dtc >= 0.0


# ============================================================
# O-Information Tests
# ============================================================


class TestOInformation:
    """Tests for O-information."""

    def test_redundant_system(self) -> None:
        """Copies of the same variable have positive O-information (redundancy)."""
        rng = np.random.default_rng(42)
        x = list(rng.choice([0, 1], size=500))
        o = o_information([x, x, x])
        assert o > 0.0

    def test_independent_variables(self) -> None:
        """Independent variables have O-information near 0."""
        rng = np.random.default_rng(42)
        x = list(rng.choice([0, 1], size=500))
        y = list(rng.choice([0, 1], size=500))
        o = o_information([x, y])
        assert abs(o) < 0.1

    def test_minimum_variables_raises(self) -> None:
        """Fewer than 2 variables raises ValueError."""
        with pytest.raises(ValueError, match="At least 2"):
            o_information([[1, 2, 3]])
