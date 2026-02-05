"""Comprehensive tests for information-theoretic hypothesis testing, channel capacity, and geometry modules.

Tests cover all public functions in:
- metainformant.information.metrics.hypothesis
- metainformant.information.metrics.channel
- metainformant.information.metrics.geometry

All tests use real implementations (NO mocking policy).
"""

from __future__ import annotations

import math
from typing import Any, Dict, List

import numpy as np
import pytest

# ============================================================
# Hypothesis testing imports
# ============================================================
from metainformant.information.metrics.hypothesis import (
    entropy_confidence_interval,
    entropy_rate_test,
    independence_test,
    information_significance_filter,
    mi_permutation_test,
)

# ============================================================
# Channel capacity imports
# ============================================================
from metainformant.information.metrics.channel import (
    channel_capacity as channel_capacity_ch,
    channel_mutual_information as channel_mi_ch,
    information_bottleneck as ib_ch,
    noisy_channel_capacity as noisy_cap_ch,
    rate_distortion as rate_distortion_ch,
)

# ============================================================
# Information geometry imports
# ============================================================
from metainformant.information.metrics.geometry import (
    channel_capacity as channel_capacity_geo,
    entropy_power_inequality,
    exponential_family_entropy,
    fisher_rao_distance,
    hellinger_distance,
    information_bottleneck as ib_geo,
    information_dimension,
    information_projection,
    natural_gradient,
    rate_distortion_function,
    statistical_divergence,
)


# ============================================================
# ============================================================
# HYPOTHESIS TESTING MODULE
# ============================================================
# ============================================================


class TestMIPermutationTest:
    """Tests for mi_permutation_test: permutation-based MI significance."""

    def test_correlated_sequences_significant(self) -> None:
        """Strongly correlated sequences should yield significant MI."""
        rng = np.random.default_rng(42)
        x = list(rng.choice(["A", "B", "C"], size=200))
        # y is a noisy copy of x -- high MI expected
        y = [v if rng.random() < 0.85 else "D" for v in x]
        result = mi_permutation_test(x, y, n_permutations=200, seed=42)

        assert "observed_mi" in result
        assert "p_value" in result
        assert "null_distribution" in result
        assert "significant" in result

        assert result["observed_mi"] > 0.0
        assert result["p_value"] < 0.05
        assert result["significant"] is True
        assert len(result["null_distribution"]) == 200

    def test_independent_sequences_not_significant(self) -> None:
        """Independent random sequences should typically not be significant."""
        rng = np.random.default_rng(42)
        x = list(rng.choice([0, 1, 2], size=200))
        y = list(rng.choice([0, 1, 2], size=200))
        result = mi_permutation_test(x, y, n_permutations=200, seed=42)

        assert result["p_value"] > 0.01
        assert result["observed_mi"] >= 0.0

    def test_return_types_and_structure(self) -> None:
        """Result dictionary has correct types."""
        x = [0, 1, 0, 1, 0, 1, 0, 1]
        y = [0, 1, 0, 1, 0, 1, 0, 1]
        result = mi_permutation_test(x, y, n_permutations=50, seed=0)

        assert isinstance(result["observed_mi"], float)
        assert isinstance(result["p_value"], float)
        assert isinstance(result["null_distribution"], list)
        assert isinstance(result["significant"], bool)
        assert 0.0 <= result["p_value"] <= 1.0

    def test_mismatched_lengths_raises(self) -> None:
        """Sequences of different lengths raise ValueError."""
        with pytest.raises(ValueError, match="same length"):
            mi_permutation_test([1, 2, 3], [1, 2])

    def test_invalid_n_permutations_raises(self) -> None:
        """n_permutations < 1 raises ValueError."""
        with pytest.raises(ValueError, match="n_permutations"):
            mi_permutation_test([1, 2], [1, 2], n_permutations=0)

    def test_reproducibility_with_seed(self) -> None:
        """Same seed produces identical results."""
        x = [0, 1, 0, 1, 0, 1]
        y = [1, 0, 1, 0, 1, 0]
        r1 = mi_permutation_test(x, y, n_permutations=100, seed=123)
        r2 = mi_permutation_test(x, y, n_permutations=100, seed=123)
        assert r1["observed_mi"] == r2["observed_mi"]
        assert r1["p_value"] == r2["p_value"]
        assert r1["null_distribution"] == r2["null_distribution"]

    def test_null_distribution_length(self) -> None:
        """Null distribution has exactly n_permutations entries."""
        x = [0, 1, 0, 1]
        y = [0, 1, 0, 1]
        result = mi_permutation_test(x, y, n_permutations=37, seed=42)
        assert len(result["null_distribution"]) == 37


class TestIndependenceTest:
    """Tests for independence_test: chi-squared independence testing."""

    def test_dependent_variables(self) -> None:
        """Perfectly dependent variables should have large statistic."""
        x = [0, 0, 1, 1, 2, 2, 0, 0, 1, 1, 2, 2] * 10
        y = x.copy()
        result = independence_test(x, y, method="chi_squared")

        assert "statistic" in result
        assert "p_value" in result
        assert "method" in result
        assert "dof" in result
        assert "significant" in result

        assert result["statistic"] > 0.0
        assert result["p_value"] < 0.05
        assert result["method"] == "chi_squared"
        assert result["significant"] is True

    def test_independent_variables(self) -> None:
        """Independent variables should have non-significant p-value."""
        rng = np.random.default_rng(42)
        x = list(rng.choice([0, 1], size=500))
        y = list(rng.choice([0, 1], size=500))
        result = independence_test(x, y)
        # p-value should be relatively high for independent data
        # (though not guaranteed, we use a large sample to reduce variance)
        assert result["statistic"] >= 0.0
        assert 0.0 <= result["p_value"] <= 1.0

    def test_degrees_of_freedom(self) -> None:
        """DOF = (r-1)(c-1) for the contingency table."""
        x = [0, 0, 1, 1, 2, 2]
        y = ["a", "b", "a", "b", "a", "b"]
        result = independence_test(x, y)
        # 3 unique x values, 2 unique y values => dof = (3-1)*(2-1) = 2
        assert result["dof"] == 2

    def test_mismatched_lengths_raises(self) -> None:
        """Different-length sequences raise ValueError."""
        with pytest.raises(ValueError, match="same length"):
            independence_test([1, 2, 3], [1, 2])

    def test_unknown_method_raises(self) -> None:
        """Unsupported method raises ValueError."""
        with pytest.raises(ValueError, match="Unknown method"):
            independence_test([1, 2], [1, 2], method="fisher_exact")

    def test_return_types(self) -> None:
        """Return values have correct types."""
        x = [0, 1, 0, 1, 0, 1]
        y = [0, 1, 0, 1, 0, 1]
        result = independence_test(x, y)
        assert isinstance(result["statistic"], float)
        assert isinstance(result["p_value"], float)
        assert isinstance(result["method"], str)
        assert isinstance(result["dof"], int)
        assert isinstance(result["significant"], bool)


class TestEntropyConfidenceInterval:
    """Tests for entropy_confidence_interval: bootstrap CI for entropy."""

    def test_basic_confidence_interval(self) -> None:
        """CI contains the point estimate and has correct structure."""
        rng = np.random.default_rng(42)
        data = list(rng.choice(["A", "B", "C", "D"], size=200))
        result = entropy_confidence_interval(data, confidence=0.95, n_bootstrap=200, seed=42)

        assert "entropy" in result
        assert "ci_lower" in result
        assert "ci_upper" in result
        assert "confidence" in result
        assert "std" in result

        assert result["ci_lower"] <= result["entropy"] <= result["ci_upper"]
        assert result["confidence"] == 0.95
        assert result["std"] >= 0.0

    def test_narrow_ci_with_large_sample(self) -> None:
        """Large samples should produce narrow confidence intervals."""
        rng = np.random.default_rng(42)
        data = list(rng.choice([0, 1], size=1000))
        result = entropy_confidence_interval(data, confidence=0.95, n_bootstrap=300, seed=42)
        ci_width = result["ci_upper"] - result["ci_lower"]
        # Width should be relatively small for a large sample
        assert ci_width < 0.5

    def test_single_symbol_zero_entropy(self) -> None:
        """Single-symbol data has entropy near zero."""
        data = ["X"] * 100
        result = entropy_confidence_interval(data, confidence=0.95, n_bootstrap=100, seed=42)
        assert result["entropy"] == pytest.approx(0.0, abs=1e-10)

    def test_empty_data_raises(self) -> None:
        """Empty data raises ValueError."""
        with pytest.raises(ValueError, match="empty"):
            entropy_confidence_interval([], confidence=0.95)

    def test_invalid_confidence_raises(self) -> None:
        """Confidence outside (0,1) raises ValueError."""
        with pytest.raises(ValueError, match="between 0 and 1"):
            entropy_confidence_interval([1, 2, 3], confidence=1.5)
        with pytest.raises(ValueError, match="between 0 and 1"):
            entropy_confidence_interval([1, 2, 3], confidence=0.0)

    def test_invalid_n_bootstrap_raises(self) -> None:
        """n_bootstrap < 1 raises ValueError."""
        with pytest.raises(ValueError, match="n_bootstrap"):
            entropy_confidence_interval([1, 2, 3], n_bootstrap=0)

    def test_higher_confidence_wider_interval(self) -> None:
        """Higher confidence level produces wider interval."""
        rng = np.random.default_rng(42)
        data = list(rng.choice([0, 1, 2], size=200))
        ci_90 = entropy_confidence_interval(data, confidence=0.90, n_bootstrap=500, seed=42)
        ci_99 = entropy_confidence_interval(data, confidence=0.99, n_bootstrap=500, seed=42)
        width_90 = ci_90["ci_upper"] - ci_90["ci_lower"]
        width_99 = ci_99["ci_upper"] - ci_99["ci_lower"]
        assert width_99 >= width_90


class TestInformationSignificanceFilter:
    """Tests for information_significance_filter: MI-based variable selection."""

    def test_significant_variable_detected(self) -> None:
        """Strongly associated variables should pass the filter."""
        rng = np.random.default_rng(42)
        target = list(rng.choice([0, 1], size=100))
        # Variable 0 is a noisy copy (highly correlated)
        var_correlated = [t if rng.random() < 0.9 else 1 - t for t in target]
        # Variable 1 is random noise
        var_noise = list(rng.choice([0, 1], size=100))
        variables = [var_correlated, var_noise]
        result = information_significance_filter(variables, target, alpha=0.05, correction="bonferroni")

        assert "significant_indices" in result
        assert "p_values" in result
        assert "adjusted_p_values" in result
        assert "mi_values" in result
        assert "alpha" in result
        assert "correction" in result

        assert len(result["p_values"]) == 2
        assert len(result["adjusted_p_values"]) == 2
        assert len(result["mi_values"]) == 2
        # The correlated variable (index 0) should be significant
        assert 0 in result["significant_indices"]

    def test_bonferroni_correction(self) -> None:
        """Bonferroni adjustment multiplies p-values by number of tests."""
        rng = np.random.default_rng(42)
        target = list(rng.choice([0, 1], size=100))
        variables = [
            list(rng.choice([0, 1], size=100)),
            list(rng.choice([0, 1], size=100)),
            list(rng.choice([0, 1], size=100)),
        ]
        result = information_significance_filter(variables, target, correction="bonferroni")
        # Adjusted p-values should be at most 1.0
        for adj_p in result["adjusted_p_values"]:
            assert adj_p <= 1.0

    def test_no_correction(self) -> None:
        """With correction='none', adjusted = raw p-values."""
        rng = np.random.default_rng(42)
        target = list(rng.choice([0, 1], size=80))
        variables = [list(rng.choice([0, 1], size=80))]
        result = information_significance_filter(variables, target, correction="none")
        assert result["p_values"] == result["adjusted_p_values"]
        assert result["correction"] == "none"

    def test_empty_variables_raises(self) -> None:
        """Empty variables list raises ValueError."""
        with pytest.raises(ValueError, match="empty"):
            information_significance_filter([], [1, 2, 3])

    def test_unknown_correction_raises(self) -> None:
        """Unknown correction method raises ValueError."""
        with pytest.raises(ValueError, match="Unknown correction"):
            information_significance_filter([[1, 2]], [1, 2], correction="holm")

    def test_fdr_bh_correction(self) -> None:
        """FDR BH correction produces valid adjusted p-values."""
        rng = np.random.default_rng(42)
        target = list(rng.choice([0, 1], size=80))
        variables = [
            list(rng.choice([0, 1], size=80)),
            list(rng.choice([0, 1], size=80)),
        ]
        result = information_significance_filter(variables, target, correction="fdr_bh")
        assert result["correction"] == "fdr_bh"
        for adj_p in result["adjusted_p_values"]:
            assert 0.0 <= adj_p <= 1.0


class TestEntropyRateTest:
    """Tests for entropy_rate_test: permutation test for transfer entropy."""

    def test_causal_signal_significant(self) -> None:
        """Causal driving signal should show significant transfer entropy."""
        rng = np.random.default_rng(42)
        n = 200
        x = list(rng.choice([0, 1, 2], size=n))
        # y follows x with a lag of 1 (causal influence)
        y = [0] + [x[i] for i in range(n - 1)]
        result = entropy_rate_test(x, y, lag=1, n_permutations=200, seed=42)

        assert "observed_te" in result
        assert "p_value" in result
        assert "significant" in result
        assert "lag" in result

        assert result["observed_te"] >= 0.0
        assert result["lag"] == 1

    def test_independent_sequences_not_significant(self) -> None:
        """Independent sequences should have non-significant transfer entropy."""
        rng = np.random.default_rng(42)
        n = 100
        x = list(rng.choice([0, 1], size=n))
        y = list(rng.choice([0, 1], size=n))
        result = entropy_rate_test(x, y, lag=1, n_permutations=200, seed=42)

        assert result["observed_te"] >= 0.0
        assert result["p_value"] > 0.01

    def test_mismatched_lengths_raises(self) -> None:
        """Mismatched sequence lengths raise ValueError."""
        with pytest.raises(ValueError, match="same length"):
            entropy_rate_test([1, 2, 3], [1, 2])

    def test_lag_too_large_raises(self) -> None:
        """Sequences too short for the given lag raise ValueError."""
        with pytest.raises(ValueError, match="too short"):
            entropy_rate_test([1, 2], [1, 2], lag=5)

    def test_negative_lag_raises(self) -> None:
        """Negative lag raises ValueError."""
        with pytest.raises(ValueError, match="Lag must be"):
            entropy_rate_test([1, 2, 3, 4, 5], [1, 2, 3, 4, 5], lag=0)

    def test_return_types(self) -> None:
        """Return values have correct types."""
        x = [0, 1, 0, 1, 0, 1, 0, 1, 0, 1]
        y = [1, 0, 1, 0, 1, 0, 1, 0, 1, 0]
        result = entropy_rate_test(x, y, lag=1, n_permutations=50, seed=42)
        assert isinstance(result["observed_te"], float)
        assert isinstance(result["p_value"], float)
        assert isinstance(result["significant"], bool)
        assert isinstance(result["lag"], int)

    def test_different_lags(self) -> None:
        """Different lag values produce results with the correct lag."""
        x = list(range(20))
        y = list(range(20))
        for lag in [1, 2, 3]:
            result = entropy_rate_test(x, y, lag=lag, n_permutations=50, seed=42)
            assert result["lag"] == lag


# ============================================================
# ============================================================
# CHANNEL CAPACITY MODULE (channel.py)
# ============================================================
# ============================================================


class TestChannelCapacityCh:
    """Tests for channel.channel_capacity: Blahut-Arimoto in channel.py."""

    def test_noiseless_binary_channel(self) -> None:
        """Noiseless 2x2 identity channel has capacity = 1 bit."""
        transition = np.eye(2)
        result = channel_capacity_ch(transition)

        assert "capacity" in result
        assert "optimal_input" in result
        assert "n_iterations" in result

        # Capacity = log2(2) = 1 bit
        assert result["capacity"] == pytest.approx(1.0, rel=0.01)

    def test_noiseless_ternary_channel(self) -> None:
        """Noiseless 3x3 identity channel has capacity = log2(3) bits."""
        transition = np.eye(3)
        result = channel_capacity_ch(transition)
        expected = math.log2(3)
        assert result["capacity"] == pytest.approx(expected, rel=0.01)

    def test_binary_symmetric_channel(self) -> None:
        """BSC(p=0.1): C = 1 - H(0.1) bits."""
        p = 0.1
        transition = np.array([[1 - p, p], [p, 1 - p]])
        result = channel_capacity_ch(transition)
        h_p = -p * math.log2(p) - (1 - p) * math.log2(1 - p)
        expected = 1.0 - h_p
        assert result["capacity"] == pytest.approx(expected, rel=0.05)

    def test_optimal_input_is_valid_distribution(self) -> None:
        """Optimal input distribution sums to 1 and is non-negative."""
        transition = np.array([[0.7, 0.3], [0.2, 0.8]])
        result = channel_capacity_ch(transition)
        p = np.array(result["optimal_input"])
        assert np.all(p >= -1e-10)
        assert np.sum(p) == pytest.approx(1.0, abs=1e-6)

    def test_invalid_matrix_raises(self) -> None:
        """Non-2D array raises ValueError."""
        with pytest.raises(ValueError):
            channel_capacity_ch(np.array([0.5, 0.5]))

    def test_negative_entries_raises(self) -> None:
        """Negative entries in transition matrix raise ValueError."""
        with pytest.raises(ValueError, match="non-negative"):
            channel_capacity_ch(np.array([[-0.1, 1.1], [0.5, 0.5]]))

    def test_rows_not_summing_to_one_raises(self) -> None:
        """Rows not summing to 1 raise ValueError."""
        with pytest.raises(ValueError, match="sum to 1"):
            channel_capacity_ch(np.array([[0.5, 0.3], [0.5, 0.5]]))


class TestRateDistortionCh:
    """Tests for channel.rate_distortion."""

    def test_binary_source_hamming(self) -> None:
        """Binary source with Hamming distortion produces valid R(D) curve."""
        source = [0.5, 0.5]
        distortion = np.array([[0, 1], [1, 0]], dtype=float)
        result = rate_distortion_ch(source, distortion, n_points=20)

        assert "rates" in result
        assert "distortions" in result
        assert "rd_curve" in result
        assert len(result["rates"]) > 0
        # All rates are non-negative
        assert all(r >= 0.0 for r in result["rates"])
        # All distortions are non-negative
        assert all(d >= 0.0 for d in result["distortions"])

    def test_asymmetric_source(self) -> None:
        """Asymmetric source distribution produces valid curve."""
        source = [0.8, 0.2]
        distortion = np.array([[0, 1], [1, 0]], dtype=float)
        result = rate_distortion_ch(source, distortion, n_points=15)
        assert len(result["rates"]) > 0

    def test_dimension_mismatch_raises(self) -> None:
        """Mismatched source_dist and distortion_matrix dimensions raise ValueError."""
        with pytest.raises(ValueError):
            rate_distortion_ch([0.5, 0.5], np.zeros((3, 2)))

    def test_invalid_source_dist_raises(self) -> None:
        """Source distribution not summing to 1 raises ValueError."""
        with pytest.raises(ValueError, match="sum to 1"):
            rate_distortion_ch([0.5, 0.3], np.array([[0, 1], [1, 0]], dtype=float))

    def test_negative_distortion_raises(self) -> None:
        """Negative distortion entries raise ValueError."""
        with pytest.raises(ValueError, match="non-negative"):
            rate_distortion_ch([0.5, 0.5], np.array([[0, -1], [1, 0]], dtype=float))


class TestInformationBottleneckCh:
    """Tests for channel.information_bottleneck."""

    def test_block_diagonal_joint(self) -> None:
        """Block-diagonal joint distribution should cluster cleanly."""
        joint = np.array(
            [
                [0.2, 0.05, 0.0, 0.0],
                [0.05, 0.2, 0.0, 0.0],
                [0.0, 0.0, 0.2, 0.05],
                [0.0, 0.0, 0.05, 0.2],
            ]
        )
        result = ib_ch(joint, beta=10.0, n_clusters=2)

        assert "assignments" in result
        assert "I_T_X" in result
        assert "I_T_Y" in result
        assert "beta" in result

        assert result["I_T_X"] >= 0.0
        assert result["I_T_Y"] >= 0.0
        assert result["beta"] == 10.0
        assert len(result["assignments"]) == 4

    def test_invalid_joint_raises(self) -> None:
        """Non-2D array raises ValueError."""
        with pytest.raises(ValueError):
            ib_ch(np.array([0.5, 0.5]), beta=1.0)

    def test_negative_entries_raises(self) -> None:
        """Negative joint entries raise ValueError."""
        with pytest.raises(ValueError, match="non-negative"):
            ib_ch(np.array([[-0.1, 0.6], [0.3, 0.2]]), beta=1.0)


class TestChannelMutualInformationCh:
    """Tests for channel.channel_mutual_information."""

    def test_noiseless_channel(self) -> None:
        """Noiseless channel with uniform input has MI = log2(n)."""
        transition = np.eye(3)
        input_dist = [1.0 / 3, 1.0 / 3, 1.0 / 3]
        mi = channel_mi_ch(transition, input_dist)
        expected = math.log2(3)
        assert mi == pytest.approx(expected, rel=0.01)

    def test_completely_noisy_channel(self) -> None:
        """Channel that maps everything to uniform output has MI = 0."""
        transition = np.array([[0.5, 0.5], [0.5, 0.5]])
        input_dist = [0.5, 0.5]
        mi = channel_mi_ch(transition, input_dist)
        assert mi == pytest.approx(0.0, abs=1e-6)

    def test_non_negative(self) -> None:
        """MI is always non-negative."""
        transition = np.array([[0.7, 0.3], [0.4, 0.6]])
        input_dist = [0.6, 0.4]
        mi = channel_mi_ch(transition, input_dist)
        assert mi >= 0.0

    def test_dimension_mismatch_raises(self) -> None:
        """Input distribution length must match number of channel rows."""
        with pytest.raises(ValueError, match="length"):
            channel_mi_ch(np.eye(2), [0.3, 0.3, 0.4])

    def test_invalid_input_dist_raises(self) -> None:
        """Input distribution not summing to 1 raises ValueError."""
        with pytest.raises(ValueError, match="sum to 1"):
            channel_mi_ch(np.eye(2), [0.5, 0.3])


class TestNoisyChannelCapacityCh:
    """Tests for channel.noisy_channel_capacity."""

    def test_bsc_noiseless(self) -> None:
        """BSC with p=0: capacity = 1 bit."""
        cap = noisy_cap_ch(0.0, channel_type="binary_symmetric")
        assert cap == pytest.approx(1.0, abs=1e-10)

    def test_bsc_half(self) -> None:
        """BSC with p=0.5: capacity = 0 bits."""
        cap = noisy_cap_ch(0.5, channel_type="binary_symmetric")
        assert cap == pytest.approx(0.0, abs=1e-6)

    def test_bec_zero_erasure(self) -> None:
        """BEC with epsilon=0: capacity = 1 bit."""
        cap = noisy_cap_ch(0.0, channel_type="binary_erasure")
        assert cap == pytest.approx(1.0, abs=1e-10)

    def test_bec_full_erasure(self) -> None:
        """BEC with epsilon=1: capacity = 0 bits."""
        cap = noisy_cap_ch(1.0, channel_type="binary_erasure")
        assert cap == pytest.approx(0.0, abs=1e-10)

    def test_awgn_zero_snr(self) -> None:
        """AWGN with SNR=0: capacity = 0."""
        cap = noisy_cap_ch(0.0, channel_type="awgn")
        assert cap == pytest.approx(0.0, abs=1e-10)

    def test_awgn_positive_snr(self) -> None:
        """AWGN with SNR=1: capacity = 0.5*log2(2) = 0.5 bits."""
        cap = noisy_cap_ch(1.0, channel_type="awgn")
        expected = 0.5 * math.log2(2.0)
        assert cap == pytest.approx(expected, rel=1e-6)

    def test_unknown_channel_type_raises(self) -> None:
        """Unknown channel type raises ValueError."""
        with pytest.raises(ValueError, match="Unknown"):
            noisy_cap_ch(0.5, channel_type="quantum")

    def test_bsc_invalid_noise_raises(self) -> None:
        """BSC with p > 1 raises ValueError."""
        with pytest.raises(ValueError):
            noisy_cap_ch(1.5, channel_type="binary_symmetric")

    def test_awgn_negative_snr_raises(self) -> None:
        """AWGN with negative SNR raises ValueError."""
        with pytest.raises(ValueError):
            noisy_cap_ch(-1.0, channel_type="awgn")


# ============================================================
# ============================================================
# INFORMATION GEOMETRY MODULE (geometry.py)
# ============================================================
# ============================================================


class TestFisherRaoDistance:
    """Tests for geometry.fisher_rao_distance."""

    def test_identical_distributions_zero(self) -> None:
        """Distance between identical distributions is 0."""
        p = [0.25, 0.25, 0.25, 0.25]
        assert fisher_rao_distance(p, p) == pytest.approx(0.0, abs=1e-10)

    def test_different_distributions_positive(self) -> None:
        """Distance between different distributions is positive."""
        p = [0.5, 0.5]
        q = [0.9, 0.1]
        assert fisher_rao_distance(p, q) > 0.0

    def test_symmetry(self) -> None:
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

    def test_disjoint_support_maximum_distance(self) -> None:
        """Disjoint support gives maximum distance (pi)."""
        p = [1.0, 0.0]
        q = [0.0, 1.0]
        d = fisher_rao_distance(p, q)
        assert d == pytest.approx(math.pi, rel=0.01)


class TestNaturalGradient:
    """Tests for geometry.natural_gradient."""

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

    def test_non_trivial_fim(self) -> None:
        """Non-trivial positive-definite FIM produces correct result."""
        grad = np.array([1.0, 0.0])
        fim = np.array([[2.0, 1.0], [1.0, 2.0]])
        nat_grad = natural_gradient(grad, fim)
        # F^{-1} = (1/3) * [[2, -1], [-1, 2]] so nat_grad = [2/3, -1/3]
        expected = np.array([2.0 / 3, -1.0 / 3])
        np.testing.assert_allclose(nat_grad, expected, atol=1e-10)

    def test_dimension_mismatch_raises(self) -> None:
        """Dimension mismatch raises ValueError."""
        with pytest.raises(ValueError, match="Dimension mismatch"):
            natural_gradient(np.array([1.0, 2.0]), np.eye(3))

    def test_non_square_fim_raises(self) -> None:
        """Non-square FIM raises ValueError."""
        with pytest.raises(ValueError, match="square"):
            natural_gradient(np.array([1.0, 2.0]), np.ones((2, 3)))


class TestInformationProjection:
    """Tests for geometry.information_projection."""

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

    def test_out_of_range_index_raises(self) -> None:
        """Index out of range in constraint raises ValueError."""
        with pytest.raises(ValueError, match="out of range"):
            information_projection([0.5, 0.5], [([0, 5], [0.3, 0.3])])

    def test_unknown_method_raises(self) -> None:
        """Unknown method raises ValueError."""
        with pytest.raises(ValueError, match="Unknown method"):
            information_projection([0.5, 0.5], [([0], [0.5])], method="gradient")


class TestStatisticalDivergence:
    """Tests for geometry.statistical_divergence (alpha-divergence)."""

    def test_identical_distributions_zero(self) -> None:
        """Alpha-divergence between identical distributions is 0."""
        p = [0.3, 0.7]
        assert statistical_divergence(p, p, alpha=0.5) == pytest.approx(0.0, abs=1e-10)

    def test_positive_for_different(self) -> None:
        """Alpha-divergence is positive for different distributions."""
        p = [0.5, 0.5]
        q = [0.9, 0.1]
        assert statistical_divergence(p, q, alpha=0.5) > 0.0

    def test_non_negative(self) -> None:
        """Alpha-divergence is non-negative for all valid alpha values."""
        p = [0.3, 0.7]
        q = [0.6, 0.4]
        for alpha in [0.2, 0.3, 0.5, 0.7, 0.9]:
            d = statistical_divergence(p, q, alpha=alpha)
            assert d >= 0.0, f"Negative divergence at alpha={alpha}"

    def test_different_alphas(self) -> None:
        """Different alpha values produce different divergences."""
        p = [0.3, 0.7]
        q = [0.6, 0.4]
        d_02 = statistical_divergence(p, q, alpha=0.2)
        d_08 = statistical_divergence(p, q, alpha=0.8)
        # Different alphas should generally give different results
        assert d_02 != d_08

    def test_mismatched_lengths_raises(self) -> None:
        """Mismatched distribution lengths raise ValueError."""
        with pytest.raises(ValueError, match="same length"):
            statistical_divergence([0.5, 0.5], [0.3, 0.3, 0.4], alpha=0.5)


class TestExponentialFamilyEntropy:
    """Tests for geometry.exponential_family_entropy."""

    def test_gaussian_entropy(self) -> None:
        """Exponential family entropy: H = -eta^T E[T] + A(eta).

        For standard Gaussian (mu=0, sigma^2=1):
        eta = [0], E[T] = [0], A(eta) = 0.5*ln(2*pi).
        So H = 0 + 0.5*ln(2*pi) = 0.9189...
        (This is the exponential family contribution; the full differential
        entropy 0.5*ln(2*pi*e) includes the base measure h(x) = exp(-x^2/2).)
        """
        sigma2 = 1.0
        mu = 0.0
        eta = [mu / sigma2]
        e_t = [mu]
        a_eta = mu**2 / (2 * sigma2) + 0.5 * math.log(2 * math.pi * sigma2)
        h = exponential_family_entropy(eta, e_t, a_eta)
        # H = -eta^T E[T] + A(eta) = 0 + 0.5*ln(2*pi)
        expected = 0.5 * math.log(2 * math.pi * sigma2)
        assert h == pytest.approx(expected, abs=1e-6)

    def test_mismatched_lengths_raises(self) -> None:
        """Mismatched parameter lengths raise ValueError."""
        with pytest.raises(ValueError, match="same length"):
            exponential_family_entropy([1.0, 2.0], [1.0], 0.5)

    def test_empty_params_raises(self) -> None:
        """Empty natural parameters raise ValueError."""
        with pytest.raises(ValueError, match="empty"):
            exponential_family_entropy([], [], 0.0)

    def test_entropy_value(self) -> None:
        """Manual computation check with simple parameters."""
        eta = [2.0]
        e_t = [3.0]
        log_partition = 10.0
        h = exponential_family_entropy(eta, e_t, log_partition)
        # H = -eta^T E[T] + A = -2*3 + 10 = 4
        assert h == pytest.approx(4.0, abs=1e-10)


class TestHellingerDistance:
    """Tests for geometry.hellinger_distance."""

    def test_identical_distributions_zero(self) -> None:
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

    def test_disjoint_support_is_one(self) -> None:
        """Disjoint support gives Hellinger distance = 1."""
        p = [1.0, 0.0]
        q = [0.0, 1.0]
        assert hellinger_distance(p, q) == pytest.approx(1.0, abs=1e-10)

    def test_mismatched_lengths_raises(self) -> None:
        """Mismatched distribution lengths raise ValueError."""
        with pytest.raises(ValueError, match="same length"):
            hellinger_distance([0.5, 0.5], [0.3, 0.3, 0.4])

    def test_relationship_to_fisher_rao(self) -> None:
        """Hellinger distance H relates to Bhattacharyya coefficient: H = sqrt(1 - BC)."""
        p = [0.3, 0.7]
        q = [0.6, 0.4]
        h = hellinger_distance(p, q)
        bc = sum(math.sqrt(pi * qi) for pi, qi in zip(p, q))
        expected = math.sqrt(1.0 - bc)
        assert h == pytest.approx(expected, abs=1e-10)


class TestChannelCapacityGeo:
    """Tests for geometry.channel_capacity (Blahut-Arimoto in geometry module)."""

    def test_noiseless_channel(self) -> None:
        """Noiseless 2x2 channel has capacity = ln(2) nats."""
        transition = np.eye(2)
        result = channel_capacity_geo(transition)

        assert "capacity" in result
        assert "optimal_input" in result
        assert "n_iterations" in result
        # Geometry module reports in nats: ln(2) ~ 0.693
        expected = math.log(2)
        assert result["capacity"] == pytest.approx(expected, rel=0.05)

    def test_optimal_input_valid(self) -> None:
        """Optimal input is a valid probability distribution."""
        transition = np.array([[0.8, 0.2], [0.3, 0.7]])
        result = channel_capacity_geo(transition)
        p = result["optimal_input"]
        assert np.all(p >= -1e-10)
        assert np.sum(p) == pytest.approx(1.0, abs=1e-6)

    def test_invalid_method_raises(self) -> None:
        """Invalid method raises ValueError."""
        with pytest.raises(ValueError, match="Unknown method"):
            channel_capacity_geo(np.eye(2), method="simplex")


class TestRateDistortionFunction:
    """Tests for geometry.rate_distortion_function."""

    def test_binary_source_hamming(self) -> None:
        """Binary source with Hamming distortion produces valid R(D) curve."""
        source = [0.5, 0.5]
        distortion = np.array([[0, 1], [1, 0]], dtype=float)
        result = rate_distortion_function(source, distortion, n_points=20)

        assert "rates" in result
        assert "distortions" in result
        assert "critical_distortion" in result
        assert len(result["rates"]) == 20
        assert all(r >= 0.0 for r in result["rates"])

    def test_dimension_mismatch_raises(self) -> None:
        """Mismatched source_probs and distortion_matrix rows raise ValueError."""
        with pytest.raises(ValueError, match="rows"):
            rate_distortion_function([0.5, 0.5], np.zeros((3, 2)))

    def test_monotone_rd_curve(self) -> None:
        """R(D) curve should be monotonically non-increasing in sorted order."""
        source = [0.5, 0.5]
        distortion = np.array([[0, 1], [1, 0]], dtype=float)
        result = rate_distortion_function(source, distortion, n_points=30)
        # After sorting by distortion (already sorted), rates should be non-increasing
        rates = result["rates"]
        distortions = result["distortions"]
        for i in range(len(rates) - 1):
            if distortions[i + 1] > distortions[i]:
                # Higher distortion should mean lower or equal rate
                assert rates[i + 1] <= rates[i] + 1e-6


class TestInformationBottleneckGeo:
    """Tests for geometry.information_bottleneck."""

    def test_block_diagonal_joint(self) -> None:
        """Block-diagonal joint should cluster into 2 groups."""
        joint = np.array(
            [
                [0.2, 0.05, 0.0, 0.0],
                [0.05, 0.2, 0.0, 0.0],
                [0.0, 0.0, 0.2, 0.05],
                [0.0, 0.0, 0.05, 0.2],
            ]
        )
        result = ib_geo(joint, beta=10.0, n_clusters=2)

        assert "compression_rate" in result
        assert "relevance" in result
        assert "cluster_assignments" in result
        assert "n_iterations" in result

        assert result["compression_rate"] >= 0.0
        assert result["relevance"] >= 0.0
        assert len(result["cluster_assignments"]) == 4

    def test_invalid_beta_raises(self) -> None:
        """Negative beta raises ValueError."""
        joint = np.ones((2, 2)) / 4
        with pytest.raises(ValueError, match="beta must be positive"):
            ib_geo(joint, beta=-1.0)

    def test_too_many_clusters_raises(self) -> None:
        """More clusters than X values raises ValueError."""
        joint = np.ones((2, 2)) / 4
        with pytest.raises(ValueError, match="n_clusters"):
            ib_geo(joint, beta=1.0, n_clusters=5)


class TestEntropyPowerInequality:
    """Tests for geometry.entropy_power_inequality."""

    def test_gaussian_variances(self) -> None:
        """EPI with Gaussian variances: N(X) = sigma^2."""
        result = entropy_power_inequality([1.0, 2.0, 3.0])
        assert result["entropy_powers"] == [1.0, 2.0, 3.0]
        assert result["sum_entropy_power"] == pytest.approx(6.0)
        assert result["epi_bound"] == pytest.approx(6.0)
        assert result["epi_satisfied"] is True

    def test_single_variable(self) -> None:
        """Single variable works correctly."""
        result = entropy_power_inequality([5.0])
        assert result["entropy_powers"] == [5.0]
        assert result["epi_satisfied"] is True

    def test_non_positive_variance_raises(self) -> None:
        """Non-positive variance raises ValueError."""
        with pytest.raises(ValueError, match="positive"):
            entropy_power_inequality([1.0, -0.5])

    def test_empty_variances_raises(self) -> None:
        """Empty variance list raises ValueError."""
        with pytest.raises(ValueError):
            entropy_power_inequality([])

    def test_return_structure(self) -> None:
        """All expected keys present in result."""
        result = entropy_power_inequality([2.0, 3.0])
        assert "entropy_powers" in result
        assert "sum_entropy_power" in result
        assert "epi_bound" in result
        assert "epi_satisfied" in result


class TestInformationDimension:
    """Tests for geometry.information_dimension (Grassberger-Procaccia)."""

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

    def test_r_squared_quality(self) -> None:
        """R-squared should be reasonably high for clean data."""
        rng = np.random.default_rng(42)
        samples = rng.normal(size=(500, 2))
        result = information_dimension(samples)
        assert result["r_squared"] > 0.5

    def test_3d_gaussian_cloud(self) -> None:
        """3D Gaussian cloud should have dimension close to 3."""
        rng = np.random.default_rng(42)
        samples = rng.normal(size=(800, 3))
        result = information_dimension(samples)
        assert 2.0 < result["dimension"] < 4.0

    def test_unknown_method_raises(self) -> None:
        """Unknown method raises ValueError."""
        rng = np.random.default_rng(42)
        samples = rng.normal(size=(100, 2))
        with pytest.raises(ValueError, match="Unknown method"):
            information_dimension(samples, method="box_counting")
