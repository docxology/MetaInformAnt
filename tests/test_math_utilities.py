"""Tests for math module utility functions.

Tests the utility functions defined in metainformant.math.__init__.py:
- correlation_coefficient
- linear_regression
- fisher_exact_test
- shannon_entropy
- jensen_shannon_divergence
"""

from __future__ import annotations

import math
import pytest
import numpy as np

from metainformant.math import (
    correlation_coefficient,
    linear_regression,
    fisher_exact_test,
    shannon_entropy,
    jensen_shannon_divergence,
)


class TestCorrelationCoefficient:
    """Test correlation_coefficient function."""

    def test_perfect_positive_correlation(self):
        """Test perfect positive correlation."""
        x = [1, 2, 3, 4, 5]
        y = [2, 4, 6, 8, 10]  # y = 2x
        result = correlation_coefficient(x, y)
        assert abs(result - 1.0) < 1e-10

    def test_perfect_negative_correlation(self):
        """Test perfect negative correlation."""
        x = [1, 2, 3, 4, 5]
        y = [10, 8, 6, 4, 2]  # y = 10 - 2x
        result = correlation_coefficient(x, y)
        assert abs(result - (-1.0)) < 1e-10

    def test_no_correlation(self):
        """Test uncorrelated data."""
        x = [1, 2, 3, 4, 5]
        y = [1, 1, 1, 1, 1]  # constant y
        result = correlation_coefficient(x, y)
        assert result == 0.0  # Returns 0.0 for zero variance

    def test_mismatched_lengths(self):
        """Test error with mismatched list lengths."""
        x = [1, 2, 3]
        y = [1, 2]
        with pytest.raises(ValueError, match="same length"):
            correlation_coefficient(x, y)

    def test_insufficient_data(self):
        """Test error with insufficient data points."""
        x = [1]
        y = [1]
        with pytest.raises(ValueError, match="at least 2 elements"):
            correlation_coefficient(x, y)

    def test_empty_lists(self):
        """Test error with empty lists."""
        with pytest.raises(ValueError, match="at least 2 elements"):
            correlation_coefficient([], [])

    def test_known_values(self):
        """Test with known correlation values."""
        x = [1, 2, 3, 4, 5]
        y = [1, 3, 5, 7, 9]  # y = 2x - 1
        result = correlation_coefficient(x, y)
        assert abs(result - 1.0) < 1e-10

    def test_large_dataset(self):
        """Test with larger dataset."""
        np.random.seed(42)
        x = np.random.normal(0, 1, 1000)
        y = 2 * x + np.random.normal(0, 0.1, 1000)  # strong correlation with noise
        result = correlation_coefficient(x.tolist(), y.tolist())
        assert result > 0.9  # Should be highly correlated


class TestLinearRegression:
    """Test linear_regression function."""

    def test_perfect_fit(self):
        """Test perfect linear fit."""
        x = [1, 2, 3, 4, 5]
        y = [2, 4, 6, 8, 10]  # y = 2x
        slope, intercept, r_squared = linear_regression(x, y)
        assert abs(slope - 2.0) < 1e-10
        assert abs(intercept - 0.0) < 1e-10
        assert abs(r_squared - 1.0) < 1e-10

    def test_offset_fit(self):
        """Test linear fit with offset."""
        x = [1, 2, 3, 4, 5]
        y = [3, 5, 7, 9, 11]  # y = 2x + 1
        slope, intercept, r_squared = linear_regression(x, y)
        assert abs(slope - 2.0) < 1e-10
        assert abs(intercept - 1.0) < 1e-10
        assert abs(r_squared - 1.0) < 1e-10

    def test_noisy_fit(self):
        """Test linear fit with noise."""
        np.random.seed(42)
        x = np.array([1, 2, 3, 4, 5])
        true_slope, true_intercept = 2.5, -1.0
        y = true_slope * x + true_intercept + np.random.normal(0, 0.1, len(x))
        slope, intercept, r_squared = linear_regression(x.tolist(), y.tolist())
        assert abs(slope - true_slope) < 0.1  # Should be close
        assert abs(intercept - true_intercept) < 0.1  # Should be close
        assert r_squared > 0.99  # Should have high R²

    def test_zero_variance_x(self):
        """Test with zero variance in x."""
        x = [2, 2, 2, 2, 2]  # constant x
        y = [1, 2, 3, 4, 5]
        slope, intercept, r_squared = linear_regression(x, y)
        # Should return fallback values
        assert slope == 0.0
        assert intercept == 0.0
        assert r_squared == 0.0

    def test_mismatched_lengths(self):
        """Test error with mismatched list lengths."""
        x = [1, 2, 3]
        y = [1, 2]
        with pytest.raises(ValueError, match="same length"):
            linear_regression(x, y)

    def test_insufficient_data(self):
        """Test error with insufficient data points."""
        x = [1]
        y = [1]
        with pytest.raises(ValueError, match="at least 2 elements"):
            linear_regression(x, y)

    def test_empty_lists(self):
        """Test error with empty lists."""
        with pytest.raises(ValueError, match="at least 2 elements"):
            linear_regression([], [])


class TestFisherExactTest:
    """Test fisher_exact_test function."""

    def test_basic_2x2_table(self):
        """Test basic 2x2 contingency table."""
        # Example from scipy docs
        table = [[10, 2], [3, 15]]
        odds_ratio, p_value = fisher_exact_test(10, 2, 3, 15)
        assert odds_ratio > 0  # Should be positive
        assert 0 <= p_value <= 1  # p-value should be in valid range

    def test_known_result(self):
        """Test with known expected result."""
        # Table: [[5, 0], [0, 5]]
        # Perfect association
        odds_ratio, p_value = fisher_exact_test(5, 0, 0, 5)
        assert math.isinf(odds_ratio) or odds_ratio > 1000  # Very large odds ratio
        assert p_value < 0.01  # Should be significant

    def test_no_association(self):
        """Test with no association expected."""
        # Balanced table with no association
        odds_ratio, p_value = fisher_exact_test(10, 10, 10, 10)
        assert abs(odds_ratio - 1.0) < 0.1  # Should be close to 1
        assert p_value > 0.05  # Should not be significant

    def test_scipy_dependency_error(self):
        """Test behavior when scipy is not available."""
        # This test would need to mock scipy import failure
        # For now, just verify the function works normally
        odds_ratio, p_value = fisher_exact_test(1, 1, 1, 1)
        assert isinstance(odds_ratio, float)
        assert isinstance(p_value, float)

    def test_edge_cases(self):
        """Test edge cases."""
        # Zero in numerator (odds ratio = 0)
        odds_ratio, p_value = fisher_exact_test(0, 1, 1, 1)
        assert odds_ratio == 0.0
        assert 0 <= p_value <= 1

        # Very small numbers
        odds_ratio, p_value = fisher_exact_test(1, 1, 1, 1)
        assert odds_ratio > 0
        assert 0 <= p_value <= 1


class TestShannonEntropy:
    """Test shannon_entropy function."""

    def test_uniform_distribution(self):
        """Test uniform distribution (maximum entropy)."""
        values = [0.25, 0.25, 0.25, 0.25]
        entropy = shannon_entropy(values)
        assert abs(entropy - 2.0) < 1e-10  # log2(4) = 2

    def test_deterministic_distribution(self):
        """Test deterministic distribution (zero entropy)."""
        values = [1.0, 0.0, 0.0]
        entropy = shannon_entropy(values)
        assert abs(entropy - 0.0) < 1e-10

    def test_binary_distribution(self):
        """Test binary distribution."""
        values = [0.5, 0.5]
        entropy = shannon_entropy(values)
        assert abs(entropy - 1.0) < 1e-10  # log2(2) = 1

    def test_normalized_automatically(self):
        """Test that function normalizes unnormalized inputs."""
        values = [2, 2, 2]  # Sums to 6, should be normalized to [1/3, 1/3, 1/3]
        entropy = shannon_entropy(values)
        expected = -3 * (1/3) * math.log2(1/3)  # Should be log2(3)
        assert abs(entropy - math.log2(3)) < 1e-10

    def test_empty_list(self):
        """Test with empty list."""
        entropy = shannon_entropy([])
        assert entropy == 0.0

    def test_single_nonzero_value(self):
        """Test with single non-zero value."""
        entropy = shannon_entropy([5.0])
        assert entropy == 0.0  # Deterministic

    def test_all_zeros(self):
        """Test with all zero values."""
        entropy = shannon_entropy([0, 0, 0])
        assert entropy == 0.0

    def test_large_distribution(self):
        """Test with larger distribution."""
        n = 100
        values = [1/n] * n  # Uniform over n possibilities
        entropy = shannon_entropy(values)
        expected = math.log2(n)
        assert abs(entropy - expected) < 1e-10


class TestJensenShannonDivergence:
    """Test jensen_shannon_divergence function."""

    def test_identical_distributions(self):
        """Test identical distributions (zero divergence)."""
        p = [0.5, 0.3, 0.2]
        q = [0.5, 0.3, 0.2]
        divergence = jensen_shannon_divergence(p, q)
        assert abs(divergence - 0.0) < 1e-10

    def test_different_distributions(self):
        """Test different distributions."""
        p = [1.0, 0.0]  # Deterministic
        q = [0.5, 0.5]  # Uniform
        divergence = jensen_shannon_divergence(p, q)
        assert divergence > 0  # Should be positive
        assert divergence <= 1.0  # Should be bounded

    def test_maximum_divergence(self):
        """Test maximum possible divergence."""
        p = [1.0, 0.0, 0.0]
        q = [0.0, 0.0, 1.0]
        divergence = jensen_shannon_divergence(p, q)
        # JS divergence is bounded by 1 in log base 2
        assert 0 <= divergence <= 1.0

    def test_symmetry(self):
        """Test symmetry property."""
        p = [0.7, 0.2, 0.1]
        q = [0.3, 0.4, 0.3]
        divergence_pq = jensen_shannon_divergence(p, q)
        divergence_qp = jensen_shannon_divergence(q, p)
        assert abs(divergence_pq - divergence_qp) < 1e-10

    def test_symmetry_only(self):
        """Test that JS divergence is symmetric."""
        p = [0.5, 0.3, 0.2]
        q = [0.3, 0.4, 0.3]

        divergence_pq = jensen_shannon_divergence(p, q)
        divergence_qp = jensen_shannon_divergence(q, p)

        assert abs(divergence_pq - divergence_qp) < 1e-10

    def test_normalization(self):
        """Test that function handles unnormalized inputs."""
        p = [3, 1, 2]  # Will be normalized to [0.5, 0.166..., 0.333...]
        q = [1, 1, 4]  # Will be normalized to [0.166..., 0.166..., 0.666...]
        divergence = jensen_shannon_divergence(p, q)
        assert 0 <= divergence <= 1.0

    def test_mismatched_lengths(self):
        """Test error with mismatched list lengths."""
        p = [0.5, 0.3, 0.2]
        q = [0.6, 0.4]
        with pytest.raises(ValueError, match="same length"):
            jensen_shannon_divergence(p, q)

    def test_empty_lists(self):
        """Test with empty lists."""
        result = jensen_shannon_divergence([], [])
        assert result == 0.0  # Returns 0.0 for empty distributions

    def test_large_distributions(self):
        """Test with larger distributions."""
        np.random.seed(42)
        p = np.random.dirichlet([1] * 10)  # Random 10-element distribution
        q = np.random.dirichlet([1] * 10)
        divergence = jensen_shannon_divergence(p.tolist(), q.tolist())
        assert 0 <= divergence <= 1.0