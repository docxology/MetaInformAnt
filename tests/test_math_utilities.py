"""Comprehensive tests for math module utility functions.

Tests cover all utility functions defined in math/__init__.py:
- correlation_coefficient
- linear_regression
- fisher_exact_test
- shannon_entropy
- jensen_shannon_divergence
"""

from __future__ import annotations

import math

import pytest

from metainformant.math import (
    correlation_coefficient,
    fisher_exact_test,
    jensen_shannon_divergence,
    linear_regression,
    shannon_entropy,
)


class TestCorrelationCoefficient:
    """Tests for correlation_coefficient function."""

    def test_perfect_positive_correlation(self):
        """Test perfect positive correlation (r = 1.0)."""
        x = [1.0, 2.0, 3.0, 4.0, 5.0]
        y = [2.0, 4.0, 6.0, 8.0, 10.0]
        r = correlation_coefficient(x, y)
        assert abs(r - 1.0) < 1e-10

    def test_perfect_negative_correlation(self):
        """Test perfect negative correlation (r = -1.0)."""
        x = [1.0, 2.0, 3.0, 4.0, 5.0]
        y = [5.0, 4.0, 3.0, 2.0, 1.0]
        r = correlation_coefficient(x, y)
        assert abs(r - (-1.0)) < 1e-10

    def test_no_correlation(self):
        """Test no correlation (r ≈ 0)."""
        x = [1.0, 2.0, 3.0, 4.0, 5.0]
        y = [1.0, 1.0, 1.0, 1.0, 1.0]  # Constant
        r = correlation_coefficient(x, y)
        assert abs(r) < 1e-10

    def test_positive_correlation(self):
        """Test positive correlation."""
        x = [1.0, 2.0, 3.0, 4.0, 5.0]
        y = [1.5, 2.5, 3.5, 4.5, 5.5]
        r = correlation_coefficient(x, y)
        assert 0.0 < r <= 1.0

    def test_different_lengths_raises_error(self):
        """Test that different length lists raise ValueError."""
        x = [1.0, 2.0, 3.0]
        y = [1.0, 2.0]
        with pytest.raises(ValueError, match="Lists must have same length"):
            correlation_coefficient(x, y)

    def test_too_short_returns_zero(self):
        """Test that lists with fewer than 2 elements return 0.0."""
        x = [1.0]
        y = [2.0]
        r = correlation_coefficient(x, y)
        assert r == 0.0

    def test_zero_variance_returns_zero(self):
        """Test that zero variance returns 0.0."""
        x = [1.0, 1.0, 1.0]
        y = [2.0, 3.0, 4.0]
        r = correlation_coefficient(x, y)
        assert r == 0.0


class TestLinearRegression:
    """Tests for linear_regression function."""

    def test_perfect_line(self):
        """Test perfect linear relationship (y = 2x + 1)."""
        x = [1.0, 2.0, 3.0, 4.0, 5.0]
        y = [3.0, 5.0, 7.0, 9.0, 11.0]
        slope, intercept, r2 = linear_regression(x, y)
        assert abs(slope - 2.0) < 1e-10
        assert abs(intercept - 1.0) < 1e-10
        assert abs(r2 - 1.0) < 1e-10

    def test_perfect_line_through_origin(self):
        """Test perfect line through origin (y = 2x)."""
        x = [1.0, 2.0, 3.0, 4.0, 5.0]
        y = [2.0, 4.0, 6.0, 8.0, 10.0]
        slope, intercept, r2 = linear_regression(x, y)
        assert abs(slope - 2.0) < 1e-10
        assert abs(intercept) < 1e-10
        assert abs(r2 - 1.0) < 1e-10

    def test_no_relationship(self):
        """Test no linear relationship."""
        x = [1.0, 2.0, 3.0, 4.0, 5.0]
        y = [1.0, 1.0, 1.0, 1.0, 1.0]  # Constant
        slope, intercept, r2 = linear_regression(x, y)
        assert abs(slope) < 1e-10
        assert abs(intercept - 1.0) < 1e-10
        assert abs(r2) < 1e-10

    def test_positive_slope(self):
        """Test positive slope."""
        x = [1.0, 2.0, 3.0, 4.0, 5.0]
        y = [1.5, 2.8, 3.2, 4.5, 5.1]
        slope, intercept, r2 = linear_regression(x, y)
        assert slope > 0.0
        assert 0.0 <= r2 <= 1.0

    def test_different_lengths_raises_error(self):
        """Test that different length lists raise ValueError."""
        x = [1.0, 2.0, 3.0]
        y = [1.0, 2.0]
        with pytest.raises(ValueError, match="Lists must have same length"):
            linear_regression(x, y)

    def test_too_short_returns_zeros(self):
        """Test that lists with fewer than 2 elements return (0, 0, 0)."""
        x = [1.0]
        y = [2.0]
        slope, intercept, r2 = linear_regression(x, y)
        assert slope == 0.0
        assert intercept == 0.0
        assert r2 == 0.0

    def test_zero_variance_x_returns_zeros(self):
        """Test that zero variance in x returns (0, 0, 0)."""
        x = [1.0, 1.0, 1.0]
        y = [2.0, 3.0, 4.0]
        slope, intercept, r2 = linear_regression(x, y)
        assert slope == 0.0
        assert intercept == 0.0
        assert r2 == 0.0


class TestFisherExactTest:
    """Tests for fisher_exact_test function."""

    def test_strong_association(self):
        """Test strong association in 2x2 table."""
        # Strong positive association
        odds, p = fisher_exact_test(10, 2, 3, 15)
        assert odds > 1.0
        assert 0.0 <= p <= 1.0

    def test_odds_ratio_calculation(self):
        """Test odds ratio calculation."""
        # Perfect association: all in diagonal
        odds, p = fisher_exact_test(10, 0, 0, 10)
        assert math.isinf(odds) or odds > 1000.0
        assert 0.0 <= p <= 1.0

    def test_no_association(self):
        """Test no association (independent)."""
        # Balanced table suggesting independence
        odds, p = fisher_exact_test(5, 5, 5, 5)
        assert abs(odds - 1.0) < 0.1  # Close to 1.0
        assert 0.0 <= p <= 1.0

    def test_requires_scipy(self):
        """Test that function requires scipy."""
        # This test will skip if scipy is not available
        try:
            import scipy.stats  # noqa: F401

            # If scipy is available, test should work
            odds, p = fisher_exact_test(10, 2, 3, 15)
            assert isinstance(odds, float)
            assert isinstance(p, float)
            assert 0.0 <= p <= 1.0
        except ImportError:
            # If scipy not available, should raise ImportError
            with pytest.raises(ImportError, match="scipy is required"):
                fisher_exact_test(10, 2, 3, 15)


class TestShannonEntropy:
    """Tests for shannon_entropy function."""

    def test_maximum_entropy_two_outcomes(self):
        """Test maximum entropy for 2 outcomes."""
        # Uniform distribution for 2 outcomes
        entropy = shannon_entropy([0.5, 0.5])
        assert abs(entropy - 1.0) < 1e-10

    def test_maximum_entropy_four_outcomes(self):
        """Test maximum entropy for 4 outcomes."""
        # Uniform distribution for 4 outcomes
        entropy = shannon_entropy([0.25, 0.25, 0.25, 0.25])
        assert abs(entropy - 2.0) < 1e-10

    def test_zero_entropy(self):
        """Test zero entropy (certainty)."""
        # Certain outcome
        entropy = shannon_entropy([1.0, 0.0])
        assert abs(entropy) < 1e-10

    def test_non_normalized_inputs(self):
        """Test that non-normalized inputs are handled."""
        # Inputs sum to 2.0, not 1.0
        entropy = shannon_entropy([1.0, 1.0])
        # Should still compute correctly (just considers positive values)
        assert entropy > 0.0

    def test_negative_values_ignored(self):
        """Test that negative values are ignored."""
        entropy1 = shannon_entropy([0.5, 0.5])
        entropy2 = shannon_entropy([0.5, 0.5, -0.1])
        assert abs(entropy1 - entropy2) < 1e-10

    def test_empty_list(self):
        """Test empty list returns 0.0."""
        entropy = shannon_entropy([])
        assert entropy == 0.0

    def test_single_value(self):
        """Test single value (certainty)."""
        entropy = shannon_entropy([1.0])
        assert abs(entropy) < 1e-10


class TestJensenShannonDivergence:
    """Tests for jensen_shannon_divergence function."""

    def test_identical_distributions(self):
        """Test that identical distributions have JS divergence = 0."""
        p = [0.5, 0.5]
        q = [0.5, 0.5]
        js = jensen_shannon_divergence(p, q)
        assert abs(js) < 1e-10

    def test_different_distributions(self):
        """Test that different distributions have JS divergence > 0."""
        p = [1.0, 0.0]
        q = [0.0, 1.0]
        js = jensen_shannon_divergence(p, q)
        assert js > 0.0
        assert js <= 1.0  # JS divergence is bounded in [0, 1]

    def test_similar_distributions(self):
        """Test similar distributions have small JS divergence."""
        p = [0.5, 0.5]
        q = [0.6, 0.4]
        js = jensen_shannon_divergence(p, q)
        assert 0.0 < js < 0.1

    def test_different_lengths_raises_error(self):
        """Test that different length distributions raise ValueError."""
        p = [0.5, 0.5]
        q = [0.33, 0.33, 0.34]
        with pytest.raises(ValueError, match="Distributions must have same length"):
            jensen_shannon_divergence(p, q)

    def test_zero_sum_distribution(self):
        """Test that zero-sum distribution returns 0.0."""
        p = [0.0, 0.0]
        q = [0.5, 0.5]
        js = jensen_shannon_divergence(p, q)
        assert js == 0.0

    def test_non_normalized_inputs(self):
        """Test that non-normalized inputs are normalized."""
        p = [1.0, 1.0]  # Sums to 2.0
        q = [1.0, 1.0]  # Sums to 2.0
        js = jensen_shannon_divergence(p, q)
        assert abs(js) < 1e-10  # Should be identical after normalization

    def test_symmetry(self):
        """Test that JS divergence is symmetric."""
        p = [0.8, 0.2]
        q = [0.3, 0.7]
        js1 = jensen_shannon_divergence(p, q)
        js2 = jensen_shannon_divergence(q, p)
        assert abs(js1 - js2) < 1e-10


class TestEdgeCases:
    """Tests for edge cases across all utility functions."""

    def test_correlation_large_arrays(self):
        """Test correlation with large arrays."""
        n = 1000
        x = list(range(n))
        y = [2 * xi + 1 for xi in x]
        r = correlation_coefficient(x, y)
        assert abs(r - 1.0) < 1e-10

    def test_linear_regression_large_arrays(self):
        """Test linear regression with large arrays."""
        n = 1000
        x = list(range(n))
        y = [2 * xi + 1 for xi in x]
        slope, intercept, r2 = linear_regression(x, y)
        assert abs(slope - 2.0) < 1e-6
        assert abs(intercept - 1.0) < 1e-6
        assert abs(r2 - 1.0) < 1e-6

    def test_entropy_large_distribution(self):
        """Test entropy with many outcomes."""
        n = 100
        uniform = [1.0 / n] * n
        entropy = shannon_entropy(uniform)
        expected = math.log2(n)
        assert abs(entropy - expected) < 1e-6

    def test_js_divergence_large_distributions(self):
        """Test JS divergence with many outcomes."""
        n = 100
        p = [1.0 / n] * n
        q = [1.0 / n] * n
        js = jensen_shannon_divergence(p, q)
        assert abs(js) < 1e-10

