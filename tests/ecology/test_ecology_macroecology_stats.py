"""Tests for untested macroecology statistical helpers.

Covers:
    - simple_linear_regression
    - chi_squared_gof
    - aic_from_gof
    - bootstrap_ci

Uses real implementations only (real-implementation policy).
"""

from __future__ import annotations

import math

import pytest

from metainformant.ecology.analysis.macroecology import (
    aic_from_gof,
    bootstrap_ci,
    chi_squared_gof,
    simple_linear_regression,
)


def _slope(x: list[float], y: list[float]) -> float:
    """Statistic callable for bootstrap tests: OLS slope only."""
    return simple_linear_regression(x, y)["slope"]


class TestSimpleLinearRegression:
    def test_perfect_fit(self) -> None:
        result = simple_linear_regression([1, 2, 3, 4], [2, 4, 6, 8])
        assert result["slope"] == pytest.approx(2.0)
        assert result["intercept"] == pytest.approx(0.0)
        assert result["r_squared"] == pytest.approx(1.0)

    def test_known_fit_with_noise(self) -> None:
        # y ~ 1 + 0.5x with one slightly perturbed point.
        x = [0.0, 1.0, 2.0, 3.0]
        y = [1.0, 1.5, 2.0, 2.6]
        result = simple_linear_regression(x, y)
        assert 0.4 < result["slope"] < 0.6
        assert 0.9 < result["intercept"] < 1.1
        assert 0.99 < result["r_squared"] <= 1.0

    def test_length_mismatch_raises(self) -> None:
        with pytest.raises(ValueError):
            simple_linear_regression([1, 2], [1, 2, 3])

    def test_single_point_raises(self) -> None:
        with pytest.raises(ValueError):
            simple_linear_regression([1], [1])

    def test_zero_variance_x_returns_flat(self) -> None:
        result = simple_linear_regression([5, 5, 5], [1, 2, 3])
        assert result["slope"] == 0.0
        assert result["intercept"] == pytest.approx(2.0)


class TestChiSquaredGof:
    def test_known_statistic(self) -> None:
        result = chi_squared_gof([10, 20, 30], [15, 15, 30])
        # (10-15)^2/15 + (20-15)^2/15 + 0 = 50/15 = 10/3
        assert result["chi_squared"] == pytest.approx(10.0 / 3.0)
        assert result["degrees_of_freedom"] == 2

    def test_perfect_fit_zero(self) -> None:
        result = chi_squared_gof([10, 20], [10, 20])
        assert result["chi_squared"] == pytest.approx(0.0)

    def test_length_mismatch_raises(self) -> None:
        with pytest.raises(ValueError):
            chi_squared_gof([1, 2], [1, 2, 3])


class TestAicFromGof:
    def test_basic_formula_with_aicc(self) -> None:
        # AIC = chi2 + 2k, plus AICc term when n > k + 1.
        assert aic_from_gof(10.0, 2, 100) == pytest.approx(14.0 + (2 * 2 * 3) / (100 - 3))

    def test_no_correction_when_n_too_small(self) -> None:
        # n == k + 1 -> no AICc term.
        assert aic_from_gof(10.0, 2, 3) == pytest.approx(14.0)

    def test_ranking_penalizes_parameters(self) -> None:
        # Same fit, more parameters => strictly worse (higher) AIC.
        assert aic_from_gof(10.0, 3, 100) > aic_from_gof(10.0, 2, 100)


class TestBootstrapCI:
    def test_deterministic_with_seed(self) -> None:
        x = [1.0, 2.0, 3.0, 4.0, 5.0]
        y = [2.1, 3.9, 6.2, 7.8, 10.1]
        ci1 = bootstrap_ci(x, y, _slope, n_bootstrap=200, seed=123)
        ci2 = bootstrap_ci(x, y, _slope, n_bootstrap=200, seed=123)
        assert ci1 == ci2

    def test_ci_brackets_true_slope(self) -> None:
        x = [float(i) for i in range(10)]
        y = [3.0 * xi + 1.0 for xi in x]
        lo, hi = bootstrap_ci(x, y, _slope, n_bootstrap=500, seed=7)
        assert math.isfinite(lo) and math.isfinite(hi)
        assert lo <= hi
        # Deterministic data: every resample slope is finite and near 3.
        assert lo < 3.0 < hi or abs(lo - 3.0) < 0.5 or abs(hi - 3.0) < 0.5

    def test_degenerate_input_returns_point_estimate(self) -> None:
        val = bootstrap_ci([1.0], [2.0], lambda xs, ys: float(ys[0]), n_bootstrap=50)
        assert val == (2.0, 2.0)
