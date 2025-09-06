from __future__ import annotations

import math

from metainformant.math import (
    correlation,
    covariance,
    delta_mean_trait,
    expectation,
    price_equation,
    relative_fitness,
    selection_differential,
    selection_gradient,
    selection_intensity,
    standard_deviation,
    variance,
    weighted_correlation,
    weighted_covariance,
    weighted_variance,
)


def test_expectation_and_variance_basic():
    vals = [0.2, 0.4, 0.1]
    assert abs(expectation(vals) - (sum(vals) / 3)) < 1e-12
    assert abs(variance(vals) - (sum((v - sum(vals) / 3) ** 2 for v in vals) / 3)) < 1e-12
    assert abs(standard_deviation(vals) - (variance(vals) ** 0.5)) < 1e-12
    # Weighted variance with weights [1,2,1]
    ws = [1.0, 2.0, 1.0]
    mu_w = sum(v * w for v, w in zip(vals, ws)) / sum(ws)
    var_w = sum(w * (v - mu_w) ** 2 for v, w in zip(vals, ws)) / sum(ws)
    assert abs(weighted_variance(vals, ws) - var_w) < 1e-12


def test_weighted_expectation():
    vals = [1.0, 2.0, 3.0]
    ws = [0.0, 1.0, 1.0]
    # weighted mean should be (2*1 + 3*1) / (1+1) = 2.5
    assert abs(expectation(vals, ws) - 2.5) < 1e-12
    # empty / zero-weight guards
    assert expectation([], []) == 0.0
    assert expectation([1.0], [0.0]) == 0.0


def test_covariance_and_correlation():
    x = [1.0, 2.0, 3.0]
    y = [2.0, 4.0, 6.0]
    cov_xy = covariance(x, y)
    # Perfect linear relation with positive slope => correlation ~ 1
    rho = correlation(x, y)
    assert cov_xy > 0
    assert abs(rho - 1.0) < 1e-12

    # Zero variance case -> correlation returns 0.0
    assert correlation([1.0, 1.0, 1.0], [2.0, 3.0, 4.0]) == 0.0

    # Weighted covariance/correlation
    ws = [1.0, 2.0, 1.0]
    cov_w = weighted_covariance(x, y, ws)
    rho_w = weighted_correlation(x, y, ws)
    assert cov_w > 0
    assert abs(rho_w - 1.0) < 1e-12


def test_selection_metrics_and_price_equation():
    w = [1.0, 1.2, 0.9]
    z = [0.2, 0.4, 0.1]
    z_prime = [0.25, 0.35, 0.15]

    # Selection differential, intensity, and gradient
    w_rel = relative_fitness(w)
    assert abs(sum(w_rel) / len(w_rel) - 1.0) < 1e-12
    S = selection_differential(w, z)
    i = selection_intensity(w, z)
    beta = selection_gradient(w, z)
    assert isinstance(S, float) and isinstance(i, float) and isinstance(beta, float)

    # Price decomposition with relative fitness normalization
    cov_term, trans_term, total = price_equation(w, z, z_prime)
    assert math.isclose(total, cov_term + trans_term, rel_tol=0, abs_tol=1e-15)
    assert math.isclose(delta_mean_trait(w, z, z_prime), total, rel_tol=0, abs_tol=1e-15)
    # If no offspring provided, transmission term should be 0
    cov2, trans2, total2 = price_equation(w, z)
    assert abs(trans2) < 1e-15
    assert math.isclose(total2, cov2, rel_tol=0, abs_tol=1e-15)
