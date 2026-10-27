"""Pure-Python special functions used when SciPy is not installed.

The metagenomics package treats SciPy as optional: statistics that need a
distribution tail call ``scipy.stats`` when it can be imported and fall back to
the exact implementations here otherwise. Both tails are expressed through the
regularized incomplete gamma function, evaluated with the series expansion where
it converges quickly and the continued fraction elsewhere, so tiny tails keep
full relative precision (a ``1 - cdf`` formulation would lose them to rounding).
"""

from __future__ import annotations

import math

_EPSILON = 1e-16
_MAX_ITERATIONS = 100_000
_TINY = 1e-300


def _lower_series(a: float, x: float, log_prefactor: float) -> float:
    """Regularized lower incomplete gamma ``P(a, x)`` by its series (best for ``x < a + 1``)."""
    term = 1.0 / a
    total = term
    denominator = a
    for _ in range(_MAX_ITERATIONS):
        denominator += 1.0
        term *= x / denominator
        total += term
        if abs(term) < abs(total) * _EPSILON:
            break
    return total * math.exp(log_prefactor)


def _upper_continued_fraction(a: float, x: float, log_prefactor: float) -> float:
    """Regularized upper incomplete gamma ``Q(a, x)`` by Lentz's continued fraction (best for ``x >= a + 1``)."""
    b = x + 1.0 - a
    c = 1.0 / _TINY
    d = 1.0 / b
    h = d
    for i in range(1, _MAX_ITERATIONS):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < _TINY:
            d = _TINY
        c = b + an / c
        if abs(c) < _TINY:
            c = _TINY
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < _EPSILON:
            break
    return math.exp(log_prefactor) * h


def regularized_gamma_p(a: float, x: float) -> float:
    """Regularized lower incomplete gamma function ``P(a, x)``.

    Args:
        a: Shape parameter, ``a > 0``.
        x: Upper integration limit, ``x >= 0``.

    Returns:
        ``P(a, x)`` in ``[0, 1]``.

    Raises:
        ValueError: If ``a`` is not positive.
    """
    if not a > 0:
        raise ValueError(f"a must be positive, got {a}")
    if x <= 0:
        return 0.0
    if math.isinf(x):
        return 1.0
    log_prefactor = -x + a * math.log(x) - math.lgamma(a)
    if x < a + 1.0:
        return min(1.0, _lower_series(a, x, log_prefactor))
    return max(0.0, 1.0 - _upper_continued_fraction(a, x, log_prefactor))


def regularized_gamma_q(a: float, x: float) -> float:
    """Regularized upper incomplete gamma function ``Q(a, x) = 1 - P(a, x)``.

    Args:
        a: Shape parameter, ``a > 0``.
        x: Lower integration limit, ``x >= 0``.

    Returns:
        ``Q(a, x)`` in ``[0, 1]``.

    Raises:
        ValueError: If ``a`` is not positive.
    """
    if not a > 0:
        raise ValueError(f"a must be positive, got {a}")
    if x <= 0:
        return 1.0
    if math.isinf(x):
        return 0.0
    log_prefactor = -x + a * math.log(x) - math.lgamma(a)
    if x < a + 1.0:
        return max(0.0, 1.0 - _lower_series(a, x, log_prefactor))
    return min(1.0, _upper_continued_fraction(a, x, log_prefactor))


def chi2_sf(x: float, df: float) -> float:
    """Chi-square survival function ``P(X >= x)`` for ``df`` degrees of freedom.

    Args:
        x: Statistic value.
        df: Degrees of freedom, ``df > 0``.

    Returns:
        Upper-tail probability.
    """
    return regularized_gamma_q(df / 2.0, x / 2.0)


def poisson_upper_tail(k: int, rate: float) -> float:
    """Poisson upper tail ``P(X >= k)`` for ``X ~ Poisson(rate)``.

    Uses the identity ``P(X >= k) = P(k, rate)`` (regularized lower incomplete gamma),
    which keeps tails far below double-precision rounding error exact.

    Args:
        k: Count threshold.
        rate: Poisson mean, ``rate >= 0``.

    Returns:
        ``P(X >= k)``; 1.0 for ``k <= 0`` and 0.0 for ``rate <= 0`` with ``k > 0``.
    """
    if k <= 0:
        return 1.0
    if rate <= 0:
        return 0.0
    return regularized_gamma_p(float(k), rate)
