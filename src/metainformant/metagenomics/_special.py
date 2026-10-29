"""Pure-Python special functions used when SciPy is not installed.

The metagenomics package treats SciPy as optional: statistics that need a
distribution tail call ``scipy.stats`` when it can be imported and fall back to
the exact implementations here otherwise. The chi-square and Poisson tails are
expressed through the regularized incomplete gamma function and Student's t
through the regularized incomplete beta function. Each is evaluated with the
series or continued fraction that converges quickly for its argument, and the
smaller of the two complementary tails is always computed directly, so tiny
tails keep full relative precision (a ``1 - cdf`` formulation would lose them to
rounding).
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


# The beta continued fraction stops once a term changes the product by no more than two units in
# the last place; an exact ``== 1.0`` test can oscillate between neighbouring doubles forever.
_BETA_TOLERANCE = 4.5e-16

# Degrees of freedom from which the Student t tail comes from a normalising transformation
_T_NORMAL_DF = 1e9


def _beta_continued_fraction(a: float, b: float, x: float) -> float:
    """Continued fraction ``h`` of the incomplete beta function (modified Lentz algorithm).

    ``I_x(a, b) = x**a * (1 - x)**b / (a * B(a, b)) * h``. The fraction converges quickly for
    ``x < (a + 1) / (a + b + 2)`` and slowly (or not at all) elsewhere.

    Args:
        a: First shape parameter, ``a > 0``.
        b: Second shape parameter, ``b > 0``.
        x: Evaluation point in ``(0, 1)``.

    Returns:
        The value of the continued fraction.

    Raises:
        ArithmeticError: If the fraction has not converged after ``_MAX_ITERATIONS`` terms.
    """
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < _TINY:
        d = _TINY
    d = 1.0 / d
    h = d
    for m in range(1, _MAX_ITERATIONS):
        m2 = 2 * m
        # Even step of the fraction
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < _TINY:
            d = _TINY
        c = 1.0 + aa / c
        if abs(c) < _TINY:
            c = _TINY
        d = 1.0 / d
        h *= d * c
        # Odd step of the fraction
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < _TINY:
            d = _TINY
        c = 1.0 + aa / c
        if abs(c) < _TINY:
            c = _TINY
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) <= _BETA_TOLERANCE:
            return h
    raise ArithmeticError(f"incomplete beta continued fraction did not converge for a={a}, b={b}, x={x}")


def _stirling_correction(z: float) -> float:
    """Stirling series ``1/(12 z) - 1/(360 z**3) + 1/(1260 z**5) - 1/(1680 z**7)`` (accurate for ``z >= 10``)."""
    z_squared = z * z
    return (1.0 / 12.0 - (1.0 / 360.0 - (1.0 / 1260.0 - 1.0 / (1680.0 * z_squared)) / z_squared) / z_squared) / z


def _log_gamma_ratio(a: float, b: float) -> float:
    """``lgamma(a + b) - lgamma(a)`` for ``a >= b > 0`` without the cancellation of the plain difference.

    For large ``a`` both log-gamma values are huge and nearly equal, so their difference loses
    about ``log10(a)`` digits. The Stirling expansion cancels the large terms analytically.
    """
    if a < 10.0:
        return math.lgamma(a + b) - math.lgamma(a)
    return (
        (a - 0.5) * math.log1p(b / a) + b * math.log(a + b) - b + _stirling_correction(a + b) - _stirling_correction(a)
    )


def _beta_tails(a: float, b: float, x: float, y: float) -> tuple[float, float]:
    """Return ``(I_x(a, b), 1 - I_x(a, b))`` with the smaller tail computed directly.

    The continued fraction is evaluated on the side where it converges quickly, which is also the
    side holding the smaller tail, so that tail keeps full relative precision. The other tail is
    its complement, which is accurate in absolute terms.

    Args:
        a: First shape parameter, ``a > 0``.
        b: Second shape parameter, ``b > 0``.
        x: Evaluation point in ``[0, 1]``.
        y: ``1 - x``, supplied separately because the caller often knows it without cancellation.

    Returns:
        ``(lower, upper)``, the regularized incomplete beta function and its complement.
    """
    if x <= 0.0:
        return 0.0, 1.0
    if y <= 0.0:
        return 1.0, 0.0
    # log(x) and log(y) lose all precision when x or y is close to 1, and the shapes multiply them
    log_x = math.log1p(-y) if x > 0.5 else math.log(x)
    log_y = math.log1p(-x) if y > 0.5 else math.log(y)
    large, small = (a, b) if a >= b else (b, a)
    log_front = _log_gamma_ratio(large, small) - math.lgamma(small) + a * log_x + b * log_y
    front = math.exp(log_front)
    if x < (a + 1.0) / (a + b + 2.0):
        lower = min(1.0, front * _beta_continued_fraction(a, b, x) / a)
        return lower, 1.0 - lower
    upper = min(1.0, front * _beta_continued_fraction(b, a, y) / b)
    return 1.0 - upper, upper


def _check_beta_arguments(a: float, b: float, x: float) -> None:
    if not (a > 0 and b > 0):
        raise ValueError(f"a and b must be positive, got a={a}, b={b}")
    if not 0.0 <= x <= 1.0:
        raise ValueError(f"x must lie in [0, 1], got {x}")


def regularized_beta_p(a: float, b: float, x: float) -> float:
    """Regularized incomplete beta function ``I_x(a, b) = P(Beta(a, b) <= x)``.

    Args:
        a: First shape parameter, ``a > 0``.
        b: Second shape parameter, ``b > 0``.
        x: Evaluation point in ``[0, 1]``.

    Returns:
        ``I_x(a, b)`` in ``[0, 1]``.

    Raises:
        ValueError: If ``a`` or ``b`` is not positive, or ``x`` is outside ``[0, 1]``.
    """
    _check_beta_arguments(a, b, x)
    return _beta_tails(a, b, x, 1.0 - x)[0]


def regularized_beta_q(a: float, b: float, x: float) -> float:
    """Complement ``1 - I_x(a, b) = I_{1-x}(b, a)``, accurate when it is tiny.

    Args:
        a: First shape parameter, ``a > 0``.
        b: Second shape parameter, ``b > 0``.
        x: Evaluation point in ``[0, 1]``.

    Returns:
        ``P(Beta(a, b) > x)`` in ``[0, 1]``.

    Raises:
        ValueError: If ``a`` or ``b`` is not positive, or ``x`` is outside ``[0, 1]``.
    """
    _check_beta_arguments(a, b, x)
    return _beta_tails(a, b, x, 1.0 - x)[1]


def student_t_two_sided_p(t: float, df: float) -> float:
    """Two-sided tail probability ``P(|T| >= |t|)`` of Student's t distribution.

    Uses the identity ``P(|T| >= |t|) = I_x(df / 2, 1 / 2)`` with ``x = df / (df + t**2)``, so the
    p-value is a regularized incomplete beta function evaluated at arguments that are exact even
    when ``|t|`` is large, and p-values far below double-precision rounding of ``1 - cdf`` stay
    accurate. This is the exact tail, not a normal approximation.

    Beyond ``df = 1e9`` the incomplete beta function loses digits (the shapes are enormous), so the
    tail is taken from the normalising transformation of Abramowitz and Stegun 26.7.8,
    ``z = t (1 - 1 / (4 df)) / sqrt(1 + t**2 / (2 df))``, whose relative error there is below 1e-9.

    Args:
        t: Test statistic (NaN gives NaN, infinities give 0.0).
        df: Degrees of freedom, ``df > 0`` (need not be an integer; ``inf`` gives the normal tail).

    Returns:
        ``P(|T| >= |t|)`` in ``[0, 1]``.

    Raises:
        ValueError: If ``df`` is not positive.
    """
    if not df > 0:
        raise ValueError(f"df must be positive, got {df}")
    if math.isnan(t):
        return math.nan
    if math.isinf(df):
        return math.erfc(abs(t) / math.sqrt(2.0))
    t_squared = t * t
    if df >= _T_NORMAL_DF and not math.isinf(t_squared):
        z = abs(t) * (1.0 - 0.25 / df) / math.sqrt(1.0 + t_squared / (2.0 * df))
        return math.erfc(z / math.sqrt(2.0))
    if math.isinf(t_squared):
        return 0.0
    denominator = df + t_squared
    return _beta_tails(df / 2.0, 0.5, df / denominator, t_squared / denominator)[0]


def student_t_sf(t: float, df: float) -> float:
    """Survival function ``P(T > t)`` of Student's t distribution with ``df`` degrees of freedom.

    Args:
        t: Test statistic (NaN gives NaN).
        df: Degrees of freedom, ``df > 0``.

    Returns:
        ``P(T > t)`` in ``[0, 1]``; the upper tail keeps full relative precision for large ``t``.

    Raises:
        ValueError: If ``df`` is not positive.
    """
    two_sided = student_t_two_sided_p(t, df)
    if math.isnan(two_sided):
        return math.nan
    half = 0.5 * two_sided
    if t < 0:
        return 1.0 - half
    return half
