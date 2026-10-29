"""Drift-diffusion models for decision making and cognitive processes.

The evidence follows ``dx = v dt + s dW`` between absorbing boundaries ``a`` apart (the lower one at 0,
the upper one at ``a``) and starts at ``z = starting_point * a``. The module provides the closed forms
for the choice probability and mean decision time of an unbiased start, the exact first-passage-time
density of either boundary for any start, and a maximum-likelihood fit built on that density.

First-passage density
---------------------
The density of absorption at one boundary at time ``t`` follows from the zero-drift density ``g`` of a
unit-width, unit-noise process by a change of measure and a rescaling of time (Navarro and Fuss 2009).
With ``near`` the start's distance to the boundary in units of ``a``, ``far = 1 - near``, ``u = s^2 t / a^2``
and ``away`` the drift (in the frame where the boundary reached is at 0) that carries the process away from it,

    f(t) = (s^2 / a^2) exp(-away * near * a / s^2 - away^2 t / (2 s^2)) g(u; near),

with ``away = v`` for the lower boundary and ``away = -v`` for the upper one. ``g`` has two exact series:

    large-time (eigenfunction): g = pi sum_{k>=1} k exp(-k^2 pi^2 u / 2) sin(k pi near)
    small-time (image method):  g = (2 pi u^3)^(-1/2) sum_{k in Z} (near + 2k) exp(-(near + 2k)^2 / (2u))

The large-time series needs few terms when ``u`` is large and the small-time series when ``u`` is small.
The module switches between them at ``u = 1/4``. Each series is truncated after a fixed number of terms
chosen so that the discarded terms change the density by less than a relative ``1e-12`` anywhere in its
range; see :func:`_series_term_counts` for the bound. Both series are summed in a form that keeps the
relative accuracy for a start close to either boundary, and the density is returned through its logarithm,
so response times far in the tails give a finite log-likelihood rather than an underflowed zero.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

# Standardised time u = s^2 t / a^2 below which the small-time series is used and from which the large-time
# series is used. At 1/4 each series needs only a handful of terms and both are well conditioned.
_SERIES_SWITCH = 0.25

# Relative accuracy guaranteed for the first-passage density by the number of series terms.
_DENSITY_RTOL = 1e-12

# Below this standardised time the density exp(-near^2 / 2u) is far below the smallest double and its
# logarithm approaches the largest one; such times are reported as an impossible (zero-density) response.
_MIN_STANDARD_TIME = 1e-300


def _series_term_counts(rtol: float) -> tuple[int, int]:
    """Numbers of series terms that bound the relative truncation error of the density by ``rtol``.

    Both counts are fixed for a given ``rtol``: each series is truncated for the worst case of its range,
    so the density is a smooth function of the parameters and the data (no term count changes under an
    optimiser's feet) and the bound holds for every response time.

    *Small-time series* (``u < 1/4``). With ``near`` the start's relative distance to the boundary reached
    and ``far = 1 - near``, the image sum is grouped in pairs that are each free of cancellation: the images
    ``2j + near`` and ``2j - near`` for ``near <= 1/2`` (leading term ``near``), and ``2j + near`` and
    ``2j + 2 - near`` for ``near > 1/2`` (all pairs positive). Bounding a pair by the integral of the
    derivative of ``x exp(-x^2 / 2u)`` over its width, relative to the leading term (``near`` in the first
    case, the first pair in the second), gives

        pair j  <=  m_j(u) = (8 j^2 / u) exp(-j (2 j - 1) / u),   j >= 1,

    which grows with ``u`` and whose consecutive ratio is at most ``4 exp(-5/u)``. The whole sum is at least
    0.41 times its leading term in the first case and 0.47 times in the second (the grid minima at
    ``u = 1/4`` are 0.945 and 0.4725), so truncating after ``K`` pairs changes the density by a relative
    ``2.44 * sum_{j>K} m_j`` at most. The first count is the smallest ``K`` for which that tail, at
    ``u = 1/4`` where ``m_j`` is largest, is below ``0.4 * rtol`` (pairs ``j = 1..K`` after the leading one).

    *Large-time series* (``u >= 1/4``). Since ``|sin(k x)| <= k |sin(x)|`` the ``k``-th term, divided by the
    first, is at most ``M_k(u) = k^2 exp(-(k^2 - 1) pi^2 u / 2)``, so the sum of all terms is at least
    ``1 - sum_{k>=2} M_k >= 0.9`` times the first one, for every start, and the relative error of the
    truncation is at most ``sum_{k>K} M_k / 0.9`` (consecutive ratios are below ``5e-3``). The second count
    is the smallest ``K`` for which that tail, at ``u = 1/4`` (where the terms are largest), is below
    ``0.9 * rtol``.

    Args:
        rtol: Target relative error of the density, within (0, 1)

    Returns:
        ``(small_time_pairs, large_time_terms)``: pairs beyond the leading one in the small-time series,
        and terms in the large-time series
    """
    if not 0.0 < rtol < 1.0:
        raise ValueError("rtol must lie strictly between 0 and 1")

    u = _SERIES_SWITCH
    ratio = 4.0 * math.exp(-5.0 / u)
    small_pairs = 0
    while True:
        j = small_pairs + 1
        tail = (8.0 * j * j / u) * math.exp(-j * (2 * j - 1) / u) / (1.0 - ratio)
        if tail <= 0.4 * rtol:
            break
        small_pairs += 1

    c = 0.5 * math.pi**2 * u
    large_terms = 1
    while True:
        k = large_terms + 1
        ratio = ((large_terms + 2) / (large_terms + 1)) ** 2 * math.exp(-(2 * large_terms + 3) * c)
        tail = k * k * math.exp(-(k * k - 1) * c) / (1.0 - ratio)
        if tail <= 0.9 * rtol:
            break
        large_terms += 1

    return small_pairs, large_terms


_SMALL_TIME_PAIRS, _LARGE_TIME_TERMS = _series_term_counts(_DENSITY_RTOL)


def _log_standard_density_small_time(u: np.ndarray, near: float, far: float, pairs: int) -> np.ndarray:
    """Log of the zero-drift density ``g(u; near)`` from the small-time image series (``u`` below 1/4).

    ``near`` and ``far = 1 - near`` are the start's distances to the boundary reached and to the other one,
    in units of the boundary separation. The image terms are summed in pairs that avoid cancellation, see
    :func:`_series_term_counts`; ``pairs`` is the number of pairs after the leading one.
    """
    scaled = u[:, None]
    if near <= far:
        # g = (2 pi u^3)^(-1/2) exp(-near^2 / 2u) * near * (1 + sum_j psi_j), the pair of images 2j +- near
        # relative to the leading term being psi_j
        j = np.arange(1, pairs + 1, dtype=float)
        psi = np.exp(-2.0 * j * (j - near) / scaled) * (
            2.0 + (2.0 * j + near) * np.expm1(-4.0 * j * near / scaled) / near
        )
        log_sum = math.log(near) + np.log1p(psi.sum(axis=1))
    else:
        # near is close to 1: the leading image and the one reflected in the far boundary nearly cancel,
        # so they form the first pair; g is the sum of the (positive) pairs 2j + near and 2j + 1 + far
        j = np.arange(0, pairs + 1, dtype=float)
        pair = np.exp(-2.0 * j * (j + near) / scaled) * (
            -2.0 * far - (2.0 * j + 1.0 + far) * np.expm1(-2.0 * far * (2.0 * j + 1.0) / scaled)
        )
        log_sum = np.log(pair.sum(axis=1))
    return np.asarray(-0.5 * (math.log(2.0 * math.pi) + 3.0 * np.log(u)) - near * near / (2.0 * u) + log_sum)


def _log_standard_density_large_time(u: np.ndarray, near: float, far: float, terms: int) -> np.ndarray:
    """Log of the zero-drift density ``g(u; near)`` from the large-time eigenfunction series (``u`` from 1/4)."""
    nearest = min(near, far)  # sin(k pi near) = (-1)^(k+1) sin(k pi far); the smaller argument keeps full accuracy
    k = np.arange(1, terms + 1, dtype=float)
    weight = k * np.sin(k * math.pi * nearest) / math.sin(math.pi * nearest)  # k sin(k x) / sin(x), |.| <= k^2
    if near > far:
        weight = weight * np.where(k % 2 == 1, 1.0, -1.0)
    c = 0.5 * math.pi**2 * u
    total = np.exp(-(k * k - 1.0) * c[:, None]) @ weight
    return np.asarray(math.log(math.pi) + math.log(math.sin(math.pi * nearest)) - c + np.log(total))


def _log_standard_density(u: np.ndarray, near: float, far: float, pairs: int, terms: int) -> np.ndarray:
    """Log of ``g(u; near)``: the small-time series below ``u = 1/4`` and the large-time series from there.

    ``-inf`` (a zero density) for ``u`` below ``1e-300``.
    """
    out = np.full_like(u, -np.inf)
    small = (u < _SERIES_SWITCH) & (u >= _MIN_STANDARD_TIME)
    large = u >= _SERIES_SWITCH
    if small.any():
        out[small] = _log_standard_density_small_time(u[small], near, far, pairs)
    if large.any():
        out[large] = _log_standard_density_large_time(u[large], near, far, terms)
    return out


def _log_boundary_density(
    decision_time: np.ndarray,
    drift_rate: float,
    boundary_separation: float,
    starting_point: float,
    noise_sd: float,
    upper: bool,
    pairs: int = _SMALL_TIME_PAIRS,
    terms: int = _LARGE_TIME_TERMS,
) -> np.ndarray:
    """Log of the density of absorption at one boundary at the (positive) ``decision_time`` values.

    Args:
        decision_time: Times since the start of the evidence accumulation, all positive
        drift_rate: Drift ``v``
        boundary_separation: Boundary separation ``a``
        starting_point: Start as a fraction of ``a`` measured from the lower boundary, within (0, 1)
        noise_sd: Diffusion standard deviation ``s``
        upper: True for the upper boundary, False for the lower one
        pairs: Pairs of image terms after the leading one in the small-time series
        terms: Terms of the large-time series

    Returns:
        ``log f(t)``, the log of the density of the joint event "absorbed at this boundary at time t"
    """
    near = 1.0 - starting_point if upper else starting_point
    far = starting_point if upper else 1.0 - starting_point
    away = -drift_rate if upper else drift_rate  # drift away from the boundary reached
    variance = noise_sd * noise_sd
    u = variance * decision_time / (boundary_separation * boundary_separation)
    with np.errstate(over="ignore", invalid="ignore"):
        return np.asarray(
            math.log(variance / (boundary_separation * boundary_separation))
            - away * near * boundary_separation / variance
            - away * away * decision_time / (2.0 * variance)
            + _log_standard_density(u, near, far, pairs, terms)
        )


def _validate_ddm_parameters(drift_rate: float, boundary: float, noise_sd: float) -> None:
    """Reject parameters outside the domain of the closed-form DDM expressions."""
    if not (math.isfinite(boundary) and boundary > 0):
        raise ValueError("Boundary must be positive and finite")
    if not (math.isfinite(noise_sd) and noise_sd > 0):
        raise ValueError("Noise standard deviation must be positive and finite")
    if not math.isfinite(drift_rate):
        raise ValueError("Drift rate must be finite")


def _validate_density_parameters(
    drift_rate: float, boundary: float, starting_point: float, non_decision_time: float, noise_sd: float
) -> None:
    """Reject parameters for which the first-passage density is undefined."""
    _validate_ddm_parameters(drift_rate, boundary, noise_sd)
    if not (math.isfinite(starting_point) and 0.0 < starting_point < 1.0):
        raise ValueError("Starting point must lie strictly between 0 and 1 (a fraction of the boundary separation)")
    if not (math.isfinite(non_decision_time) and non_decision_time >= 0.0):
        raise ValueError("Non-decision time must be non-negative and finite")


def ddm_mean_decision_time(drift_rate: float, boundary: float, noise_sd: float = 1.0) -> float:
    """Calculate the mean decision time of a drift-diffusion model.

    The evidence follows dx = v dt + s dW between absorbing boundaries a apart and starts midway
    between them (an unbiased start). The standard closed form (Palmer, Huk and Shadlen 2005;
    Bogacz et al. 2006, with thresholds at a/2 on either side of the start) is

        E[T] = (a / (2 v)) * tanh(v a / (2 s^2)),

    which tends to a^2 / (4 s^2) as v tends to 0. It is even in v, and with an unbiased start
    correct and error responses share this mean. Non-decision time is not included.

    Args:
        drift_rate: Drift rate (v); only its magnitude matters
        boundary: Boundary separation (a), positive
        noise_sd: Standard deviation of the diffusion noise (s), positive

    Returns:
        Mean decision time

    Raises:
        ValueError: If boundary or noise_sd is not positive, or drift_rate is not finite
    """
    _validate_ddm_parameters(drift_rate, boundary, noise_sd)

    # E[T] = a^2 / (4 s^2) * tanh(x) / x with x = |v| a / (2 s^2); tanh(x) / x tends to 1 as x -> 0.
    x = abs(drift_rate) * boundary / (2.0 * noise_sd**2)
    shape = math.tanh(x) / x if x > 0 else 1.0
    return boundary**2 / (4.0 * noise_sd**2) * shape


def ddm_analytic_accuracy(drift_rate: float, boundary: float, noise_sd: float = 1.0) -> float:
    """Calculate the analytic accuracy of a drift-diffusion model.

    For the diffusion of :func:`ddm_mean_decision_time` started midway between the boundaries,
    the probability of absorption at the upper boundary (the correct response when the drift is
    positive) is the standard closed form (Palmer, Huk and Shadlen 2005; Bogacz et al. 2006)

        P(correct) = 1 / (1 + exp(-v a / s^2)).

    Args:
        drift_rate: Drift rate (v)
        boundary: Boundary separation (a), positive
        noise_sd: Standard deviation of the diffusion noise (s), positive

    Returns:
        Probability of a correct (upper-boundary) decision

    Raises:
        ValueError: If boundary or noise_sd is not positive, or drift_rate is not finite
    """
    _validate_ddm_parameters(drift_rate, boundary, noise_sd)

    return float(expit(drift_rate * boundary / noise_sd**2))


def ddm_first_passage_density(
    rt: Any,
    choice: Any,
    drift_rate: float,
    boundary_separation: float,
    starting_point: float = 0.5,
    non_decision_time: float = 0.0,
    noise_sd: float = 1.0,
    *,
    log: bool = False,
) -> Any:
    """Density of a response at a given time through a given boundary (Wiener first passage time).

    This is the joint density of "the evidence reaches boundary ``choice`` at decision time
    ``rt - non_decision_time``" for dx = v dt + s dW between absorbing boundaries at 0 and ``a``,
    started at ``starting_point * a``. It is zero for ``rt <= non_decision_time``, and integrated over
    time it gives the probability of that choice: ``expit(v a / s^2)`` for the upper boundary and an
    unbiased start (see :func:`ddm_analytic_accuracy`), and in general
    ``(1 - exp(-2 v z / s^2)) / (1 - exp(-2 v a / s^2))`` with ``z = starting_point * a`` (``z / a`` when
    ``v = 0``). Both boundaries are covered: the lower-boundary density at drift ``v`` and start ``w``
    is the upper-boundary density at drift ``-v`` and start ``1 - w``.

    The density is evaluated with the small-time and large-time series of Navarro and Fuss (2009)
    (see the module documentation), truncated so that the relative error stays below ``1e-12``.

    Args:
        rt: Response time(s), a number or an array
        choice: Boundary reached, 1 for the upper and 2 for the lower; a number or an array that
            broadcasts against ``rt``
        drift_rate: Drift rate (v)
        boundary_separation: Distance between the boundaries (a), positive
        starting_point: Start as a fraction of ``boundary_separation`` from the lower boundary (z/a),
            strictly between 0 and 1
        non_decision_time: Non-decision time component (t0), non-negative
        noise_sd: Standard deviation of the diffusion noise (s), positive. Drift, boundary separation
            and noise are only identified up to a common factor, so the density is unchanged when
            ``v``, ``a`` and ``s`` are multiplied by the same constant
        log: Return the natural logarithm of the density, which stays finite far in the tails where the
            density itself underflows

    Returns:
        The density (or its logarithm) as a float for scalar input and an array otherwise; 0 (or
        ``-inf`` on the log scale) where ``rt <= non_decision_time``

    Raises:
        ValueError: If ``choice`` is not 1 or 2, ``rt`` is not finite, the boundary separation or noise
            is not positive and finite, the drift is not finite, the starting point is not strictly
            inside (0, 1) or the non-decision time is negative or not finite
    """
    _validate_density_parameters(drift_rate, boundary_separation, starting_point, non_decision_time, noise_sd)

    response_time, response_choice = np.broadcast_arrays(np.asarray(rt, dtype=float), np.asarray(choice))
    if not np.all(np.isin(response_choice, (1, 2))):
        raise ValueError("Choice must be 1 or 2")
    if not np.all(np.isfinite(response_time)):
        raise ValueError("Response times must be finite")

    decision_time = response_time - non_decision_time
    log_density = np.full(decision_time.shape, -np.inf)
    for code, upper in ((1, True), (2, False)):
        selected = (response_choice == code) & (decision_time > 0.0)
        if selected.any():
            log_density[selected] = _log_boundary_density(
                decision_time[selected], drift_rate, boundary_separation, starting_point, noise_sd, upper
            )

    result = log_density if log else np.exp(log_density)
    return float(result) if result.ndim == 0 else result


def ddm_log_likelihood(
    rt: float,
    choice: int,
    drift_rate: float,
    boundary_separation: float,
    starting_point: float = 0.5,
    non_decision_time: float = 0.0,
    noise_sd: float = 1.0,
) -> float:
    """Calculate the exact log-likelihood of a single DDM trial.

    The log-likelihood is the logarithm of the Wiener first-passage-time density
    (:func:`ddm_first_passage_density`) of the observed boundary at the observed time, for the diffusion
    dx = v dt + s dW between absorbing boundaries at 0 and ``a``, started at ``starting_point * a``. Both
    boundaries are modelled, so a response against the sign of the drift has a finite log-likelihood.

    Args:
        rt: Response time
        choice: Choice (1 or 2, where 1 is upper boundary)
        drift_rate: Drift rate (v)
        boundary_separation: Distance between boundaries (a), positive
        starting_point: Starting point as fraction of boundary separation (z/a), strictly between 0 and 1
        non_decision_time: Non-decision time component (t0), non-negative
        noise_sd: Standard deviation of the diffusion noise (s), positive; the default 1 is the usual
            scaling convention, under which drift and boundary separation are in units of ``s``

    Returns:
        Log-likelihood of the trial; ``-inf`` (an impossible trial) when ``rt <= non_decision_time``

    Raises:
        ValueError: If ``choice`` is not 1 or 2, ``rt`` is not finite, or a parameter is outside its range
            (see :func:`ddm_first_passage_density`)
    """
    return float(
        ddm_first_passage_density(
            rt, choice, drift_rate, boundary_separation, starting_point, non_decision_time, noise_sd, log=True
        )
    )


_PARAMETER_NAMES = ("drift_rate", "boundary_separation", "starting_point", "non_decision_time")
_DEFAULT_FIT_BOUNDS: dict[str, tuple[float, float]] = {
    "drift_rate": (0.01, 2.0),
    "boundary_separation": (0.5, 3.0),
    "starting_point": (0.2, 0.8),
    "non_decision_time": (0.0, 0.5),
}
# The likelihood is zero once the non-decision time reaches the fastest response, so the optimiser's upper
# bound for it stays this fraction below that response.
_NON_DECISION_CEILING = 1.0 - 1e-6
_FIT_STARTS = 3
_FIT_RESTARTS = 3


def fit_ddm_parameters(
    rt_data: list[float],
    choice_data: list[int],
    bounds: dict[str, tuple[float, float]] | None = None,
    *,
    noise_sd: float = 1.0,
) -> dict[str, float]:
    """Fit DDM parameters to data by maximum likelihood.

    Maximises the exact Wiener first-passage log-likelihood (:func:`ddm_first_passage_density`) summed
    over trials with respect to the drift, the boundary separation, the starting point and the
    non-decision time, within box bounds. The diffusion noise is not fitted: drift, boundary separation
    and noise are identified only up to a common scale, so ``noise_sd`` fixes it (and with it the units
    of the fitted drift and boundary separation).

    The optimiser (L-BFGS-B) is started from the three best points of a fixed grid over the bounds, so a
    given dataset always gives the same fit and no random numbers are drawn. The likelihood is zero once
    the non-decision time reaches the fastest response, so the upper bound of ``non_decision_time`` is
    lowered to just below that response whenever the supplied bound is higher.

    Args:
        rt_data: Response times, positive and finite
        choice_data: Choices (1 for the upper boundary, 2 for the lower), one per response time
        bounds: Box bounds ``(low, high)`` by parameter name (``drift_rate``, ``boundary_separation``,
            ``starting_point``, ``non_decision_time``); a parameter left out keeps its default
            (``drift_rate`` 0.01 to 2.0, ``boundary_separation`` 0.5 to 3.0, ``starting_point`` 0.2 to
            0.8, ``non_decision_time`` 0 to 0.5). The drift bounds may be negative; ``low == high``
            fixes a parameter
        noise_sd: Fixed diffusion standard deviation, positive and finite (default 1)

    Returns:
        Dictionary with the fitted ``drift_rate``, ``boundary_separation``, ``starting_point`` and
        ``non_decision_time`` and the maximised total ``log_likelihood``

    Raises:
        ValueError: If the data are empty, differ in length, hold a non-positive or non-finite response
            time or a choice other than 1 or 2; if ``noise_sd`` is not positive and finite; if a bound
            names an unknown parameter, is reversed, or lies outside the parameter's domain (boundary
            separation above 0, starting point strictly inside (0, 1), non-decision time at least 0); or
            if no non-decision time in its bounds lies below the fastest response
        RuntimeError: If the optimiser does not converge. There is no fallback to default parameters: a
            fit that did not converge is not returned
    """
    if len(rt_data) != len(choice_data):
        raise ValueError("rt_data and choice_data must have same length")

    if not len(rt_data):
        raise ValueError("No data provided")

    response_times = np.asarray(rt_data, dtype=float)
    choices = np.asarray(choice_data)
    if response_times.ndim != 1 or choices.ndim != 1:
        raise ValueError("rt_data and choice_data must be one-dimensional")
    if not np.all(np.isfinite(response_times)) or np.any(response_times <= 0.0):
        raise ValueError("Response times must be positive and finite")
    if not np.all(np.isin(choices, (1, 2))):
        raise ValueError("Choice must be 1 or 2")
    if not (math.isfinite(noise_sd) and noise_sd > 0):
        raise ValueError("Noise standard deviation must be positive and finite")

    limits = _resolve_fit_bounds(bounds, float(response_times.min()))
    upper_times = response_times[choices == 1]
    lower_times = response_times[choices == 2]
    n_trials = response_times.size

    def log_likelihood(params: np.ndarray) -> float:
        drift, boundary, start, non_decision = (float(p) for p in params)
        total = 0.0
        for times, upper in ((upper_times, True), (lower_times, False)):
            if times.size:
                total += float(
                    _log_boundary_density(times - non_decision, drift, boundary, start, noise_sd, upper).sum()
                )
        return total

    def objective(params: np.ndarray) -> float:
        # Mean negative log-likelihood per trial: unit-sized values keep finite-difference gradients accurate
        value = -log_likelihood(params) / n_trials
        if not math.isfinite(value):
            raise RuntimeError(f"The DDM log-likelihood is not finite at parameters {list(map(float, params))}")
        return value

    starts = _fit_starts(objective, limits)
    best = None
    for start in starts:
        result = minimize(objective, start, method="L-BFGS-B", bounds=limits, options={"maxiter": 1000})
        if best is None or result.fun < best.fun:
            best = result
    assert best is not None  # _fit_starts always returns at least one start

    for _ in range(_FIT_RESTARTS):
        if best.success:
            break
        # L-BFGS-B sometimes stops on a line-search failure next to the optimum; a restart rebuilds its Hessian
        best = minimize(objective, best.x, method="L-BFGS-B", bounds=limits, options={"maxiter": 1000})
    if not best.success:
        raise RuntimeError(f"DDM fit did not converge: {best.message}")

    fitted = {name: float(value) for name, value in zip(_PARAMETER_NAMES, best.x)}
    fitted["log_likelihood"] = log_likelihood(best.x)
    return fitted


def _resolve_fit_bounds(
    bounds: dict[str, tuple[float, float]] | None, fastest_response: float
) -> list[tuple[float, float]]:
    """Merge ``bounds`` with the defaults, validate them and cap the non-decision time below ``fastest_response``."""
    supplied = dict(bounds) if bounds is not None else {}
    unknown = sorted(set(supplied) - set(_PARAMETER_NAMES))
    if unknown:
        raise ValueError(f"Unknown parameter(s) in bounds: {', '.join(unknown)}; expected {list(_PARAMETER_NAMES)}")

    limits: list[tuple[float, float]] = []
    for name in _PARAMETER_NAMES:
        low, high = (float(x) for x in supplied.get(name, _DEFAULT_FIT_BOUNDS[name]))
        if not (math.isfinite(low) and math.isfinite(high)) or low > high:
            raise ValueError(f"Bounds for {name} must be finite and ordered (low <= high), got ({low}, {high})")
        if name == "boundary_separation" and low <= 0.0:
            raise ValueError("Bounds for boundary_separation must be positive")
        if name == "starting_point" and not (0.0 < low and high < 1.0):
            raise ValueError("Bounds for starting_point must lie strictly inside (0, 1)")
        if name == "non_decision_time":
            if low < 0.0:
                raise ValueError("Bounds for non_decision_time must be non-negative")
            high = min(high, fastest_response * _NON_DECISION_CEILING)
            if low > high:
                raise ValueError(
                    f"No non_decision_time within its bounds lies below the fastest response time ({fastest_response})"
                )
        limits.append((low, high))
    return limits


def _fit_starts(objective: Any, limits: list[tuple[float, float]]) -> list[np.ndarray]:
    """The best few points of a fixed grid over ``limits``, as starting points for the optimiser."""
    axes = []
    for (low, high), points in zip(limits, (5, 5, 3, 3)):
        axes.append(np.unique(np.linspace(low, high, points)))
    grid = np.array(np.meshgrid(*axes, indexing="ij")).reshape(len(axes), -1).T
    values = np.array([objective(point) for point in grid])
    order = np.argsort(values, kind="stable")[:_FIT_STARTS]
    return [grid[i] for i in order]


def island_model_update(p: float, m: float, pm: float) -> float:
    """Update allele frequency in an island model of migration.

    In an island model, each subpopulation receives migrants from a mainland
    population with frequency pm at rate m: the new frequency is the weighted
    average (1 - m) * p + m * pm. Genetic drift is not applied.

    Args:
        p: Current allele frequency in subpopulation, within [0, 1]
        m: Migration rate (proportion of population replaced by migrants), within [0, 1]
        pm: Allele frequency in mainland population, within [0, 1]

    Returns:
        New allele frequency after migration (no drift)

    Raises:
        ValueError: If p, m or pm is outside [0, 1]
    """
    if not (0 <= p <= 1):
        raise ValueError("Allele frequency must be between 0 and 1")
    if not (0 <= pm <= 1):
        raise ValueError("Migrant allele frequency must be between 0 and 1")
    if not (0 <= m <= 1):
        raise ValueError("Migration rate must be between 0 and 1")

    return (1 - m) * p + m * pm
