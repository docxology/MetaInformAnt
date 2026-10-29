"""Tests for the exact Wiener first-passage-time density and the maximum-likelihood fit built on it.

REAL IMPLEMENTATION POLICY: the real functions run against independent oracles, never mocks or patches:

- closed forms for the boundary probabilities (scale function) and for the Laplace transform of the
  absorption-time density (solution of the backward equation), integrated numerically with SciPy;
- the inverse Gaussian (Wald) density that the Wiener density approaches when the far boundary is remote;
- the leading large-time asymptote of the eigenfunction series;
- an Euler-Maruyama simulation of the diffusion with the Broadie-Glasserman-Kou barrier correction, which
  shares no code with the density;
- sampling error: a planted parameter set is recovered within four standard errors taken from the observed
  information of the likelihood, and a replicate study checks that the estimates are unbiased.

``ddm_log_likelihood`` and ``fit_ddm_parameters`` used to approximate the first-passage time by a normal
distribution that ignored the opposite boundary and the noise, gave a response against the sign of the drift
a likelihood of ``-inf`` and returned invented default parameters when no fit was found.
"""

from __future__ import annotations

import math
from typing import Callable

import numpy as np
import pytest

from metainformant.math.decision_theory import ddm

integrate = pytest.importorskip("scipy.integrate")
stats = pytest.importorskip("scipy.stats")

# (drift v, boundary separation a, starting point w = z / a, noise s)
CASES = [
    (1.0, 2.0, 0.5, 1.0),
    (0.7, 1.3, 0.3, 1.0),
    (-0.9, 1.7, 0.6, 1.0),
    (0.0, 1.0, 0.4, 1.0),
    (1.5, 0.8, 0.5, 0.5),
    (0.4, 3.0, 0.8, 2.0),
    (2.5, 1.0, 0.1, 1.0),
]
UPPER, LOWER = 1, 2


def _noise_kwargs(noise: float) -> dict[str, float]:
    return {} if noise == 1.0 else {"noise_sd": noise}


def _density(choice: int, drift: float, boundary: float, start: float, noise: float) -> Callable[[float], float]:
    """The density of a response at time ``t`` as a scalar function, through the pre-existing scalar API."""

    def f(t: float) -> float:
        return math.exp(ddm.ddm_log_likelihood(t, choice, drift, boundary, start, 0.0, **_noise_kwargs(noise)))

    return f


_EDGES = (0.0, 0.02, 0.2, 1.0, 5.0, 25.0, 200.0)


def _integral(
    function: Callable[[float], float], weight: Callable[[float], float] = lambda t: 1.0
) -> tuple[float, float]:
    """Integral of ``function(t) * weight(t)`` over (0, inf), cut at fixed times: (value, summed error estimate)."""
    total = error = 0.0
    for low, high in zip(_EDGES[:-1], _EDGES[1:]):
        value, err = integrate.quad(lambda t: function(t) * weight(t), low, high, limit=200, epsabs=1e-13, epsrel=1e-13)
        total, error = total + value, error + err
    value, err = integrate.quad(lambda t: function(t) * weight(t), _EDGES[-1], np.inf, limit=200)
    return total + value, error + err


_DECADES = tuple(10.0**k for k in range(-18, 3))


def _integral_by_decades(function: Callable[[float], float]) -> tuple[float, float]:
    """Integral over (0, inf) cut at every power of ten, which resolves the spike of a start next to a boundary."""
    total = error = 0.0
    for low, high in zip((0.0,) + _DECADES[:-1], _DECADES):
        value, err = integrate.quad(function, low, high, limit=200, epsabs=1e-17, epsrel=1e-13)
        total, error = total + value, error + err
    value, err = integrate.quad(function, _DECADES[-1], np.inf, limit=200)
    return total + value, error + err


def _boundary_probabilities(drift: float, boundary: float, start: float, noise: float) -> tuple[float, float]:
    """(P(upper), P(lower)) from the scale function, each written without cancellation so a tiny one stays exact."""
    if drift == 0.0:
        return start, 1.0 - start
    z, k = start * boundary, 2.0 * drift / noise**2
    upper = math.expm1(-k * z) / math.expm1(-k * boundary)
    lower = math.exp(-k * z) * math.expm1(-k * (boundary - z)) / math.expm1(-k * boundary)
    return upper, lower


def _upper_probability(drift: float, boundary: float, start: float, noise: float) -> float:
    """P(absorbed at the upper boundary): the scale function of dx = v dt + s dW, started at start * boundary."""
    if drift == 0.0:
        return start
    return math.expm1(-2.0 * drift * start * boundary / noise**2) / math.expm1(-2.0 * drift * boundary / noise**2)


def _laplace_transform(upper: bool, lam: float, drift: float, boundary: float, start: float, noise: float) -> float:
    """E[exp(-lam T); absorbed at this boundary], the solution of (s^2/2) phi'' + v phi' = lam phi.

    The solution is 1 at the boundary in question and 0 at the other one.
    """
    gamma = math.sqrt(drift**2 + 2.0 * noise**2 * lam) / noise**2
    z = start * boundary
    distance = z if upper else boundary - z  # to the boundary whose transform is wanted, measured from the other one
    ratio = distance / boundary if gamma == 0.0 else math.sinh(gamma * distance) / math.sinh(gamma * boundary)
    if upper:
        return math.exp(drift * (boundary - z) / noise**2) * ratio
    return math.exp(-drift * z / noise**2) * ratio


class TestBoundaryProbabilities:
    """The density of each boundary integrates to the probability of absorption there."""

    @pytest.mark.parametrize("drift,boundary,start,noise", CASES)
    def test_density_integrates_to_the_scale_function_probability(
        self, drift: float, boundary: float, start: float, noise: float
    ) -> None:
        upper, upper_error = _integral(_density(UPPER, drift, boundary, start, noise))
        lower, lower_error = _integral(_density(LOWER, drift, boundary, start, noise))
        expected = _upper_probability(drift, boundary, start, noise)
        # The series are accurate to a relative 1e-12; quad reports its own absolute error
        assert upper == pytest.approx(expected, abs=max(1e-10, 10 * upper_error))
        assert lower == pytest.approx(1.0 - expected, abs=max(1e-10, 10 * lower_error))
        assert upper + lower == pytest.approx(1.0, abs=1e-9)

    @pytest.mark.parametrize(
        "drift,boundary,noise", [(1.0, 2.0, 1.0), (0.5, 1.5, 0.7), (-1.2, 0.9, 1.3), (3.0, 1.0, 2.0)]
    )
    def test_unbiased_start_gives_expit_of_drift_times_boundary_over_noise_squared(
        self, drift: float, boundary: float, noise: float
    ) -> None:
        upper, error = _integral(_density(UPPER, drift, boundary, 0.5, noise))
        expit = 1.0 / (1.0 + math.exp(-drift * boundary / noise**2))
        assert upper == pytest.approx(expit, abs=max(1e-10, 10 * error))
        assert upper == pytest.approx(ddm.ddm_analytic_accuracy(drift, boundary, noise), abs=max(1e-10, 10 * error))

    @pytest.mark.parametrize("drift", [0.0, 0.8, -0.8])
    @pytest.mark.parametrize("start", [1e-6, 1e-3, 0.999, 1.0 - 1e-6])
    def test_a_start_next_to_a_boundary_keeps_the_probability_of_the_far_one_to_full_relative_accuracy(
        self, drift: float, start: float
    ) -> None:
        """The unlikely boundary can have probability 1e-7; the regrouped sums give it to 1e-9 relative accuracy."""
        upper, lower = _boundary_probabilities(drift, 1.3, start, 1.0)
        for choice, expected in ((UPPER, upper), (LOWER, lower)):
            value, error = _integral_by_decades(_density(choice, drift, 1.3, start, 1.0))
            # series relative error 1e-12 of the mass; quad reports its own absolute error
            assert value == pytest.approx(expected, rel=1e-9, abs=max(1e-12, 10 * error)), (choice, value, expected)

    def test_a_start_near_either_boundary_sends_the_responses_there(self) -> None:
        near_upper, _ = _integral(_density(UPPER, 0.0, 1.0, 0.95, 1.0))
        near_lower, _ = _integral(_density(LOWER, 0.0, 1.0, 0.05, 1.0))
        assert near_upper == pytest.approx(0.95, abs=1e-9)
        assert near_lower == pytest.approx(0.95, abs=1e-9)


class TestAbsorptionTimeDistribution:
    """The shape of the density, not only its mass: Laplace transform and mean decision time."""

    @pytest.mark.parametrize("drift,boundary,start,noise", CASES)
    @pytest.mark.parametrize("lam", [0.0, 0.3, 1.5, 6.0])
    def test_laplace_transform_matches_the_backward_equation_solution(
        self, drift: float, boundary: float, start: float, noise: float, lam: float
    ) -> None:
        for choice, upper in ((UPPER, True), (LOWER, False)):
            value, error = _integral(_density(choice, drift, boundary, start, noise), lambda t: math.exp(-lam * t))
            expected = _laplace_transform(upper, lam, drift, boundary, start, noise)
            assert value == pytest.approx(expected, abs=max(1e-10, 10 * error))

    @pytest.mark.parametrize(
        "drift,boundary,noise", [(1.0, 2.0, 1.0), (0.3, 1.4, 0.8), (-0.8, 1.1, 1.0), (0.0, 1.0, 1.5)]
    )
    def test_mean_decision_time_of_an_unbiased_start_matches_the_closed_form(
        self, drift: float, boundary: float, noise: float
    ) -> None:
        total = error = 0.0
        for choice in (UPPER, LOWER):
            value, err = _integral(_density(choice, drift, boundary, 0.5, noise), lambda t: t)
            total, error = total + value, error + err
        assert total == pytest.approx(ddm.ddm_mean_decision_time(drift, boundary, noise), abs=max(1e-9, 10 * error))

    def test_far_boundary_limit_is_the_inverse_gaussian_density(self) -> None:
        """With the other boundary 40 units away only the nearby one matters: the Wald (inverse Gaussian) density."""
        boundary, distance, drift, noise = 40.0, 1.2, 0.6, 1.5
        times = np.linspace(0.05, 6.0, 40)
        near_lower = ddm.ddm_first_passage_density(times, LOWER, drift, boundary, distance / boundary, noise_sd=noise)
        wald_lower = (
            distance
            / (noise * np.sqrt(2 * np.pi * times**3))
            * np.exp(-((distance + drift * times) ** 2) / (2 * noise**2 * times))
        )
        # The neglected image of the far boundary is of relative size exp(-200) here, so only rounding is left
        np.testing.assert_allclose(near_lower, wald_lower, rtol=1e-10)

        near_upper = ddm.ddm_first_passage_density(
            times, UPPER, drift, boundary, 1.0 - distance / boundary, noise_sd=noise
        )
        wald_upper = (
            distance
            / (noise * np.sqrt(2 * np.pi * times**3))
            * np.exp(-((distance - drift * times) ** 2) / (2 * noise**2 * times))
        )
        np.testing.assert_allclose(near_upper, wald_upper, rtol=1e-10)

    def test_slow_responses_follow_the_leading_eigenfunction_and_stay_finite_on_the_log_scale(self) -> None:
        """Late on only the first eigenfunction is left, decaying at drift^2 / 2 + pi^2 / (2 a^2) (units of s)."""
        drift, boundary, start, noise = 0.5, 2.0, 0.35, 1.0
        times = np.array([40.0, 400.0, 2000.0])
        for choice, away, near in ((LOWER, drift, start), (UPPER, -drift, 1.0 - start)):
            expected = (
                math.log(math.pi * noise**2 / boundary**2 * math.sin(math.pi * near))
                - away * near * boundary / noise**2
                - (away**2 / (2 * noise**2) + math.pi**2 * noise**2 / (2 * boundary**2)) * times
            )
            log_density = ddm.ddm_first_passage_density(times, choice, drift, boundary, start, noise_sd=noise, log=True)
            np.testing.assert_allclose(log_density, expected, rtol=1e-12)
            # At t = 2000 the density is below the smallest double (about exp(-2700)) but its logarithm is finite
            assert ddm.ddm_first_passage_density(2000.0, choice, drift, boundary, start) == 0.0
            assert ddm.ddm_log_likelihood(2000.0, choice, drift, boundary, start) == pytest.approx(
                expected[-1], rel=1e-12
            )

    def test_very_fast_responses_follow_the_first_image(self) -> None:
        """For t much shorter than the boundary distance squared only the direct path counts: the Wald log-density."""
        drift, boundary, start, noise = 0.8, 2.0, 0.5, 1.0
        for t in (1e-3, 4e-3, 2e-2):
            distance = start * boundary
            wald = math.log(distance / (noise * math.sqrt(2 * math.pi * t**3))) - (distance + drift * t) ** 2 / (
                2 * noise**2 * t
            )
            assert ddm.ddm_log_likelihood(t, LOWER, drift, boundary, start) == pytest.approx(wald, rel=1e-12)
        assert math.isfinite(ddm.ddm_log_likelihood(1e-6, LOWER, drift, boundary, start))


class TestSymmetriesAndShifts:
    def test_lower_boundary_is_the_mirror_image_of_the_upper_boundary(self) -> None:
        times = np.linspace(0.02, 4.0, 50)
        for drift, boundary, start, noise in CASES:
            upper = ddm.ddm_first_passage_density(times, UPPER, drift, boundary, start, noise_sd=noise)
            mirrored = ddm.ddm_first_passage_density(times, LOWER, -drift, boundary, 1.0 - start, noise_sd=noise)
            np.testing.assert_allclose(upper, mirrored, rtol=1e-12)

    def test_both_boundaries_have_one_density_when_the_diffusion_is_symmetric(self) -> None:
        times = np.linspace(0.02, 4.0, 50)
        upper = ddm.ddm_first_passage_density(times, UPPER, 0.0, 1.4, 0.5)
        lower = ddm.ddm_first_passage_density(times, LOWER, 0.0, 1.4, 0.5)
        np.testing.assert_allclose(upper, lower, rtol=1e-12)

    @pytest.mark.parametrize("factor", [0.1, 3.0, 40.0])
    def test_density_is_invariant_to_a_common_scale_of_drift_boundary_and_noise(self, factor: float) -> None:
        times = np.linspace(0.05, 5.0, 30)
        for choice in (UPPER, LOWER):
            base = ddm.ddm_first_passage_density(times, choice, 0.8, 1.5, 0.4, noise_sd=0.9)
            scaled = ddm.ddm_first_passage_density(
                times, choice, 0.8 * factor, 1.5 * factor, 0.4, noise_sd=0.9 * factor
            )
            np.testing.assert_allclose(scaled, base, rtol=1e-11)

    def test_noise_changes_the_density_and_defaults_to_one(self) -> None:
        assert ddm.ddm_log_likelihood(0.8, UPPER, 1.0, 2.0) == ddm.ddm_log_likelihood(
            0.8, UPPER, 1.0, 2.0, noise_sd=1.0
        )
        assert ddm.ddm_log_likelihood(0.8, UPPER, 1.0, 2.0, noise_sd=2.0) != pytest.approx(
            ddm.ddm_log_likelihood(0.8, UPPER, 1.0, 2.0), abs=0.05
        )

    def test_non_decision_time_shifts_the_density(self) -> None:
        for choice in (UPPER, LOWER):
            shifted = ddm.ddm_first_passage_density([0.5, 0.9, 2.0], choice, 0.9, 1.6, 0.45, non_decision_time=0.3)
            unshifted = ddm.ddm_first_passage_density([0.2, 0.6, 1.7], choice, 0.9, 1.6, 0.45)
            np.testing.assert_allclose(shifted, unshifted, rtol=1e-13)

    def test_no_density_before_the_non_decision_time(self) -> None:
        assert ddm.ddm_log_likelihood(0.3, UPPER, 1.0, 2.0, non_decision_time=0.3) == float("-inf")
        assert ddm.ddm_log_likelihood(0.1, LOWER, 1.0, 2.0, non_decision_time=0.3) == float("-inf")
        values = ddm.ddm_first_passage_density([0.1, 0.3, 0.5], UPPER, 1.0, 2.0, non_decision_time=0.3)
        assert values[0] == 0.0 and values[1] == 0.0 and values[2] > 0.0

    def test_a_response_against_the_drift_is_possible(self) -> None:
        """The opposite boundary is modelled, so error responses have a finite likelihood (it used to be -inf)."""
        assert math.isfinite(ddm.ddm_log_likelihood(1.5, LOWER, 1.0, 2.0))
        assert math.isfinite(ddm.ddm_log_likelihood(1.5, UPPER, -1.0, 2.0))
        assert math.isfinite(ddm.ddm_log_likelihood(1.5, UPPER, 0.0, 2.0))
        # ... and is rarer than the response with the drift
        assert ddm.ddm_log_likelihood(1.5, LOWER, 1.0, 2.0) < ddm.ddm_log_likelihood(1.5, UPPER, 1.0, 2.0)

    def test_scalar_and_array_forms_agree(self) -> None:
        times = np.array([0.15, 0.4, 0.9, 1.8])
        choices = np.array([1, 2, 2, 1])
        array = ddm.ddm_first_passage_density(times, choices, 0.7, 1.4, 0.55, 0.1, log=True)
        scalar = [ddm.ddm_log_likelihood(float(t), int(c), 0.7, 1.4, 0.55, 0.1) for t, c in zip(times, choices)]
        np.testing.assert_allclose(array, scalar, rtol=1e-14)
        assert isinstance(ddm.ddm_first_passage_density(0.5, 1, 0.7, 1.4), float)
        assert ddm.ddm_first_passage_density(times, 1, 0.7, 1.4).shape == (4,)
        assert ddm.ddm_first_passage_density(0.5, choices, 0.7, 1.4).shape == (4,)


def _reference_log_density(t: float, drift: float, boundary: float, start: float, noise: float, upper: bool) -> float:
    """Plain Navarro-Fuss evaluation: many terms and none of the module's regrouping (moderate starting points)."""
    near = 1.0 - start if upper else start
    away = -drift if upper else drift
    u = noise**2 * t / boundary**2
    if u < 0.25:
        k = np.arange(-60, 61)
        x = near + 2 * k
        g = float(np.sum(x * np.exp(-(x**2) / (2 * u)))) / math.sqrt(2 * math.pi * u**3)
    else:
        k = np.arange(1, 80)
        g = math.pi * float(np.sum(k * np.exp(-(k**2) * math.pi**2 * u / 2) * np.sin(k * math.pi * near)))
    return (
        math.log(noise**2 / boundary**2)
        - away * near * boundary / noise**2
        - away**2 * t / (2 * noise**2)
        + math.log(g)
    )


class TestSeries:
    def test_density_matches_a_plain_many_term_evaluation(self) -> None:
        for drift, boundary, start, noise in CASES:
            for t in (0.01, 0.05, 0.2, 0.6, 1.0, 2.5, 8.0):
                for choice, upper in ((UPPER, True), (LOWER, False)):
                    expected = _reference_log_density(t, drift, boundary, start, noise, upper)
                    value = ddm.ddm_log_likelihood(t, choice, drift, boundary, start, 0.0, **_noise_kwargs(noise))
                    # Truncation is below 1e-12 (relative) in both evaluations; 1e-9 leaves room for the rounding
                    # of the plain sum, whose pairs of images cancel by up to a factor of about 10 here
                    assert value == pytest.approx(expected, abs=1e-9), (drift, boundary, start, noise, t, choice)

    @pytest.mark.parametrize("far", [1e-9, 1e-6, 1e-3, 0.05, 0.3, 0.5, 0.7, 0.95, 1 - 1e-3, 1 - 1e-6, 1 - 1e-9])
    def test_small_and_large_time_series_agree_where_both_converge(self, far: float) -> None:
        """The two series are independent derivations; they must meet, also for a start next to either boundary.

        ``near = 1 - far`` is given with the exactly known distance ``far`` to the other boundary, which is how the
        module passes it, so the check covers the regrouped small-time sums at both extremes.
        """
        near = 1.0 - far
        u = np.linspace(0.12, 0.6, 25)
        small = ddm._log_standard_density_small_time(u, near, far, pairs=10)
        large = ddm._log_standard_density_large_time(u, near, far, terms=60)
        # Over this range both sums have more terms than needed, and neither loses more than about one digit to
        # cancellation (the conditioning factors in _series_term_counts), so they agree to a few 1e-14
        np.testing.assert_allclose(small, large, rtol=0, atol=1e-11)

    def test_one_boundary_is_much_closer_than_the_other_at_extreme_starts(self) -> None:
        """A start 1e-9 from the lower boundary sends nearly everything there; the log-density is regular, not NaN."""
        for start in (1e-9, 1e-6, 1.0 - 1e-6, 1.0 - 1e-9):
            for t in (0.001, 0.05, 0.4, 3.0):
                for choice in (UPPER, LOWER):
                    assert math.isfinite(ddm.ddm_log_likelihood(t, choice, 0.5, 1.5, start))

    @pytest.mark.parametrize("rtol", [1e-2, 1e-5, 1e-8, 1e-12])
    def test_truncation_error_stays_within_the_stated_bound(self, rtol: float) -> None:
        pairs, terms = ddm._series_term_counts(rtol)
        worst = 0.0
        for far in (1e-6, 0.01, 0.2, 0.5, 0.8, 0.99, 1 - 1e-6):
            near = 1.0 - far
            small_u = np.array([1e-3, 0.02, 0.1, 0.2, 0.2499])
            exact = ddm._log_standard_density_small_time(small_u, near, far, pairs=12)
            truncated = ddm._log_standard_density_small_time(small_u, near, far, pairs=pairs)
            worst = max(worst, float(np.max(np.abs(np.expm1(truncated - exact)))))
            large_u = np.array([0.25, 0.3, 0.5, 1.0, 3.0])
            exact = ddm._log_standard_density_large_time(large_u, near, far, terms=40)
            truncated = ddm._log_standard_density_large_time(large_u, near, far, terms=terms)
            worst = max(worst, float(np.max(np.abs(np.expm1(truncated - exact)))))
        # The bound is on the truncation; leave room for double-precision rounding of the sums themselves
        assert worst <= rtol + 1e-13

    def test_the_bound_needs_the_terms_it_asks_for(self) -> None:
        """One term fewer than the 1e-12 bound asks for really exceeds 1e-12 at the edges of the two ranges.

        A start next to the boundary reached (near = 1e-3) is where the neglected terms weigh most: their sines
        are all close to k times the first one, the largest ratio the bound has to cover.
        """
        pairs, terms = ddm._series_term_counts(1e-12)
        near, far = 1e-3, 1.0 - 1e-3
        small = np.array([0.2499])
        large = np.array([0.25])
        full_small = ddm._log_standard_density_small_time(small, near, far, pairs=12)
        full_large = ddm._log_standard_density_large_time(large, near, far, terms=40)
        short_small = ddm._log_standard_density_small_time(small, near, far, pairs=pairs - 1)
        short_large = ddm._log_standard_density_large_time(large, near, far, terms=terms - 1)
        assert abs(math.expm1(float(short_small[0] - full_small[0]))) > 1e-12
        assert abs(math.expm1(float(short_large[0] - full_large[0]))) > 1e-12

    def test_term_counts_grow_as_the_tolerance_tightens(self) -> None:
        counts = [ddm._series_term_counts(r) for r in (1e-2, 1e-4, 1e-8, 1e-12, 1e-15)]
        assert all(a[0] <= b[0] and a[1] <= b[1] for a, b in zip(counts, counts[1:]))
        assert counts[0] != counts[3]
        assert counts[3] == (ddm._SMALL_TIME_PAIRS, ddm._LARGE_TIME_TERMS)

    @pytest.mark.parametrize("rtol", [0.0, 1.0, -1e-3, 2.0])
    def test_tolerance_must_lie_strictly_between_zero_and_one(self, rtol: float) -> None:
        with pytest.raises(ValueError, match="rtol"):
            ddm._series_term_counts(rtol)


class TestValidation:
    @pytest.mark.parametrize("choice", [0, 3, -1, 1.5])
    def test_choice_must_be_one_or_two(self, choice: float) -> None:
        with pytest.raises(ValueError, match="Choice must be 1 or 2"):
            ddm.ddm_log_likelihood(1.0, choice, 1.0, 2.0)  # type: ignore[arg-type]

    def test_choice_is_checked_even_for_an_impossible_response_time(self) -> None:
        with pytest.raises(ValueError, match="Choice must be 1 or 2"):
            ddm.ddm_log_likelihood(0.1, 3, 1.0, 2.0, non_decision_time=0.5)

    @pytest.mark.parametrize("rt", [float("nan"), float("inf"), float("-inf")])
    def test_response_time_must_be_finite(self, rt: float) -> None:
        with pytest.raises(ValueError, match="finite"):
            ddm.ddm_log_likelihood(rt, 1, 1.0, 2.0)

    @pytest.mark.parametrize(
        "kwargs",
        [
            {"boundary_separation": 0.0},
            {"boundary_separation": -1.0},
            {"boundary_separation": float("nan")},
            {"drift_rate": float("nan")},
            {"drift_rate": float("inf")},
            {"noise_sd": 0.0},
            {"noise_sd": -1.0},
            {"starting_point": 0.0},
            {"starting_point": 1.0},
            {"starting_point": -0.2},
            {"starting_point": 1.2},
            {"starting_point": float("nan")},
            {"non_decision_time": -0.1},
            {"non_decision_time": float("nan")},
        ],
    )
    def test_parameters_outside_their_domain_raise(self, kwargs: dict) -> None:
        params = {"rt": 1.0, "choice": 1, "drift_rate": 1.0, "boundary_separation": 2.0} | kwargs
        with pytest.raises(ValueError):
            ddm.ddm_log_likelihood(**params)

    def test_trial_and_density_functions_validate_alike(self) -> None:
        with pytest.raises(ValueError):
            ddm.ddm_first_passage_density([0.5, 1.0], [1, 3], 1.0, 2.0)
        with pytest.raises(ValueError):
            ddm.ddm_first_passage_density([0.5, float("nan")], 1, 1.0, 2.0)


def _simulate_trials(
    drift: float, boundary: float, start: float, noise: float, n_trials: int, dt: float, seed: int
) -> tuple[np.ndarray, np.ndarray]:
    """Euler-Maruyama paths from ``start * boundary`` until they leave (0, boundary): (times, hit_upper).

    Stepping misses excursions between steps, so each barrier is pulled in by 0.5826 s sqrt(dt) (Broadie,
    Glasserman and Kou 1997), which makes the stepped walk match the continuous process to O(dt). The
    simulation shares no code with the density.
    """
    rng = np.random.default_rng(seed)
    shift = 0.5826 * noise * math.sqrt(dt)
    lower, upper = shift, boundary - shift
    step_sd = noise * math.sqrt(dt)
    position = np.full(n_trials, start * boundary)
    elapsed = np.zeros(n_trials)
    hit_upper = np.zeros(n_trials, dtype=bool)
    active = np.arange(n_trials)
    while active.size:
        position[active] += drift * dt + step_sd * rng.standard_normal(active.size)
        elapsed[active] += dt
        above = position[active] >= upper
        below = position[active] <= lower
        hit_upper[active[above]] = True
        active = active[~(above | below)]
    return elapsed, hit_upper


class TestAgainstSimulation:
    @pytest.mark.parametrize(
        "drift,boundary,start,noise", [(0.8, 1.5, 0.5, 1.0), (1.2, 2.0, 0.35, 1.0), (-0.6, 1.2, 0.6, 0.8)]
    )
    def test_simulated_trials_follow_the_density(
        self, drift: float, boundary: float, start: float, noise: float
    ) -> None:
        """Chi-square test of (choice, time-bin) counts against the cell probabilities integrated from the density."""
        n_trials, dt = 40000, 1e-3
        times, hit_upper = _simulate_trials(drift, boundary, start, noise, n_trials, dt, seed=17)

        # Bin edges come from an independent pilot run, so they do not depend on the data under test
        pilot, _ = _simulate_trials(drift, boundary, start, noise, 20000, dt, seed=5)
        edges = np.concatenate([[0.0], np.quantile(pilot, np.linspace(0.08, 0.92, 11)), [np.inf]])

        observed, expected = [], []
        for choice, selected in ((UPPER, hit_upper), (LOWER, ~hit_upper)):
            counts, _ = np.histogram(times[selected], bins=edges)
            observed.extend(counts)
            density = _density(choice, drift, boundary, start, noise)
            for low, high in zip(edges[:-1], edges[1:]):
                mass, _ = integrate.quad(density, low, high if np.isfinite(high) else 200.0, limit=200)
                expected.append(n_trials * mass)
        assert min(expected) > 100
        assert sum(expected) == pytest.approx(n_trials, rel=1e-6)
        expected = list(np.array(expected) * (n_trials / sum(expected)))
        assert stats.chisquare(observed, expected).pvalue > 1e-3


def _negative_log_likelihood(theta: np.ndarray, times: np.ndarray, choices: np.ndarray, noise: float = 1.0) -> float:
    return -float(np.sum(ddm.ddm_first_passage_density(times, choices, *theta, noise, log=True)))


def _standard_errors(theta: np.ndarray, times: np.ndarray, choices: np.ndarray, noise: float = 1.0) -> np.ndarray:
    """Asymptotic standard errors: the inverse of the observed information (central-difference Hessian of -log L)."""
    n = len(theta)
    step = 1e-4 * (1.0 + np.abs(theta))
    hessian = np.zeros((n, n))
    centre = _negative_log_likelihood(theta, times, choices, noise)
    for i in range(n):
        for j in range(i, n):
            ei, ej = np.zeros(n), np.zeros(n)
            ei[i], ej[j] = step[i], step[j]
            if i == j:
                value = (
                    _negative_log_likelihood(theta + ei, times, choices, noise)
                    - 2.0 * centre
                    + _negative_log_likelihood(theta - ei, times, choices, noise)
                ) / step[i] ** 2
            else:
                value = (
                    _negative_log_likelihood(theta + ei + ej, times, choices, noise)
                    - _negative_log_likelihood(theta + ei - ej, times, choices, noise)
                    - _negative_log_likelihood(theta - ei + ej, times, choices, noise)
                    + _negative_log_likelihood(theta - ei - ej, times, choices, noise)
                ) / (4.0 * step[i] * step[j])
            hessian[i, j] = hessian[j, i] = value
    return np.sqrt(np.diag(np.linalg.inv(hessian)))


PARAMETERS = ("drift_rate", "boundary_separation", "starting_point", "non_decision_time")


def _fit(times: np.ndarray, choices: np.ndarray, **kwargs: object) -> dict[str, float]:
    return ddm.fit_ddm_parameters(list(times), [int(c) for c in choices], **kwargs)  # type: ignore[arg-type]


def _as_array(fit: dict[str, float]) -> np.ndarray:
    return np.array([fit[name] for name in PARAMETERS])


def _simulate_data(truth: tuple[float, float, float, float], n_trials: int, seed: int, noise: float = 1.0):
    drift, boundary, start, non_decision = truth
    times, hit_upper = _simulate_trials(drift, boundary, start, noise, n_trials, 1e-3, seed)
    return times + non_decision, np.where(hit_upper, UPPER, LOWER)


class TestFit:
    @pytest.mark.parametrize(
        "truth,seed",
        [
            ((1.0, 1.6, 0.45, 0.25), 1),
            ((0.8, 1.2, 0.50, 0.30), 2),
            ((1.4, 2.0, 0.60, 0.20), 3),
            ((0.5, 1.0, 0.35, 0.15), 4),
        ],
    )
    def test_planted_parameters_are_recovered_within_four_standard_errors(
        self, truth: tuple[float, float, float, float], seed: int
    ) -> None:
        times, choices = _simulate_data(truth, 3000, seed)
        fit = _fit(times, choices)
        estimate = _as_array(fit)
        errors = _standard_errors(estimate, times, choices)
        z = (estimate - np.array(truth)) / errors
        assert np.all(np.abs(z) < 4.0), dict(zip(PARAMETERS, z))
        # The error bars themselves are informative: a few percent of the planted values
        assert np.all(errors < 0.25 * np.array(truth))

    @pytest.mark.parametrize("truth,seed", [((1.0, 1.6, 0.45, 0.25), 11), ((0.7, 1.1, 0.55, 0.35), 12)])
    def test_the_fit_is_at_least_as_likely_as_the_planted_parameters(
        self, truth: tuple[float, float, float, float], seed: int
    ) -> None:
        """Deterministic: the planted point is inside the bounds, so no maximum lies below its log-likelihood."""
        times, choices = _simulate_data(truth, 1500, seed)
        fit = _fit(times, choices)
        planted = -_negative_log_likelihood(np.array(truth), times, choices)
        assert fit["log_likelihood"] >= planted - 1e-8

    def test_reported_log_likelihood_is_the_sum_over_trials(self) -> None:
        times, choices = _simulate_data((1.0, 1.6, 0.45, 0.25), 800, 21)
        fit = _fit(times, choices)
        total = sum(
            ddm.ddm_log_likelihood(
                float(t),
                int(c),
                fit["drift_rate"],
                fit["boundary_separation"],
                fit["starting_point"],
                fit["non_decision_time"],
            )
            for t, c in zip(times, choices)
        )
        assert fit["log_likelihood"] == pytest.approx(total, rel=1e-12)

    def test_fit_sits_at_a_maximum(self) -> None:
        """Moving any parameter by a small step lowers the log-likelihood (all four lie inside their bounds here)."""
        times, choices = _simulate_data((1.0, 1.6, 0.45, 0.25), 1500, 22)
        fit = _fit(times, choices)
        estimate = _as_array(fit)
        errors = _standard_errors(estimate, times, choices)
        best = -_negative_log_likelihood(estimate, times, choices)
        for i in range(4):
            for sign in (-1.0, 1.0):
                moved = estimate.copy()
                moved[i] += sign * 0.5 * errors[i]
                assert -_negative_log_likelihood(moved, times, choices) < best

    def test_replicate_study_shows_no_bias(self) -> None:
        """24 datasets of 1500 trials: each parameter's standardised errors have mean 0 and standard deviation 1."""
        truth = (1.0, 1.6, 0.45, 0.25)
        replicates = 24
        z = np.empty((replicates, 4))
        for r in range(replicates):
            times, choices = _simulate_data(truth, 1500, 100 + r)
            estimate = _as_array(_fit(times, choices))
            z[r] = (estimate - np.array(truth)) / _standard_errors(estimate, times, choices)
        # The mean of 24 standard normals has standard error 0.2; allow 4 of them
        assert np.all(np.abs(z.mean(axis=0)) < 4 / math.sqrt(replicates)), z.mean(axis=0)
        assert np.all((z.std(axis=0, ddof=1) > 0.55) & (z.std(axis=0, ddof=1) < 1.6)), z.std(axis=0, ddof=1)

    def test_negative_drift_is_recovered_with_wide_drift_bounds(self) -> None:
        truth = (-0.7, 1.4, 0.55, 0.2)
        times, choices = _simulate_data(truth, 3000, 31)
        assert (choices == LOWER).mean() > 0.6
        fit = _fit(times, choices, bounds={"drift_rate": (-2.0, 2.0), "starting_point": (0.1, 0.9)})
        estimate = _as_array(fit)
        z = (estimate - np.array(truth)) / _standard_errors(estimate, times, choices)
        assert np.all(np.abs(z) < 4.0), dict(zip(PARAMETERS, z))

    def test_noise_sd_fixes_the_scale_of_drift_and_boundary(self) -> None:
        """``noise_sd=0.5`` on the same data gives the s = 1 fit with the drift and the boundary halved."""
        times, choices = _simulate_data((1.0, 1.6, 0.45, 0.25), 1500, 41)
        unit = _fit(times, choices)
        half = _fit(
            times,
            choices,
            noise_sd=0.5,
            bounds={"drift_rate": (0.005, 1.0), "boundary_separation": (0.25, 1.5)},
        )
        assert half["drift_rate"] == pytest.approx(0.5 * unit["drift_rate"], rel=2e-3)
        assert half["boundary_separation"] == pytest.approx(0.5 * unit["boundary_separation"], rel=2e-3)
        assert half["starting_point"] == pytest.approx(unit["starting_point"], abs=2e-3)
        assert half["non_decision_time"] == pytest.approx(unit["non_decision_time"], abs=2e-3)
        assert half["log_likelihood"] == pytest.approx(unit["log_likelihood"], abs=1e-4)

    def test_fit_is_deterministic(self) -> None:
        times, choices = _simulate_data((1.0, 1.6, 0.45, 0.25), 600, 51)
        assert _fit(times, choices) == _fit(times, choices)

    def test_missing_bounds_keep_their_defaults(self) -> None:
        times, choices = _simulate_data((1.0, 1.6, 0.45, 0.25), 600, 52)
        partial = _fit(times, choices, bounds={"drift_rate": (0.01, 2.0)})
        assert partial == _fit(times, choices)

    def test_a_parameter_can_be_fixed_with_equal_bounds(self) -> None:
        times, choices = _simulate_data((1.0, 1.6, 0.45, 0.25), 800, 53)
        fit = _fit(times, choices, bounds={"starting_point": (0.5, 0.5)})
        assert fit["starting_point"] == 0.5

    def test_the_non_decision_time_stays_below_the_fastest_response(self) -> None:
        """Default bounds allow 0.5 s, but the fastest response below is 0.2 s: the likelihood is zero past it."""
        times, choices = _simulate_data((1.4, 1.0, 0.5, 0.05), 500, 54)
        assert times.min() < 0.5
        fit = _fit(times, choices)
        assert 0.0 <= fit["non_decision_time"] < times.min()
        assert math.isfinite(fit["log_likelihood"])

    def test_impossible_bounds_raise_instead_of_returning_defaults(self) -> None:
        """No fit exists when every allowed non-decision time exceeds the fastest response (defaults were returned)."""
        with pytest.raises(ValueError, match="fastest response"):
            ddm.fit_ddm_parameters([0.1, 0.2, 0.3], [1, 2, 1], bounds={"non_decision_time": (0.5, 0.8)})

    @pytest.mark.parametrize(
        "args,message",
        [
            (([1.0, 2.0], [1]), "same length"),
            (([], []), "No data"),
            (([0.5, 0.6], [1, 3]), "Choice must be 1 or 2"),
            (([0.5, -0.6], [1, 2]), "positive and finite"),
            (([0.5, float("nan")], [1, 2]), "positive and finite"),
            (([0.5, 0.0], [1, 2]), "positive and finite"),
        ],
    )
    def test_data_are_validated(self, args: tuple[list, list], message: str) -> None:
        with pytest.raises(ValueError, match=message):
            ddm.fit_ddm_parameters(*args)

    @pytest.mark.parametrize(
        "bounds,message",
        [
            ({"drift": (0.0, 1.0)}, "Unknown parameter"),
            ({"drift_rate": (1.0, 0.5)}, "ordered"),
            ({"boundary_separation": (0.0, 2.0)}, "boundary_separation must be positive"),
            ({"starting_point": (0.0, 0.5)}, "strictly inside"),
            ({"starting_point": (0.5, 1.0)}, "strictly inside"),
            ({"non_decision_time": (-0.1, 0.3)}, "non-negative"),
            ({"drift_rate": (float("nan"), 1.0)}, "finite"),
        ],
    )
    def test_bounds_are_validated(self, bounds: dict, message: str) -> None:
        with pytest.raises(ValueError, match=message):
            ddm.fit_ddm_parameters([0.5, 0.7, 0.9], [1, 2, 1], bounds=bounds)

    @pytest.mark.parametrize("noise", [0.0, -1.0, float("nan"), float("inf")])
    def test_noise_must_be_positive_and_finite(self, noise: float) -> None:
        with pytest.raises(ValueError, match="Noise"):
            ddm.fit_ddm_parameters([0.5, 0.7, 0.9], [1, 2, 1], noise_sd=noise)

    def test_data_that_do_not_come_from_a_diffusion_still_give_a_bounded_maximum(self) -> None:
        """The pre-existing smoke-test inputs (exponential response times, random choices) fit within the bounds."""
        rng = np.random.default_rng(42)
        times = rng.exponential(0.5, size=20) + 0.3
        choices = [1 if rng.random() > 0.3 else 2 for _ in range(20)]
        fit = ddm.fit_ddm_parameters(list(times), choices)
        assert 0.01 <= fit["drift_rate"] <= 2.0
        assert 0.5 <= fit["boundary_separation"] <= 3.0
        assert 0.2 <= fit["starting_point"] <= 0.8
        assert 0.0 <= fit["non_decision_time"] < times.min()
        assert math.isfinite(fit["log_likelihood"])


class TestIslandModel:
    """``island_model_update`` promised drift in its docstring and accepted any frequency or rate."""

    def test_weighted_average_of_the_two_frequencies(self) -> None:
        assert ddm.island_model_update(0.2, 0.25, 0.8) == pytest.approx(0.75 * 0.2 + 0.25 * 0.8)

    @pytest.mark.parametrize(
        "p,m,pm", [(-0.1, 0.1, 0.5), (1.1, 0.1, 0.5), (0.5, -0.1, 0.5), (0.5, 1.1, 0.5), (0.5, 0.1, 1.2)]
    )
    def test_values_outside_the_unit_interval_raise(self, p: float, m: float, pm: float) -> None:
        with pytest.raises(ValueError):
            ddm.island_model_update(p, m, pm)


class TestExtremeResponseTimes:
    def test_an_absurdly_fast_response_has_zero_density_without_warnings(self) -> None:
        """Times below 1e-300 (in units of a^2 / s^2) underflow every double: reported as -inf, not NaN."""
        assert ddm.ddm_log_likelihood(1e-310, UPPER, 1.0, 2.0) == float("-inf")
        assert ddm.ddm_first_passage_density(1e-310, LOWER, 1.0, 2.0) == 0.0
        values = ddm.ddm_first_passage_density([1e-310, 0.5], UPPER, 1.0, 2.0, log=True)
        assert values[0] == float("-inf") and math.isfinite(values[1])
