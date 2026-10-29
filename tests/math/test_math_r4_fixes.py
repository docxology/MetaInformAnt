"""Regression tests for math defects fixed in the round-4 audit.

REAL IMPLEMENTATION POLICY: every test runs the real functions and checks them against independent oracles
(hand computations, closed forms, exact enumeration, the standard library and a small coalescent
simulation); nothing is mocked or patched.

Covered here:

- ``effective_size_from_family_size_variance`` computed ``k / (V - 1)`` in list mode, which is not an
  estimator of the effective size. It now returns the standard ``(4N - 2) / (Vk + 2)``.
- ``permutation_test`` answered an empty group with a p-value of 1.0 and could return a p-value of 0.
- ``compare_population_statistic`` answered a zero standard error with a p-value of 1.0.
- ``conjugate_normal`` invented a variance of 1.0 for one observation and of 1e-10 for identical ones.
- ``conjugate_beta_binomial`` called a normal approximation the 95% credible interval.
- ``compute_dic`` reported the mean deviance as the "deviance at the posterior mean".
- ``metropolis_hastings`` accepted NaN log-posteriors as improvements.
- ``ld_coefficients``, ``r_squared`` and ``multilevel_selection_decomposition`` died with a bare
  ``AssertionError`` on valid input types.
- ``hardy_weinberg_genotype_freqs`` answered an impossible frequency with zeros.
- the egt ``replicator_derivative`` returned zeros whenever the mean fitness was 0, and ``replicator_step``
  invented a uniform population when every frequency was zero.
- ``sdt_metrics`` invented a trial count of 1 for an empty trial class; ``d_prime`` clipped rates above 1.
- ``expected_time_to_mrca`` answered a non-positive effective size with 0.0; ``watterson_theta`` and
  ``expected_total_branch_length`` replaced the harmonic number by ``log(n) + gamma``.
- the statistics ``expected_pairwise_diversity`` returned ``theta * n / (n - 1) * H_{n-1}`` instead of theta,
  and ``expected_coalescent_waiting_times`` ignored its effective size.
- ``detect_outliers`` divided by a zero standard deviation; ``haldane_d_to_c`` accepted negative distances.
- ``island_model_update`` of the decision-theory module and the ``covariance`` docstring promised more
  than the code does.
"""

from __future__ import annotations

import doctest
import itertools
import math
import os
import statistics as std
import subprocess
import sys
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pytest

import metainformant
from metainformant.math.bayesian.inference import (
    compute_dic,
    conjugate_beta_binomial,
    conjugate_normal,
    metropolis_hastings,
)
from metainformant.math.core.utilities import covariance, r_squared
from metainformant.math.evolutionary_dynamics import egt
from metainformant.math.perception.signal_detection import criterion_c, d_prime, sdt_metrics
from metainformant.math.population_genetics import coalescent, statistics
from metainformant.math.population_genetics.core import hardy_weinberg_genotype_freqs
from metainformant.math.population_genetics.ld import haldane_d_to_c, ld_coefficients
from metainformant.math.population_genetics.selection import multilevel_selection_decomposition
from metainformant.math.population_genetics.statistics import (
    compare_population_statistic,
    detect_outliers,
    effective_size_from_family_size_variance,
    permutation_test,
)

scipy_stats = pytest.importorskip("scipy.stats")


def _docstring_examples_pass(function: Callable[..., Any], **extra_globals: Any) -> bool:
    """Run the ``>>>`` examples of ``function``'s docstring and report whether every one of them holds."""
    finder = doctest.DocTestFinder(recurse=False)
    runner = doctest.DocTestRunner(verbose=False, optionflags=doctest.ELLIPSIS)
    tests = finder.find(function, function.__name__, globs={function.__name__: function, **extra_globals})
    for test in tests:
        runner.run(test)
    return runner.failures == 0 and runner.tries > 0


# ---------------------------------------------------------------------------
# effective_size_from_family_size_variance: the standard (4N - 2) / (Vk + 2)
# ---------------------------------------------------------------------------


class TestFamilySizeVarianceEstimator:
    @pytest.mark.parametrize(
        "counts,expected",
        [
            # N = 4, every parent leaves 2 offspring: Vk = 0, Ne = (16 - 2) / 2
            ([2, 2, 2, 2], 7.0),
            # N = 5, deviations -2..2 give a sum of squares of 10, Vk = 10 / 5 = 2: Ne = 18 / 4
            ([0, 1, 2, 3, 4], 4.5),
            # N = 2, mean 2, Vk = (1 + 1) / 2 = 1: Ne = 6 / 3
            ([1, 3], 2.0),
            # N = 2, Vk = (4 + 4) / 2 = 4: Ne = 6 / 6
            ([0, 4], 1.0),
            # N = 6, mean 2, sum of squares 4, Vk = 2/3: Ne = 22 / (2/3 + 2) = 8.25
            ([1, 2, 3, 2, 1, 3], 8.25),
            # N = 1, a single parent: Vk = 0, Ne = 2 / 2
            ([3], 1.0),
            # N = 4, mean 2, sum of squares 2, Vk = 0.5: Ne = 14 / 2.5
            ([1, 3, 2, 2], 5.6),
        ],
    )
    def test_hand_computed_values(self, counts: list[int], expected: float) -> None:
        assert effective_size_from_family_size_variance(counts) == pytest.approx(expected, rel=1e-12)

    def test_the_old_expression_k_over_v_minus_one_is_gone(self) -> None:
        """k / (V - 1) gave inf for these counts (V <= 1); the standard estimator is finite, above N for equal sizes."""
        assert effective_size_from_family_size_variance([1, 3, 2, 2]) != float("inf")
        assert effective_size_from_family_size_variance([2, 2, 2, 2]) > 4  # equal family sizes: Ne = 2N - 1 > N

    def test_variance_is_the_population_variance_of_the_parents(self) -> None:
        """[0, 1, 2, 3, 4]: the population variance 2.0 gives Ne = 4.5; the sample variance 2.5 would give 4.0."""
        counts = [0, 1, 2, 3, 4]
        assert np.var(counts, ddof=0) == 2.0 and np.var(counts, ddof=1) == 2.5
        assert effective_size_from_family_size_variance(counts) == pytest.approx(18 / 4)
        assert effective_size_from_family_size_variance(counts) != pytest.approx(18 / 4.5)

    @pytest.mark.parametrize("counts", [[0, 1, 2, 3, 4], [1, 2, 3, 2, 1, 3], [5, 0, 1, 1, 2, 3, 0, 4], [0.5, 2.5, 3.0]])
    def test_list_form_is_the_scalar_form_at_the_list_length_and_population_variance(self, counts: list[float]) -> None:
        scalar = effective_size_from_family_size_variance(float(len(counts)), Vk=float(np.var(counts)))
        assert effective_size_from_family_size_variance(counts) == pytest.approx(scalar, rel=1e-13)

    @pytest.mark.parametrize("counts", [[0, 1, 2, 3, 4], [1, 2, 3, 2, 1, 3]])
    def test_agrees_with_the_crow_kimura_general_formula_when_the_mean_is_two(self, counts: list[int]) -> None:
        """(N k - 1) / (k - 1 + Vk / k) at the mean family size k = 2 is (4N - 2) / (Vk + 2)."""
        n, k, v = len(counts), float(np.mean(counts)), float(np.var(counts))
        assert k == 2.0
        assert effective_size_from_family_size_variance(counts) == pytest.approx((n * k - 1) / (k - 1 + v / k))

    def test_only_the_number_and_the_spread_of_the_counts_matter(self) -> None:
        """The formula assumes a mean family size of 2; a shifted or reordered list has the same N and Vk."""
        base = effective_size_from_family_size_variance([1, 3, 2, 2])
        assert effective_size_from_family_size_variance([2, 4, 3, 3]) == pytest.approx(base)
        assert effective_size_from_family_size_variance([3, 2, 1, 2]) == pytest.approx(base)

    def test_more_variance_in_family_size_lowers_the_effective_size(self) -> None:
        even = effective_size_from_family_size_variance([2] * 10)
        skewed = effective_size_from_family_size_variance([0, 0, 0, 0, 0, 0, 5, 5, 5, 5])
        assert skewed < even

    def test_input_types(self) -> None:
        assert effective_size_from_family_size_variance((1, 3, 2, 2)) == pytest.approx(5.6)
        assert effective_size_from_family_size_variance(np.array([1, 3, 2, 2])) == pytest.approx(5.6)
        assert effective_size_from_family_size_variance([]) == 0.0

    @pytest.mark.parametrize("counts", [[1, -1, 2], [1, float("nan")], [1, float("inf")], [[1, 2], [3, 4]]])
    def test_impossible_counts_raise(self, counts: list) -> None:
        with pytest.raises(ValueError):
            effective_size_from_family_size_variance(counts)

    @pytest.mark.parametrize("parents", [0.0, 0.5, -3, float("nan")])
    def test_fewer_than_one_parent_raises(self, parents: float) -> None:
        """(4N - 2) / (Vk + 2) is negative for N below 1/2, which is not a population size."""
        with pytest.raises(ValueError):
            effective_size_from_family_size_variance(parents, Vk=1.0)


# ---------------------------------------------------------------------------
# permutation_test and compare_population_statistic: no invented p-values
# ---------------------------------------------------------------------------


def _exact_permutation_p_value(first: list[float], second: list[float]) -> float:
    """Exact two-sided p-value of the difference of means by enumerating every split of the pooled values."""
    pooled = first + second
    n1 = len(first)
    observed = abs(np.mean(first) - np.mean(second))
    hits = total = 0
    for chosen in itertools.combinations(range(len(pooled)), n1):
        group = [pooled[i] for i in chosen]
        rest = [pooled[i] for i in range(len(pooled)) if i not in chosen]
        hits += abs(np.mean(group) - np.mean(rest)) >= observed - 1e-12
        total += 1
    return hits / total


class TestPermutationTest:
    @pytest.mark.parametrize("first,second", [([], [1.0, 2.0]), ([1.0, 2.0], []), ([], [])])
    def test_an_empty_group_has_no_p_value(self, first: list[float], second: list[float]) -> None:
        """The old answer for an empty group was a p-value of 1.0 and a statistic of 0.0."""
        with pytest.raises(ValueError, match="at least one value"):
            permutation_test(first, second, n_permutations=100)

    @pytest.mark.parametrize("n_permutations", [0, -5])
    def test_at_least_one_permutation_is_needed(self, n_permutations: int) -> None:
        with pytest.raises(ValueError, match="n_permutations"):
            permutation_test([1.0, 2.0], [3.0, 4.0], n_permutations=n_permutations)

    def test_a_monte_carlo_p_value_is_never_zero(self) -> None:
        """No permutation of well separated groups reaches the observed difference, yet p is 1 / (n + 1), not 0."""
        group1 = [float(x) for x in range(1, 11)]
        group2 = [float(x) for x in range(100, 111)]
        result = permutation_test(group1, group2, n_permutations=500, random_state=1)
        assert result["p_value"] == pytest.approx(1 / 501)

    @pytest.mark.parametrize(
        "first,second",
        [
            ([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]),
            ([1.0, 2.0, 3.0, 10.0], [4.0, 5.0, 6.0]),
            ([2.0, 4.0], [3.0, 5.0, 4.0, 6.0]),
        ],
    )
    def test_estimates_the_exact_permutation_p_value(self, first: list[float], second: list[float]) -> None:
        n_permutations = 20000
        exact = _exact_permutation_p_value(first, second)
        estimate = permutation_test(first, second, n_permutations=n_permutations, random_state=7)["p_value"]
        standard_error = math.sqrt(exact * (1 - exact) / n_permutations) + 1.0 / (n_permutations + 1)
        assert abs(estimate - exact) < 4 * standard_error

    def test_is_reproducible_with_a_seed(self) -> None:
        a = permutation_test([1.0, 2.5, 3.0], [2.0, 6.0, 7.5], n_permutations=300, random_state=11)
        b = permutation_test([1.0, 2.5, 3.0], [2.0, 6.0, 7.5], n_permutations=300, random_state=11)
        assert a == b


class TestComparePopulationStatistic:
    def test_zero_standard_error_is_an_error_not_a_p_value_of_one(self) -> None:
        """Two populations with no variance used to give t = 0 and p = 1.0 however far apart their statistics were."""
        with pytest.raises(ValueError, match="zero"):
            compare_population_statistic({"pi": 0.9, "variance": 0.0}, {"pi": 0.1, "variance": 0.0}, "pi")

    def test_negative_variance_is_rejected(self) -> None:
        with pytest.raises(ValueError, match="non-negative"):
            compare_population_statistic({"pi": 0.9, "variance": -0.1}, {"pi": 0.1, "variance": 0.2}, "pi")

    def test_p_value_is_the_normal_tail_of_the_documented_standard_error(self) -> None:
        """SE = sqrt(mean variance * 2 / 10) (ten observations assumed per population), two-sided normal p-value."""
        result = compare_population_statistic({"pi": 0.5, "variance": 0.2}, {"pi": 0.1, "variance": 0.4}, "pi")
        se = math.sqrt(0.3 * 2 / 10)
        t = 0.4 / se
        assert result["test_statistic"] == pytest.approx(t)
        assert result["p_value"] == pytest.approx(math.erfc(t / math.sqrt(2)), rel=1e-12)

    def test_one_zero_variance_is_fine(self) -> None:
        result = compare_population_statistic({"pi": 0.5, "variance": 0.0}, {"pi": 0.1, "variance": 0.2}, "pi")
        assert 0.0 <= result["p_value"] <= 1.0


# ---------------------------------------------------------------------------
# bayesian.inference
# ---------------------------------------------------------------------------


class TestConjugateNormalInventsNoVariance:
    def test_one_observation_has_no_sample_variance(self) -> None:
        """A variance of 1.0 used to be substituted for it."""
        with pytest.raises(ValueError, match="known_var"):
            conjugate_normal([4.0])

    def test_identical_observations_have_a_zero_sample_variance(self) -> None:
        """A variance of 1e-10 used to be substituted for it, giving a posterior sd of 1e-5."""
        with pytest.raises(ValueError, match="known_var"):
            conjugate_normal([4.0, 4.0, 4.0])

    @pytest.mark.parametrize("known_var", [0.0, -1.0, float("nan")])
    def test_a_known_variance_must_be_positive(self, known_var: float) -> None:
        with pytest.raises(ValueError, match="known_var"):
            conjugate_normal([1.0, 2.0, 3.0], known_var=known_var)

    def test_posterior_matches_the_hand_computation(self) -> None:
        # prior N(0, 4): precision 0.25; data [1, 2, 3], known variance 2: n / sigma^2 = 1.5 and mean 2
        result = conjugate_normal([1.0, 2.0, 3.0], prior_mean=0.0, prior_var=4.0, known_var=2.0)
        assert result["posterior_var"] == pytest.approx(1 / 1.75)
        assert result["posterior_mean"] == pytest.approx(1.5 * 2.0 / 1.75)

    def test_plug_in_variance_is_the_sample_variance(self) -> None:
        plug_in = conjugate_normal([1.0, 2.0, 3.0], prior_var=9.0)
        known = conjugate_normal([1.0, 2.0, 3.0], prior_var=9.0, known_var=1.0)  # sample variance of 1, 2, 3 is 1
        assert plug_in == known


class TestBetaBinomialInterval:
    def test_interval_is_exact_for_a_posterior_with_a_closed_form_quantile(self) -> None:
        """0 of 5 successes, uniform prior: Beta(1, 6) has CDF 1 - (1 - x)^6, so its q-quantile is 1 - (1 - q)^(1/6)."""
        low, high = conjugate_beta_binomial(0, 5)["credible_interval_95"]
        assert low == pytest.approx(1 - 0.975 ** (1 / 6), rel=1e-9)
        assert high == pytest.approx(1 - 0.025 ** (1 / 6), rel=1e-9)
        # the old normal approximation put the limit at 0.385 and clipped the lower limit to exactly 0
        assert low > 0.0 and high > 0.45

    def test_interval_is_exact_for_the_other_closed_form(self) -> None:
        """4 successes in 4 trials, uniform prior: Beta(5, 1) has CDF x^5, so its quantile q is q^(1/5)."""
        low, high = conjugate_beta_binomial(4, 4)["credible_interval_95"]
        assert low == pytest.approx(0.025 ** (1 / 5), rel=1e-9)
        assert high == pytest.approx(0.975 ** (1 / 5), rel=1e-9)

    def test_interval_holds_95_percent_of_the_posterior_mass(self) -> None:
        result = conjugate_beta_binomial(13, 40, prior_alpha=2.0, prior_beta=3.0)
        low, high = result["credible_interval_95"]
        posterior = scipy_stats.beta(result["posterior_alpha"], result["posterior_beta"])
        assert posterior.cdf(high) - posterior.cdf(low) == pytest.approx(0.95, abs=1e-9)
        assert posterior.cdf(low) == pytest.approx(0.025, abs=1e-9)

    @pytest.mark.parametrize("alpha,beta", [(0.0, 1.0), (1.0, 0.0), (-1.0, 2.0)])
    def test_prior_parameters_must_be_positive(self, alpha: float, beta: float) -> None:
        with pytest.raises(ValueError, match="positive"):
            conjugate_beta_binomial(3, 10, prior_alpha=alpha, prior_beta=beta)


class TestComputeDic:
    LOG_LIKELIHOODS = [-10.0, -12.0, -11.0, -9.0]  # deviances 20, 24, 22, 18: mean 21, population variance 5
    PARAMETERS = [[0.1], [0.2], [0.3], [0.4]]

    def test_deviance_at_the_mean_is_not_the_mean_deviance(self) -> None:
        """Without the likelihood function the deviance at the posterior mean is unknown (it was the mean deviance)."""
        result = compute_dic(self.LOG_LIKELIHOODS, self.PARAMETERS)
        assert result["mean_deviance"] == pytest.approx(21.0)
        assert result["deviance_at_mean"] is None
        assert result["effective_parameters"] == pytest.approx(2.5)  # var(D) / 2
        assert result["dic"] == pytest.approx(23.5)

    def test_classic_p_d_when_the_log_likelihood_at_the_mean_is_given(self) -> None:
        result = compute_dic(self.LOG_LIKELIHOODS, self.PARAMETERS, log_likelihood_at_mean=-9.5)
        assert result["deviance_at_mean"] == pytest.approx(19.0)
        assert result["effective_parameters"] == pytest.approx(21.0 - 19.0)
        assert result["dic"] == pytest.approx(21.0 + 2.0)


class TestMetropolisHastingsGuards:
    @staticmethod
    def _log_posterior_nan_outside(params: list[float]) -> float:
        return -0.5 * params[0] ** 2 if abs(params[0]) < 1.0 else float("nan")

    def test_a_nan_log_posterior_is_an_error_not_an_improvement(self) -> None:
        """NaN used to pass the acceptance test and be accepted, so the chain wandered off through undefined values."""
        with pytest.raises(ValueError, match="nan"):
            metropolis_hastings(
                self._log_posterior_nan_outside, [0.0], n_samples=200, burn_in=10, proposal_scale=1.0, random_state=1
            )

    @pytest.mark.parametrize("initial_value", [float("nan"), float("-inf"), float("inf")])
    def test_the_initial_point_needs_a_finite_log_posterior(self, initial_value: float) -> None:
        with pytest.raises(ValueError, match="finite"):
            metropolis_hastings(lambda params: initial_value, [0.0], n_samples=50, burn_in=5, random_state=1)

    def test_zero_density_proposals_are_rejected_not_errors(self) -> None:
        def inside_the_interval(params: list[float]) -> float:
            return -0.5 * params[0] ** 2 if abs(params[0]) <= 1.0 else float("-inf")

        chain = metropolis_hastings(
            inside_the_interval, [0.0], n_samples=3000, burn_in=100, proposal_scale=0.8, random_state=2
        )
        assert max(abs(s[0]) for s in chain["samples"]) <= 1.0
        assert 0.0 < chain["acceptance_rate"] < 1.0

    @pytest.mark.parametrize("scale", [0.0, -1.0, float("nan")])
    def test_proposal_scale_must_be_positive(self, scale: float) -> None:
        with pytest.raises(ValueError, match="proposal_scale"):
            metropolis_hastings(lambda p: -0.5 * p[0] ** 2, [0.0], n_samples=50, burn_in=5, proposal_scale=scale)

    def test_burn_in_must_not_be_negative(self) -> None:
        """A negative burn-in used to keep only the last |burn_in| draws (a slice from the end)."""
        with pytest.raises(ValueError, match="burn_in"):
            metropolis_hastings(lambda p: -0.5 * p[0] ** 2, [0.0], n_samples=50, burn_in=-5)


# ---------------------------------------------------------------------------
# Bare assertions on valid input types
# ---------------------------------------------------------------------------


class TestInputTypesAreNotAssertions:
    def test_ld_frequency_mode_accepts_numpy_scalars_and_ints(self) -> None:
        """``assert isinstance(pA, float)`` rejected np.float32 (and ints) with an AssertionError."""
        d, d_prime_value = ld_coefficients(np.float32(0.5), 0.5, 0.5, 0.5, 0.3)
        assert d == pytest.approx(0.05)
        assert d_prime_value == pytest.approx(0.05 / 0.25)
        assert ld_coefficients(1, 0.0, 1.0, 0.0, 1.0) == (0.0, 0)

    def test_ld_genotype_mode_accepts_tuples_and_arrays(self) -> None:
        rows = [[1, 1], [1, 0], [0, 1], [0, 0], [1, 1], [0, 0]]
        expected = ld_coefficients(rows)
        assert ld_coefficients(tuple(map(tuple, rows))) == expected
        assert ld_coefficients(np.array(rows)) == expected

    def test_ld_unphased_genotypes_as_an_array(self) -> None:
        rows = [[[0, 1], [0, 1]], [[0, 0], [0, 1]], [[1, 1], [0, 0]], [[0, 1], [1, 1]], [[1, 1], [1, 1]]]
        assert ld_coefficients(np.array(rows), phased=False) == ld_coefficients(rows, phased=False)

    @pytest.mark.parametrize("bad", ["0.5", None, {"a": 1}])
    def test_ld_wrong_types_raise_type_error(self, bad: Any) -> None:
        with pytest.raises(TypeError):
            ld_coefficients(bad)
        with pytest.raises(TypeError):
            ld_coefficients(bad, 0.5, 0.5, 0.5, 0.3)

    def test_r_squared_regression_mode_accepts_arrays_and_tuples(self) -> None:
        """``assert isinstance(x, list)`` rejected the arrays that ``linear_regression`` itself accepts."""
        x, y = [1.0, 2.0, 3.0, 4.0], [2.1, 3.9, 6.2, 7.8]
        expected = r_squared(x, y)
        assert 0.99 < expected < 1.0
        assert r_squared(np.array(x), np.array(y)) == pytest.approx(expected)
        assert r_squared(tuple(x), tuple(y)) == pytest.approx(expected)

    def test_r_squared_ld_mode_accepts_numpy_scalars_and_ints(self) -> None:
        # D = 0.3 - 0.25; r^2 = D^2 / (pA pa pB pb)
        assert r_squared(np.float32(0.5), np.float32(0.5), 0.5, 0.5, 0.3) == pytest.approx(0.05**2 / 0.0625)
        assert r_squared(1, 0, 1, 0, 1) == 0.0

    @pytest.mark.parametrize(
        "args", [("a", "b"), (3.0, [1.0, 2.0]), ([1.0, 2.0], 3.0)], ids=["strings", "number-list", "list-number"]
    )
    def test_r_squared_wrong_types_raise_type_error(self, args: tuple) -> None:
        with pytest.raises(TypeError):
            r_squared(*args)

    def test_multilevel_decomposition_requires_its_variances(self) -> None:
        with pytest.raises(ValueError, match="required"):
            multilevel_selection_decomposition(group_size=5)
        with pytest.raises(ValueError, match="required"):
            multilevel_selection_decomposition(group_size=5, individual_trait_variance=1.0)
        result = multilevel_selection_decomposition(
            group_size=4, individual_trait_variance=2.0, group_trait_variance=0.5
        )
        assert result["total_variance"] == pytest.approx(2.0 + 0.5 * 4)

    def test_haldane_rejects_a_negative_distance(self) -> None:
        with pytest.raises(ValueError, match="negative"):
            haldane_d_to_c(-0.1)
        assert haldane_d_to_c(0.0) == 0.0
        assert haldane_d_to_c(0.5) == pytest.approx(0.5 * (1 - math.exp(-1.0)))


# ---------------------------------------------------------------------------
# Placeholder values for impossible input
# ---------------------------------------------------------------------------


class TestHardyWeinberg:
    @pytest.mark.parametrize("frequency", [-0.1, 1.5, float("nan"), float("inf")])
    def test_a_frequency_outside_the_unit_interval_raises(self, frequency: float) -> None:
        """Zeros summing to 0 used to be returned, which is not a genotype distribution."""
        with pytest.raises(ValueError, match="between 0 and 1"):
            hardy_weinberg_genotype_freqs(frequency)

    @pytest.mark.parametrize("p", [0.0, 0.25, 0.6, 1.0])
    def test_valid_frequencies_give_p_squared_2pq_q_squared(self, p: float) -> None:
        aa, het, hom = hardy_weinberg_genotype_freqs(p)
        assert (aa, het) == (pytest.approx(p * p), pytest.approx(2 * p * (1 - p)))
        assert hom == pytest.approx((1 - p) ** 2)
        assert aa + het + hom == pytest.approx(1.0)


class TestReplicatorDynamics:
    def test_derivative_with_zero_mean_fitness_is_x_times_f(self) -> None:
        """dx_i/dt = x_i (f_i - mean f) has no division by the mean; zeros used to be returned at mean fitness 0."""
        assert egt.replicator_derivative([1.0, -1.0], [0.5, 0.5]) == [0.5, -0.5]
        # mean fitness 0.25 * 3 + 0.75 * (-1) = 0
        assert egt.replicator_derivative([3.0, -1.0], [0.25, 0.75]) == [0.75, -0.75]

    def test_derivatives_of_a_distribution_sum_to_zero_and_follow_the_definition(self) -> None:
        rng = np.random.default_rng(3)
        for _ in range(20):
            frequencies = rng.dirichlet(np.ones(5))
            fitnesses = rng.normal(size=5)
            derivatives = egt.replicator_derivative(list(fitnesses), list(frequencies))
            mean = float(np.dot(fitnesses, frequencies))
            assert derivatives == pytest.approx(list(frequencies * (fitnesses - mean)), abs=1e-15)
            assert sum(derivatives) == pytest.approx(0.0, abs=1e-14)

    def test_step_follows_the_euler_update(self) -> None:
        # dx = (0.5, -0.5) at zero mean fitness: x' = (0.5, 0.5) + 0.1 * (0.5, -0.5) = (0.55, 0.45)
        assert egt.replicator_step([1.0, -1.0], [0.5, 0.5], time_step=0.1) == pytest.approx([0.55, 0.45])

    @pytest.mark.parametrize("frequencies", [[0.0, 0.0], [-0.1, 1.1]])
    def test_impossible_frequencies_raise_instead_of_becoming_uniform(self, frequencies: list[float]) -> None:
        with pytest.raises(ValueError):
            egt.replicator_step([1.0, 2.0], frequencies, time_step=0.1)

    def test_a_step_that_zeroes_every_frequency_raises(self) -> None:
        with pytest.raises(ValueError, match="every frequency"):
            egt.replicator_step([-10.0, -10.0], [0.1, 0.1], time_step=1.0)


class TestSignalDetection:
    def test_empty_trial_classes_leave_the_rates_undefined(self) -> None:
        """One trial used to be assumed for an empty class, giving a hit or false-alarm rate of 0 and a made-up d'."""
        with pytest.raises(ValueError, match="No noise trials"):
            sdt_metrics(hits=8, misses=2, false_alarms=0, correct_rejections=0)
        with pytest.raises(ValueError, match="No signal trials"):
            sdt_metrics(hits=0, misses=0, false_alarms=2, correct_rejections=8)

    def test_negative_counts_raise(self) -> None:
        with pytest.raises(ValueError, match="non-negative"):
            sdt_metrics(hits=-1, misses=5, false_alarms=2, correct_rejections=8)

    @pytest.mark.parametrize("rate", [-0.1, 1.5, float("nan")])
    def test_rates_must_be_proportions(self, rate: float) -> None:
        """Rates above 1 or below 0 used to be clipped silently to a finite d'."""
        with pytest.raises(ValueError, match="proportion"):
            d_prime(rate, 0.2)
        with pytest.raises(ValueError, match="proportion"):
            d_prime(0.8, rate)
        with pytest.raises(ValueError, match="proportion"):
            criterion_c(rate, 0.2)

    def test_d_prime_is_the_difference_of_normal_quantiles(self) -> None:
        normal = std.NormalDist()
        for hit, false_alarm in [(0.8, 0.2), (0.95, 0.4), (0.6, 0.55)]:
            assert d_prime(hit, false_alarm) == pytest.approx(
                normal.inv_cdf(hit) - normal.inv_cdf(false_alarm), abs=1e-9
            )
            assert criterion_c(hit, false_alarm) == pytest.approx(
                -(normal.inv_cdf(hit) + normal.inv_cdf(false_alarm)) / 2, abs=1e-9
            )

    def test_exact_zero_and_one_rates_are_still_clipped_to_finite_values(self) -> None:
        assert math.isfinite(d_prime(1.0, 0.0))
        assert math.isfinite(criterion_c(1.0, 0.0))

    def test_metrics_from_counts(self) -> None:
        result = sdt_metrics(hits=40, misses=10, false_alarms=10, correct_rejections=40)
        assert result["hit_rate"] == pytest.approx(0.8)
        assert result["false_alarm_rate"] == pytest.approx(0.2)
        assert result["d_prime"] == pytest.approx(2 * std.NormalDist().inv_cdf(0.8), abs=1e-9)
        assert result["c"] == pytest.approx(0.0, abs=1e-9)


# ---------------------------------------------------------------------------
# Coalescent and population-genetic expectations
# ---------------------------------------------------------------------------


def _harmonic(m: int) -> float:
    return math.fsum(1.0 / i for i in range(1, m + 1))


def _simulated_mean_pairwise_difference(
    n_samples: int, theta: float, replicates: int, seed: int
) -> tuple[float, float]:
    """Mean number of pairwise differences over simulated Kingman genealogies: (mean, standard error).

    Time is in units of 2N generations, where k lineages coalesce at rate k (k - 1) / 2 and a branch of length
    t carries Poisson(theta t / 2) mutations, so two sequences whose MRCA is T back differ at
    Poisson(theta T) sites. The simulation shares no code with the module.
    """
    rng = np.random.default_rng(seed)
    pair_count = n_samples * (n_samples - 1) // 2
    means = np.empty(replicates)
    for r in range(replicates):
        groups = [[i] for i in range(n_samples)]
        elapsed = 0.0
        pair_times = []
        while len(groups) > 1:
            k = len(groups)
            elapsed += rng.exponential(2.0 / (k * (k - 1)))
            i, j = sorted(rng.choice(k, 2, replace=False))
            pair_times.extend([elapsed] * (len(groups[i]) * len(groups[j])))
            groups[i] = groups[i] + groups[j]
            del groups[j]
        assert len(pair_times) == pair_count
        means[r] = rng.poisson(theta * np.array(pair_times)).mean()
    return float(means.mean()), float(means.std(ddof=1) / math.sqrt(replicates))


class TestCoalescentExpectations:
    @pytest.mark.parametrize("n_samples", [2, 5, 12])
    def test_expected_pairwise_diversity_is_theta_whatever_the_sample_size(self, n_samples: int) -> None:
        """Every pair has the same expected number of differences (Tajima 1983); the old value grew with n."""
        theta = 2.0
        simulated, standard_error = _simulated_mean_pairwise_difference(n_samples, theta, 3000, seed=n_samples)
        assert abs(simulated - theta) < 4 * standard_error
        assert statistics.expected_pairwise_diversity(n_samples, theta) == pytest.approx(theta)

    def test_a_single_sequence_has_no_pairs(self) -> None:
        assert statistics.expected_pairwise_diversity(1, 5.0) == 0.0

    def test_expected_segregating_sites_still_grows_with_the_sample(self) -> None:
        assert statistics.expected_segregating_sites(5, 2.0) == pytest.approx(2.0 * _harmonic(4))

    def test_waiting_times_use_the_effective_size(self) -> None:
        """4 Ne / (k (k - 1)) generations for k = 4, 3, 2 lineages; the effective size used to be ignored."""
        assert statistics.expected_coalescent_waiting_times(4, 100.0) == pytest.approx([400 / 12, 400 / 6, 400 / 2])
        assert statistics.expected_coalescent_waiting_times(4, 200.0) == pytest.approx([800 / 12, 800 / 6, 800 / 2])

    def test_waiting_times_add_up_to_the_expected_time_to_the_mrca(self) -> None:
        for n_samples in (2, 7, 30):
            total = sum(statistics.expected_coalescent_waiting_times(n_samples, 250.0))
            assert total == pytest.approx(coalescent.expected_time_to_mrca(n_samples, 250.0))
            assert total == pytest.approx(4 * 250.0 * (1 - 1 / n_samples))

    @pytest.mark.parametrize("n_samples,effective_size", [(1, 100.0), (5, 0.0), (5, -3.0)])
    def test_waiting_times_validate(self, n_samples: int, effective_size: float) -> None:
        with pytest.raises(ValueError):
            statistics.expected_coalescent_waiting_times(n_samples, effective_size)

    @pytest.mark.parametrize("effective_size", [0.0, -1.0])
    def test_time_to_mrca_needs_a_positive_effective_size(self, effective_size: float) -> None:
        with pytest.raises(ValueError, match="positive"):
            coalescent.expected_time_to_mrca(5, effective_size)

    @pytest.mark.parametrize("n_samples", [2, 50, 101, 102, 5000])
    def test_total_branch_length_uses_the_exact_harmonic_number(self, n_samples: int) -> None:
        """Beyond 101 lineages 4 Ne (log(n - 1) + gamma) replaced 4 Ne H_{n-1}, an error of about 1 / (2 n)."""
        assert coalescent.expected_total_branch_length(n_samples, 250.0) == pytest.approx(
            4 * 250.0 * _harmonic(n_samples - 1), rel=1e-14
        )

    @pytest.mark.parametrize("effective_size", [0.0, -5.0])
    def test_total_branch_length_validates_the_effective_size(self, effective_size: float) -> None:
        with pytest.raises(ValueError, match="positive"):
            coalescent.expected_total_branch_length(10, effective_size)

    @pytest.mark.parametrize("n_samples", [2, 10, 100, 101, 1000, 20000])
    def test_watterson_theta_divides_by_the_exact_harmonic_number(self, n_samples: int) -> None:
        assert coalescent.watterson_theta(37, n_samples, sequence_length=1000) == pytest.approx(
            37 / (_harmonic(n_samples - 1) * 1000), rel=1e-14
        )

    def test_watterson_theta_needs_a_positive_sequence_length(self) -> None:
        with pytest.raises(ValueError, match="length"):
            coalescent.watterson_theta(10, 5, sequence_length=0)

    def test_watterson_theta_without_scipy_is_exact_not_approximated(self, tmp_path: Path) -> None:
        """Without scipy the estimator fell back to log(n - 1) + gamma for n > 100 (relative error near 5e-4).

        scipy is blocked after the package is imported, in a fresh interpreter pinned to the very package this
        process imported.
        """
        code = (
            "import sys, math\n"
            "from metainformant.math.population_genetics.coalescent import watterson_theta\n"
            "sys.modules['scipy'] = None\n"
            "sys.modules['scipy.special'] = None\n"
            "exact = 50 / (math.fsum(1.0 / i for i in range(1, 1000)) * 100)\n"
            "print(abs(watterson_theta(50, 1000, sequence_length=100) - exact) / exact)\n"
        )
        package_root = str(Path(metainformant.__file__).resolve().parent.parent)
        env = {**os.environ, "PYTHONDONTWRITEBYTECODE": "1", "PYTHONPATH": package_root}
        proc = subprocess.run(
            [sys.executable, "-c", code], capture_output=True, text=True, cwd=tmp_path, env=env, timeout=120
        )
        assert proc.returncode == 0, proc.stderr
        assert float(proc.stdout) < 1e-14


# ---------------------------------------------------------------------------
# Outliers, docstring promises
# ---------------------------------------------------------------------------


class TestOutliersAndDocstrings:
    def test_identical_values_have_no_z_score_outliers(self) -> None:
        """The z-scores of identical values divided 0 by 0 (a RuntimeWarning, an error in this suite)."""
        result = detect_outliers([2.0, 2.0, 2.0, 2.0], method="zscore", threshold=2.0)
        assert result["outlier_indices"] == [] and result["outliers"] == []
        assert result["total_values"] == 4

    def test_covariance_is_the_population_covariance_and_its_docstring_says_so(self) -> None:
        x, y = [1, 2, 3, 4, 5], [2, 4, 6, 8, 10]
        assert covariance(x, y) == 4.0  # sum of (x - 3)(y - 6) is 20, over n = 5
        assert covariance(x, y) == pytest.approx(np.cov(x, y, ddof=1)[0, 1] * 4 / 5)
        assert _docstring_examples_pass(covariance)

    def test_simulation_docstring_example_runs(self) -> None:
        from metainformant.dna.variation import mutations

        assert _docstring_examples_pass(mutations.simulate_sequence_evolution, random=__import__("random"))
