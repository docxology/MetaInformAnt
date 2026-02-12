"""Tests for Bayesian inference: MH MCMC, ABC, Bayes factor, conjugate priors, DIC, WAIC.

NO MOCKING POLICY: All tests use real implementations.
"""
from __future__ import annotations

import math

import pytest

from metainformant.math.bayesian.inference import (
    abc_rejection,
    compute_bayes_factor,
    compute_dic,
    compute_waic,
    conjugate_beta_binomial,
    conjugate_normal,
    metropolis_hastings,
)


# ---------------------------------------------------------------------------
# metropolis_hastings
# ---------------------------------------------------------------------------


class TestMetropolisHastings:
    def test_basic_sampling(self):
        def log_post(params):
            return -0.5 * params[0] ** 2  # Standard normal

        result = metropolis_hastings(log_post, [0.0], n_samples=2000, burn_in=500, proposal_scale=1.0)
        assert "samples" in result
        assert len(result["samples"]) == 1500
        assert result["n_params"] == 1

    def test_acceptance_rate_reasonable(self):
        def log_post(params):
            return -0.5 * params[0] ** 2

        result = metropolis_hastings(log_post, [0.0], n_samples=5000, burn_in=1000, proposal_scale=0.5)
        assert 0.05 < result["acceptance_rate"] < 0.99

    def test_burn_in_exceeds_samples_raises(self):
        def log_post(params):
            return 0.0

        with pytest.raises(ValueError, match="burn_in"):
            metropolis_hastings(log_post, [0.0], n_samples=100, burn_in=200)

    def test_multivariate_sampling(self):
        def log_post(params):
            return -0.5 * (params[0] ** 2 + params[1] ** 2)

        result = metropolis_hastings(log_post, [0.0, 0.0], n_samples=3000, burn_in=500, proposal_scale=0.5)
        assert result["n_params"] == 2
        assert len(result["samples"][0]) == 2

    def test_effective_sample_size_positive(self):
        def log_post(params):
            return -0.5 * params[0] ** 2

        result = metropolis_hastings(log_post, [0.0], n_samples=3000, burn_in=500)
        assert all(e > 0 for e in result["effective_sample_size"])

    def test_log_posteriors_returned(self):
        def log_post(params):
            return -0.5 * params[0] ** 2

        result = metropolis_hastings(log_post, [0.0], n_samples=500, burn_in=100)
        assert len(result["log_posteriors"]) == 400
        assert all(isinstance(lp, float) for lp in result["log_posteriors"])


# ---------------------------------------------------------------------------
# abc_rejection
# ---------------------------------------------------------------------------


class TestABCRejection:
    def test_basic_rejection(self):
        import random

        random.seed(42)

        def simulator(params):
            return [params[0] + random.gauss(0, 0.01)]

        def prior_sampler():
            return [random.uniform(-1, 1)]

        result = abc_rejection(simulator, [0.5], prior_sampler, n_simulations=5000, tolerance=0.1)
        assert "accepted_params" in result
        assert "acceptance_rate" in result
        assert result["acceptance_rate"] >= 0.0

    def test_manhattan_distance(self):
        import random

        random.seed(42)

        def simulator(params):
            return [params[0]]

        def prior_sampler():
            return [random.uniform(0, 1)]

        result = abc_rejection(simulator, [0.5], prior_sampler, n_simulations=2000, tolerance=0.2, distance="manhattan")
        assert isinstance(result["distances"], list)

    def test_posterior_summary_keys(self):
        import random

        random.seed(42)

        def simulator(params):
            return [params[0] + random.gauss(0, 0.01)]

        def prior_sampler():
            return [random.uniform(0, 1)]

        result = abc_rejection(simulator, [0.5], prior_sampler, n_simulations=5000, tolerance=0.2)
        if result["accepted_params"]:
            assert "param_0" in result["posterior_summary"]
            summary = result["posterior_summary"]["param_0"]
            assert "mean" in summary
            assert "std" in summary
            assert "median" in summary


# ---------------------------------------------------------------------------
# compute_bayes_factor
# ---------------------------------------------------------------------------


class TestComputeBayesFactor:
    def test_equal_models(self):
        result = compute_bayes_factor(0.0, 0.0)
        assert result["bf"] == pytest.approx(1.0)
        assert result["log_bf"] == pytest.approx(0.0)
        assert result["favored_model"] == 1

    def test_model1_favored(self):
        result = compute_bayes_factor(-10.0, -20.0)
        assert result["log_bf"] == pytest.approx(10.0)
        assert result["favored_model"] == 1

    def test_model2_favored(self):
        result = compute_bayes_factor(-30.0, -10.0)
        assert result["log_bf"] < 0
        assert result["favored_model"] == 2

    def test_extreme_overflow_protection(self):
        result = compute_bayes_factor(1000.0, 0.0)
        assert result["bf"] == float("inf")

    def test_interpretation_present(self):
        result = compute_bayes_factor(-5.0, -15.0)
        assert "interpretation" in result
        assert "favoring Model" in result["interpretation"]


# ---------------------------------------------------------------------------
# conjugate_beta_binomial
# ---------------------------------------------------------------------------


class TestConjugateBetaBinomial:
    def test_uniform_prior(self):
        result = conjugate_beta_binomial(successes=7, trials=10)
        assert result["posterior_alpha"] == 8.0
        assert result["posterior_beta"] == 4.0
        assert 0.0 < result["posterior_mean"] < 1.0

    def test_informative_prior(self):
        result = conjugate_beta_binomial(successes=5, trials=10, prior_alpha=10.0, prior_beta=10.0)
        assert result["posterior_alpha"] == 15.0
        assert result["posterior_beta"] == 15.0
        assert result["posterior_mean"] == pytest.approx(0.5)

    def test_credible_interval_within_bounds(self):
        result = conjugate_beta_binomial(successes=50, trials=100)
        ci = result["credible_interval_95"]
        assert 0.0 <= ci[0] < ci[1] <= 1.0

    def test_map_estimate(self):
        result = conjugate_beta_binomial(successes=8, trials=10, prior_alpha=2.0, prior_beta=2.0)
        assert 0.0 < result["map_estimate"] < 1.0

    def test_invalid_successes_raises(self):
        with pytest.raises(ValueError):
            conjugate_beta_binomial(successes=11, trials=10)

    def test_negative_values_raise(self):
        with pytest.raises(ValueError):
            conjugate_beta_binomial(successes=-1, trials=10)


# ---------------------------------------------------------------------------
# conjugate_normal
# ---------------------------------------------------------------------------


class TestConjugateNormal:
    def test_basic_update(self):
        data = [5.0, 6.0, 7.0, 5.5, 6.5]
        result = conjugate_normal(data)
        assert "posterior_mean" in result
        assert "posterior_var" in result
        assert result["n_observations"] == 5

    def test_strong_prior(self):
        data = [10.0, 11.0, 12.0]
        result = conjugate_normal(data, prior_mean=0.0, prior_var=0.001)
        # Strong prior should keep posterior mean close to 0
        assert abs(result["posterior_mean"]) < 5.0

    def test_known_variance(self):
        data = [5.0, 6.0, 7.0]
        result = conjugate_normal(data, known_var=1.0)
        assert result["data_var"] == 1.0

    def test_empty_data_raises(self):
        with pytest.raises(ValueError, match="at least one"):
            conjugate_normal([])

    def test_negative_prior_var_raises(self):
        with pytest.raises(ValueError, match="positive"):
            conjugate_normal([1.0], prior_var=-1.0)

    def test_credible_interval(self):
        data = [10.0] * 50
        result = conjugate_normal(data, prior_mean=10.0, prior_var=100.0)
        ci = result["credible_interval_95"]
        assert ci[0] <= result["posterior_mean"] <= ci[1]


# ---------------------------------------------------------------------------
# compute_dic
# ---------------------------------------------------------------------------


class TestComputeDIC:
    def test_basic_dic(self):
        log_likelihoods = [-10.0, -11.0, -9.5, -10.5, -10.2]
        parameters = [[0.5], [0.6], [0.4], [0.55], [0.52]]
        result = compute_dic(log_likelihoods, parameters)
        assert "dic" in result
        assert "effective_parameters" in result
        assert "mean_deviance" in result

    def test_mismatched_lengths_raises(self):
        with pytest.raises(ValueError):
            compute_dic([-10.0, -11.0], [[0.5]])

    def test_empty_input_raises(self):
        with pytest.raises(ValueError, match="at least one"):
            compute_dic([], [])


# ---------------------------------------------------------------------------
# compute_waic
# ---------------------------------------------------------------------------


class TestComputeWAIC:
    def test_basic_waic(self):
        # 3 samples, 4 observations
        pointwise = [
            [-1.0, -2.0, -1.5, -1.2],
            [-1.1, -2.1, -1.4, -1.3],
            [-0.9, -1.9, -1.6, -1.1],
        ]
        result = compute_waic(pointwise)
        assert "waic" in result
        assert "p_waic" in result
        assert "lppd" in result
        assert len(result["pointwise_waic"]) == 4

    def test_empty_raises(self):
        with pytest.raises(ValueError):
            compute_waic([])

    def test_ragged_raises(self):
        with pytest.raises(ValueError, match="observations"):
            compute_waic([[-1.0, -2.0], [-1.0]])
