"""Bayesian inference methods for biological and statistical modelling.

Provides Markov Chain Monte Carlo sampling (Metropolis-Hastings), Approximate
Bayesian Computation (rejection sampling), Bayes factor computation, conjugate
prior analysis (Beta-Binomial, Normal-Normal), and model comparison criteria
(DIC, WAIC).

All algorithms are pure Python implementations with optional NumPy
acceleration for matrix operations.
"""

from __future__ import annotations

import math
import random
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _mean(values: list[float]) -> float:
    """Arithmetic mean."""
    return sum(values) / len(values) if values else 0.0


def _variance(values: list[float], ddof: int = 1) -> float:
    """Sample variance."""
    n = len(values)
    if n < 2:
        return 0.0
    mu = _mean(values)
    return sum((x - mu) ** 2 for x in values) / (n - ddof)


def _std(values: list[float], ddof: int = 1) -> float:
    """Sample standard deviation."""
    return math.sqrt(max(0.0, _variance(values, ddof)))


def _quantile(sorted_data: list[float], q: float) -> float:
    """Compute quantile from sorted data."""
    n = len(sorted_data)
    if n == 0:
        return 0.0
    idx = q * (n - 1)
    low = int(math.floor(idx))
    high = min(int(math.ceil(idx)), n - 1)
    if low == high:
        return sorted_data[low]
    frac = idx - low
    return sorted_data[low] * (1 - frac) + sorted_data[high] * frac


def _euclidean_distance(a: list[float], b: list[float]) -> float:
    """Euclidean distance between two vectors."""
    return math.sqrt(sum((ai - bi) ** 2 for ai, bi in zip(a, b)))


def _effective_sample_size(samples: list[float]) -> float:
    """Estimate effective sample size from autocorrelation.

    Uses the initial monotone sequence estimator (truncated at first
    negative autocorrelation).
    """
    n = len(samples)
    if n < 10:
        return float(n)

    mu = _mean(samples)
    var = _variance(samples, ddof=0)
    if var == 0:
        return float(n)

    # Compute autocorrelations until they become negative
    total_acf = 0.0
    for lag in range(1, n // 2):
        acf = sum((samples[i] - mu) * (samples[i + lag] - mu) for i in range(n - lag)) / ((n - lag) * var)
        if acf < 0:
            break
        total_acf += acf

    ess = n / (1.0 + 2.0 * total_acf)
    return max(1.0, ess)


# ---------------------------------------------------------------------------
# Metropolis-Hastings MCMC
# ---------------------------------------------------------------------------


def _resolve_rng(random_state: int | random.Random | None) -> Any:
    """Random source for the samplers in this module.

    ``None`` keeps the long-standing behaviour of drawing from (and advancing) Python's global
    generator, so ``random.seed`` still makes a run reproducible. An ``int`` seeds a private
    ``random.Random`` and a ``random.Random`` is used as given; in both cases the global generator
    is neither read nor advanced, so the run is reproducible whatever the rest of the program does.
    """
    if random_state is None:
        return random
    if isinstance(random_state, random.Random):
        return random_state
    return random.Random(random_state)


def metropolis_hastings(
    log_posterior: Any,
    initial: list[float],
    n_samples: int = 10000,
    proposal_scale: float = 0.1,
    burn_in: int = 1000,
    *,
    random_state: int | random.Random | None = None,
) -> dict:
    """MCMC sampling using the Metropolis-Hastings algorithm.

    Draws samples from an unnormalised posterior distribution using a
    symmetric Gaussian random walk proposal.

    Args:
        log_posterior: Callable that takes a list of floats (parameter vector)
            and returns the log-posterior density (up to a normalising
            constant). Must accept ``list[float]`` and return ``float``.
        initial: Starting parameter values.
        n_samples: Total number of samples to draw (including burn-in).
        proposal_scale: Standard deviation of the Gaussian proposal
            distribution for each parameter.
        burn_in: Number of initial samples to discard as burn-in.
        random_state: Seed (an ``int``) or a ``random.Random`` for the proposal
            and accept/reject draws. The same seed always gives the same chain,
            independent of Python's global generator. ``None`` (default) draws
            from the global generator instead, which ``random.seed`` controls.

    Returns:
        Dictionary with keys:
            - samples: 2D list of shape ``(n_retained, n_params)`` after
              burn-in.
            - acceptance_rate: Fraction of proposals accepted.
            - log_posteriors: List of log-posterior values for retained samples.
            - effective_sample_size: Estimated ESS for each parameter.
            - n_params: Number of parameters.

    Raises:
        ValueError: If burn_in >= n_samples, burn_in is negative, proposal_scale is not positive, or
            ``log_posterior`` is not finite at ``initial`` or returns NaN or +inf at a proposal (-inf is
            a legitimate zero density and the proposal is rejected). A NaN used to be accepted as if it
            were an improvement, silently turning the chain into an unguided random walk.
    """
    if burn_in >= n_samples:
        raise ValueError(f"burn_in ({burn_in}) must be less than n_samples ({n_samples})")
    if burn_in < 0:
        raise ValueError(f"burn_in must be non-negative, got {burn_in}")
    if not proposal_scale > 0:
        raise ValueError(f"proposal_scale must be positive, got {proposal_scale}")

    rng = _resolve_rng(random_state)

    n_params = len(initial)
    current = list(initial)
    current_lp = float(log_posterior(current))
    if not math.isfinite(current_lp):
        raise ValueError(f"log_posterior must be finite at the initial point, got {current_lp} at {current}")

    all_samples: list[list[float]] = []
    all_lp: list[float] = []
    n_accepted = 0

    for step in range(n_samples):
        # Propose new state
        proposal = [current[i] + rng.gauss(0, proposal_scale) for i in range(n_params)]

        proposal_lp = float(log_posterior(proposal))
        if math.isnan(proposal_lp) or proposal_lp == math.inf:
            raise ValueError(f"log_posterior returned {proposal_lp} at the proposal {proposal}")

        # Accept/reject
        log_alpha = proposal_lp - current_lp
        if log_alpha >= 0 or rng.random() < math.exp(min(0.0, log_alpha)):
            current = proposal
            current_lp = proposal_lp
            n_accepted += 1

        all_samples.append(list(current))
        all_lp.append(current_lp)

    # Discard burn-in
    retained = all_samples[burn_in:]
    retained_lp = all_lp[burn_in:]

    # Compute ESS per parameter
    ess = []
    for p in range(n_params):
        param_chain = [s[p] for s in retained]
        ess.append(_effective_sample_size(param_chain))

    acceptance_rate = n_accepted / n_samples

    logger.info(
        "MH sampling: %d samples (burn-in %d), acceptance=%.3f, ESS=%s",
        len(retained),
        burn_in,
        acceptance_rate,
        [f"{e:.0f}" for e in ess],
    )

    return {
        "samples": retained,
        "acceptance_rate": acceptance_rate,
        "log_posteriors": retained_lp,
        "effective_sample_size": ess,
        "n_params": n_params,
    }


# ---------------------------------------------------------------------------
# Approximate Bayesian Computation
# ---------------------------------------------------------------------------


def abc_rejection(
    simulator: Any,
    observed_summary: list[float],
    prior_sampler: Any,
    n_simulations: int = 10000,
    tolerance: float = 0.01,
    distance: str = "euclidean",
) -> dict:
    """Approximate Bayesian Computation using rejection sampling.

    Samples parameters from the prior, simulates data, computes summary
    statistics, and accepts parameters whose summaries are within
    ``tolerance`` of the observed summaries.

    This function draws no random numbers itself: ``prior_sampler`` and
    ``simulator`` own all of the randomness, so seed those (for example with
    a private ``random.Random`` or ``numpy`` ``Generator``) to make a run
    reproducible.

    Args:
        simulator: Callable that takes a parameter vector and returns a list
            of summary statistics. ``simulator(params) -> list[float]``.
        observed_summary: List of observed summary statistics.
        prior_sampler: Callable that returns a random parameter draw from the
            prior. ``prior_sampler() -> list[float]``.
        n_simulations: Number of simulation rounds.
        tolerance: Distance threshold for acceptance.
        distance: Distance metric: ``"euclidean"`` (default) or ``"manhattan"``.

    Returns:
        Dictionary with keys:
            - accepted_params: List of accepted parameter vectors.
            - acceptance_rate: Fraction of simulations accepted.
            - posterior_summary: Dict with mean, std, and quantiles per param.
            - distances: List of distances for accepted simulations.
    """
    if distance == "manhattan":

        def dist_fn(a: list[float], b: list[float]) -> float:
            return sum(abs(ai - bi) for ai, bi in zip(a, b))

    else:
        dist_fn = _euclidean_distance

    accepted_params: list[list[float]] = []
    accepted_dists: list[float] = []

    for sim_idx in range(n_simulations):
        params = prior_sampler()
        sim_summary = simulator(params)
        d = dist_fn(sim_summary, observed_summary)

        if d <= tolerance:
            accepted_params.append(list(params))
            accepted_dists.append(d)

    acceptance_rate = len(accepted_params) / n_simulations if n_simulations > 0 else 0.0

    # Posterior summary
    posterior_summary: dict[str, dict[str, float]] = {}
    if accepted_params:
        n_params = len(accepted_params[0])
        for p in range(n_params):
            vals = sorted(ap[p] for ap in accepted_params)
            posterior_summary[f"param_{p}"] = {
                "mean": _mean(vals),
                "std": _std(vals) if len(vals) > 1 else 0.0,
                "median": _quantile(vals, 0.5),
                "q025": _quantile(vals, 0.025),
                "q975": _quantile(vals, 0.975),
            }

    logger.info(
        "ABC rejection: %d/%d accepted (rate=%.4f, tol=%.4f)",
        len(accepted_params),
        n_simulations,
        acceptance_rate,
        tolerance,
    )

    return {
        "accepted_params": accepted_params,
        "acceptance_rate": acceptance_rate,
        "posterior_summary": posterior_summary,
        "distances": accepted_dists,
    }


# ---------------------------------------------------------------------------
# Bayes Factor
# ---------------------------------------------------------------------------


def compute_bayes_factor(
    log_marginal_1: float,
    log_marginal_2: float,
) -> dict:
    """Compute Bayes factor and provide interpretation.

    The Bayes factor BF_12 = p(data|M1) / p(data|M2). Here we compute
    it from log marginal likelihoods.

    Args:
        log_marginal_1: Log marginal likelihood of model 1.
        log_marginal_2: Log marginal likelihood of model 2.

    Returns:
        Dictionary with keys:
            - bf: Bayes factor (M1 vs M2).
            - log_bf: Log Bayes factor.
            - interpretation: Jeffreys' scale interpretation.
            - favored_model: 1 or 2.
    """
    log_bf = log_marginal_1 - log_marginal_2

    # Prevent overflow
    if log_bf > 700:
        bf = float("inf")
    elif log_bf < -700:
        bf = 0.0
    else:
        bf = math.exp(log_bf)

    # Jeffreys' scale interpretation
    abs_log_bf = abs(log_bf)
    if abs_log_bf < math.log(3):
        interpretation = "Anecdotal evidence"
    elif abs_log_bf < math.log(10):
        interpretation = "Moderate evidence"
    elif abs_log_bf < math.log(30):
        interpretation = "Strong evidence"
    elif abs_log_bf < math.log(100):
        interpretation = "Very strong evidence"
    else:
        interpretation = "Decisive evidence"

    favored = 1 if log_bf >= 0 else 2
    interpretation += f" favoring Model {favored}"

    logger.info(
        "Bayes factor: BF=%.4g, log(BF)=%.4f, %s",
        bf,
        log_bf,
        interpretation,
    )

    return {
        "bf": bf,
        "log_bf": log_bf,
        "interpretation": interpretation,
        "favored_model": favored,
    }


# ---------------------------------------------------------------------------
# Conjugate analysis: Beta-Binomial
# ---------------------------------------------------------------------------


def conjugate_beta_binomial(
    successes: int,
    trials: int,
    prior_alpha: float = 1.0,
    prior_beta: float = 1.0,
) -> dict:
    """Beta-Binomial conjugate posterior update.

    Given a Beta(alpha, beta) prior and Binomial(n, k) data, the posterior
    is Beta(alpha + k, beta + n - k).

    Args:
        successes: Number of successes observed.
        trials: Total number of trials.
        prior_alpha: Alpha parameter of the Beta prior (default 1 = uniform).
        prior_beta: Beta parameter of the Beta prior (default 1 = uniform).

    Returns:
        Dictionary with keys:
            - posterior_alpha: Updated alpha.
            - posterior_beta: Updated beta.
            - posterior_mean: Posterior mean.
            - credible_interval_95: Tuple of (lower, upper), the equal-tailed 95% credible interval
              from the exact quantiles (2.5% and 97.5%) of the Beta posterior.
            - prior_mean: Prior mean.
            - map_estimate: Maximum a posteriori estimate (the posterior mean when the posterior has no
              interior mode, that is when posterior_alpha or posterior_beta is at most 1).

    Raises:
        ValueError: If successes > trials, a count is negative, or a prior parameter is not positive.
    """
    if successes < 0 or trials < 0:
        raise ValueError("successes and trials must be non-negative")
    if successes > trials:
        raise ValueError(f"successes ({successes}) > trials ({trials})")
    if not (prior_alpha > 0 and prior_beta > 0):
        raise ValueError(f"prior_alpha and prior_beta must be positive, got {prior_alpha} and {prior_beta}")

    post_alpha = prior_alpha + successes
    post_beta = prior_beta + (trials - successes)

    prior_mean = prior_alpha / (prior_alpha + prior_beta)
    post_mean = post_alpha / (post_alpha + post_beta)

    # MAP estimate
    if post_alpha > 1 and post_beta > 1:
        map_est = (post_alpha - 1) / (post_alpha + post_beta - 2)
    else:
        map_est = post_mean

    # Equal-tailed 95% credible interval from the exact Beta quantiles (a normal approximation is poor for
    # small counts and for posteriors piled against 0 or 1)
    from scipy.stats import beta as beta_distribution

    ci_lower = float(beta_distribution.ppf(0.025, post_alpha, post_beta))
    ci_upper = float(beta_distribution.ppf(0.975, post_alpha, post_beta))

    logger.info(
        "Beta-Binomial update: %d/%d successes, posterior=Beta(%.2f, %.2f), " "mean=%.4f, CI=[%.4f, %.4f]",
        successes,
        trials,
        post_alpha,
        post_beta,
        post_mean,
        ci_lower,
        ci_upper,
    )

    return {
        "posterior_alpha": post_alpha,
        "posterior_beta": post_beta,
        "posterior_mean": post_mean,
        "credible_interval_95": (ci_lower, ci_upper),
        "prior_mean": prior_mean,
        "map_estimate": map_est,
    }


# ---------------------------------------------------------------------------
# Conjugate analysis: Normal-Normal
# ---------------------------------------------------------------------------


def conjugate_normal(
    data: list[float],
    prior_mean: float = 0.0,
    prior_var: float = 100.0,
    known_var: float | None = None,
) -> dict:
    """Normal-Normal conjugate posterior update.

    Given a Normal prior on the mean mu ~ N(mu_0, tau_0^2) and Normal
    likelihood y_i ~ N(mu, sigma^2), the posterior is Normal with updated
    mean and variance.

    Args:
        data: Observed data values.
        prior_mean: Prior mean (mu_0).
        prior_var: Prior variance (tau_0^2).
        known_var: Known data variance (sigma^2), positive. If ``None``, the
            sample variance of the data is used as a plug-in value.

    Returns:
        Dictionary with keys:
            - posterior_mean: Updated posterior mean.
            - posterior_var: Updated posterior variance.
            - credible_interval_95: Tuple of (lower, upper) 95% CI.
            - n_observations: Number of data points.
            - data_mean: Sample mean.
            - data_var: Sample (or known) variance.

    Raises:
        ValueError: If data is empty, prior_var <= 0, or known_var is not positive. Without
            ``known_var`` the data variance must be estimable: a single observation has no sample
            variance and identical observations give zero, and no stand-in value (1.0 or 1e-10 used
            to be substituted) is invented for them.
    """
    if not data:
        raise ValueError("Data must contain at least one observation")
    if prior_var <= 0:
        raise ValueError(f"prior_var must be positive, got {prior_var}")

    n = len(data)
    data_mean = sum(data) / n

    if known_var is None:
        if n < 2:
            raise ValueError("A single observation has no sample variance: pass known_var")
        data_var = sum((x - data_mean) ** 2 for x in data) / (n - 1)
        if data_var <= 0:
            raise ValueError("The observations are all equal, so the sample variance is zero: pass known_var")
    else:
        if not known_var > 0:
            raise ValueError(f"known_var must be positive, got {known_var}")
        data_var = known_var

    # Posterior precision = prior precision + data precision
    prior_prec = 1.0 / prior_var
    data_prec = n / data_var
    post_prec = prior_prec + data_prec
    post_var = 1.0 / post_prec

    # Posterior mean = weighted average
    post_mean = (prior_prec * prior_mean + data_prec * data_mean) / post_prec

    # 95% credible interval
    post_sd = math.sqrt(post_var)
    ci_lower = post_mean - 1.96 * post_sd
    ci_upper = post_mean + 1.96 * post_sd

    logger.info(
        "Normal conjugate update: n=%d, prior=(%.2f, %.2f), " "posterior=(%.4f, %.6f), CI=[%.4f, %.4f]",
        n,
        prior_mean,
        prior_var,
        post_mean,
        post_var,
        ci_lower,
        ci_upper,
    )

    return {
        "posterior_mean": post_mean,
        "posterior_var": post_var,
        "credible_interval_95": (ci_lower, ci_upper),
        "n_observations": n,
        "data_mean": data_mean,
        "data_var": data_var,
    }


# ---------------------------------------------------------------------------
# Deviance Information Criterion (DIC)
# ---------------------------------------------------------------------------


def compute_dic(
    log_likelihoods: list[float],
    parameters: list[list[float]],
    log_likelihood_at_mean: float | None = None,
) -> dict:
    """Compute the Deviance Information Criterion.

    DIC = D_bar + p_D, where D_bar is the mean deviance of the samples and p_D is the effective number
    of parameters. The two estimators of p_D need different inputs:

    - With ``log_likelihood_at_mean`` (the log-likelihood evaluated at the posterior mean of the
      parameters, which only the caller's model can compute), p_D = D_bar - D(theta_bar)
      (Spiegelhalter et al. 2002).
    - Without it, p_D = var(D) / 2, with var the population variance of the sampled deviances
      (Gelman et al. 2004), which needs no further model evaluation. The deviance at the mean is then
      unavailable, and is reported as ``None`` rather than as the mean deviance.

    Args:
        log_likelihoods: List of log-likelihood values at each MCMC sample.
        parameters: List of parameter vectors at each MCMC sample; only checked to hold one vector per
            sample (the likelihood function is needed to use them).
        log_likelihood_at_mean: Log-likelihood at the posterior mean parameters, or None.

    Returns:
        Dictionary with keys:
            - dic: DIC value.
            - effective_parameters: Estimated effective number of parameters (p_D).
            - mean_deviance: Mean deviance D_bar.
            - deviance_at_mean: -2 * ``log_likelihood_at_mean``, or None when it was not supplied.

    Raises:
        ValueError: If inputs have different lengths or are empty.
    """
    if len(log_likelihoods) != len(parameters):
        raise ValueError(
            f"log_likelihoods ({len(log_likelihoods)}) and parameters " f"({len(parameters)}) must have same length"
        )
    if not log_likelihoods:
        raise ValueError("Need at least one sample")

    # Deviance = -2 * log_likelihood
    deviances = [-2.0 * ll for ll in log_likelihoods]
    d_bar = _mean(deviances)

    d_at_mean: float | None
    if log_likelihood_at_mean is not None:
        # p_D = D_bar - D(theta_bar)
        d_at_mean = -2.0 * log_likelihood_at_mean
        p_d = d_bar - d_at_mean
    else:
        # Without the likelihood function the deviance at theta_bar is unknown: p_D = 0.5 * var(deviance)
        d_at_mean = None
        p_d = 0.5 * _variance(deviances, ddof=0) if len(deviances) > 1 else 0.0

    dic = d_bar + p_d

    logger.info(
        "DIC = %.2f (mean_deviance=%.2f, p_D=%.2f)",
        dic,
        d_bar,
        p_d,
    )

    return {
        "dic": dic,
        "effective_parameters": p_d,
        "mean_deviance": d_bar,
        "deviance_at_mean": d_at_mean,
    }


# ---------------------------------------------------------------------------
# Widely Applicable Information Criterion (WAIC)
# ---------------------------------------------------------------------------


def compute_waic(
    pointwise_log_lik: list[list[float]],
) -> dict:
    """Compute the Widely Applicable Information Criterion.

    WAIC = -2 * (lppd - p_waic), where lppd is the log pointwise
    predictive density and p_waic is the effective number of parameters.

    Args:
        pointwise_log_lik: 2D list of shape ``(n_samples, n_observations)``.
            Element ``[s][i]`` is the log-likelihood of observation *i*
            evaluated at posterior sample *s*.

    Returns:
        Dictionary with keys:
            - waic: WAIC value.
            - p_waic: Effective number of parameters.
            - lppd: Log pointwise predictive density.
            - pointwise_waic: Per-observation WAIC contributions.

    Raises:
        ValueError: If input is empty or ragged.
    """
    if not pointwise_log_lik:
        raise ValueError("Need at least one MCMC sample")

    n_samples = len(pointwise_log_lik)
    n_obs = len(pointwise_log_lik[0])

    for s in range(n_samples):
        if len(pointwise_log_lik[s]) != n_obs:
            raise ValueError(f"Sample {s} has {len(pointwise_log_lik[s])} observations, " f"expected {n_obs}")

    # lppd = sum_i log(mean_s exp(log_lik[s][i]))
    # p_waic = sum_i var_s(log_lik[s][i])
    lppd = 0.0
    p_waic = 0.0
    pointwise_waic: list[float] = []

    for i in range(n_obs):
        ll_i = [pointwise_log_lik[s][i] for s in range(n_samples)]

        # Log-sum-exp for numerical stability
        max_ll = max(ll_i)
        log_mean_exp = max_ll + math.log(sum(math.exp(ll - max_ll) for ll in ll_i) / n_samples)
        lppd_i = log_mean_exp

        # Variance of log-likelihoods
        var_ll = _variance(ll_i, ddof=1) if n_samples > 1 else 0.0

        lppd += lppd_i
        p_waic += var_ll
        pointwise_waic.append(-2.0 * (lppd_i - var_ll))

    waic = -2.0 * (lppd - p_waic)

    logger.info(
        "WAIC = %.2f (lppd=%.2f, p_waic=%.2f, n_obs=%d)",
        waic,
        lppd,
        p_waic,
        n_obs,
    )

    return {
        "waic": waic,
        "p_waic": p_waic,
        "lppd": lppd,
        "pointwise_waic": pointwise_waic,
    }
