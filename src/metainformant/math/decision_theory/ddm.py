"""Drift-diffusion models for decision making and cognitive processes.

This module provides mathematical models for drift-diffusion processes,
commonly used in cognitive psychology and neuroscience to model decision
making under uncertainty.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import scipy.stats as stats
from scipy.special import expit

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def _validate_ddm_parameters(drift_rate: float, boundary: float, noise_sd: float) -> None:
    """Reject parameters outside the domain of the closed-form DDM expressions."""
    if not (math.isfinite(boundary) and boundary > 0):
        raise ValueError("Boundary must be positive and finite")
    if not (math.isfinite(noise_sd) and noise_sd > 0):
        raise ValueError("Noise standard deviation must be positive and finite")
    if not math.isfinite(drift_rate):
        raise ValueError("Drift rate must be finite")


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


def ddm_log_likelihood(
    rt: float,
    choice: int,
    drift_rate: float,
    boundary_separation: float,
    starting_point: float = 0.5,
    non_decision_time: float = 0.0,
) -> float:
    """Calculate log-likelihood for a single DDM trial.

    Args:
        rt: Response time
        choice: Choice (1 or 2, where 1 is upper boundary)
        drift_rate: Drift rate (v)
        boundary_separation: Distance between boundaries (a)
        starting_point: Starting point as fraction of boundary separation (z/a)
        non_decision_time: Non-decision time component (t0)

    Returns:
        Log-likelihood of the trial
    """
    if rt <= non_decision_time:
        return float("-inf")  # Impossible trial

    if choice not in [1, 2]:
        raise ValueError("Choice must be 1 or 2")

    decision_time = rt - non_decision_time
    z = starting_point * boundary_separation

    # For choice 1 (upper boundary), drift_rate is positive
    # For choice 2 (lower boundary), we flip the sign
    if choice == 2:
        drift_rate = -drift_rate
        z = boundary_separation - z

    # Approximate log-likelihood using normal approximation
    # This is a simplified version - full DDM likelihood requires numerical integration

    # Mean time to boundary
    if drift_rate > 0:
        mean_time = (boundary_separation - z) / drift_rate
        # Variance approximation
        var_time = (boundary_separation - z) * boundary_separation / (3 * (drift_rate**2))
    else:
        # If drift_rate is negative, it's impossible to reach upper boundary
        return float("-inf")

    if var_time <= 0:
        return float("-inf")

    # Log-likelihood assuming normal distribution
    try:
        log_lik = stats.norm.logpdf(decision_time, loc=mean_time, scale=math.sqrt(var_time))
        return float(log_lik)
    except (ValueError, ZeroDivisionError):
        return float("-inf")


def fit_ddm_parameters(
    rt_data: list[float], choice_data: list[int], bounds: dict[str, tuple[float, float]] | None = None
) -> dict[str, float]:
    """Fit DDM parameters to data using maximum likelihood.

    Args:
        rt_data: List of response times
        choice_data: List of choices (1 or 2)
        bounds: Parameter bounds for fitting

    Returns:
        Fitted parameters
    """
    if len(rt_data) != len(choice_data):
        raise ValueError("rt_data and choice_data must have same length")

    if not rt_data:
        raise ValueError("No data provided")

    # Default bounds
    if bounds is None:
        bounds = {
            "drift_rate": (0.01, 2.0),
            "boundary_separation": (0.5, 3.0),
            "starting_point": (0.2, 0.8),
            "non_decision_time": (0.0, 0.5),
        }

    # Try to use scipy.optimize for better fitting
    try:
        from scipy.optimize import minimize

        # Negative log-likelihood function for minimization
        def neg_log_lik(params: np.ndarray) -> float:
            v, a, z, t0 = params
            # Pentalize invalid parameters
            if a <= 0 or not (0 < z < 1) or t0 < 0:
                return 1e10

            total_log_lik = 0.0
            for rt, choice in zip(rt_data, choice_data):
                # Check for impossible RTs (t0 > rt) to avoid infinite penalties blowing up optimization
                if rt <= t0:
                    return 1e10

                ll = ddm_log_likelihood(rt, choice, v, a, z, t0)
                if ll == float("-inf"):
                    return 1e10
                total_log_lik += ll
            return -total_log_lik

        # Initial guess (middle of bounds)
        initial_guess = [
            np.mean(bounds["drift_rate"]),
            np.mean(bounds["boundary_separation"]),
            np.mean(bounds["starting_point"]),
            np.mean(bounds["non_decision_time"]),
        ]

        # Bounds for optimization
        # Note: L-BFGS-B handles bounds well
        opt_bounds = [
            bounds["drift_rate"],
            bounds["boundary_separation"],
            bounds["starting_point"],
            bounds["non_decision_time"],
        ]

        result = minimize(neg_log_lik, initial_guess, method="L-BFGS-B", bounds=opt_bounds)

        if result.success:
            best_params = {
                "drift_rate": result.x[0],
                "boundary_separation": result.x[1],
                "starting_point": result.x[2],
                "non_decision_time": result.x[3],
                "log_likelihood": -result.fun,
            }
            return best_params

    except ImportError:
        logger.warning("scipy.optimize not available, falling back to grid search")

    # Simple grid search (Fallback)
    best_params_grid: dict[str, Any] | None = None
    best_likelihood = float("-inf")

    # Coarse grid search
    drift_rates = np.linspace(bounds["drift_rate"][0], bounds["drift_rate"][1], 10)
    boundary_separations = np.linspace(bounds["boundary_separation"][0], bounds["boundary_separation"][1], 10)
    starting_points = np.linspace(bounds["starting_point"][0], bounds["starting_point"][1], 5)
    non_decision_times = np.linspace(bounds["non_decision_time"][0], bounds["non_decision_time"][1], 5)

    for v in drift_rates:
        for a in boundary_separations:
            for z in starting_points:
                for t0 in non_decision_times:
                    total_log_lik = 0.0
                    for rt, choice in zip(rt_data, choice_data):
                        log_lik = ddm_log_likelihood(rt, choice, v, a, z, t0)
                        total_log_lik += log_lik

                    if total_log_lik > best_likelihood:
                        best_likelihood = total_log_lik
                        best_params_grid = {
                            "drift_rate": v,
                            "boundary_separation": a,
                            "starting_point": z,
                            "non_decision_time": t0,
                            "log_likelihood": best_likelihood,
                        }

    if best_params_grid is None:
        # Fallback to reasonable defaults
        best_params_grid = {
            "drift_rate": 0.5,
            "boundary_separation": 1.0,
            "starting_point": 0.5,
            "non_decision_time": 0.2,
            "log_likelihood": float("-inf"),
        }

    return best_params_grid


def island_model_update(p: float, m: float, pm: float) -> float:
    """Update allele frequency in an island model of migration.

    In an island model, each subpopulation receives migrants from a mainland
    population with frequency pm at rate m, and then experiences drift.

    Args:
        p: Current allele frequency in subpopulation
        m: Migration rate (proportion of population replaced by migrants)
        pm: Allele frequency in mainland population

    Returns:
        New allele frequency after migration and drift
    """
    # Migration: weighted average between current and mainland frequency
    p_migrated = (1 - m) * p + m * pm

    # In this simple model, we return the migrated frequency
    # (drift would be applied separately if needed)
    return p_migrated
