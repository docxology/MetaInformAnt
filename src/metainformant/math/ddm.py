from __future__ import annotations

import math


def ddm_analytic_accuracy(drift_rate: float, boundary: float, noise_sd: float = 1.0) -> float:
    """Calculate closed-form accuracy for symmetric drift-diffusion model (DDM).
    
    Computes the probability of correct choice in a two-alternative forced choice
    decision-making task. The DDM models decision-making as a stochastic process
    where evidence accumulates until a threshold is reached.
    
    Args:
        drift_rate: Mean rate of evidence accumulation (v)
        boundary: Decision threshold (a)
        noise_sd: Standard deviation of noise (sigma), default 1.0
        
    Returns:
        Probability of correct choice in [0, 1]. Returns 0.5 for invalid inputs.
        Formula: P(correct) = 1 / (1 + exp(-2*v*a/sigma^2)) when starting at 0.
        
    Examples:
        >>> ddm_analytic_accuracy(drift_rate=0.5, boundary=1.0, noise_sd=1.0)
        0.731...
        >>> ddm_analytic_accuracy(drift_rate=0.0, boundary=1.0, noise_sd=1.0)
        0.5
        
    References:
        Ratcliff, R., & McKoon, G. (2008). The diffusion decision model: theory
        and data for two-choice decision tasks. Neural computation, 20(4), 873-922.
    """
    if noise_sd <= 0 or boundary <= 0:
        return 0.5
    k = -2.0 * drift_rate * boundary / (noise_sd**2)
    try:
        return 1.0 / (1.0 + math.exp(k))
    except OverflowError:
        return 1.0 if k < 0 else 0.0


def ddm_mean_decision_time(drift_rate: float, boundary: float, noise_sd: float = 1.0) -> float:
    """Calculate approximate mean decision time for symmetric drift-diffusion model.
    
    Estimates the expected time until a decision boundary is reached in a
    two-alternative forced choice task.
    
    Args:
        drift_rate: Mean rate of evidence accumulation (v)
        boundary: Decision threshold (a)
        noise_sd: Standard deviation of noise (sigma), default 1.0
        
    Returns:
        Expected decision time. Returns 0.0 for invalid inputs.
        Formula: E[T] ≈ (a/v) * tanh(a*v/sigma^2) for v != 0.
        For small v (drift_rate < 1e-8), uses diffusion-limited approximation a^2/sigma^2.
        
    Examples:
        >>> ddm_mean_decision_time(drift_rate=0.5, boundary=1.0, noise_sd=1.0)
        1.928...
        >>> ddm_mean_decision_time(drift_rate=0.0, boundary=1.0, noise_sd=1.0)
        1.0
        
    References:
        Ratcliff, R., & McKoon, G. (2008). The diffusion decision model: theory
        and data for two-choice decision tasks. Neural computation, 20(4), 873-922.
    """
    if noise_sd <= 0 or boundary <= 0:
        return 0.0
    if abs(drift_rate) < 1e-8:
        # diffusion-limited; scale with a^2 / sigma^2
        return (boundary**2) / (noise_sd**2)
    x = (boundary * drift_rate) / (noise_sd**2)
    return (boundary / drift_rate) * math.tanh(x)
