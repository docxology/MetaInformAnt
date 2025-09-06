from __future__ import annotations

import math


def ddm_analytic_accuracy(drift_rate: float, boundary: float, noise_sd: float = 1.0) -> float:
    """Closed-form accuracy for symmetric DDM with unbiased start.

    P(correct) = 1 / (1 + exp(-2*v*a/sigma^2)) when starting at 0.
    """
    if noise_sd <= 0 or boundary <= 0:
        return 0.5
    k = -2.0 * drift_rate * boundary / (noise_sd**2)
    try:
        return 1.0 / (1.0 + math.exp(k))
    except OverflowError:
        return 1.0 if k < 0 else 0.0


def ddm_mean_decision_time(drift_rate: float, boundary: float, noise_sd: float = 1.0) -> float:
    """Approximate mean decision time for symmetric DDM.

    E[T] ≈ (a/v) * tanh(a*v/sigma^2) for v != 0; for small v, use series.
    """
    if noise_sd <= 0 or boundary <= 0:
        return 0.0
    if abs(drift_rate) < 1e-8:
        # diffusion-limited; scale with a^2 / sigma^2
        return (boundary**2) / (noise_sd**2)
    x = (boundary * drift_rate) / (noise_sd**2)
    return (boundary / drift_rate) * math.tanh(x)
