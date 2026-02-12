"""Information geometry measures for statistical manifolds.

This module implements foundational information geometry functions including
Fisher-Rao geodesic distance, natural gradient computation, information
projection (m-projection via iterative scaling), alpha-divergences,
exponential family entropy, Hellinger distance, channel capacity
(Blahut-Arimoto), rate-distortion theory, the information bottleneck,
entropy power inequality, and information dimension
(Grassberger-Procaccia). These measures operate on the statistical manifold
of probability distributions and arise naturally in the analysis of
evolutionary fitness landscapes, population genetics parameter estimation,
neural coding optimality, and fractal structure in point-cloud omics data.

.. deprecated::
    This module is a backward-compatibility shim.  The implementation has
    been split into :mod:`fisher_rao` (distances, natural gradient, entropy,
    Hellinger) and :mod:`information_projection` (projection, divergences,
    channel capacity, rate-distortion, bottleneck, EPI, dimension).
    Import directly from those modules for new code.
"""

from __future__ import annotations

# Re-export everything from the split modules so that existing imports
# like ``from metainformant.information.metrics.advanced.geometry import fisher_rao_distance``
# continue to work.

from .fisher_rao import (
    _validate_distribution,
    exponential_family_entropy,
    fisher_rao_distance,
    hellinger_distance,
    natural_gradient,
)
from .information_projection import (
    _normalize_distribution,
    _safe_log,
    _validate_transition_matrix,
    channel_capacity,
    entropy_power_inequality,
    information_bottleneck,
    information_dimension,
    information_projection,
    rate_distortion_function,
    statistical_divergence,
)

__all__ = [
    # Fisher-Rao module
    "fisher_rao_distance",
    "natural_gradient",
    "exponential_family_entropy",
    "hellinger_distance",
    # Information projection module
    "information_projection",
    "statistical_divergence",
    "channel_capacity",
    "rate_distortion_function",
    "information_bottleneck",
    "entropy_power_inequality",
    "information_dimension",
]
