"""Perception and Psychophysics submodule.

This module provides tools for modeling sensory perception and decision making,
including Psychophysics (Weber-Fechner, Stevens) and Signal Detection Theory.
"""

from .psychophysics import fechner_law, michelson_contrast, stevens_power_law, weber_contrast
from .signal_detection import criterion_c, d_prime, likelihood_ratio_beta, sdt_metrics

__all__ = [
    "weber_contrast",
    "michelson_contrast",
    "fechner_law",
    "stevens_power_law",
    "d_prime",
    "criterion_c",
    "likelihood_ratio_beta",
    "sdt_metrics",
]
