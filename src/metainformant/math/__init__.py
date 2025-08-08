"""Mathematical and theoretical biology utilities.

This subpackage provides quantitative biology primitives used across domains:
- Price equation and selection decomposition
- Kin and multilevel selection toy calculators
- Drift–Diffusion models (DDM)

These are deliberately lightweight, dependency-minimal, and intended for
composition with `core` utilities and `dna`/`rna`/`protein` domain modules.
"""

from .price import price_equation, covariance, expectation
from .selection import kin_selection_response, multilevel_selection_decomposition
from .ddm import ddm_analytic_accuracy, ddm_mean_decision_time

__all__ = [
    "price_equation",
    "covariance",
    "expectation",
    "kin_selection_response",
    "multilevel_selection_decomposition",
    "ddm_analytic_accuracy",
    "ddm_mean_decision_time",
]


