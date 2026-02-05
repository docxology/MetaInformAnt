"""Morphological phenotype analysis module.

Provides tools for morphometric measurements, shape analysis,
allometric regressions, and cross-specimen comparisons.
"""

from .measurement import Measurement
from .profile import (
    MorphometricProfile,
    coefficient_of_variation,
    allometric_regression,
    compare_profiles,
    summary_statistics,
)

__all__ = [
    "Measurement",
    "MorphometricProfile",
    "coefficient_of_variation",
    "allometric_regression",
    "compare_profiles",
    "summary_statistics",
]
