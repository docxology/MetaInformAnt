"""Machine learning utilities for biological data analysis in METAINFORMANT.

This module provides machine learning tools tailored for biological applications,
including classification, regression, feature selection, dimensionality reduction,
and model validation for genomics, transcriptomics, and other biological data.
"""

from __future__ import annotations

# Import all machine learning submodules
from . import (
    classification,
    dimensionality,
    features,
    regression,
    validation,
)

# Optional imports with graceful fallbacks
try:
    from . import dimensionality
except ImportError:
    dimensionality = None

# Type checking imports
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
    # Core ML tasks
    "classification",
    "regression",

    # Feature engineering
    "features",
    "dimensionality",

    # Model evaluation
    "validation",
]








