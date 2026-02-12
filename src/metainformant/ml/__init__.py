"""Machine learning module for METAINFORMANT.

This module provides machine learning capabilities for bioinformatics analysis,
including classification, regression, dimensionality reduction, feature
selection/engineering, and local LLM inference.
"""

from __future__ import annotations

from . import automl, evaluation, features, interpretability, llm, models

__all__ = [
    "automl",
    "evaluation",
    "features",
    "interpretability",
    "llm",
    "models",
]
