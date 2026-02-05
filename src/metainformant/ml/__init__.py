"""Machine learning module for METAINFORMANT.

This module provides machine learning capabilities for bioinformatics analysis,
including classification, regression, dimensionality reduction, feature
selection/engineering, and local LLM inference.
"""

from __future__ import annotations

# Import subpackages
from . import automl, evaluation, features, interpretability, llm, models
from .evaluation import (
    validation,
)
from .features import (
    dimensionality,
    features,
)

# Import modules from subpackages for backward compatibility
from .models import (
    classification,
    regression,
)

# Optional imports with graceful fallbacks
try:
    from .models import classification
except ImportError:
    classification = None

try:
    from .models import regression
except ImportError:
    regression = None

try:
    from .features import dimensionality
except ImportError:
    dimensionality = None

# Interpretability subpackage
try:
    from .interpretability.explainers import (
        compute_attention_weights,
        compute_lime_explanation,
        compute_permutation_importance,
        compute_shap_values_kernel,
        feature_interaction,
        partial_dependence,
    )
    from .interpretability.feature_selection import (
        boruta_selection,
        mutual_information_selection,
        recursive_elimination,
        stability_selection,
    )
except ImportError:
    compute_permutation_importance = None  # type: ignore[assignment]
    compute_shap_values_kernel = None  # type: ignore[assignment]
    compute_lime_explanation = None  # type: ignore[assignment]
    feature_interaction = None  # type: ignore[assignment]
    partial_dependence = None  # type: ignore[assignment]
    compute_attention_weights = None  # type: ignore[assignment]
    boruta_selection = None  # type: ignore[assignment]
    recursive_elimination = None  # type: ignore[assignment]
    stability_selection = None  # type: ignore[assignment]
    mutual_information_selection = None  # type: ignore[assignment]

# AutoML subpackage
try:
    from .automl.optimization import (
        auto_preprocess,
        bayesian_optimization,
        grid_search,
        model_selection,
        random_search,
    )
except ImportError:
    random_search = None  # type: ignore[assignment]
    bayesian_optimization = None  # type: ignore[assignment]
    grid_search = None  # type: ignore[assignment]
    model_selection = None  # type: ignore[assignment]
    auto_preprocess = None  # type: ignore[assignment]

# Type checking imports
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
    # Subpackages
    "evaluation",
    "features",
    "models",
    "llm",
    "interpretability",
    "automl",
    # Models
    "classification",
    "regression",
    # Feature engineering
    "features",
    "dimensionality",
    # Validation and metrics
    "validation",
    # Interpretability
    "compute_permutation_importance",
    "compute_shap_values_kernel",
    "compute_lime_explanation",
    "feature_interaction",
    "partial_dependence",
    "compute_attention_weights",
    # Feature selection (interpretability)
    "boruta_selection",
    "recursive_elimination",
    "stability_selection",
    "mutual_information_selection",
    # AutoML
    "random_search",
    "bayesian_optimization",
    "grid_search",
    "model_selection",
    "auto_preprocess",
]
