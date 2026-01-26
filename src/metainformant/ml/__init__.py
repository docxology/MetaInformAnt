"""Machine learning module for METAINFORMANT.

This module provides machine learning capabilities for bioinformatics analysis,
including classification, regression, dimensionality reduction, feature
selection/engineering, and local LLM inference.
"""

from __future__ import annotations

# Import subpackages
from . import evaluation
from . import features
from . import models
from . import llm

# Import modules from subpackages for backward compatibility
from .models import (
    classification,
    regression,
)
from .features import (
    dimensionality,
    features,
)
from .evaluation import (
    validation,
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
    # Models
    "classification",
    "regression",
    # Feature engineering
    "features",
    "dimensionality",
    # Validation and metrics
    "validation",
]
