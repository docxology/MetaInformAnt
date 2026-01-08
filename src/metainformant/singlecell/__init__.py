"""Single-cell analysis module for METAINFORMANT.

This module provides tools for single-cell data analysis, including
dimensionality reduction, clustering, trajectory inference, and
multi-sample integration.
"""

from __future__ import annotations

# Import subpackages
from . import analysis
from . import data
from . import visualization

# Import modules from subpackages for backward compatibility
from .analysis import (
    clustering,
    dimensionality,
    trajectory,
)
from .data import (
    integration,
    preprocessing,
)
from .visualization import visualization as visualization_module

# Optional imports with graceful fallbacks
try:
    from .analysis import clustering
except ImportError:
    clustering = None

try:
    from .analysis import trajectory
except ImportError:
    trajectory = None

try:
    from .data import integration
except ImportError:
    integration = None

# Type checking imports
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
    # Subpackages
    "analysis",
    "data",
    "visualization",

    # Core analysis
    "clustering",
    "dimensionality",
    "trajectory",

    # Data processing
    "integration",
    "preprocessing",

    # Visualization
    "visualization_module",
]



