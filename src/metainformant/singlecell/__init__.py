"""Single-cell RNA-seq analysis module for METAINFORMANT.

This module provides comprehensive tools for single-cell transcriptomics analysis,
including preprocessing, dimensionality reduction, clustering, trajectory inference,
differential expression analysis, and visualization.
"""

from __future__ import annotations

# Import all single-cell analysis submodules
from . import (
    clustering,
    dimensionality,
    integration,
    preprocessing,
    trajectory,
    visualization,
)

# Optional imports with graceful fallbacks
try:
    from . import trajectory
except ImportError:
    trajectory = None

try:
    from . import integration
except ImportError:
    integration = None

try:
    from . import visualization
except ImportError:
    visualization = None

# Type checking imports
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
    # Core single-cell analysis
    "preprocessing",
    "dimensionality",
    "clustering",

    # Advanced analysis
    "trajectory",
    "integration",

    # Visualization
    "visualization",
]







