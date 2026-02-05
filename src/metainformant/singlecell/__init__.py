"""Single-cell analysis module for METAINFORMANT.

This module provides tools for single-cell data analysis, including
dimensionality reduction, clustering, trajectory inference, multi-sample
integration, cell type annotation, differential expression, and RNA
velocity estimation.

Subpackages:
    analysis: Clustering, dimensionality reduction, trajectory inference
    data: Preprocessing, multi-sample integration
    visualization: Single-cell visualization tools
    celltyping: Cell type annotation and label transfer
    differential: Differential expression analysis
    velocity: RNA velocity estimation and analysis
"""

from __future__ import annotations

# Import subpackages
from . import analysis
from . import data
from . import visualization

# Import new subpackages
from . import celltyping
from . import differential
from . import velocity

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

# Import new modules for convenience
from .celltyping import annotation
from .differential import expression
from .velocity import rna_velocity

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
    "celltyping",
    "differential",
    "velocity",
    # Core analysis
    "clustering",
    "dimensionality",
    "trajectory",
    # Data processing
    "integration",
    "preprocessing",
    # Visualization
    "visualization_module",
    # Cell typing
    "annotation",
    # Differential expression
    "expression",
    # RNA velocity
    "rna_velocity",
]
