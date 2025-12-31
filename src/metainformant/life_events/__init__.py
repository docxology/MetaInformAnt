"""Life course event analysis and temporal modeling module for METAINFORMANT.

This module provides tools for analyzing life course trajectories, event sequence modeling,
temporal pattern recognition, and outcome prediction from longitudinal data.
"""

from __future__ import annotations

# Import all life events analysis submodules
from . import (
    config,
    embeddings,
    events,
    interpretability,
    models,
    utils,
    visualization,
    workflow,
)

# Optional imports with graceful fallbacks
try:
    from . import interpretability
except ImportError:
    interpretability = None

try:
    from . import visualization
except ImportError:
    visualization = None

# Type checking imports
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
    # Core data structures and processing
    "events",
    "utils",

    # Modeling and analysis
    "embeddings",
    "models",
    "workflow",

    # Configuration and visualization
    "config",
    "visualization",

    # Advanced analysis
    "interpretability",
]


