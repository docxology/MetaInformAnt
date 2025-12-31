"""Information theory analysis module for METAINFORMANT.

This module provides comprehensive information-theoretic methods for biological data,
including syntactic information (Shannon entropy, mutual information),
semantic information (information content, semantic similarity),
continuous information theory, and analysis workflows.
"""

from __future__ import annotations

# Import all information theory submodules
from . import (
    analysis,
    continuous,
    estimation,
    integration,
    networks,
    semantic,
    syntactic,
    visualization,
    workflows,
)

# Optional imports with graceful fallbacks
try:
    from . import integration
except ImportError:
    integration = None

try:
    from . import networks
except ImportError:
    networks = None

try:
    from . import visualization
except ImportError:
    visualization = None

# Type checking imports
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
    # Core information measures
    "syntactic",
    "semantic",
    "continuous",

    # Analysis and estimation
    "analysis",
    "estimation",
    "workflows",

    # Specialized applications
    "integration",
    "networks",
    "visualization",
]
