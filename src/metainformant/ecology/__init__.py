"""Ecology and community analysis module for METAINFORMANT.

This module provides tools for ecological analysis, including community
structure, diversity metrics, and ecological visualization.
"""

from __future__ import annotations

# Import subpackages
from . import analysis
from . import visualization

# Import modules from subpackages for backward compatibility
from .analysis import (
    community,
)
from .visualization import visualization as visualization_module

# Optional imports with graceful fallbacks
try:
    from .analysis import community
except ImportError:
    community = None

# Type checking imports
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
    # Subpackages
    "analysis",
    "visualization",
]
