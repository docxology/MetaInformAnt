"""Multi-omics integration module for METAINFORMANT.

This module provides tools for integrating multiple omics layers, including
genomics, transcriptomics, epigenomics, and proteomics.
"""

from __future__ import annotations

# Import subpackages
from . import analysis
from . import visualization

# Import modules from subpackages for backward compatibility
from .analysis import (
    integration,
)
from .visualization import visualization as visualization_module

# Optional imports with graceful fallbacks
try:
    from .analysis import integration
except ImportError:
    integration = None

# Type checking imports
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
    # Subpackages
    "analysis",
    "visualization",
]

