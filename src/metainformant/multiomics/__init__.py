"""Multi-omic data integration module for METAINFORMANT.

This module provides tools for integrating and analyzing multiple omics data types,
including cross-platform harmonization, joint dimensionality reduction,
correlation analysis, and systems-level biological interpretation.
"""

from __future__ import annotations

# Import all multiomics submodules
from . import integration

# Type checking imports
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
    # Core integration functionality
    "integration",
]


