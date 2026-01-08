"""Phenotype analysis module for METAINFORMANT.

This module provides tools for analyzing phenotypic data, including
trait extraction from text (AntWiki), life course analysis, and
morphological data visualization.
"""

from __future__ import annotations

# Import subpackages
from . import analysis
from . import data
from . import visualization

# Import modules from subpackages for backward compatibility
from .data import (
    antwiki,
    scraper,
)
from .analysis import (
    life_course,
)
from .visualization import visualization as visualization_module

# Optional imports with graceful fallbacks
try:
    from .data import scraper
except ImportError:
    scraper = None

# Type checking imports
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
    # Subpackages
    "analysis",
    "data",
    "visualization",

    # Data sources
    "antwiki",
    "scraper",

    # Analysis
    "life_course",

    # Visualization
    "visualization_module",
]
