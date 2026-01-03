"""Phenotypic trait analysis and curation module for METAINFORMANT.

This module provides tools for phenotypic data analysis, trait curation,
life course analysis, and integration with biological databases like AntWiki.
"""

from __future__ import annotations

# Import all phenotype analysis submodules
from . import (
    antwiki,
    life_course,
    scraper,
    visualization,
)

# Direct imports of commonly used classes and functions
from .life_course import analyze_life_course

# Optional imports with graceful fallbacks
try:
    from . import antwiki
except ImportError:
    antwiki = None

try:
    from . import scraper
except ImportError:
    scraper = None

# Type checking imports
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
    # Core phenotype analysis
    "life_course",
    "analyze_life_course",

    # Data sources and curation
    "antwiki",
    "scraper",

    # Visualization
    "visualization",
]



