"""Visualization and plotting utilities module for METAINFORMANT.

This module provides comprehensive plotting and visualization capabilities
for biological data, including statistical plots, genomic visualizations,
network graphs, animations, and publication-quality figure generation.
"""

from __future__ import annotations

# Import all visualization submodules
from . import (
    amalgkit_visualization,
    animations,
    basic,
    dimred,
    expression,
    genomics,
    information,
    multidim,
    networks,
    quality,
    specialized,
    statistical,
    timeseries,
    trees,
)

# Optional imports with graceful fallbacks
try:
    from . import amalgkit_visualization
except ImportError:
    amalgkit_visualization = None

try:
    from metainformant.protein import visualization as protein
except ImportError:
    protein = None

try:
    from metainformant.epigenome import visualization as epigenome
except ImportError:
    epigenome = None

try:
    from metainformant.ontology import visualization as ontology
except ImportError:
    ontology = None

try:
    from metainformant.phenotype import visualization as phenotype
except ImportError:
    phenotype = None

try:
    from metainformant.ecology import visualization as ecology
except ImportError:
    ecology = None

try:
    from metainformant.math import visualization as math
except ImportError:
    math = None

try:
    from metainformant.multiomics import visualization as multiomics
except ImportError:
    multiomics = None

# Type checking imports
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
    # Core plotting functionality
    "basic",

    # Statistical visualizations
    "statistical",
    "quality",

    # Biological data visualizations
    "genomics",
    "expression",
    "networks",
    "trees",

    # Advanced plotting
    "dimred",
    "multidim",
    "timeseries",
    "animations",
    "specialized",

    # Information theory visualizations
    "information",

    # Module integrations
    "amalgkit_visualization",

    # Domain-specific visualizations
    "protein",
    "epigenome",
    "ontology",
    "phenotype",
    "ecology",
    "math",
    "multiomics",
]






