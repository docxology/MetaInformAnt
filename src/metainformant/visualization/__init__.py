"""Visualization and plotting utilities module for METAINFORMANT.

This module provides comprehensive plotting and visualization capabilities
for biological data, including statistical plots, genomic visualizations,
genomic networks, animations, and publication-quality figure generation.
"""

from __future__ import annotations

# Import subpackages
from . import analysis
from . import genomics
from . import plots

# Import modules from subpackages for backward compatibility
from .plots.animations import *
from .plots.basic import *
from .plots.general import *
from .plots.multidim import *
from .plots.specialized import *

# Re-export submodules for explicit imports
from .plots import (
    animations,
    basic,
    general,
    multidim,
    specialized,
)
from .genomics import (
    expression,
    genomics as genomics_module,
    networks,
    trees,
)
from .analysis import (
    dimred,
    information,
    quality,
    statistical,
    timeseries,
)

# Optional imports with graceful fallbacks
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
    # Subpackages
    "analysis",
    "genomics",
    "plots",

    # Core plotting functionality
    "basic",
    "plots",

    # Statistical visualizations
    "statistical",
    "quality",

    # Biological data visualizations
    "genomics_module",
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

    # Domain-specific visualizations
    "protein",
    "epigenome",
    "ontology",
    "phenotype",
    "ecology",
    "math",
    "multiomics",
]






