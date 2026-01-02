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
    export,
    expression,
    genomics,
    gwas_integration,
    information,
    information_integration,
    interactive,
    layout,
    life_events_integration,
    multidim,
    networks,
    plots,
    quality,
    singlecell_integration,
    statistical,
    style,
    timeseries,
    trees,
)

# Optional imports with graceful fallbacks
try:
    from . import amalgkit_visualization
except ImportError:
    amalgkit_visualization = None

try:
    from . import interactive
except ImportError:
    interactive = None

try:
    from . import gwas_integration
except ImportError:
    gwas_integration = None

try:
    from . import information_integration
except ImportError:
    information_integration = None

try:
    from . import life_events_integration
except ImportError:
    life_events_integration = None

try:
    from . import singlecell_integration
except ImportError:
    singlecell_integration = None

# Type checking imports
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
    # Core plotting functionality
    "plots",
    "basic",
    "style",
    "export",

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

    # Information theory visualizations
    "information",
    "information_integration",

    # Module integrations
    "amalgkit_visualization",
    "gwas_integration",
    "life_events_integration",
    "singlecell_integration",

    # Interactive and layout
    "interactive",
    "layout",
]



