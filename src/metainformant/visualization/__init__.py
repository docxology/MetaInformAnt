"""Visualization and plotting utilities module for METAINFORMANT.

This module provides comprehensive plotting and visualization capabilities
for biological data, including statistical plots, genomic visualizations,
genomic networks, animations, publication-quality themes, interactive plots,
and composite dashboards.
"""

from __future__ import annotations

# Core utilities
# Import subpackages
from . import analysis, composite, genomics, interactive, interactive_dashboards, palettes, plots, themes
from .analysis import (  # Quality submodules
    dimred,
    information,
    quality_assessment,
    quality_omics,
    quality_sequencing,
    statistical,
    timeseries,
)

# Composite dashboards
from .composite import (
    genomic_overview,
    multi_panel,
    qc_summary,
)
from .genomics import (
    expression,
)
from .genomics import genomics as genomics_module
from .genomics import (
    networks,
    trees,
)
from .genomics.trees import (
    circular_tree_plot,
    plot_phylo_tree,
    unrooted_tree_plot,
)

# Interactive plots
from .interactive import (
    interactive_heatmap,
    interactive_manhattan,
    interactive_scatter,
    interactive_volcano,
)

# Interactive dashboards
from .interactive_dashboards.dashboards import (
    create_dashboard,
    create_genome_browser_track,
    create_interactive_heatmap,
    create_interactive_scatter,
    create_interactive_volcano,
    export_to_html,
)
from .palettes import (
    CHROMOSOME_COLORS,
    IBM_COLORBLIND,
    TOL_BRIGHT,
    WONG,
    alternating_pair,
    categorical,
    chromosome_palette,
    expression_gradient,
    heatmap_cmap,
    significance_color,
    significance_palette,
)

# Re-export submodules for explicit imports
from .plots import (
    animations,
    basic,
    general,
    multidim,
    specialized,
)

# Import modules from subpackages for backward compatibility
from .plots.animations import *
from .plots.basic import *
from .plots.general import *
from .plots.multidim import *
from .plots.specialized import *

# Theme and palette utilities
from .themes import (
    apply_theme,
    get_theme,
    list_themes,
    register_theme,
    reset_theme,
    theme,
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
    # Core modules
    "themes",
    "palettes",
    "composite",
    "interactive",
    # Core plotting functionality
    "basic",
    "general",
    "multidim",
    "specialized",
    "animations",
    # Statistical visualizations
    "statistical",
    # Biological data visualizations
    "genomics_module",
    "expression",
    "networks",
    "trees",
    "plot_phylo_tree",
    "circular_tree_plot",
    "unrooted_tree_plot",
    # Advanced plotting
    "dimred",
    "timeseries",
    # Information theory visualizations
    "information",
    # Quality submodules
    "quality_sequencing",
    "quality_omics",
    "quality_assessment",
    # Theme utilities
    "list_themes",
    "get_theme",
    "apply_theme",
    "reset_theme",
    "theme",
    "register_theme",
    # Palette utilities
    "chromosome_palette",
    "categorical",
    "expression_gradient",
    "significance_palette",
    "significance_color",
    "heatmap_cmap",
    "alternating_pair",
    "CHROMOSOME_COLORS",
    "WONG",
    "TOL_BRIGHT",
    "IBM_COLORBLIND",
    # Composite dashboards
    "multi_panel",
    "genomic_overview",
    "qc_summary",
    # Interactive plots
    "interactive_scatter",
    "interactive_heatmap",
    "interactive_volcano",
    "interactive_manhattan",
    # Interactive dashboards
    "interactive_dashboards",
    "create_interactive_scatter",
    "create_interactive_heatmap",
    "create_genome_browser_track",
    "create_interactive_volcano",
    "export_to_html",
    "create_dashboard",
    # Domain-specific visualizations
    "protein",
    "epigenome",
    "ontology",
    "phenotype",
    "ecology",
    "math",
    "multiomics",
]
