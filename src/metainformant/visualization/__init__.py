"""Visualization and plotting utilities module for METAINFORMANT.

This module provides comprehensive plotting and visualization capabilities
for biological data, including statistical plots, genomic visualizations,
genomic networks, animations, publication-quality themes, interactive plots,
and composite dashboards.
"""

from __future__ import annotations

# Import subpackages
from . import analysis
from . import genomics
from . import interactive_dashboards
from . import plots

# Core utilities
from . import themes
from . import palettes
from . import composite
from . import interactive

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
from .genomics.trees import (
    plot_phylo_tree,
    circular_tree_plot,
    unrooted_tree_plot,
)
from .analysis import (
    dimred,
    information,
    statistical,
    timeseries,
    # Quality submodules
    quality_sequencing,
    quality_omics,
    quality_assessment,
)

# Theme and palette utilities
from .themes import (
    list_themes,
    get_theme,
    apply_theme,
    reset_theme,
    theme,
    register_theme,
)
from .palettes import (
    chromosome_palette,
    categorical,
    expression_gradient,
    significance_palette,
    significance_color,
    heatmap_cmap,
    alternating_pair,
    CHROMOSOME_COLORS,
    WONG,
    TOL_BRIGHT,
    IBM_COLORBLIND,
)

# Composite dashboards
from .composite import (
    multi_panel,
    genomic_overview,
    qc_summary,
)

# Interactive plots
from .interactive import (
    interactive_scatter,
    interactive_heatmap,
    interactive_volcano,
    interactive_manhattan,
)

# Interactive dashboards
from .interactive_dashboards.dashboards import (
    create_interactive_scatter,
    create_interactive_heatmap,
    create_genome_browser_track,
    create_interactive_volcano,
    export_to_html,
    create_dashboard,
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
