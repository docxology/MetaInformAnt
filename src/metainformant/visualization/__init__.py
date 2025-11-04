"""Visualization and animation utilities for METAINFORMANT.

Provides consistent, dependency-light wrappers around matplotlib/seaborn and
integrations with domain modules (e.g., phylogenetic trees).

Public API returns matplotlib objects, enabling composition in scripts and
notebooks, and keeping I/O concerns at the call site.
"""

from .animations import animate_time_series
from .plots import (
    box_plot,
    correlation_heatmap,
    enrichment_plot,
    expression_heatmap,
    heatmap,
    histogram,
    lineplot,
    manhattan_plot,
    network_plot,
    pairplot_dataframe,
    pca_plot,
    qq_plot,
    scatter_plot,
    violin_plot,
    volcano_plot,
)
from .trees import plot_phylo_tree

# Backward-compat alias expected by tests
plot_tree = plot_phylo_tree

__all__ = [
    # Basic plots
    "lineplot",
    "scatter_plot",
    "histogram",
    "box_plot",
    "violin_plot",
    # Advanced plots
    "heatmap",
    "correlation_heatmap",
    "pairplot_dataframe",
    # Domain-specific plots
    "volcano_plot",
    "expression_heatmap",
    "pca_plot",
    "manhattan_plot",
    "qq_plot",
    "enrichment_plot",
    "network_plot",
    # Phylogenetic trees
    "plot_tree",
    "plot_phylo_tree",
    # Animation
    "animate_time_series",
]
