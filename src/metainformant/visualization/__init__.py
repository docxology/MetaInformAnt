"""Visualization and animation utilities for METAINFORMANT.

Provides consistent, dependency-light wrappers around matplotlib/seaborn and
integrations with domain modules (e.g., phylogenetic trees).

Public API returns matplotlib objects, enabling composition in scripts and
notebooks, and keeping I/O concerns at the call site.
"""

from .plots import heatmap, lineplot, pairplot_dataframe
from .animations import animate_time_series
from .trees import plot_phylo_tree

# Backward-compat alias expected by tests
plot_tree = plot_phylo_tree

__all__ = [
    "lineplot",
    "heatmap",
    "pairplot_dataframe",
    "animate_time_series",
    "plot_tree",
    "plot_phylo_tree",
]




