"""Unified plotting and animation utilities.

Lightweight wrappers around matplotlib/seaborn and Biopython Phylo to provide a
consistent API for common visualizations used across domains.
"""

from .plots import lineplot, heatmap
from .animations import animate_time_series
from .trees import plot_tree

__all__ = [
    "lineplot",
    "heatmap",
    "animate_time_series",
    "plot_tree",
]

"""Visualization and animation utilities for METAINFORMANT.

Provides consistent, dependency-light wrappers around matplotlib/seaborn and
integrations with domain modules (e.g., phylogenetic trees).

Public API returns matplotlib objects, enabling composition in scripts and
notebooks, and keeping I/O concerns at the call site.
"""

from .plots import lineplot, heatmap, pairplot_dataframe
from .animations import animate_time_series
from .trees import plot_phylo_tree

__all__ = [
    "lineplot",
    "heatmap",
    "pairplot_dataframe",
    "animate_time_series",
    "plot_phylo_tree",
]




