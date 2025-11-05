"""Backward compatibility wrapper for plots module.

This module maintains backward compatibility by re-exporting functions
from the new modular structure. New code should import directly from
the specific modules (basic, statistical, genomics, etc.).
"""

from __future__ import annotations

# Re-export from new modules for backward compatibility
from .basic import heatmap, lineplot, scatter_plot
from .expression import enrichment_plot, expression_heatmap
from .genomics import manhattan_plot, volcano_plot
from .dimred import pca_plot
from .networks import network_plot
from .statistical import (
    box_plot,
    correlation_heatmap,
    histogram,
    qq_plot,
    violin_plot,
)
from .multidim import pairplot_dataframe

__all__ = [
    "lineplot",
    "heatmap",
    "pairplot_dataframe",
    "scatter_plot",
    "histogram",
    "box_plot",
    "violin_plot",
    "correlation_heatmap",
    "volcano_plot",
    "expression_heatmap",
    "pca_plot",
    "manhattan_plot",
    "qq_plot",
    "enrichment_plot",
    "network_plot",
]
