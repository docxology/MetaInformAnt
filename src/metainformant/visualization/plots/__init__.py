"""Visualization plots subpackage.

This module provides access to various plotting functions for data visualization.
"""

from __future__ import annotations

# Basic plots
from .basic import (
    lineplot,
    scatter_plot,
    heatmap,
    bar_plot,
    pie_chart,
    area_plot,
    step_plot,
)

# Specialized plots
from .specialized import (
    plot_venn_diagram,
    plot_sankey_diagram,
    plot_chord_diagram,
    plot_alluvial_diagram,
    plot_circular_barplot,
    plot_network_circular_layout,
    plot_upset_plot,
)

# Multi-dimensional plots
from .multidim import (
    plot_pairwise_relationships,
    plot_parallel_coordinates,
    plot_radar_chart,
    plot_3d_scatter,
)

# General plots
from .general import (
    expression_heatmap,
    pca_plot,
    correlation_heatmap,
    qq_plot,
    volcano_plot,
    manhattan_plot,
)

# Animations (optional)
try:
    from .animations import (
        create_animation,
        animate_time_series,
    )
except ImportError:
    pass

__all__ = [
    # Basic plots
    "lineplot",
    "scatter_plot",
    "heatmap",
    "bar_plot",
    "pie_chart",
    "area_plot",
    "step_plot",
    # Specialized plots
    "plot_venn_diagram",
    "plot_sankey_diagram",
    "plot_chord_diagram",
    "plot_alluvial_diagram",
    "plot_circular_barplot",
    "plot_network_circular_layout",
    "plot_upset_plot",
    # Multi-dimensional plots
    "plot_pairwise_relationships",
    "plot_parallel_coordinates",
    "plot_radar_chart",
    "plot_3d_scatter",
    # General plots
    "expression_heatmap",
    "pca_plot",
    "correlation_heatmap",
    "qq_plot",
    "volcano_plot",
    "manhattan_plot",
]
