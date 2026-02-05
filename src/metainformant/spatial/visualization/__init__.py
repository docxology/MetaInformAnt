"""Spatial visualization module for METAINFORMANT.

Provides publication-quality plotting functions for spatial transcriptomics
data, including spatial scatter plots, tissue overlays, expression maps,
and analysis result visualizations.
"""

from __future__ import annotations

from . import plots

from .plots import (
    plot_cell_type_map,
    plot_deconvolution_pie,
    plot_gene_expression_map,
    plot_interaction_heatmap,
    plot_neighborhood_graph,
    plot_spatial_autocorrelation,
    plot_spatial_scatter,
    plot_tissue_overlay,
)

__all__ = [
    "plots",
    "plot_cell_type_map",
    "plot_deconvolution_pie",
    "plot_gene_expression_map",
    "plot_interaction_heatmap",
    "plot_neighborhood_graph",
    "plot_spatial_autocorrelation",
    "plot_spatial_scatter",
    "plot_tissue_overlay",
]
