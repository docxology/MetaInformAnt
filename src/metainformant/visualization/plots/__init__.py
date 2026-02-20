"""Visualization plots subpackage.

This module provides access to various plotting functions for data visualization."""
from __future__ import annotations

from . import animations, basic, general, multidim, specialized

# Re-export commonly used plot functions for direct access
from .basic import scatter_plot
from .general import (
    correlation_heatmap,
    expression_heatmap,
    manhattan_plot,
    pca_plot,
    qq_plot,
    volcano_plot,
)

__all__ = [
    'animations', 'basic', 'general', 'multidim', 'specialized',
    'scatter_plot', 'correlation_heatmap', 'expression_heatmap',
    'manhattan_plot', 'pca_plot', 'qq_plot', 'volcano_plot',
]
