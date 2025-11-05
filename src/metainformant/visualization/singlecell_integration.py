"""Single-cell visualization integration.

This module provides unified access to single-cell visualization functions
from the singlecell module, integrating them into the central visualization package.
"""

from __future__ import annotations

# Re-export single-cell visualization functions
try:
    from ..singlecell.visualization import (
        plot_embedding,
        plot_gene_expression,
        plot_pca,
        plot_qc_metrics,
        plot_qc_scatter,
    )

    SINGLECELL_VISUALIZATION_AVAILABLE = True
except ImportError:
    SINGLECELL_VISUALIZATION_AVAILABLE = False
    # Define placeholder functions
    def _not_available(*args, **kwargs):
        raise ImportError("Single-cell visualization functions not available. Install singlecell module dependencies.")
    
    plot_qc_metrics = _not_available
    plot_qc_scatter = _not_available
    plot_embedding = _not_available
    plot_pca = _not_available
    plot_gene_expression = _not_available

__all__ = [
    "plot_qc_metrics",
    "plot_qc_scatter",
    "plot_embedding",
    "plot_pca",
    "plot_gene_expression",
    "SINGLECELL_VISUALIZATION_AVAILABLE",
]

