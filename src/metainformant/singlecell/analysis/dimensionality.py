"""Dimensionality reduction methods for single-cell data.

This module provides functions for reducing the dimensionality of single-cell
expression data using various algorithms including PCA, t-SNE, UMAP, and
diffusion maps. These methods help visualize and analyze high-dimensional
single-cell data.

.. deprecated::
    This module is a backward-compatibility shim.  The implementation has
    been split into :mod:`pca_methods` (linear methods) and
    :mod:`nonlinear_methods` (nonlinear methods).  Import directly from
    those modules for new code.
"""

from __future__ import annotations

# Re-export everything from the split modules so that existing imports
# like ``from metainformant.singlecell.analysis.dimensionality import pca_reduction``
# continue to work.

from .pca_methods import (
    compute_pca,
    factor_analysis_reduction,
    ica_reduction,
    pca_reduction,
    run_pca,
    select_hvgs,
)
from .nonlinear_methods import (
    SKLEARN_AVAILABLE,
    compute_diffusion_map,
    compute_dimensionality_metrics,
    compute_neighbors,
    compute_tsne,
    compute_umap,
    diffusion_map_reduction,
    mds_reduction,
    run_tsne,
    run_umap,
    tsne_reduction,
    umap_reduction,
)

__all__ = [
    # Linear methods (pca_methods)
    "pca_reduction",
    "ica_reduction",
    "factor_analysis_reduction",
    "run_pca",
    "compute_pca",
    "select_hvgs",
    # Nonlinear methods (nonlinear_methods)
    "tsne_reduction",
    "umap_reduction",
    "diffusion_map_reduction",
    "mds_reduction",
    "compute_dimensionality_metrics",
    "run_tsne",
    "run_umap",
    "compute_tsne",
    "compute_umap",
    "compute_neighbors",
    "compute_diffusion_map",
]
