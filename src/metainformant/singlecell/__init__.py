"""Single-cell genomics analysis toolkit.

This module provides comprehensive tools for single-cell RNA sequencing (scRNA-seq)
and other single-cell omics data analysis, following best practices from Scanpy,
Seurat, and other leading single-cell analysis frameworks.

Key Features:
- Data loading and preprocessing
- Quality control and filtering
- Normalization and scaling
- Dimensionality reduction (PCA, UMAP, t-SNE)
- Clustering and cell type identification  
- Differential expression analysis
- Trajectory inference and pseudotime analysis
- Batch effect correction
- Multi-modal integration

Real Implementation Policy:
All functions perform actual computations without mocking. External dependencies
like scanpy or anndata are used when available, with graceful degradation when not.
"""

from __future__ import annotations

# Core functionality imports
from .preprocessing import (
    load_count_matrix,
    calculate_qc_metrics,
    filter_cells,
    filter_genes,
    normalize_counts,
    log_transform,
    scale_data,
)

from .dimensionality import (
    compute_pca,
    compute_umap,
    compute_tsne,
    select_hvgs,
    compute_neighbors,
)

from .clustering import (
    leiden_clustering,
    louvain_clustering,
    hierarchical_clustering,
    find_marker_genes,
)

from .visualization import (
    plot_qc_metrics,
    plot_dimensionality_reduction,
    plot_gene_expression,
    plot_clusters,
)

from .trajectory import (
    compute_pseudotime,
    trajectory_analysis,
    lineage_analysis,
)

from .integration import (
    batch_correction,
    integrate_datasets,
    harmony_integration,
)

__all__ = [
    # Preprocessing
    "load_count_matrix",
    "calculate_qc_metrics", 
    "filter_cells",
    "filter_genes",
    "normalize_counts",
    "log_transform",
    "scale_data",
    # Dimensionality reduction
    "compute_pca",
    "compute_umap", 
    "compute_tsne",
    "select_hvgs",
    "compute_neighbors",
    # Clustering
    "leiden_clustering",
    "louvain_clustering",
    "hierarchical_clustering",
    "find_marker_genes",
    # Visualization
    "plot_qc_metrics",
    "plot_dimensionality_reduction",
    "plot_gene_expression", 
    "plot_clusters",
    # Trajectory
    "compute_pseudotime",
    "trajectory_analysis",
    "lineage_analysis",
    # Integration
    "batch_correction",
    "integrate_datasets",
    "harmony_integration",
]
