"""Visualization and animation utilities for METAINFORMANT.

Provides consistent, dependency-light wrappers around matplotlib/seaborn and
integrations with domain modules (e.g., phylogenetic trees).

Public API returns matplotlib objects, enabling composition in scripts and
notebooks, and keeping I/O concerns at the call site.
"""

# Basic plots
from .basic import (
    area_plot,
    bar_plot,
    heatmap,
    lineplot,
    pie_chart,
    scatter_plot,
    step_plot,
)

# Statistical plots
from .statistical import (
    box_plot,
    correlation_heatmap,
    density_plot,
    histogram,
    leverage_plot,
    precision_recall_curve,
    qq_plot,
    residual_plot,
    ridge_plot,
    roc_curve,
    violin_plot,
)

# Genomics plots
from .genomics import (
    chromosome_ideogram,
    circular_manhattan_plot,
    coverage_plot,
    manhattan_plot,
    regional_plot,
    variant_plot,
    volcano_plot,
)

# Expression plots
from .expression import (
    differential_expression_plot,
    enrichment_plot,
    expression_heatmap,
    gene_expression_plot,
    log_fold_change_plot,
)

# Dimensionality reduction
from .dimred import (
    biplot,
    pca_loadings_plot,
    pca_plot,
    pca_scree_plot,
    tsne_plot,
    umap_plot,
)

# Network plots
from .networks import (
    circular_network_plot,
    community_network_plot,
    force_directed_plot,
    hierarchical_network_plot,
    network_plot,
)

# Time series
from .timeseries import (
    autocorrelation_plot,
    forecast_plot,
    seasonal_decomposition_plot,
    time_series_plot,
    trend_plot,
)

# Multi-dimensional plots
from .multidim import (
    pairplot_dataframe,
    parallel_coordinates_plot,
    radar_chart,
    scatter_3d,
    splom_plot,
)

# Quality control
from .quality import (
    adapter_content_plot,
    per_base_quality_plot,
    qc_metrics_plot,
    quality_score_plot,
    sequence_length_distribution,
)

# Information theory
from .information import (
    entropy_plot,
    information_network_plot,
    information_profile_plot,
    mutual_information_plot,
    renyi_spectrum_plot,
)

# Trees
from .trees import (
    circular_tree_plot,
    plot_phylo_tree,
    tree_annotation_plot,
    tree_comparison_plot,
    unrooted_tree_plot,
)

# Animations
from .animations import (
    animate_clustering,
    animate_evolution,
    animate_network,
    animate_time_series,
    animate_trajectory,
)

# Backward-compat alias expected by tests
plot_tree = plot_phylo_tree

__all__ = [
    # Basic plots
    "lineplot",
    "scatter_plot",
    "bar_plot",
    "pie_chart",
    "area_plot",
    "step_plot",
    "heatmap",
    # Statistical plots
    "histogram",
    "box_plot",
    "violin_plot",
    "qq_plot",
    "correlation_heatmap",
    "density_plot",
    "ridge_plot",
    "roc_curve",
    "precision_recall_curve",
    "residual_plot",
    "leverage_plot",
    # Genomics plots
    "manhattan_plot",
    "volcano_plot",
    "regional_plot",
    "circular_manhattan_plot",
    "chromosome_ideogram",
    "coverage_plot",
    "variant_plot",
    # Expression plots
    "expression_heatmap",
    "enrichment_plot",
    "gene_expression_plot",
    "differential_expression_plot",
    "log_fold_change_plot",
    # Dimensionality reduction
    "pca_plot",
    "umap_plot",
    "tsne_plot",
    "pca_scree_plot",
    "pca_loadings_plot",
    "biplot",
    # Network plots
    "network_plot",
    "circular_network_plot",
    "hierarchical_network_plot",
    "force_directed_plot",
    "community_network_plot",
    # Time series
    "time_series_plot",
    "autocorrelation_plot",
    "seasonal_decomposition_plot",
    "forecast_plot",
    "trend_plot",
    # Multi-dimensional
    "pairplot_dataframe",
    "parallel_coordinates_plot",
    "radar_chart",
    "splom_plot",
    "scatter_3d",
    # Quality control
    "qc_metrics_plot",
    "quality_score_plot",
    "per_base_quality_plot",
    "adapter_content_plot",
    "sequence_length_distribution",
    # Information theory
    "entropy_plot",
    "mutual_information_plot",
    "information_profile_plot",
    "renyi_spectrum_plot",
    "information_network_plot",
    # Phylogenetic trees
    "plot_tree",
    "plot_phylo_tree",
    "circular_tree_plot",
    "unrooted_tree_plot",
    "tree_comparison_plot",
    "tree_annotation_plot",
    # Animation
    "animate_time_series",
    "animate_evolution",
    "animate_clustering",
    "animate_network",
    "animate_trajectory",
]
