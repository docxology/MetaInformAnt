"""Analysis visualization subpackage.

Provides dimensionality reduction, statistical, time-series, quality control,
and information-theoretic visualization functions.
"""

from __future__ import annotations

# Dimensionality reduction
from .dimred import (
    plot_pca,
    plot_umap,
    plot_tsne,
    plot_pca_loadings,
    biplot,
)

# Statistical visualizations
from .statistical import (
    histogram,
    box_plot,
    violin_plot,
    qq_plot,
    correlation_heatmap,
    density_plot,
    ridge_plot,
    roc_curve,
    precision_recall_curve,
    residual_plot,
    leverage_plot,
)

# Time series
from .timeseries import (
    plot_time_series,
    plot_autocorrelation,
    plot_seasonal_decomposition,
    plot_forecast,
    plot_trend_analysis,
)

# Information theory
from .information import (
    plot_entropy_profile,
    plot_mutual_information_matrix,
    plot_renyi_spectra,
    plot_information_landscape,
    plot_information_network,
)

# Quality control - sequencing
from .quality_sequencing import (
    plot_quality_metrics,
    plot_adapter_content,
    plot_gc_distribution,
    plot_length_distribution,
    plot_per_base_quality_boxplot,
    plot_sequence_duplication_levels,
    plot_overrepresented_sequences,
    plot_kmer_profiles,
)

# Quality control - omics
from .quality_omics import (
    plot_vcf_quality_metrics,
    plot_singlecell_qc_metrics,
    plot_protein_structure_quality,
    plot_multiomics_quality_overview,
)

# Quality control - assessment
from .quality_assessment import (
    plot_coverage_uniformity,
    plot_error_profiles,
    plot_batch_effects_qc,
    plot_data_integrity_metrics,
)

# Backward compatibility re-export from original quality.py
# (keeping the old module importable but directing to submodules)
from . import quality_sequencing as quality

__all__ = [
    # Submodules
    "quality",
    # Dimensionality reduction
    "plot_pca",
    "plot_umap",
    "plot_tsne",
    "plot_pca_loadings",
    "biplot",
    # Statistical
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
    # Time series
    "plot_time_series",
    "plot_autocorrelation",
    "plot_seasonal_decomposition",
    "plot_forecast",
    "plot_trend_analysis",
    # Information
    "plot_entropy_profile",
    "plot_mutual_information_matrix",
    "plot_renyi_spectra",
    "plot_information_landscape",
    "plot_information_network",
    # Quality - sequencing
    "plot_quality_metrics",
    "plot_adapter_content",
    "plot_gc_distribution",
    "plot_length_distribution",
    "plot_per_base_quality_boxplot",
    "plot_sequence_duplication_levels",
    "plot_overrepresented_sequences",
    "plot_kmer_profiles",
    # Quality - omics
    "plot_vcf_quality_metrics",
    "plot_singlecell_qc_metrics",
    "plot_protein_structure_quality",
    "plot_multiomics_quality_overview",
    # Quality - assessment
    "plot_coverage_uniformity",
    "plot_error_profiles",
    "plot_batch_effects_qc",
    "plot_data_integrity_metrics",
]
