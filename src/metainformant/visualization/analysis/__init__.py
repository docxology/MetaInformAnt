"""Analysis visualization subpackage.

Provides dimensionality reduction, statistical, time-series, quality control,
and information-theoretic visualization functions."""
from __future__ import annotations

from . import dimred, information, quality, quality_assessment, quality_omics, quality_sequencing, statistical, timeseries

# Re-export key functions for direct access via `analysis.plot_pca()` etc.
from .dimred import biplot, plot_pca, plot_tsne, plot_umap
from .information import plot_entropy_profile, plot_mutual_information_matrix
from .quality_assessment import plot_coverage_uniformity, plot_error_profiles
from .quality_omics import plot_singlecell_qc_metrics, plot_vcf_quality_metrics
from .quality_sequencing import plot_gc_distribution, plot_kmer_profiles, plot_quality_metrics
from .statistical import box_plot, correlation_heatmap, histogram, qq_plot, roc_curve, violin_plot
from .timeseries import plot_autocorrelation, plot_forecast, plot_time_series

__all__ = [
    'dimred', 'information', 'quality', 'quality_assessment',
    'quality_omics', 'quality_sequencing', 'statistical', 'timeseries',
    'plot_pca', 'plot_umap', 'plot_tsne', 'biplot',
    'histogram', 'box_plot', 'violin_plot', 'qq_plot',
    'correlation_heatmap', 'roc_curve',
    'plot_time_series', 'plot_autocorrelation', 'plot_forecast',
    'plot_entropy_profile', 'plot_mutual_information_matrix',
    'plot_quality_metrics', 'plot_gc_distribution', 'plot_kmer_profiles',
    'plot_vcf_quality_metrics', 'plot_singlecell_qc_metrics',
    'plot_coverage_uniformity', 'plot_error_profiles',
]
