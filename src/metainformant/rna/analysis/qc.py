"""RNA-seq quality control module.

Comprehensive quality control metrics for RNA-seq count data including
sample-level and gene-level statistics, library complexity analysis,
batch effect detection, and bias assessment.

This module provides REAL implementations using numpy, scipy, and pandas.
No mocking, no placeholder data.

This module re-exports all public symbols from :mod:`qc_metrics` and
:mod:`qc_filtering` for backward compatibility.
"""

from __future__ import annotations

from metainformant.rna.analysis.qc_metrics import (
    classify_expression_level,
    compute_correlation_matrix,
    compute_gene_metrics,
    compute_sample_metrics,
    compute_saturation_curve,
    detect_outlier_samples,
    estimate_library_complexity,
)
from metainformant.rna.analysis.qc_filtering import (
    detect_batch_effects,
    detect_gc_bias,
    detect_length_bias,
    generate_qc_report,
)

__all__ = [
    "compute_sample_metrics",
    "detect_outlier_samples",
    "compute_gene_metrics",
    "classify_expression_level",
    "estimate_library_complexity",
    "compute_saturation_curve",
    "compute_correlation_matrix",
    "detect_batch_effects",
    "detect_gc_bias",
    "detect_length_bias",
    "generate_qc_report",
]
