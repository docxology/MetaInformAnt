"""Differential expression analysis module for RNA-seq data.

This module provides comprehensive tools for normalizing count matrices,
performing differential expression analysis, dimensionality reduction,
and preparing data for visualization. All implementations are pure Python
using numpy, scipy, and pandas - no R dependency required.

Main Functions:
    Count Normalization:
        - normalize_counts: Apply various normalization methods (CPM, TPM, RPKM, etc.)
        - estimate_size_factors: DESeq2-style geometric mean normalization

    Differential Expression:
        - differential_expression: Run DE analysis with multiple methods
        - adjust_pvalues: Multiple testing correction (BH, Bonferroni)

    Dimensionality Reduction:
        - pca_analysis: Principal component analysis with loadings
        - compute_sample_distances: Distance matrix computation

    Gene Filtering:
        - filter_low_expression: Remove lowly expressed genes
        - get_highly_variable_genes: Identify most variable genes

    Visualization Helpers:
        - prepare_volcano_data: Add regulation status for volcano plots
        - prepare_ma_data: Prepare data for MA plots

Example:
    >>> from metainformant.rna.analysis import expression
    >>> import pandas as pd
    >>> counts = pd.DataFrame(...)  # genes x samples count matrix
    >>> conditions = ["control", "control", "treatment", "treatment"]
    >>> # Normalize counts
    >>> normalized = expression.normalize_counts(counts, method="cpm")
    >>> # Run differential expression
    >>> de_results = expression.differential_expression(
    ...     counts, conditions, method="deseq2_like"
    ... )
    >>> # Prepare for volcano plot
    >>> volcano_df = expression.prepare_volcano_data(de_results, fc_threshold=1.0)

.. deprecated::
    This module is a backward-compatibility shim.  The implementation has
    been split into :mod:`expression_core` (normalization, filtering) and
    :mod:`expression_analysis` (differential expression, PCA, visualization).
    Import directly from those modules for new code.
"""

from __future__ import annotations

# Re-export everything from the split modules so that existing imports
# like ``from metainformant.rna.analysis.expression import normalize_counts``
# continue to work.

from .expression_core import (
    NormalizationMethod,
    VarianceMethod,
    estimate_size_factors,
    filter_low_expression,
    get_highly_variable_genes,
    normalize_counts,
)
from .expression_analysis import (
    DEMethod,
    DistanceMethod,
    PValueMethod,
    adjust_pvalues,
    compute_sample_distances,
    differential_expression,
    pca_analysis,
    prepare_ma_data,
    prepare_volcano_data,
)

__all__ = [
    # Type aliases
    "NormalizationMethod",
    "DEMethod",
    "PValueMethod",
    "DistanceMethod",
    "VarianceMethod",
    # Normalization (expression_core)
    "normalize_counts",
    "estimate_size_factors",
    "filter_low_expression",
    "get_highly_variable_genes",
    # Differential expression (expression_analysis)
    "differential_expression",
    "adjust_pvalues",
    "pca_analysis",
    "compute_sample_distances",
    "prepare_volcano_data",
    "prepare_ma_data",
]
