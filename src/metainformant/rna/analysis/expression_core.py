"""Count matrix normalization, size factor estimation, and gene filtering.

This module provides tools for normalizing RNA-seq count matrices using
various methods (CPM, TPM, RPKM, quantile, median ratio), estimating
DESeq2-style size factors, filtering lowly expressed genes, and selecting
highly variable genes.

All implementations are pure Python using numpy, scipy, and pandas.
"""

from __future__ import annotations

from typing import Dict, List, Literal, Optional, Union

import numpy as np
import pandas as pd

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


# =============================================================================
# Type Definitions
# =============================================================================

NormalizationMethod = Literal["cpm", "tpm", "rpkm", "log2cpm", "quantile", "median_ratio"]
VarianceMethod = Literal["cv", "variance"]


# =============================================================================
# Count Matrix Normalization
# =============================================================================


def normalize_counts(
    counts_df: pd.DataFrame,
    method: NormalizationMethod = "cpm",
    gene_lengths: Optional[pd.Series] = None,
) -> pd.DataFrame:
    """Normalize a count matrix using the specified method.

    Applies standard RNA-seq normalization methods to transform raw counts
    into comparable expression values across samples.

    Args:
        counts_df: Raw count matrix with genes as rows and samples as columns.
            All values must be non-negative integers or floats.
        method: Normalization method to apply:
            - "cpm": Counts per million - simple library size normalization
            - "tpm": Transcripts per million - accounts for gene length and library size
            - "rpkm": Reads per kilobase per million - for single-end reads
            - "log2cpm": Log2-transformed CPM with pseudocount of 1
            - "quantile": Quantile normalization across samples
            - "median_ratio": DESeq2-style size factor normalization
        gene_lengths: Gene lengths in base pairs, required for "tpm" and "rpkm".
            Index must match counts_df row index.

    Returns:
        Normalized expression matrix with same dimensions as input.

    Raises:
        ValueError: If gene_lengths is required but not provided, or if
            counts contain negative values, or if method is unknown.

    Example:
        >>> counts = pd.DataFrame({
        ...     "sample1": [100, 200, 50],
        ...     "sample2": [120, 180, 60]
        ... }, index=["gene1", "gene2", "gene3"])
        >>> normalized = normalize_counts(counts, method="cpm")
    """
    if counts_df.empty:
        logger.warning("Empty count matrix provided, returning empty DataFrame")
        return counts_df.copy()

    # Validate counts are non-negative
    if (counts_df.values < 0).any():
        raise ValueError("Count matrix contains negative values")

    # Convert to float for calculations
    counts = counts_df.astype(float)

    if method == "cpm":
        return _normalize_cpm(counts)
    elif method == "tpm":
        if gene_lengths is None:
            raise ValueError("gene_lengths required for TPM normalization")
        return _normalize_tpm(counts, gene_lengths)
    elif method == "rpkm":
        if gene_lengths is None:
            raise ValueError("gene_lengths required for RPKM normalization")
        return _normalize_rpkm(counts, gene_lengths)
    elif method == "log2cpm":
        cpm = _normalize_cpm(counts)
        return np.log2(cpm + 1)
    elif method == "quantile":
        return _normalize_quantile(counts)
    elif method == "median_ratio":
        return _normalize_median_ratio(counts)
    else:
        raise ValueError(
            f"Unknown normalization method: {method}. "
            f"Valid methods: cpm, tpm, rpkm, log2cpm, quantile, median_ratio"
        )


def _normalize_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    """Normalize counts to counts per million (CPM).

    Args:
        counts: Raw count matrix (genes x samples).

    Returns:
        CPM-normalized expression matrix.
    """
    library_sizes = counts.sum(axis=0)
    # Avoid division by zero
    library_sizes = library_sizes.replace(0, 1)
    cpm = counts.div(library_sizes, axis=1) * 1e6
    return cpm


def _normalize_tpm(counts: pd.DataFrame, gene_lengths: pd.Series) -> pd.DataFrame:
    """Normalize counts to transcripts per million (TPM).

    TPM first normalizes by gene length (reads per kilobase), then by
    library size, ensuring values sum to 1 million per sample.

    Args:
        counts: Raw count matrix (genes x samples).
        gene_lengths: Gene lengths in base pairs.

    Returns:
        TPM-normalized expression matrix.
    """
    # Align gene lengths with counts index
    lengths = gene_lengths.reindex(counts.index)
    if lengths.isna().any():
        missing = lengths.isna().sum()
        logger.warning(f"{missing} genes missing length information, using median length")
        median_length = lengths.median()
        if pd.isna(median_length):
            median_length = 1000  # Default fallback
        lengths = lengths.fillna(median_length)

    # Avoid division by zero
    lengths = lengths.replace(0, 1)

    # Calculate reads per kilobase (RPK)
    rpk = counts.div(lengths / 1000, axis=0)

    # Normalize to per million
    rpk_sum = rpk.sum(axis=0)
    rpk_sum = rpk_sum.replace(0, 1)
    tpm = rpk.div(rpk_sum, axis=1) * 1e6

    return tpm


def _normalize_rpkm(counts: pd.DataFrame, gene_lengths: pd.Series) -> pd.DataFrame:
    """Normalize counts to RPKM (reads per kilobase per million).

    Args:
        counts: Raw count matrix (genes x samples).
        gene_lengths: Gene lengths in base pairs.

    Returns:
        RPKM-normalized expression matrix.
    """
    # Align gene lengths with counts index
    lengths = gene_lengths.reindex(counts.index)
    if lengths.isna().any():
        missing = lengths.isna().sum()
        logger.warning(f"{missing} genes missing length information, using median length")
        median_length = lengths.median()
        if pd.isna(median_length):
            median_length = 1000
        lengths = lengths.fillna(median_length)

    # Avoid division by zero
    lengths = lengths.replace(0, 1)

    # Calculate library sizes
    library_sizes = counts.sum(axis=0)
    library_sizes = library_sizes.replace(0, 1)

    # RPKM = (reads * 10^9) / (gene_length * total_reads)
    rpkm = counts.div(lengths / 1000, axis=0).div(library_sizes, axis=1) * 1e6

    return rpkm


def _normalize_quantile(counts: pd.DataFrame) -> pd.DataFrame:
    """Apply quantile normalization across samples.

    Forces samples to have identical distributions by replacing values
    with the mean of values at each rank position.

    Args:
        counts: Raw count matrix (genes x samples).

    Returns:
        Quantile-normalized expression matrix.
    """
    if counts.empty:
        return counts.copy()

    # Rank values within each sample
    ranked = counts.rank(axis=0, method="average")

    # Sort each column and compute mean across samples at each position
    sorted_counts = np.sort(counts.values, axis=0)
    mean_distribution = sorted_counts.mean(axis=1)

    # Create mapping from rank to normalized value
    n_genes = len(counts)
    rank_to_value = pd.Series(mean_distribution, index=np.arange(1, n_genes + 1))

    # Apply normalization
    normalized = pd.DataFrame(index=counts.index, columns=counts.columns, dtype=float)

    for col in counts.columns:
        # Map ranks to normalized values (handle ties by averaging)
        col_ranks = ranked[col]
        # For fractional ranks (from ties), interpolate
        floor_ranks = np.floor(col_ranks).astype(int).clip(1, n_genes)
        ceil_ranks = np.ceil(col_ranks).astype(int).clip(1, n_genes)
        frac = col_ranks - floor_ranks

        floor_values = rank_to_value.iloc[floor_ranks - 1].values
        ceil_values = rank_to_value.iloc[ceil_ranks - 1].values

        normalized[col] = floor_values * (1 - frac) + ceil_values * frac

    return normalized


def _normalize_median_ratio(counts: pd.DataFrame) -> pd.DataFrame:
    """Apply DESeq2-style median ratio normalization.

    Estimates size factors using the median of ratios method and
    divides counts by these factors.

    Args:
        counts: Raw count matrix (genes x samples).

    Returns:
        Size factor-normalized expression matrix.
    """
    size_factors = estimate_size_factors(counts)
    # Avoid division by zero
    size_factors = size_factors.replace(0, 1)
    normalized = counts.div(size_factors, axis=1)
    return normalized


def estimate_size_factors(counts_df: pd.DataFrame) -> pd.Series:
    """Estimate size factors using DESeq2's median-of-ratios method.

    Computes normalization factors based on the geometric mean of
    expression across samples, following the DESeq2 methodology.

    Args:
        counts_df: Raw count matrix with genes as rows and samples as columns.

    Returns:
        Size factors for each sample. Values > 1 indicate samples with
        higher sequencing depth than average.

    Raises:
        ValueError: If counts contain negative values.

    Example:
        >>> counts = pd.DataFrame({
        ...     "sample1": [100, 200, 50],
        ...     "sample2": [200, 400, 100]
        ... })
        >>> sf = estimate_size_factors(counts)
        >>> # sample2 has ~2x coverage, so sf["sample2"] ~ 2
    """
    if counts_df.empty:
        return pd.Series(dtype=float)

    if (counts_df.values < 0).any():
        raise ValueError("Count matrix contains negative values")

    counts = counts_df.astype(float)

    # Calculate geometric mean of each gene (excluding zeros)
    # Use log-transform for numerical stability
    with np.errstate(divide="ignore"):
        log_counts = np.log(counts)
        log_counts = log_counts.replace(-np.inf, np.nan)

    # Geometric mean per gene (across samples)
    log_geo_mean = log_counts.mean(axis=1)

    # Filter genes with valid geometric means (non-zero in all samples)
    valid_genes = ~log_geo_mean.isna()
    if valid_genes.sum() == 0:
        logger.warning("No genes with non-zero counts in all samples, using library size normalization")
        library_sizes = counts.sum(axis=0)
        median_lib = library_sizes.median()
        if median_lib == 0:
            return pd.Series(1.0, index=counts.columns)
        return library_sizes / median_lib

    # Calculate log-ratios to geometric mean for valid genes
    log_ratios = log_counts.loc[valid_genes].sub(log_geo_mean[valid_genes], axis=0)

    # Size factor is median of ratios (in exp space)
    size_factors = np.exp(log_ratios.median(axis=0))

    return size_factors


# =============================================================================
# Gene Filtering
# =============================================================================


def filter_low_expression(
    counts_df: pd.DataFrame,
    min_count: int = 10,
    min_samples: int = 2,
) -> pd.DataFrame:
    """Filter out lowly expressed genes.

    Removes genes that don't meet minimum expression thresholds across
    samples, reducing noise and multiple testing burden.

    Args:
        counts_df: Raw count matrix with genes as rows and samples as columns.
        min_count: Minimum count threshold for a gene to be considered expressed.
        min_samples: Minimum number of samples where gene must be expressed.

    Returns:
        Filtered count matrix with lowly expressed genes removed.

    Example:
        >>> counts = pd.DataFrame({
        ...     "s1": [100, 5, 200], "s2": [90, 3, 180],
        ...     "s3": [110, 1, 190]
        ... }, index=["gene1", "gene2", "gene3"])
        >>> filtered = filter_low_expression(counts, min_count=10, min_samples=2)
        >>> # gene2 removed (never reaches 10 counts)
    """
    if counts_df.empty:
        return counts_df.copy()

    # Count samples where each gene exceeds threshold
    expressed_samples = (counts_df >= min_count).sum(axis=1)

    # Keep genes expressed in at least min_samples
    keep_genes = expressed_samples >= min_samples

    n_removed = (~keep_genes).sum()
    logger.info(
        f"Removed {n_removed}/{len(counts_df)} genes with low expression "
        f"(< {min_count} counts in < {min_samples} samples)"
    )

    return counts_df.loc[keep_genes]


def get_highly_variable_genes(
    expression_df: pd.DataFrame,
    n_top: int = 2000,
    method: VarianceMethod = "cv",
) -> List[str]:
    """Identify highly variable genes.

    Selects genes with highest variability across samples, useful for
    dimensionality reduction and clustering.

    Args:
        expression_df: Expression matrix with genes as rows and samples as columns.
            Should be normalized (e.g., log-transformed CPM).
        n_top: Number of top variable genes to return.
        method: Variability measure:
            - "cv": Coefficient of variation (sd/mean)
            - "variance": Raw variance

    Returns:
        List of gene names sorted by variability (most variable first).

    Example:
        >>> normalized = normalize_counts(counts, method="log2cpm")
        >>> hvg = get_highly_variable_genes(normalized, n_top=1000)
    """
    if expression_df.empty:
        return []

    if method == "cv":
        # Coefficient of variation: std / mean
        means = expression_df.mean(axis=1)
        stds = expression_df.std(axis=1)

        # Avoid division by zero
        means_safe = means.replace(0, 1e-10)
        variability = stds / means_safe

    elif method == "variance":
        variability = expression_df.var(axis=1)

    else:
        raise ValueError(f"Unknown variance method: {method}. Valid: cv, variance")

    # Sort by variability and select top genes
    sorted_genes = variability.sort_values(ascending=False)

    n_top = min(n_top, len(sorted_genes))
    top_genes = sorted_genes.head(n_top).index.tolist()

    logger.info(f"Selected {len(top_genes)} highly variable genes using {method}")

    return top_genes
