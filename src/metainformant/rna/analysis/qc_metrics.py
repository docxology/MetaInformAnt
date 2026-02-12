"""RNA-seq quality control metric computation.

Provides sample-level and gene-level quality metrics, library complexity
estimation, and saturation curve analysis for RNA-seq count data.

This module provides REAL implementations using numpy, scipy, and pandas.
No mocking, no placeholder data.
"""

from __future__ import annotations

from typing import Any, Literal, Optional, Sequence

import numpy as np
import pandas as pd
from scipy import stats

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


# =============================================================================
# Sample Quality Metrics
# =============================================================================


def compute_sample_metrics(counts_df: pd.DataFrame) -> pd.DataFrame:
    """Compute per-sample quality metrics from count data.

    Calculates comprehensive statistics for each sample (column) in the
    expression count matrix, useful for identifying low-quality samples.

    Args:
        counts_df: DataFrame with genes as rows and samples as columns.
            Values should be raw or normalized counts (non-negative).

    Returns:
        DataFrame with samples as rows and metrics as columns:
            - total_counts: Sum of all counts in the sample
            - detected_genes: Number of genes with count > 0
            - median_expression: Median expression of detected genes
            - mean_expression: Mean expression across all genes
            - pct_zero: Percentage of genes with zero counts
            - cv: Coefficient of variation (std/mean) across all genes

    Raises:
        ValueError: If counts_df is empty or contains negative values.

    Examples:
        >>> counts = pd.DataFrame({
        ...     'sample1': [10, 20, 0, 5],
        ...     'sample2': [15, 25, 3, 8]
        ... }, index=['geneA', 'geneB', 'geneC', 'geneD'])
        >>> metrics = compute_sample_metrics(counts)
        >>> metrics.columns.tolist()
        ['total_counts', 'detected_genes', 'median_expression', 'mean_expression', 'pct_zero', 'cv']
    """
    if counts_df.empty:
        raise ValueError("counts_df cannot be empty")

    if (counts_df < 0).any().any():
        raise ValueError("counts_df cannot contain negative values")

    n_genes = len(counts_df)
    results = []

    for sample in counts_df.columns:
        sample_counts = counts_df[sample].values.astype(np.float64)

        total_counts = float(np.sum(sample_counts))
        detected_genes = int(np.sum(sample_counts > 0))
        detected_values = sample_counts[sample_counts > 0]

        median_expression = float(np.median(detected_values)) if len(detected_values) > 0 else 0.0
        mean_expression = float(np.mean(sample_counts))
        pct_zero = 100.0 * (n_genes - detected_genes) / n_genes if n_genes > 0 else 0.0

        # Coefficient of variation: std/mean (avoid division by zero)
        std_val = float(np.std(sample_counts))
        cv = std_val / mean_expression if mean_expression > 0 else 0.0

        results.append(
            {
                "sample": sample,
                "total_counts": total_counts,
                "detected_genes": detected_genes,
                "median_expression": median_expression,
                "mean_expression": mean_expression,
                "pct_zero": pct_zero,
                "cv": cv,
            }
        )

    result_df = pd.DataFrame(results)
    result_df.set_index("sample", inplace=True)
    return result_df


def detect_outlier_samples(
    counts_df: pd.DataFrame,
    method: Literal["mad", "isolation_forest", "pca_distance"] = "mad",
    threshold: float = 3.0,
) -> list[str]:
    """Identify outlier samples using various detection methods.

    Detects samples that deviate significantly from the population based
    on expression patterns or summary statistics.

    Args:
        counts_df: DataFrame with genes as rows and samples as columns.
        method: Outlier detection method:
            - "mad": Median Absolute Deviation on total counts and detected genes
            - "isolation_forest": Isolation Forest on sample feature vectors
            - "pca_distance": Distance from centroid in PCA space
        threshold: Detection threshold:
            - For "mad": number of MADs from median (default 3.0)
            - For "isolation_forest": contamination proportion (default 0.1 used)
            - For "pca_distance": number of standard deviations (default 3.0)

    Returns:
        List of sample names identified as outliers.

    Raises:
        ValueError: If method is not recognized or counts_df is empty.

    Examples:
        >>> outliers = detect_outlier_samples(counts_df, method="mad", threshold=3.0)
        >>> print(f"Found {len(outliers)} outlier samples")
    """
    if counts_df.empty:
        raise ValueError("counts_df cannot be empty")

    if method not in ("mad", "isolation_forest", "pca_distance"):
        raise ValueError(f"Unknown method: {method}. Must be 'mad', 'isolation_forest', or 'pca_distance'")

    samples = counts_df.columns.tolist()
    outliers: list[str] = []

    if method == "mad":
        # MAD-based outlier detection using multiple metrics
        metrics = compute_sample_metrics(counts_df)
        total_counts = metrics["total_counts"].values
        detected_genes = metrics["detected_genes"].values

        # MAD calculation for total counts
        median_tc = np.median(total_counts)
        mad_tc = np.median(np.abs(total_counts - median_tc))
        # Scale factor for consistency with normal distribution
        mad_tc_scaled = 1.4826 * mad_tc if mad_tc > 0 else 1.0

        # MAD calculation for detected genes
        median_dg = np.median(detected_genes)
        mad_dg = np.median(np.abs(detected_genes - median_dg))
        mad_dg_scaled = 1.4826 * mad_dg if mad_dg > 0 else 1.0

        for i, sample in enumerate(samples):
            tc_zscore = abs(total_counts[i] - median_tc) / mad_tc_scaled if mad_tc_scaled > 0 else 0
            dg_zscore = abs(detected_genes[i] - median_dg) / mad_dg_scaled if mad_dg_scaled > 0 else 0

            # Sample is outlier if either metric exceeds threshold
            if tc_zscore > threshold or dg_zscore > threshold:
                outliers.append(sample)

    elif method == "isolation_forest":
        # Isolation Forest requires sklearn
        try:
            from sklearn.ensemble import IsolationForest
        except ImportError:
            raise ImportError(
                "sklearn is required for isolation_forest method. Install with: uv pip install scikit-learn"
            )

        # Prepare feature matrix: samples as rows, gene expressions as features
        # For efficiency with many genes, use log-transformed summary stats
        metrics = compute_sample_metrics(counts_df)
        feature_cols = ["total_counts", "detected_genes", "median_expression", "mean_expression", "pct_zero", "cv"]
        X = metrics[feature_cols].values

        # Handle potential inf/nan values
        X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)

        # Standardize features
        X_mean = np.mean(X, axis=0)
        X_std = np.std(X, axis=0)
        X_std[X_std == 0] = 1.0  # Avoid division by zero
        X_scaled = (X - X_mean) / X_std

        # Fit isolation forest
        # contamination is the expected proportion of outliers
        contamination = min(threshold, 0.5) if threshold < 1.0 else 0.1
        clf = IsolationForest(contamination=contamination, random_state=42, n_estimators=100)
        predictions = clf.fit_predict(X_scaled)

        # Outliers are marked as -1
        for i, sample in enumerate(samples):
            if predictions[i] == -1:
                outliers.append(sample)

    elif method == "pca_distance":
        # PCA-based outlier detection
        try:
            from sklearn.decomposition import PCA
        except ImportError:
            raise ImportError("sklearn is required for pca_distance method. Install with: uv pip install scikit-learn")

        # Transpose: samples as rows, genes as columns
        X = counts_df.T.values.astype(np.float64)

        # Log-transform for better distribution (add pseudocount)
        X_log = np.log1p(X)

        # Handle potential inf/nan
        X_log = np.nan_to_num(X_log, nan=0.0, posinf=0.0, neginf=0.0)

        # Standardize
        X_mean = np.mean(X_log, axis=0)
        X_std = np.std(X_log, axis=0)
        X_std[X_std == 0] = 1.0
        X_scaled = (X_log - X_mean) / X_std

        # Fit PCA - use min of n_samples-1 and n_features for n_components
        n_samples, n_features = X_scaled.shape
        n_components = min(min(n_samples, n_features) - 1, 10)  # Cap at 10 components
        n_components = max(n_components, 2)  # At least 2 components

        if n_samples < 3:
            logger.warning("Too few samples for reliable PCA-based outlier detection")
            return []

        pca = PCA(n_components=n_components)
        X_pca = pca.fit_transform(X_scaled)

        # Calculate distance from centroid
        centroid = np.mean(X_pca, axis=0)
        distances = np.sqrt(np.sum((X_pca - centroid) ** 2, axis=1))

        # Outliers based on z-score of distances
        dist_mean = np.mean(distances)
        dist_std = np.std(distances)
        if dist_std > 0:
            dist_zscores = (distances - dist_mean) / dist_std
            for i, sample in enumerate(samples):
                if dist_zscores[i] > threshold:
                    outliers.append(sample)

    logger.info(f"Detected {len(outliers)} outlier samples using {method} method")
    return outliers


# =============================================================================
# Gene Quality Metrics
# =============================================================================


def compute_gene_metrics(counts_df: pd.DataFrame) -> pd.DataFrame:
    """Compute per-gene quality metrics from count data.

    Calculates comprehensive statistics for each gene (row) in the
    expression count matrix, useful for gene filtering decisions.

    Args:
        counts_df: DataFrame with genes as rows and samples as columns.
            Values should be raw or normalized counts (non-negative).

    Returns:
        DataFrame with genes as rows and metrics as columns:
            - mean_expression: Mean expression across all samples
            - variance: Variance of expression across samples
            - cv: Coefficient of variation (std/mean)
            - pct_zero: Percentage of samples with zero counts
            - n_samples_detected: Number of samples where gene is detected (count > 0)
            - dispersion_estimate: Estimated dispersion (variance/mean, approximates negative binomial dispersion)

    Raises:
        ValueError: If counts_df is empty or contains negative values.

    Examples:
        >>> gene_metrics = compute_gene_metrics(counts_df)
        >>> lowly_expressed = gene_metrics[gene_metrics['mean_expression'] < 1.0]
    """
    if counts_df.empty:
        raise ValueError("counts_df cannot be empty")

    if (counts_df < 0).any().any():
        raise ValueError("counts_df cannot contain negative values")

    n_samples = len(counts_df.columns)
    results = []

    for gene in counts_df.index:
        gene_counts = counts_df.loc[gene].values.astype(np.float64)

        mean_expr = float(np.mean(gene_counts))
        variance = float(np.var(gene_counts, ddof=1)) if n_samples > 1 else 0.0
        std_val = float(np.std(gene_counts, ddof=1)) if n_samples > 1 else 0.0

        # Coefficient of variation
        cv = std_val / mean_expr if mean_expr > 0 else 0.0

        # Detection statistics
        n_detected = int(np.sum(gene_counts > 0))
        pct_zero = 100.0 * (n_samples - n_detected) / n_samples if n_samples > 0 else 0.0

        # Dispersion estimate (variance/mean) - approximation for negative binomial
        # This is the "index of dispersion" or Fano factor
        dispersion = variance / mean_expr if mean_expr > 0 else 0.0

        results.append(
            {
                "gene": gene,
                "mean_expression": mean_expr,
                "variance": variance,
                "cv": cv,
                "pct_zero": pct_zero,
                "n_samples_detected": n_detected,
                "dispersion_estimate": dispersion,
            }
        )

    result_df = pd.DataFrame(results)
    result_df.set_index("gene", inplace=True)
    return result_df


def classify_expression_level(
    counts_df: pd.DataFrame,
    thresholds: Optional[dict[str, float]] = None,
) -> pd.DataFrame:
    """Classify genes into expression level categories.

    Assigns each gene to a category (low/medium/high) based on mean
    expression across samples.

    Args:
        counts_df: DataFrame with genes as rows and samples as columns.
        thresholds: Dictionary with threshold values:
            - "low_high": boundary between low and medium (default: 1.0)
            - "medium_high": boundary between medium and high (default: 10.0)
            If None, uses default thresholds based on count data distribution.

    Returns:
        DataFrame with genes as rows and columns:
            - mean_expression: Mean expression value
            - expression_level: Category ("low", "medium", "high")
            - log2_mean: Log2(mean_expression + 1)

    Raises:
        ValueError: If counts_df is empty.

    Examples:
        >>> classified = classify_expression_level(counts_df)
        >>> print(classified['expression_level'].value_counts())
    """
    if counts_df.empty:
        raise ValueError("counts_df cannot be empty")

    # Default thresholds if not provided
    if thresholds is None:
        # Use adaptive thresholds based on data distribution
        all_means = counts_df.mean(axis=1)
        non_zero_means = all_means[all_means > 0]

        if len(non_zero_means) > 0:
            q25 = float(np.percentile(non_zero_means, 25))
            q75 = float(np.percentile(non_zero_means, 75))
            thresholds = {"low_high": q25, "medium_high": q75}
        else:
            thresholds = {"low_high": 1.0, "medium_high": 10.0}

    low_high = thresholds.get("low_high", 1.0)
    medium_high = thresholds.get("medium_high", 10.0)

    results = []
    for gene in counts_df.index:
        mean_expr = float(counts_df.loc[gene].mean())
        log2_mean = float(np.log2(mean_expr + 1))

        if mean_expr < low_high:
            level = "low"
        elif mean_expr < medium_high:
            level = "medium"
        else:
            level = "high"

        results.append(
            {
                "gene": gene,
                "mean_expression": mean_expr,
                "expression_level": level,
                "log2_mean": log2_mean,
            }
        )

    result_df = pd.DataFrame(results)
    result_df.set_index("gene", inplace=True)
    return result_df


# =============================================================================
# Library Complexity
# =============================================================================


def estimate_library_complexity(counts_df: pd.DataFrame) -> pd.DataFrame:
    """Estimate library complexity metrics for each sample.

    Calculates diversity indices that measure how evenly reads are
    distributed across genes, indicating library complexity.

    Args:
        counts_df: DataFrame with genes as rows and samples as columns.
            Values should be raw counts (non-negative integers).

    Returns:
        DataFrame with samples as rows and columns:
            - shannon_entropy: Shannon diversity index (higher = more diverse)
            - simpson_diversity: Simpson's diversity index (1 - sum(p_i^2))
            - effective_gene_count: Effective number of genes (exp(shannon_entropy))
            - max_possible_entropy: Maximum possible entropy for this gene count
            - normalized_entropy: Shannon entropy / max possible entropy (0-1 scale)

    Raises:
        ValueError: If counts_df is empty or contains negative values.

    Examples:
        >>> complexity = estimate_library_complexity(counts_df)
        >>> low_complexity = complexity[complexity['normalized_entropy'] < 0.5]
    """
    if counts_df.empty:
        raise ValueError("counts_df cannot be empty")

    if (counts_df < 0).any().any():
        raise ValueError("counts_df cannot contain negative values")

    n_genes = len(counts_df)
    max_entropy = np.log(n_genes) if n_genes > 1 else 1.0

    results = []
    for sample in counts_df.columns:
        sample_counts = counts_df[sample].values.astype(np.float64)
        total = np.sum(sample_counts)

        if total == 0:
            results.append(
                {
                    "sample": sample,
                    "shannon_entropy": 0.0,
                    "simpson_diversity": 0.0,
                    "effective_gene_count": 1.0,
                    "max_possible_entropy": max_entropy,
                    "normalized_entropy": 0.0,
                }
            )
            continue

        # Calculate proportions
        proportions = sample_counts / total
        # Filter out zeros to avoid log(0)
        non_zero_props = proportions[proportions > 0]

        # Shannon entropy: H = -sum(p_i * log(p_i))
        shannon = float(-np.sum(non_zero_props * np.log(non_zero_props)))

        # Simpson's diversity: D = 1 - sum(p_i^2)
        simpson = float(1.0 - np.sum(proportions**2))

        # Effective gene count (Hill number q=1): exp(H)
        effective_genes = float(np.exp(shannon))

        # Normalized entropy
        normalized = shannon / max_entropy if max_entropy > 0 else 0.0

        results.append(
            {
                "sample": sample,
                "shannon_entropy": shannon,
                "simpson_diversity": simpson,
                "effective_gene_count": effective_genes,
                "max_possible_entropy": max_entropy,
                "normalized_entropy": normalized,
            }
        )

    result_df = pd.DataFrame(results)
    result_df.set_index("sample", inplace=True)
    return result_df


def compute_saturation_curve(
    counts_df: pd.DataFrame,
    fractions: Optional[Sequence[float]] = None,
    n_iterations: int = 10,
    random_state: Optional[int] = 42,
) -> pd.DataFrame:
    """Compute saturation curves by subsampling reads.

    Estimates how many genes would be detected at various sequencing depths
    by randomly subsampling counts. Useful for assessing whether samples
    have sufficient sequencing depth.

    Args:
        counts_df: DataFrame with genes as rows and samples as columns.
            Values should be raw counts (non-negative integers).
        fractions: Sequence of fractions (0-1) to subsample.
            Default: [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        n_iterations: Number of random subsampling iterations for variance estimation.
        random_state: Random seed for reproducibility.

    Returns:
        DataFrame with columns:
            - sample: Sample name
            - fraction: Subsampling fraction
            - detected_genes_mean: Mean detected genes across iterations
            - detected_genes_std: Standard deviation across iterations
            - total_counts: Total counts at this fraction

    Raises:
        ValueError: If counts_df is empty or fractions are invalid.

    Examples:
        >>> saturation = compute_saturation_curve(counts_df)
        >>> # Plot saturation curves
        >>> for sample in saturation['sample'].unique():
        ...     data = saturation[saturation['sample'] == sample]
        ...     plt.plot(data['fraction'], data['detected_genes_mean'], label=sample)
    """
    if counts_df.empty:
        raise ValueError("counts_df cannot be empty")

    if fractions is None:
        fractions = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    fractions = list(fractions)
    for f in fractions:
        if not 0 < f <= 1.0:
            raise ValueError(f"Fractions must be in (0, 1], got {f}")

    rng = np.random.default_rng(random_state)
    results = []

    for sample in counts_df.columns:
        sample_counts = counts_df[sample].values.astype(np.int64)
        total_counts = int(np.sum(sample_counts))

        if total_counts == 0:
            for frac in fractions:
                results.append(
                    {
                        "sample": sample,
                        "fraction": frac,
                        "detected_genes_mean": 0.0,
                        "detected_genes_std": 0.0,
                        "total_counts": 0,
                    }
                )
            continue

        # Build read pool (gene indices for each read)
        # For efficiency, we approximate by multinomial sampling
        for frac in fractions:
            subsample_size = int(total_counts * frac)
            detected_counts = []

            for _ in range(n_iterations):
                # Multinomial subsampling
                probabilities = sample_counts / total_counts
                subsampled = rng.multinomial(subsample_size, probabilities)
                detected = int(np.sum(subsampled > 0))
                detected_counts.append(detected)

            results.append(
                {
                    "sample": sample,
                    "fraction": frac,
                    "detected_genes_mean": float(np.mean(detected_counts)),
                    "detected_genes_std": float(np.std(detected_counts)),
                    "total_counts": subsample_size,
                }
            )

    return pd.DataFrame(results)


def compute_correlation_matrix(
    expression_df: pd.DataFrame,
    method: Literal["pearson", "spearman"] = "pearson",
) -> pd.DataFrame:
    """Compute sample-sample correlation matrix.

    Calculates pairwise correlations between all samples based on their
    expression profiles. Useful for identifying similar samples and detecting outliers.

    Args:
        expression_df: DataFrame with genes as rows and samples as columns.
            Should be normalized expression values.
        method: Correlation method:
            - "pearson": Pearson correlation coefficient (linear relationship)
            - "spearman": Spearman rank correlation (monotonic relationship)

    Returns:
        DataFrame with samples as both rows and columns, containing correlation values.
        Values range from -1 (perfect negative correlation) to 1 (perfect positive).

    Raises:
        ValueError: If expression_df is empty or method is invalid.

    Examples:
        >>> corr_matrix = compute_correlation_matrix(expression_df, method="spearman")
        >>> # Find samples with low correlation to others
        >>> mean_corr = corr_matrix.mean(axis=1)
        >>> outliers = mean_corr[mean_corr < 0.8].index.tolist()
    """
    if expression_df.empty:
        raise ValueError("expression_df cannot be empty")

    if method not in ("pearson", "spearman"):
        raise ValueError(f"Unknown method: {method}. Must be 'pearson' or 'spearman'")

    # pandas corr() computes column-column correlation
    if method == "pearson":
        corr_matrix = expression_df.corr(method="pearson")
    else:
        corr_matrix = expression_df.corr(method="spearman")

    return corr_matrix
