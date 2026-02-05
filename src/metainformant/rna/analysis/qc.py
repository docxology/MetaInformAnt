"""RNA-seq quality control module.

Comprehensive quality control metrics for RNA-seq count data including
sample-level and gene-level statistics, library complexity analysis,
batch effect detection, and bias assessment.

This module provides REAL implementations using numpy, scipy, and pandas.
No mocking, no placeholder data.
"""

from __future__ import annotations

from typing import Any, Literal, Optional, Sequence

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import pdist, squareform

from metainformant.core import logging

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


# =============================================================================
# Batch Effect Detection
# =============================================================================


def detect_batch_effects(
    expression_df: pd.DataFrame,
    batch_labels: pd.Series,
    method: Literal["kruskal", "silhouette", "pvca"] = "kruskal",
) -> dict[str, Any]:
    """Detect batch effects in expression data.

    Tests for systematic differences between batches that may confound
    biological signal.

    Args:
        expression_df: DataFrame with genes as rows and samples as columns.
            Should be normalized expression values (e.g., log-transformed).
        batch_labels: Series with sample names as index and batch IDs as values.
            Must contain labels for all samples in expression_df.
        method: Detection method:
            - "kruskal": Kruskal-Wallis test per gene, returns fraction of significant genes
            - "silhouette": Silhouette score in PCA space (higher = more separation by batch)
            - "pvca": Principal Variance Component Analysis approximation

    Returns:
        Dictionary with method-specific results:
            - method: The detection method used
            - batch_effect_detected: Boolean indicating if batch effect is significant
            - Additional keys depend on method

        For "kruskal":
            - pct_significant_genes: Percentage of genes with p < 0.05
            - n_significant_genes: Number of significant genes
            - median_pvalue: Median p-value across all genes

        For "silhouette":
            - silhouette_score: Mean silhouette score (-1 to 1)
            - batch_separation: "strong" if > 0.3, "moderate" if > 0.1, else "weak"

        For "pvca":
            - variance_explained_batch: Proportion of variance explained by batch
            - variance_explained_residual: Proportion unexplained

    Raises:
        ValueError: If expression_df is empty, method is invalid, or batch_labels don't match.

    Examples:
        >>> batch_info = pd.Series({'sample1': 'A', 'sample2': 'A', 'sample3': 'B', 'sample4': 'B'})
        >>> result = detect_batch_effects(expression_df, batch_info, method="kruskal")
        >>> if result['batch_effect_detected']:
        ...     print("Warning: batch effects detected!")
    """
    if expression_df.empty:
        raise ValueError("expression_df cannot be empty")

    if method not in ("kruskal", "silhouette", "pvca"):
        raise ValueError(f"Unknown method: {method}. Must be 'kruskal', 'silhouette', or 'pvca'")

    samples = expression_df.columns.tolist()
    missing_labels = [s for s in samples if s not in batch_labels.index]
    if missing_labels:
        raise ValueError(f"Missing batch labels for samples: {missing_labels[:5]}...")

    # Align batch labels with expression data
    batch_aligned = batch_labels.loc[samples]
    unique_batches = batch_aligned.unique()

    if len(unique_batches) < 2:
        return {
            "method": method,
            "batch_effect_detected": False,
            "message": "Only one batch present, cannot detect batch effects",
        }

    if method == "kruskal":
        # Kruskal-Wallis test for each gene
        pvalues = []
        for gene in expression_df.index:
            gene_values = expression_df.loc[gene]
            groups = [gene_values[batch_aligned == batch].values for batch in unique_batches]
            # Filter out empty groups
            groups = [g for g in groups if len(g) > 0]
            if len(groups) >= 2 and all(len(g) >= 1 for g in groups):
                try:
                    stat, pval = stats.kruskal(*groups)
                    pvalues.append(pval)
                except ValueError:
                    # All values might be identical
                    pvalues.append(1.0)

        if not pvalues:
            return {
                "method": method,
                "batch_effect_detected": False,
                "message": "Could not perform Kruskal-Wallis test",
            }

        pvalues = np.array(pvalues)
        n_significant = int(np.sum(pvalues < 0.05))
        pct_significant = 100.0 * n_significant / len(pvalues)

        # Batch effect detected if >10% of genes show significant differences
        batch_detected = pct_significant > 10.0

        return {
            "method": method,
            "batch_effect_detected": batch_detected,
            "pct_significant_genes": pct_significant,
            "n_significant_genes": n_significant,
            "n_genes_tested": len(pvalues),
            "median_pvalue": float(np.median(pvalues)),
        }

    elif method == "silhouette":
        try:
            from sklearn.decomposition import PCA
            from sklearn.metrics import silhouette_score
            from sklearn.preprocessing import LabelEncoder
        except ImportError:
            raise ImportError("sklearn is required for silhouette method. Install with: uv pip install scikit-learn")

        # Transpose: samples as rows
        X = expression_df.T.values.astype(np.float64)
        X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)

        # PCA for dimensionality reduction
        n_samples = X.shape[0]
        n_components = min(n_samples - 1, 10, X.shape[1])
        if n_components < 2:
            return {
                "method": method,
                "batch_effect_detected": False,
                "message": "Not enough samples for PCA",
            }

        pca = PCA(n_components=n_components)
        X_pca = pca.fit_transform(X)

        # Encode batch labels
        le = LabelEncoder()
        batch_encoded = le.fit_transform(batch_aligned.values)

        # Calculate silhouette score
        try:
            sil_score = float(silhouette_score(X_pca, batch_encoded))
        except ValueError:
            sil_score = 0.0

        if sil_score > 0.3:
            separation = "strong"
            batch_detected = True
        elif sil_score > 0.1:
            separation = "moderate"
            batch_detected = True
        else:
            separation = "weak"
            batch_detected = False

        return {
            "method": method,
            "batch_effect_detected": batch_detected,
            "silhouette_score": sil_score,
            "batch_separation": separation,
            "n_components_used": n_components,
            "variance_explained": float(np.sum(pca.explained_variance_ratio_)),
        }

    elif method == "pvca":
        # Principal Variance Component Analysis approximation
        # Simplified version: compute variance explained by batch vs residual

        # Transpose: samples as rows
        X = expression_df.T.values.astype(np.float64)
        X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)

        # Total variance
        total_var = float(np.var(X))
        if total_var == 0:
            return {
                "method": method,
                "batch_effect_detected": False,
                "message": "No variance in expression data",
            }

        # Compute batch means
        batch_means = {}
        for batch in unique_batches:
            mask = batch_aligned == batch
            batch_means[batch] = X[mask.values].mean(axis=0)

        # Between-batch variance: variance of batch means weighted by batch size
        batch_sizes = batch_aligned.value_counts()
        weighted_batch_means = np.zeros(X.shape[1])
        for batch in unique_batches:
            weighted_batch_means += (batch_sizes[batch] / len(samples)) * batch_means[batch]

        between_var = 0.0
        for batch in unique_batches:
            weight = batch_sizes[batch] / len(samples)
            between_var += weight * np.mean((batch_means[batch] - weighted_batch_means) ** 2)

        variance_batch = between_var / total_var if total_var > 0 else 0.0
        variance_residual = 1.0 - variance_batch

        # Batch effect detected if batch explains >10% of variance
        batch_detected = variance_batch > 0.1

        return {
            "method": method,
            "batch_effect_detected": batch_detected,
            "variance_explained_batch": float(variance_batch),
            "variance_explained_residual": float(variance_residual),
            "n_batches": len(unique_batches),
            "batch_sizes": {str(k): int(v) for k, v in batch_sizes.items()},
        }

    # Should not reach here
    return {"method": method, "batch_effect_detected": False}


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


# =============================================================================
# GC and Length Bias
# =============================================================================


def detect_gc_bias(
    expression_df: pd.DataFrame,
    gc_content: pd.Series,
) -> dict[str, Any]:
    """Detect GC content bias in expression data.

    Tests whether expression levels correlate with gene GC content, which
    can indicate technical bias from library preparation or sequencing.

    Args:
        expression_df: DataFrame with genes as rows and samples as columns.
            Should be normalized expression values.
        gc_content: Series with gene names as index and GC content (0-1) as values.
            Must contain values for genes in expression_df.

    Returns:
        Dictionary with:
            - overall_correlation: Mean correlation across samples
            - overall_pvalue: P-value for overall correlation test
            - bias_detected: True if significant correlation detected
            - per_sample: Dict mapping sample names to correlation results:
                - correlation: Pearson correlation coefficient
                - pvalue: P-value for correlation test

    Raises:
        ValueError: If expression_df is empty or gc_content doesn't match genes.

    Examples:
        >>> gc_bias = detect_gc_bias(expression_df, gene_gc_content)
        >>> if gc_bias['bias_detected']:
        ...     print("GC bias detected, consider correction")
    """
    if expression_df.empty:
        raise ValueError("expression_df cannot be empty")

    genes = expression_df.index.tolist()
    common_genes = [g for g in genes if g in gc_content.index]

    if len(common_genes) < 10:
        return {
            "overall_correlation": 0.0,
            "overall_pvalue": 1.0,
            "bias_detected": False,
            "message": f"Too few genes with GC content data ({len(common_genes)})",
            "per_sample": {},
        }

    # Subset to common genes
    expr_subset = expression_df.loc[common_genes]
    gc_subset = gc_content.loc[common_genes].values.astype(np.float64)

    per_sample_results = {}
    correlations = []

    for sample in expr_subset.columns:
        sample_expr = expr_subset[sample].values.astype(np.float64)

        # Filter out zero expression genes for correlation
        mask = sample_expr > 0
        if np.sum(mask) < 10:
            per_sample_results[sample] = {"correlation": 0.0, "pvalue": 1.0}
            continue

        expr_filtered = sample_expr[mask]
        gc_filtered = gc_subset[mask]

        # Pearson correlation
        corr, pval = stats.pearsonr(expr_filtered, gc_filtered)
        per_sample_results[sample] = {"correlation": float(corr), "pvalue": float(pval)}
        correlations.append(corr)

    if not correlations:
        return {
            "overall_correlation": 0.0,
            "overall_pvalue": 1.0,
            "bias_detected": False,
            "message": "Could not compute correlations",
            "per_sample": per_sample_results,
        }

    overall_corr = float(np.mean(correlations))
    # Test if overall correlation is significantly different from zero
    # Using one-sample t-test on correlations
    if len(correlations) >= 2:
        t_stat, overall_pval = stats.ttest_1samp(correlations, 0)
        overall_pval = float(overall_pval)
    else:
        overall_pval = 1.0

    # Bias detected if significant correlation with |r| > 0.1
    bias_detected = overall_pval < 0.05 and abs(overall_corr) > 0.1

    return {
        "overall_correlation": overall_corr,
        "overall_pvalue": overall_pval,
        "bias_detected": bias_detected,
        "n_genes_tested": len(common_genes),
        "per_sample": per_sample_results,
    }


def detect_length_bias(
    expression_df: pd.DataFrame,
    gene_lengths: pd.Series,
) -> dict[str, Any]:
    """Detect gene length bias in expression data.

    Tests whether expression levels correlate with gene length, which
    can indicate technical bias (longer genes may have higher counts
    simply due to more mapping opportunities).

    Args:
        expression_df: DataFrame with genes as rows and samples as columns.
            Should be normalized expression values.
        gene_lengths: Series with gene names as index and lengths (bp) as values.
            Must contain values for genes in expression_df.

    Returns:
        Dictionary with:
            - overall_correlation: Mean correlation across samples
            - overall_pvalue: P-value for overall correlation test
            - bias_detected: True if significant correlation detected
            - per_sample: Dict mapping sample names to correlation results:
                - correlation: Pearson correlation coefficient
                - pvalue: P-value for correlation test

    Raises:
        ValueError: If expression_df is empty or gene_lengths doesn't match genes.

    Examples:
        >>> length_bias = detect_length_bias(expression_df, transcript_lengths)
        >>> if length_bias['bias_detected']:
        ...     print("Length bias detected, consider TPM/RPKM normalization")
    """
    if expression_df.empty:
        raise ValueError("expression_df cannot be empty")

    genes = expression_df.index.tolist()
    common_genes = [g for g in genes if g in gene_lengths.index]

    if len(common_genes) < 10:
        return {
            "overall_correlation": 0.0,
            "overall_pvalue": 1.0,
            "bias_detected": False,
            "message": f"Too few genes with length data ({len(common_genes)})",
            "per_sample": {},
        }

    # Subset to common genes
    expr_subset = expression_df.loc[common_genes]
    length_subset = gene_lengths.loc[common_genes].values.astype(np.float64)

    # Log-transform lengths for better correlation assessment
    log_lengths = np.log10(length_subset + 1)

    per_sample_results = {}
    correlations = []

    for sample in expr_subset.columns:
        sample_expr = expr_subset[sample].values.astype(np.float64)

        # Filter out zero expression genes for correlation
        mask = sample_expr > 0
        if np.sum(mask) < 10:
            per_sample_results[sample] = {"correlation": 0.0, "pvalue": 1.0}
            continue

        expr_filtered = sample_expr[mask]
        log_expr = np.log10(expr_filtered + 1)
        lengths_filtered = log_lengths[mask]

        # Pearson correlation on log-transformed values
        corr, pval = stats.pearsonr(log_expr, lengths_filtered)
        per_sample_results[sample] = {"correlation": float(corr), "pvalue": float(pval)}
        correlations.append(corr)

    if not correlations:
        return {
            "overall_correlation": 0.0,
            "overall_pvalue": 1.0,
            "bias_detected": False,
            "message": "Could not compute correlations",
            "per_sample": per_sample_results,
        }

    overall_corr = float(np.mean(correlations))
    # Test if overall correlation is significantly different from zero
    if len(correlations) >= 2:
        t_stat, overall_pval = stats.ttest_1samp(correlations, 0)
        overall_pval = float(overall_pval)
    else:
        overall_pval = 1.0

    # Bias detected if significant correlation with |r| > 0.1
    bias_detected = overall_pval < 0.05 and abs(overall_corr) > 0.1

    return {
        "overall_correlation": overall_corr,
        "overall_pvalue": overall_pval,
        "bias_detected": bias_detected,
        "n_genes_tested": len(common_genes),
        "per_sample": per_sample_results,
    }


# =============================================================================
# Summary Report
# =============================================================================


def generate_qc_report(
    counts_df: pd.DataFrame,
    sample_metadata: Optional[pd.DataFrame] = None,
    gene_metadata: Optional[pd.DataFrame] = None,
) -> dict[str, Any]:
    """Generate comprehensive QC report for RNA-seq data.

    Combines all QC metrics into a single summary report suitable for
    JSON serialization and downstream analysis.

    Args:
        counts_df: DataFrame with genes as rows and samples as columns.
            Values should be raw counts (non-negative integers).
        sample_metadata: Optional DataFrame with sample annotations.
            If provided and contains 'batch' column, batch effects will be tested.
            Index should be sample names.
        gene_metadata: Optional DataFrame with gene annotations.
            If provided and contains 'gc_content' and/or 'length' columns,
            bias detection will be performed. Index should be gene names.

    Returns:
        Nested dictionary with QC results:
            - summary: High-level statistics
                - n_samples: Number of samples
                - n_genes: Number of genes
                - total_counts: Sum of all counts
                - median_counts_per_sample: Median library size
            - sample_metrics: Per-sample quality metrics (from compute_sample_metrics)
            - gene_metrics: Per-gene quality metrics (from compute_gene_metrics)
            - library_complexity: Library diversity indices (from estimate_library_complexity)
            - outlier_samples: List of detected outlier samples
            - expression_classification: Gene expression level distribution
            - batch_effects: Batch effect detection results (if batch info provided)
            - gc_bias: GC content bias results (if gc_content provided)
            - length_bias: Gene length bias results (if length provided)
            - correlation_stats: Sample correlation statistics
            - warnings: List of QC warnings/issues detected

    Raises:
        ValueError: If counts_df is empty.

    Examples:
        >>> report = generate_qc_report(counts_df, sample_metadata=sample_info)
        >>> with open('qc_report.json', 'w') as f:
        ...     json.dump(report, f, indent=2)
    """
    if counts_df.empty:
        raise ValueError("counts_df cannot be empty")

    logger.info(f"Generating QC report for {len(counts_df.columns)} samples and {len(counts_df)} genes")

    warnings: list[str] = []
    report: dict[str, Any] = {}

    # =================
    # Summary statistics
    # =================
    total_counts = float(counts_df.values.sum())
    counts_per_sample = counts_df.sum(axis=0)

    report["summary"] = {
        "n_samples": len(counts_df.columns),
        "n_genes": len(counts_df),
        "total_counts": total_counts,
        "median_counts_per_sample": float(counts_per_sample.median()),
        "mean_counts_per_sample": float(counts_per_sample.mean()),
        "min_counts_per_sample": float(counts_per_sample.min()),
        "max_counts_per_sample": float(counts_per_sample.max()),
    }

    # Check for low-count samples
    low_count_threshold = 1_000_000  # 1M reads typically minimum for RNA-seq
    low_count_samples = counts_per_sample[counts_per_sample < low_count_threshold].index.tolist()
    if low_count_samples:
        warnings.append(f"{len(low_count_samples)} samples have < {low_count_threshold:,} total counts")

    # =================
    # Sample metrics
    # =================
    logger.info("Computing sample quality metrics...")
    sample_metrics = compute_sample_metrics(counts_df)
    report["sample_metrics"] = sample_metrics.to_dict(orient="index")

    # Check for high zero percentage
    high_zero_samples = sample_metrics[sample_metrics["pct_zero"] > 80].index.tolist()
    if high_zero_samples:
        warnings.append(f"{len(high_zero_samples)} samples have >80% zero counts")

    # =================
    # Gene metrics
    # =================
    logger.info("Computing gene quality metrics...")
    gene_metrics = compute_gene_metrics(counts_df)
    # Store summary statistics instead of full per-gene data (can be very large)
    report["gene_metrics"] = {
        "mean_expression": {
            "mean": float(gene_metrics["mean_expression"].mean()),
            "median": float(gene_metrics["mean_expression"].median()),
            "std": float(gene_metrics["mean_expression"].std()),
        },
        "pct_zero": {
            "mean": float(gene_metrics["pct_zero"].mean()),
            "median": float(gene_metrics["pct_zero"].median()),
        },
        "n_genes_detected_in_all": int((gene_metrics["pct_zero"] == 0).sum()),
        "n_genes_detected_in_none": int((gene_metrics["n_samples_detected"] == 0).sum()),
        "n_genes_high_variance": int((gene_metrics["cv"] > 1.0).sum()),
    }

    # Warn about genes detected in no samples
    n_undetected = int((gene_metrics["n_samples_detected"] == 0).sum())
    if n_undetected > 0:
        pct_undetected = 100.0 * n_undetected / len(gene_metrics)
        if pct_undetected > 10:
            warnings.append(f"{n_undetected} genes ({pct_undetected:.1f}%) detected in no samples")

    # =================
    # Library complexity
    # =================
    logger.info("Estimating library complexity...")
    complexity = estimate_library_complexity(counts_df)
    report["library_complexity"] = complexity.to_dict(orient="index")

    # Check for low complexity samples
    low_complexity = complexity[complexity["normalized_entropy"] < 0.5].index.tolist()
    if low_complexity:
        warnings.append(f"{len(low_complexity)} samples have low library complexity (normalized entropy < 0.5)")

    # =================
    # Outlier detection
    # =================
    logger.info("Detecting outlier samples...")
    try:
        outliers_mad = detect_outlier_samples(counts_df, method="mad", threshold=3.0)
    except Exception as e:
        logger.warning(f"MAD outlier detection failed: {e}")
        outliers_mad = []

    report["outlier_samples"] = {
        "method": "mad",
        "threshold": 3.0,
        "outliers": outliers_mad,
        "n_outliers": len(outliers_mad),
    }

    if outliers_mad:
        warnings.append(f"{len(outliers_mad)} outlier samples detected")

    # =================
    # Expression classification
    # =================
    logger.info("Classifying gene expression levels...")
    classification = classify_expression_level(counts_df)
    level_counts = classification["expression_level"].value_counts().to_dict()
    report["expression_classification"] = {
        "low": level_counts.get("low", 0),
        "medium": level_counts.get("medium", 0),
        "high": level_counts.get("high", 0),
    }

    # =================
    # Sample correlation
    # =================
    logger.info("Computing sample correlations...")
    corr_matrix = compute_correlation_matrix(counts_df, method="spearman")
    # Get off-diagonal correlations (copy to avoid read-only array issue)
    corr_values = corr_matrix.to_numpy(copy=True)
    np.fill_diagonal(corr_values, np.nan)
    corr_matrix = pd.DataFrame(corr_values, index=corr_matrix.index, columns=corr_matrix.columns)
    mean_corr = corr_matrix.mean(axis=1)
    report["correlation_stats"] = {
        "method": "spearman",
        "mean_pairwise_correlation": float(np.nanmean(corr_matrix.values)),
        "min_pairwise_correlation": float(np.nanmin(corr_matrix.values)),
        "samples_with_low_correlation": mean_corr[mean_corr < 0.8].index.tolist(),
    }

    low_corr_samples = mean_corr[mean_corr < 0.7].index.tolist()
    if low_corr_samples:
        warnings.append(f"{len(low_corr_samples)} samples have low mean correlation (<0.7) with other samples")

    # =================
    # Batch effects (if batch info provided)
    # =================
    if sample_metadata is not None and "batch" in sample_metadata.columns:
        logger.info("Testing for batch effects...")
        batch_labels = sample_metadata["batch"]
        # Ensure batch_labels is a Series with sample names as index
        if not isinstance(batch_labels.index, pd.Index) or not all(s in batch_labels.index for s in counts_df.columns):
            # Try to align
            common_samples = [s for s in counts_df.columns if s in batch_labels.index]
            if len(common_samples) >= 2:
                batch_labels = batch_labels.loc[common_samples]
                counts_subset = counts_df[common_samples]
                # Log-transform for batch effect testing
                log_counts = np.log1p(counts_subset)
                batch_result = detect_batch_effects(log_counts, batch_labels, method="kruskal")
                report["batch_effects"] = batch_result
                if batch_result.get("batch_effect_detected", False):
                    warnings.append("Significant batch effects detected")
            else:
                report["batch_effects"] = {"message": "Not enough samples with batch labels"}
        else:
            log_counts = np.log1p(counts_df)
            batch_result = detect_batch_effects(log_counts, batch_labels, method="kruskal")
            report["batch_effects"] = batch_result
            if batch_result.get("batch_effect_detected", False):
                warnings.append("Significant batch effects detected")
    else:
        report["batch_effects"] = {"message": "No batch information provided"}

    # =================
    # GC bias (if gene metadata provided)
    # =================
    if gene_metadata is not None and "gc_content" in gene_metadata.columns:
        logger.info("Testing for GC content bias...")
        gc_content = gene_metadata["gc_content"]
        log_counts = np.log1p(counts_df)
        gc_result = detect_gc_bias(log_counts, gc_content)
        # Simplify per_sample for report
        gc_summary = {k: v for k, v in gc_result.items() if k != "per_sample"}
        gc_summary["per_sample_summary"] = {
            "mean_correlation": float(np.mean([v["correlation"] for v in gc_result.get("per_sample", {}).values()])),
            "n_samples_with_bias": sum(1 for v in gc_result.get("per_sample", {}).values() if v["pvalue"] < 0.05),
        }
        report["gc_bias"] = gc_summary
        if gc_result.get("bias_detected", False):
            warnings.append("GC content bias detected")
    else:
        report["gc_bias"] = {"message": "No GC content information provided"}

    # =================
    # Length bias (if gene metadata provided)
    # =================
    if gene_metadata is not None and "length" in gene_metadata.columns:
        logger.info("Testing for gene length bias...")
        gene_lengths = gene_metadata["length"]
        log_counts = np.log1p(counts_df)
        length_result = detect_length_bias(log_counts, gene_lengths)
        # Simplify per_sample for report
        length_summary = {k: v for k, v in length_result.items() if k != "per_sample"}
        length_summary["per_sample_summary"] = {
            "mean_correlation": float(
                np.mean([v["correlation"] for v in length_result.get("per_sample", {}).values()])
            ),
            "n_samples_with_bias": sum(1 for v in length_result.get("per_sample", {}).values() if v["pvalue"] < 0.05),
        }
        report["length_bias"] = length_summary
        if length_result.get("bias_detected", False):
            warnings.append("Gene length bias detected")
    else:
        report["length_bias"] = {"message": "No gene length information provided"}

    # =================
    # Warnings summary
    # =================
    report["warnings"] = warnings
    report["n_warnings"] = len(warnings)

    if warnings:
        logger.warning(f"QC report generated with {len(warnings)} warnings")
        for w in warnings:
            logger.warning(f"  - {w}")
    else:
        logger.info("QC report generated with no warnings")

    return report
