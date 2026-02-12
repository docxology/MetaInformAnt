"""Differential expression analysis, PCA, and visualization data preparation.

This module provides tools for performing differential expression analysis
between conditions, multiple testing correction, PCA on expression data,
computing sample distance matrices, and preparing data for volcano and
MA plots.

All implementations are pure Python using numpy, scipy, and pandas.
"""

from __future__ import annotations

from typing import Dict, List, Literal, Optional, Tuple, Union

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import gammaln

from metainformant.core.utils import logging

from .expression_core import estimate_size_factors

logger = logging.get_logger(__name__)


# =============================================================================
# Type Definitions
# =============================================================================

DEMethod = Literal["deseq2_like", "ttest", "wilcoxon"]
PValueMethod = Literal["bh", "bonferroni", "fdr"]
DistanceMethod = Literal["euclidean", "correlation", "cosine"]


# =============================================================================
# Differential Expression Analysis
# =============================================================================


def differential_expression(
    counts_df: pd.DataFrame,
    conditions: Union[List[str], pd.Series],
    method: DEMethod = "deseq2_like",
    reference: Optional[str] = None,
    **kwargs,
) -> pd.DataFrame:
    """Perform differential expression analysis between conditions.

    Compares gene expression between two conditions using the specified
    statistical method. Returns log2 fold changes, p-values, and adjusted
    p-values for each gene.

    Args:
        counts_df: Raw count matrix with genes as rows and samples as columns.
        conditions: Condition labels for each sample, matching column order.
            Must contain exactly two unique conditions.
        method: Statistical method for DE analysis:
            - "deseq2_like": Negative binomial model with Wald test
            - "ttest": Student's t-test on normalized counts
            - "wilcoxon": Wilcoxon rank-sum (Mann-Whitney U) test
        reference: Reference condition for fold change calculation.
            If None, uses the first condition alphabetically.
        **kwargs: Additional method-specific parameters:
            - min_count (int): Minimum total count for gene inclusion (default: 10)
            - pvalue_method (str): P-value adjustment method (default: "bh")

    Returns:
        DataFrame with columns:
            - gene: Gene identifier
            - log2_fold_change: Log2 fold change (treatment vs reference)
            - p_value: Raw p-value from statistical test
            - adjusted_p_value: Multiple testing-adjusted p-value
            - base_mean: Mean expression across all samples
            - stat: Test statistic (t-stat, Wald stat, or U stat)

    Raises:
        ValueError: If conditions has != 2 unique values, or samples don't
            match count columns.

    Example:
        >>> counts = pd.DataFrame({
        ...     "ctrl1": [100, 50], "ctrl2": [110, 45],
        ...     "treat1": [200, 25], "treat2": [180, 30]
        ... }, index=["gene1", "gene2"])
        >>> conditions = ["control", "control", "treatment", "treatment"]
        >>> de = differential_expression(counts, conditions)
    """
    if counts_df.empty:
        logger.warning("Empty count matrix provided")
        return pd.DataFrame(columns=["gene", "log2_fold_change", "p_value", "adjusted_p_value", "base_mean", "stat"])

    # Convert conditions to series
    if isinstance(conditions, list):
        conditions = pd.Series(conditions, index=counts_df.columns)

    # Validate conditions
    unique_conditions = conditions.unique()
    if len(unique_conditions) != 2:
        raise ValueError(f"Expected exactly 2 conditions, got {len(unique_conditions)}: {unique_conditions}")

    if len(conditions) != len(counts_df.columns):
        raise ValueError(f"Conditions length ({len(conditions)}) doesn't match samples ({len(counts_df.columns)})")

    # Determine reference and treatment conditions
    if reference is None:
        reference = sorted(unique_conditions)[0]

    treatment = [c for c in unique_conditions if c != reference][0]
    logger.info(f"Comparing {treatment} vs {reference} (reference)")

    # Get sample indices for each condition
    ref_samples = conditions[conditions == reference].index.tolist()
    treat_samples = conditions[conditions == treatment].index.tolist()

    if len(ref_samples) < 2 or len(treat_samples) < 2:
        logger.warning(f"Small sample sizes: {len(ref_samples)} reference, {len(treat_samples)} treatment")

    # Filter low-expression genes
    min_count = kwargs.get("min_count", 10)
    gene_totals = counts_df.sum(axis=1)
    valid_genes = gene_totals >= min_count
    filtered_counts = counts_df.loc[valid_genes]

    logger.info(f"Analyzing {valid_genes.sum()}/{len(counts_df)} genes (min_count={min_count})")

    # Run differential expression
    if method == "deseq2_like":
        results = _de_deseq2_like(filtered_counts, ref_samples, treat_samples)
    elif method == "ttest":
        results = _de_ttest(filtered_counts, ref_samples, treat_samples)
    elif method == "wilcoxon":
        results = _de_wilcoxon(filtered_counts, ref_samples, treat_samples)
    else:
        raise ValueError(f"Unknown DE method: {method}. Valid methods: deseq2_like, ttest, wilcoxon")

    # Adjust p-values
    pvalue_method = kwargs.get("pvalue_method", "bh")
    results["adjusted_p_value"] = adjust_pvalues(results["p_value"].values, method=pvalue_method)

    # Sort by adjusted p-value
    results = results.sort_values("adjusted_p_value")

    return results


def _de_deseq2_like(
    counts: pd.DataFrame,
    ref_samples: List[str],
    treat_samples: List[str],
) -> pd.DataFrame:
    """Perform DESeq2-like analysis using negative binomial model.

    Args:
        counts: Filtered count matrix.
        ref_samples: Reference condition sample names.
        treat_samples: Treatment condition sample names.

    Returns:
        DataFrame with DE results (without adjusted p-values).
    """
    results = []

    # Normalize counts for fold change calculation
    size_factors = estimate_size_factors(counts)
    normalized = counts.div(size_factors, axis=1)

    for gene in counts.index:
        ref_counts = counts.loc[gene, ref_samples].values.astype(float)
        treat_counts = counts.loc[gene, treat_samples].values.astype(float)

        # Calculate base mean
        base_mean = normalized.loc[gene].mean()

        # Calculate log2 fold change
        ref_mean = ref_counts.mean() + 0.5  # Pseudocount
        treat_mean = treat_counts.mean() + 0.5
        log2fc = np.log2(treat_mean / ref_mean)

        # Negative binomial test
        log2fc_nb, pvalue, _ = _negative_binomial_test(ref_counts, treat_counts)

        # Use the NB-derived fold change if valid
        if not np.isnan(log2fc_nb):
            log2fc = log2fc_nb

        # Wald statistic approximation
        # Standard error from negative binomial model
        all_counts = np.concatenate([ref_counts, treat_counts])
        if all_counts.var() > all_counts.mean():
            # Overdispersed - use NB variance estimate
            dispersion = _estimate_dispersion(all_counts)
            variance = all_counts.mean() + dispersion * all_counts.mean() ** 2
        else:
            variance = all_counts.mean()  # Poisson-like

        se = np.sqrt(variance / len(all_counts)) if variance > 0 else 1.0
        wald_stat = log2fc / se if se > 0 else 0.0

        results.append(
            {
                "gene": gene,
                "log2_fold_change": log2fc,
                "p_value": pvalue,
                "base_mean": base_mean,
                "stat": wald_stat,
            }
        )

    return pd.DataFrame(results)


def _de_ttest(
    counts: pd.DataFrame,
    ref_samples: List[str],
    treat_samples: List[str],
) -> pd.DataFrame:
    """Perform differential expression using t-tests.

    Args:
        counts: Filtered count matrix.
        ref_samples: Reference condition sample names.
        treat_samples: Treatment condition sample names.

    Returns:
        DataFrame with DE results.
    """
    # Log-transform for t-test
    log_counts = np.log2(counts + 1)

    results = []
    for gene in counts.index:
        ref_vals = log_counts.loc[gene, ref_samples].values
        treat_vals = log_counts.loc[gene, treat_samples].values

        # Calculate fold change from original counts
        ref_mean = counts.loc[gene, ref_samples].mean() + 0.5
        treat_mean = counts.loc[gene, treat_samples].mean() + 0.5
        log2fc = np.log2(treat_mean / ref_mean)

        # T-test
        t_stat, pvalue = stats.ttest_ind(treat_vals, ref_vals, equal_var=False)

        # Handle edge cases
        if np.isnan(pvalue):
            pvalue = 1.0
            t_stat = 0.0

        results.append(
            {
                "gene": gene,
                "log2_fold_change": log2fc,
                "p_value": pvalue,
                "base_mean": counts.loc[gene].mean(),
                "stat": t_stat,
            }
        )

    return pd.DataFrame(results)


def _de_wilcoxon(
    counts: pd.DataFrame,
    ref_samples: List[str],
    treat_samples: List[str],
) -> pd.DataFrame:
    """Perform differential expression using Wilcoxon rank-sum test.

    Args:
        counts: Filtered count matrix.
        ref_samples: Reference condition sample names.
        treat_samples: Treatment condition sample names.

    Returns:
        DataFrame with DE results.
    """
    results = []

    for gene in counts.index:
        ref_vals = counts.loc[gene, ref_samples].values.astype(float)
        treat_vals = counts.loc[gene, treat_samples].values.astype(float)

        # Calculate fold change
        ref_mean = ref_vals.mean() + 0.5
        treat_mean = treat_vals.mean() + 0.5
        log2fc = np.log2(treat_mean / ref_mean)

        # Wilcoxon rank-sum (Mann-Whitney U) test
        try:
            u_stat, pvalue = stats.mannwhitneyu(treat_vals, ref_vals, alternative="two-sided")
        except ValueError:
            # All values identical
            u_stat = 0.0
            pvalue = 1.0

        if np.isnan(pvalue):
            pvalue = 1.0

        results.append(
            {
                "gene": gene,
                "log2_fold_change": log2fc,
                "p_value": pvalue,
                "base_mean": counts.loc[gene].mean(),
                "stat": u_stat,
            }
        )

    return pd.DataFrame(results)


def _negative_binomial_test(
    counts_a: np.ndarray,
    counts_b: np.ndarray,
) -> Tuple[float, float, float]:
    """Perform per-gene negative binomial test.

    Uses maximum likelihood estimation to fit negative binomial parameters
    and computes a likelihood ratio test for differential expression.

    Args:
        counts_a: Counts from condition A (reference).
        counts_b: Counts from condition B (treatment).

    Returns:
        Tuple of (log2_fold_change, p_value, base_mean).
    """
    counts_a = np.asarray(counts_a, dtype=float)
    counts_b = np.asarray(counts_b, dtype=float)

    # Calculate means with pseudocount
    mean_a = counts_a.mean() + 0.5
    mean_b = counts_b.mean() + 0.5
    log2fc = np.log2(mean_b / mean_a)

    # Combined data for null model
    all_counts = np.concatenate([counts_a, counts_b])
    base_mean = all_counts.mean()

    # Estimate dispersions
    dispersion_pooled = _estimate_dispersion(all_counts)
    dispersion_a = _estimate_dispersion(counts_a)
    dispersion_b = _estimate_dispersion(counts_b)

    # Log-likelihood for negative binomial
    def nb_loglik(counts: np.ndarray, mu: float, dispersion: float) -> float:
        """Compute negative binomial log-likelihood."""
        if mu <= 0 or dispersion <= 0:
            return -np.inf

        r = 1.0 / dispersion  # Size parameter
        loglik = 0.0
        for k in counts:
            if k < 0:
                return -np.inf
            # NB log-likelihood: log(C(k+r-1, k)) + k*log(p) + r*log(1-p)
            # where p = mu/(mu + r)
            p = mu / (mu + r)
            loglik += gammaln(k + r) - gammaln(k + 1) - gammaln(r)
            loglik += k * np.log(p + 1e-10) + r * np.log(1 - p + 1e-10)
        return loglik

    # Null model: same mean for both groups
    ll_null = nb_loglik(all_counts, base_mean, dispersion_pooled)

    # Alternative model: different means
    ll_alt = nb_loglik(counts_a, mean_a, dispersion_a) + nb_loglik(counts_b, mean_b, dispersion_b)

    # Likelihood ratio test (chi-squared with 1 df)
    lr_stat = 2 * (ll_alt - ll_null)

    if lr_stat < 0 or np.isnan(lr_stat):
        # Model fitting issue, fall back to simple test
        pvalue = 1.0
    else:
        pvalue = stats.chi2.sf(lr_stat, df=1)

    return log2fc, pvalue, base_mean


def _estimate_dispersion(counts: np.ndarray) -> float:
    """Estimate negative binomial dispersion parameter.

    Uses method of moments estimation with shrinkage toward a
    reasonable default for small samples.

    Args:
        counts: Array of count values.

    Returns:
        Estimated dispersion parameter (always positive).
    """
    counts = np.asarray(counts, dtype=float)

    if len(counts) < 2:
        return 0.1  # Default dispersion

    mean_val = counts.mean()
    var_val = counts.var(ddof=1)

    if mean_val <= 0:
        return 0.1

    # Method of moments: var = mu + alpha * mu^2
    # alpha = (var - mu) / mu^2
    dispersion = (var_val - mean_val) / (mean_val**2) if mean_val > 0 else 0.1

    # Shrink toward prior and ensure positive
    prior_dispersion = 0.1
    shrinkage = 0.5  # Weight toward prior for small samples
    dispersion = shrinkage * prior_dispersion + (1 - shrinkage) * dispersion

    return max(dispersion, 1e-6)  # Ensure positive


def adjust_pvalues(
    pvalues: np.ndarray,
    method: PValueMethod = "bh",
) -> np.ndarray:
    """Adjust p-values for multiple testing.

    Args:
        pvalues: Array of raw p-values.
        method: Adjustment method:
            - "bh" or "fdr": Benjamini-Hochberg FDR correction
            - "bonferroni": Bonferroni correction

    Returns:
        Array of adjusted p-values, same length as input.

    Raises:
        ValueError: If method is unknown.

    Example:
        >>> pvals = np.array([0.01, 0.04, 0.03, 0.05])
        >>> adj = adjust_pvalues(pvals, method="bh")
    """
    pvalues = np.asarray(pvalues, dtype=float)
    n = len(pvalues)

    if n == 0:
        return np.array([])

    # Handle NaN values
    nan_mask = np.isnan(pvalues)
    valid_pvals = pvalues.copy()
    valid_pvals[nan_mask] = 1.0

    if method in ("bh", "fdr"):
        # Benjamini-Hochberg procedure
        sorted_idx = np.argsort(valid_pvals)
        sorted_pvals = valid_pvals[sorted_idx]

        # Calculate adjusted p-values: p_adj = p * n / rank
        ranks = np.arange(1, n + 1)
        adjusted = sorted_pvals * n / ranks

        # Ensure monotonicity (cumulative minimum from the end)
        adjusted_monotonic = np.minimum.accumulate(adjusted[::-1])[::-1]

        # Cap at 1.0
        adjusted_monotonic = np.clip(adjusted_monotonic, 0, 1)

        # Return to original order
        result = np.empty(n)
        result[sorted_idx] = adjusted_monotonic

    elif method == "bonferroni":
        result = np.clip(valid_pvals * n, 0, 1)

    else:
        raise ValueError(f"Unknown p-value adjustment method: {method}. Valid: bh, fdr, bonferroni")

    # Restore NaN positions
    result[nan_mask] = np.nan

    return result


# =============================================================================
# Dimensionality Reduction
# =============================================================================


def pca_analysis(
    expression_df: pd.DataFrame,
    n_components: int = 2,
    scale: bool = True,
) -> Dict:
    """Perform PCA on expression data.

    Reduces dimensionality of expression data for visualization and
    exploration, returning transformed coordinates, variance explained,
    and gene loadings.

    Args:
        expression_df: Expression matrix with genes as rows and samples as columns.
            Should be normalized (e.g., log-transformed CPM).
        n_components: Number of principal components to compute.
        scale: Whether to standardize features (zero mean, unit variance).

    Returns:
        Dictionary with keys:
            - "transformed": DataFrame of PC coordinates (samples x components)
            - "explained_variance_ratio": Array of variance explained per PC
            - "loadings": DataFrame of gene loadings (genes x components)
            - "components": Principal component vectors (components x genes)

    Raises:
        ValueError: If n_components > min(n_samples, n_genes).

    Example:
        >>> normalized = normalize_counts(counts, method="log2cpm")
        >>> pca_result = pca_analysis(normalized.T, n_components=3)
        >>> pc_coords = pca_result["transformed"]
    """
    if expression_df.empty:
        return {
            "transformed": pd.DataFrame(),
            "explained_variance_ratio": np.array([]),
            "loadings": pd.DataFrame(),
            "components": np.array([]),
        }

    # Transpose so samples are rows (standard PCA input)
    # Input: genes x samples -> samples x genes
    X = expression_df.T.values.astype(float)
    n_samples, n_features = X.shape

    max_components = min(n_samples, n_features)
    if n_components > max_components:
        logger.warning(f"Reducing n_components from {n_components} to {max_components}")
        n_components = max_components

    # Handle missing values
    if np.isnan(X).any():
        logger.warning("Missing values detected, imputing with column means")
        col_means = np.nanmean(X, axis=0)
        nan_idx = np.where(np.isnan(X))
        X[nan_idx] = np.take(col_means, nan_idx[1])

    # Center data
    X_mean = X.mean(axis=0)
    X_centered = X - X_mean

    # Scale if requested
    if scale:
        X_std = X.std(axis=0)
        X_std[X_std == 0] = 1  # Avoid division by zero
        X_scaled = X_centered / X_std
    else:
        X_scaled = X_centered

    # Compute SVD
    try:
        U, S, Vt = np.linalg.svd(X_scaled, full_matrices=False)
    except np.linalg.LinAlgError:
        logger.error("SVD did not converge")
        return {
            "transformed": pd.DataFrame(),
            "explained_variance_ratio": np.array([]),
            "loadings": pd.DataFrame(),
            "components": np.array([]),
        }

    # Select components
    U = U[:, :n_components]
    S = S[:n_components]
    Vt = Vt[:n_components, :]

    # Transformed coordinates (PC scores)
    transformed = U * S

    # Explained variance ratio
    total_variance = (X_scaled**2).sum()
    explained_variance = S**2 / (n_samples - 1)
    explained_variance_ratio = explained_variance / (total_variance / (n_samples - 1))

    # Gene loadings (correlation of genes with PCs)
    components = Vt  # Principal directions
    loadings = Vt.T * S / np.sqrt(n_samples - 1)  # Scaled loadings

    # Create DataFrames
    pc_names = [f"PC{i+1}" for i in range(n_components)]

    transformed_df = pd.DataFrame(
        transformed,
        index=expression_df.columns,  # Sample names
        columns=pc_names,
    )

    loadings_df = pd.DataFrame(
        loadings,
        index=expression_df.index,  # Gene names
        columns=pc_names,
    )

    return {
        "transformed": transformed_df,
        "explained_variance_ratio": explained_variance_ratio,
        "loadings": loadings_df,
        "components": components,
    }


def compute_sample_distances(
    expression_df: pd.DataFrame,
    method: DistanceMethod = "euclidean",
) -> pd.DataFrame:
    """Compute pairwise distances between samples.

    Args:
        expression_df: Expression matrix with genes as rows and samples as columns.
            Should be normalized (e.g., log-transformed).
        method: Distance metric:
            - "euclidean": Euclidean distance
            - "correlation": 1 - Pearson correlation
            - "cosine": Cosine distance (1 - cosine similarity)

    Returns:
        Square DataFrame of pairwise distances (samples x samples).

    Example:
        >>> log_counts = normalize_counts(counts, method="log2cpm")
        >>> dist_matrix = compute_sample_distances(log_counts, method="correlation")
    """
    if expression_df.empty:
        return pd.DataFrame()

    # Transpose to samples x genes
    X = expression_df.T.values.astype(float)
    samples = expression_df.columns

    n_samples = X.shape[0]
    distances = np.zeros((n_samples, n_samples))

    if method == "euclidean":
        for i in range(n_samples):
            for j in range(i + 1, n_samples):
                dist = np.sqrt(np.sum((X[i] - X[j]) ** 2))
                distances[i, j] = dist
                distances[j, i] = dist

    elif method == "correlation":
        # 1 - Pearson correlation
        for i in range(n_samples):
            for j in range(i + 1, n_samples):
                corr, _ = stats.pearsonr(X[i], X[j])
                if np.isnan(corr):
                    corr = 0
                dist = 1 - corr
                distances[i, j] = dist
                distances[j, i] = dist

    elif method == "cosine":
        # Cosine distance = 1 - cosine similarity
        for i in range(n_samples):
            for j in range(i + 1, n_samples):
                norm_i = np.linalg.norm(X[i])
                norm_j = np.linalg.norm(X[j])
                if norm_i == 0 or norm_j == 0:
                    dist = 1.0
                else:
                    cos_sim = np.dot(X[i], X[j]) / (norm_i * norm_j)
                    dist = 1 - cos_sim
                distances[i, j] = dist
                distances[j, i] = dist

    else:
        raise ValueError(f"Unknown distance method: {method}. Valid: euclidean, correlation, cosine")

    return pd.DataFrame(distances, index=samples, columns=samples)


# =============================================================================
# Visualization Data Preparation
# =============================================================================


def prepare_volcano_data(
    de_results: pd.DataFrame,
    fc_threshold: float = 1.0,
    pvalue_threshold: float = 0.05,
    use_adjusted: bool = True,
) -> pd.DataFrame:
    """Prepare differential expression results for volcano plot.

    Adds a "regulation" column indicating whether each gene is significantly
    upregulated, downregulated, or not significant.

    Args:
        de_results: DataFrame from differential_expression() with columns:
            gene, log2_fold_change, p_value, adjusted_p_value
        fc_threshold: Minimum absolute log2 fold change for significance.
        pvalue_threshold: Maximum p-value for significance.
        use_adjusted: Whether to use adjusted_p_value (True) or p_value (False).

    Returns:
        DataFrame with additional columns:
            - regulation: "up", "down", or "ns" (not significant)
            - neg_log10_pvalue: -log10(p-value) for y-axis plotting

    Example:
        >>> de = differential_expression(counts, conditions)
        >>> volcano_df = prepare_volcano_data(de, fc_threshold=1.0, pvalue_threshold=0.05)
        >>> # Use for plotting: x=log2_fold_change, y=neg_log10_pvalue, color=regulation
    """
    if de_results.empty:
        result = de_results.copy()
        result["regulation"] = pd.Series(dtype=str)
        result["neg_log10_pvalue"] = pd.Series(dtype=float)
        return result

    result = de_results.copy()

    # Select p-value column
    pval_col = "adjusted_p_value" if use_adjusted else "p_value"
    if pval_col not in result.columns:
        pval_col = "p_value"

    # Calculate -log10(pvalue)
    pvals = result[pval_col].values.astype(float)
    # Handle zeros and very small values
    pvals = np.clip(pvals, 1e-300, 1)
    result["neg_log10_pvalue"] = -np.log10(pvals)

    # Determine regulation status
    log2fc = result["log2_fold_change"].values
    pvals = result[pval_col].values

    regulation = np.full(len(result), "ns", dtype=object)

    # Upregulated: positive fold change, significant
    up_mask = (log2fc >= fc_threshold) & (pvals <= pvalue_threshold)
    regulation[up_mask] = "up"

    # Downregulated: negative fold change, significant
    down_mask = (log2fc <= -fc_threshold) & (pvals <= pvalue_threshold)
    regulation[down_mask] = "down"

    result["regulation"] = regulation

    # Log summary
    n_up = (regulation == "up").sum()
    n_down = (regulation == "down").sum()
    n_ns = (regulation == "ns").sum()
    logger.info(
        f"Volcano plot data: {n_up} up, {n_down} down, {n_ns} not significant "
        f"(|log2FC| >= {fc_threshold}, p <= {pvalue_threshold})"
    )

    return result


def prepare_ma_data(de_results: pd.DataFrame) -> pd.DataFrame:
    """Prepare differential expression results for MA plot.

    An MA plot shows log fold change (M) vs average expression (A),
    useful for identifying expression-dependent biases.

    Args:
        de_results: DataFrame from differential_expression() with columns:
            gene, log2_fold_change, base_mean

    Returns:
        DataFrame with additional columns:
            - A: Average expression (log2(base_mean + 1))
            - M: Log2 fold change (same as log2_fold_change)

    Example:
        >>> de = differential_expression(counts, conditions)
        >>> ma_df = prepare_ma_data(de)
        >>> # Use for plotting: x=A, y=M
    """
    if de_results.empty:
        result = de_results.copy()
        result["A"] = pd.Series(dtype=float)
        result["M"] = pd.Series(dtype=float)
        return result

    result = de_results.copy()

    # A = average expression (log scale)
    base_mean = result["base_mean"].values.astype(float)
    result["A"] = np.log2(base_mean + 1)

    # M = log fold change
    result["M"] = result["log2_fold_change"]

    return result
