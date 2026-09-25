"""Differential expression analysis, PCA, and visualization data preparation.

This module provides tools for performing differential expression analysis
between conditions, multiple testing correction, PCA on expression data,
computing sample distance matrices, and preparing data for volcano and
MA plots.

All implementations are pure Python using numpy, scipy, and pandas.
"""

from __future__ import annotations

import warnings

from typing import Any, Dict, List, Literal, Optional, Tuple, Union

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

DE_RESULT_COLUMNS = [
    "gene",
    "log2_fold_change",
    "p_value",
    "adjusted_p_value",
    "base_mean",
    "stat",
]


def _empty_de_results() -> pd.DataFrame:
    """Return an empty differential-expression result with the public schema."""
    return pd.DataFrame(columns=DE_RESULT_COLUMNS)


def _align_conditions_to_counts(
    counts_df: pd.DataFrame, conditions: Union[List[str], pd.Series]
) -> pd.Series:
    """Validate and align condition labels to count-matrix columns."""
    if len(conditions) != len(counts_df.columns):
        raise ValueError(
            f"Conditions length ({len(conditions)}) doesn't match samples ({len(counts_df.columns)})"
        )

    if isinstance(conditions, list):
        aligned = pd.Series(conditions, index=counts_df.columns)
    elif isinstance(conditions, pd.Series):
        aligned = conditions.copy()
        if aligned.index.equals(counts_df.columns):
            aligned = aligned.loc[counts_df.columns]
        elif isinstance(aligned.index, pd.RangeIndex) and aligned.index.equals(
            pd.RangeIndex(len(counts_df.columns))
        ):
            aligned.index = counts_df.columns
        else:
            raise ValueError(
                "Conditions Series index must match count matrix columns or use a positional RangeIndex"
            )
    else:
        raise TypeError("conditions must be a list or pandas Series")

    if aligned.isna().any():
        raise ValueError("Conditions contain missing values")

    return aligned


# =============================================================================
# Differential Expression Analysis
# =============================================================================


def differential_expression(
    counts_df: pd.DataFrame,
    conditions: Union[List[str], pd.Series],
    method: DEMethod = "deseq2_like",
    reference: Optional[str] = None,
    **kwargs: Any,
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
            - "deseq2_like": Negative binomial likelihood ratio test with
              log(size factor) offsets (a Wald statistic is also reported)
            - "ttest": Welch's t-test on log2 size-factor-normalized counts
            - "wilcoxon": Wilcoxon rank-sum (Mann-Whitney U) test
        reference: Reference condition for fold change calculation.
            If None, uses the first condition alphabetically.
        **kwargs: Additional method-specific parameters:
            - min_count (int): Minimum total count for gene inclusion (default: 10)
            - pvalue_method (str): P-value adjustment method (default: "bh")

    Returns:
        DataFrame with columns:
            - gene: Gene identifier
            - log2_fold_change: Log2 fold change (treatment vs reference),
              computed from size-factor-normalized group means
            - p_value: Raw p-value from statistical test
            - adjusted_p_value: Multiple testing-adjusted p-value
            - base_mean: Size-factor-normalized mean expression across all
              samples
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
        return _empty_de_results()

    conditions = _align_conditions_to_counts(counts_df, conditions)

    # Validate conditions
    unique_conditions = conditions.unique()
    if len(unique_conditions) != 2:
        raise ValueError(
            f"Expected exactly 2 conditions, got {len(unique_conditions)}: {unique_conditions}"
        )

    # Determine reference and treatment conditions
    if reference is None:
        reference = sorted(unique_conditions)[0]
    elif reference not in unique_conditions:
        raise ValueError(
            f"Reference condition '{reference}' is not present in conditions: {list(unique_conditions)}"
        )

    treatment = [c for c in unique_conditions if c != reference][0]
    logger.info(f"Comparing {treatment} vs {reference} (reference)")

    # Get sample indices for each condition
    ref_samples = conditions[conditions == reference].index.tolist()
    treat_samples = conditions[conditions == treatment].index.tolist()

    if len(ref_samples) < 2 or len(treat_samples) < 2:
        logger.warning(
            f"Small sample sizes: {len(ref_samples)} reference, {len(treat_samples)} treatment"
        )

    # Filter low-expression genes
    min_count = kwargs.get("min_count", 10)
    gene_totals = counts_df.sum(axis=1)
    valid_genes = gene_totals >= min_count
    filtered_counts = counts_df.loc[valid_genes]

    logger.info(
        f"Analyzing {valid_genes.sum()}/{len(counts_df)} genes (min_count={min_count})"
    )

    if filtered_counts.empty:
        return _empty_de_results()

    # Run differential expression
    if method == "deseq2_like":
        results = _de_deseq2_like(filtered_counts, ref_samples, treat_samples)
    elif method == "ttest":
        results = _de_ttest(filtered_counts, ref_samples, treat_samples)
    elif method == "wilcoxon":
        results = _de_wilcoxon(filtered_counts, ref_samples, treat_samples)
    else:
        raise ValueError(
            f"Unknown DE method: {method}. Valid methods: deseq2_like, ttest, wilcoxon"
        )

    if results.empty:
        return _empty_de_results()

    # Adjust p-values
    pvalue_method = kwargs.get("pvalue_method", "bh")
    results["adjusted_p_value"] = adjust_pvalues(
        results["p_value"].values, method=pvalue_method
    )

    # Sort by adjusted p-value
    results = results.sort_values("adjusted_p_value")

    return results


def _de_deseq2_like(
    counts: pd.DataFrame,
    ref_samples: List[str],
    treat_samples: List[str],
) -> pd.DataFrame:
    """Perform DESeq2-like analysis using negative binomial model.

    The negative binomial likelihood is fitted on the raw counts with
    log(size factor) as an offset, so library-depth differences between
    samples are absorbed before testing. Fold changes, base means, and
    the Wald statistic are reported on the size-factor-normalized scale.

    Args:
        counts: Filtered count matrix.
        ref_samples: Reference condition sample names.
        treat_samples: Treatment condition sample names.

    Returns:
        DataFrame with DE results (without adjusted p-values).
    """
    results = []

    # Size factors drive the NB offset and all normalized-scale summaries
    size_factors = estimate_size_factors(counts)
    normalized = counts.div(size_factors, axis=1)
    ref_sf = size_factors[ref_samples].to_numpy(dtype=float)
    treat_sf = size_factors[treat_samples].to_numpy(dtype=float)

    for gene in counts.index:
        ref_counts = counts.loc[gene, ref_samples].values.astype(float)
        treat_counts = counts.loc[gene, treat_samples].values.astype(float)
        ref_norm = normalized.loc[gene, ref_samples].to_numpy(dtype=float)
        treat_norm = normalized.loc[gene, treat_samples].to_numpy(dtype=float)

        # Base mean on the size-factor-normalized scale
        base_mean = normalized.loc[gene].mean()

        # Fallback log2 fold change from normalized group means
        ref_mean = ref_norm.mean() + 0.5  # Pseudocount
        treat_mean = treat_norm.mean() + 0.5
        log2fc = float(np.log2(treat_mean / ref_mean))

        # Negative binomial test with log(size factor) offsets
        log2fc_nb, pvalue, _ = _negative_binomial_test(
            ref_counts, treat_counts, size_factors_a=ref_sf, size_factors_b=treat_sf
        )

        # Use the NB-derived fold change if valid
        if not np.isnan(log2fc_nb):
            log2fc = log2fc_nb

        # Wald statistic on the log2 scale via the delta method: the SE of a
        # log2 fold change is derived from the offset-adjusted NB variance
        # (Var(y_j / sf_j) = m / sf_j + dispersion * m^2) propagated through
        # the log transform, so the statistic is unit-consistent (z-like)
        # with the normalized log2 fold change. The p-value comes from the
        # NB test above, never from this statistic.
        all_norm = np.concatenate([ref_norm, treat_norm])
        dispersion = (
            _estimate_dispersion(all_norm) if all_norm.var() > all_norm.mean() else 0.0
        )

        def _se_log2_term(sf_values: "np.ndarray", mean_norm: float) -> float:
            mean = mean_norm + 0.5  # match the log2fc pseudocount
            total_variance = float(np.sum(mean / sf_values + dispersion * mean**2))
            n = len(sf_values)
            return float(total_variance / (n * mean) ** 2)

        se_log2 = float(
            np.sqrt(
                _se_log2_term(ref_sf, ref_norm.mean())
                + _se_log2_term(treat_sf, treat_norm.mean())
            )
            / np.log(2)
        )
        wald_stat = float(log2fc / se_log2) if se_log2 > 0 else 0.0

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

    Welch's t-test on log2 size-factor-normalized counts; fold changes
    and base means are also reported on the normalized scale.

    Args:
        counts: Filtered count matrix.
        ref_samples: Reference condition sample names.
        treat_samples: Treatment condition sample names.

    Returns:
        DataFrame with DE results.
    """
    # Size-factor normalize, then log-transform for the t-test
    size_factors = estimate_size_factors(counts)
    normalized = counts.div(size_factors, axis=1)
    log_counts = np.log2(normalized + 1)

    results = []
    for gene in counts.index:
        ref_vals = log_counts.loc[gene, ref_samples].values
        treat_vals = log_counts.loc[gene, treat_samples].values

        # Fold change from normalized group means
        ref_mean = normalized.loc[gene, ref_samples].mean() + 0.5
        treat_mean = normalized.loc[gene, treat_samples].mean() + 0.5
        log2fc = np.log2(treat_mean / ref_mean)

        # T-test. Zero-variance groups (fully tied normalized values, which
        # depth-shifted designs produce) are handled analytically: scipy's
        # moment calculation would warn on catastrophic cancellation, and the
        # deterministic answers are exact anyway.
        treat_arr = np.asarray(treat_vals, dtype=float)
        ref_arr = np.asarray(ref_vals, dtype=float)
        var_treat = float(treat_arr.var(ddof=1)) if treat_arr.size > 1 else 0.0
        var_ref = float(ref_arr.var(ddof=1)) if ref_arr.size > 1 else 0.0
        if var_treat == 0.0 and var_ref == 0.0:
            t_stat = 0.0
            pvalue = 1.0 if float(treat_arr.mean()) == float(ref_arr.mean()) else 0.0
        elif var_treat == 0.0 or var_ref == 0.0:
            # Exactly one constant group: the Welch SE collapses to the other
            # group's variance and df collapses to n_other - 1.
            se = float(np.sqrt(var_treat / treat_arr.size + var_ref / ref_arr.size))
            t_stat = (float(treat_arr.mean()) - float(ref_arr.mean())) / se
            df = (treat_arr.size - 1) if var_treat == 0.0 else (ref_arr.size - 1)
            pvalue = float(2.0 * stats.t.sf(abs(t_stat), df))
        else:
            # scipy warns on near-tied inputs (catastrophic cancellation in
            # moment calculations); the test result is still the best
            # available estimate, so keep it and don't let the warning
            # escalate under strict warning filters.
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                t_stat, pvalue = stats.ttest_ind(treat_arr, ref_arr, equal_var=False)
            if np.isnan(pvalue):
                pvalue = 1.0
                t_stat = 0.0

        results.append(
            {
                "gene": gene,
                "log2_fold_change": log2fc,
                "p_value": pvalue,
                "base_mean": normalized.loc[gene].mean(),
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

    The rank-sum test runs on the raw counts; fold changes and base means
    are reported on the size-factor-normalized scale.

    Args:
        counts: Filtered count matrix.
        ref_samples: Reference condition sample names.
        treat_samples: Treatment condition sample names.

    Returns:
        DataFrame with DE results.
    """
    size_factors = estimate_size_factors(counts)
    normalized = counts.div(size_factors, axis=1)

    results = []

    for gene in counts.index:
        ref_vals = counts.loc[gene, ref_samples].values.astype(float)
        treat_vals = counts.loc[gene, treat_samples].values.astype(float)

        # Fold change from normalized group means
        ref_mean = normalized.loc[gene, ref_samples].mean() + 0.5
        treat_mean = normalized.loc[gene, treat_samples].mean() + 0.5
        log2fc = np.log2(treat_mean / ref_mean)

        # Wilcoxon rank-sum (Mann-Whitney U) test
        try:
            u_stat, pvalue = stats.mannwhitneyu(
                treat_vals, ref_vals, alternative="two-sided"
            )
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
                "base_mean": normalized.loc[gene].mean(),
                "stat": u_stat,
            }
        )

    return pd.DataFrame(results)


def _negative_binomial_test(
    counts_a: np.ndarray,
    counts_b: np.ndarray,
    size_factors_a: Optional[np.ndarray] = None,
    size_factors_b: Optional[np.ndarray] = None,
) -> Tuple[float, float, float]:
    """Perform per-gene negative binomial offset likelihood ratio test.

    Fits negative binomial means on the raw counts with log(size factor)
    as an offset (exposure), so library-depth differences between samples
    are absorbed before testing. The likelihood ratio test compares a null
    model with a single shared normalized mean against an alternative with
    one normalized mean per group; both models share the pooled dispersion
    estimate, so the alternative adds exactly one mean parameter and the
    statistic is referred to chi-square with 1 df.

    Fold change and base mean are reported on the size-factor-normalized
    scale (raw counts divided by their sample size factors).

    Args:
        counts_a: Counts from condition A (reference).
        counts_b: Counts from condition B (treatment).
        size_factors_a: Size factors for condition A samples
            (defaults to all 1).
        size_factors_b: Size factors for condition B samples
            (defaults to all 1).

    Returns:
        Tuple of (log2_fold_change, p_value, base_mean) on the
        size-factor-normalized scale.
    """
    counts_a = np.asarray(counts_a, dtype=float)
    counts_b = np.asarray(counts_b, dtype=float)
    sf_a = (
        np.ones(counts_a.size)
        if size_factors_a is None
        else np.asarray(size_factors_a, dtype=float)
    )
    sf_b = (
        np.ones(counts_b.size)
        if size_factors_b is None
        else np.asarray(size_factors_b, dtype=float)
    )

    # Offset-adjusted (normalized) values drive all fitted means
    norm_a = counts_a / sf_a
    norm_b = counts_b / sf_b
    all_norm = np.concatenate([norm_a, norm_b])
    base_mean = float(all_norm.mean())

    # Group means on the normalized scale, with pseudocount so that
    # all-zero groups stay finite
    mean_a = float(norm_a.mean()) + 0.5
    mean_b = float(norm_b.mean()) + 0.5
    log2fc = float(np.log2(mean_b / mean_a))

    # Pooled dispersion from the offset-adjusted values, shared by the
    # null and alternative models
    dispersion = _estimate_dispersion(all_norm)

    # Log-likelihood for negative binomial with per-observation exposure
    def nb_loglik(
        counts: np.ndarray, sf_values: np.ndarray, mean_norm: float, dispersion: float
    ) -> float:
        """Compute negative binomial log-likelihood for offset means."""
        if mean_norm <= 0 or dispersion <= 0:
            return float("-inf")

        r = 1.0 / dispersion  # Size parameter
        loglik = 0.0
        for k, sf in zip(counts, sf_values):
            if k < 0:
                return float("-inf")
            mu = sf * mean_norm
            # NB log-likelihood: log(C(k+r-1, k)) + k*log(p) + r*log(1-p)
            # where p = mu/(mu + r)
            p = mu / (mu + r)
            loglik += gammaln(k + r) - gammaln(k + 1) - gammaln(r)
            loglik += k * np.log(p + 1e-10) + r * np.log(1 - p + 1e-10)
        return loglik

    # Null model: one shared normalized mean for both groups
    all_counts = np.concatenate([counts_a, counts_b])
    all_sf = np.concatenate([sf_a, sf_b])
    mean_null = float(all_norm.mean()) + 0.5
    ll_null = nb_loglik(all_counts, all_sf, mean_null, dispersion)

    # Alternative model: different normalized means per group, same
    # dispersion (single extra parameter => chi-square with 1 df)
    ll_alt = nb_loglik(counts_a, sf_a, mean_a, dispersion) + nb_loglik(
        counts_b, sf_b, mean_b, dispersion
    )

    # Likelihood ratio test
    lr_stat = 2 * (ll_alt - ll_null)

    if lr_stat < 0 or np.isnan(lr_stat):
        # Model fitting issue, fall back to no evidence of differential expression
        pvalue = 1.0
    else:
        pvalue = stats.chi2.sf(lr_stat, df=1)

    return log2fc, pvalue, base_mean


def _estimate_dispersion(counts: np.ndarray) -> float:
    """Estimate negative binomial dispersion parameter.

    Uses method of moments estimation with shrinkage toward a
    reasonable prior; the shrinkage weight decays with sample size
    (half-weight at n=2, near-zero at realistic replicate counts), so
    small samples lean on the prior while larger samples trust the
    method-of-moments estimate.

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

    # Shrink toward the prior with a weight that decays as the sample
    # grows: w = 2 / (2 + n) gives half-weight at n=2 and concentrates on
    # the method-of-moments estimate as replicate counts increase.
    prior_dispersion = 0.1
    shrinkage = 2.0 / (2.0 + len(counts))
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
        raise ValueError(
            f"Unknown p-value adjustment method: {method}. Valid: bh, fdr, bonferroni"
        )

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
) -> Dict[str, Any]:
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
        >>> pca_result = pca_analysis(normalized, n_components=3)
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
    pc_names = [f"PC{i + 1}" for i in range(n_components)]

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
        raise ValueError(
            f"Unknown distance method: {method}. Valid: euclidean, correlation, cosine"
        )

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
