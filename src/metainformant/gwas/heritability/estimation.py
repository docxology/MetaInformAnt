"""Heritability estimation methods for GWAS summary statistics.

Implements LD Score Regression (LDSC), partitioned heritability by functional
annotation, cross-trait genetic correlation, Haseman-Elston regression,
simplified GREML (Genomic REML), and liability-scale transformation for
case-control studies.

References:
    - Bulik-Sullivan et al. (2015) Nature Genetics 47:291-295 (LDSC)
    - Finucane et al. (2015) Nature Genetics 47:1228-1235 (partitioned LDSC)
    - Bulik-Sullivan et al. (2015) Nature Genetics 47:1236-1241 (genetic correlation)
    - Haseman & Elston (1972) Behavior Genetics 2:3-19 (HE regression)
    - Yang et al. (2011) Nature Genetics 43:519-525 (GREML)
    - Lee et al. (2012) Genetic Epidemiology 36:214-224 (liability scale)
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Tuple

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None  # type: ignore[assignment]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def estimate_h2_ldsc(
    chi2_stats: list[float],
    ld_scores: list[float],
    n_samples: int,
    intercept: bool = True,
) -> dict:
    """Estimate SNP heritability using LD Score Regression.

    Regresses chi-squared statistics on LD scores using weighted least
    squares. The slope estimates h2 = slope * M / N, where M is the
    number of SNPs and N is sample size. The intercept captures
    confounding bias (population stratification, cryptic relatedness).

    Args:
        chi2_stats: Chi-squared statistics for each SNP (z^2 values).
        ld_scores: LD scores for each SNP (sum of r^2 with all other SNPs).
        n_samples: GWAS sample size.
        intercept: Whether to fit an intercept (default True). If False,
            the intercept is constrained to 1.0 (no confounding).

    Returns:
        Dictionary with keys:
            - status: "success" or "error"
            - h2: Estimated SNP heritability
            - h2_se: Standard error of h2
            - intercept: Estimated intercept (1.0 if not fitted)
            - intercept_se: Standard error of intercept
            - lambda_gc: Genomic inflation factor (median chi2 / 0.4549)
            - mean_chi2: Mean chi-squared statistic
    """
    if not chi2_stats or not ld_scores:
        return {"status": "error", "message": "Chi2 stats and LD scores are required"}

    if len(chi2_stats) != len(ld_scores):
        return {
            "status": "error",
            "message": (f"Length mismatch: {len(chi2_stats)} chi2 stats vs " f"{len(ld_scores)} LD scores"),
        }

    m = len(chi2_stats)  # Number of SNPs
    n = n_samples

    if n < 1:
        return {"status": "error", "message": f"n_samples must be positive, got {n}"}

    if m < 3:
        return {"status": "error", "message": f"Need at least 3 SNPs, got {m}"}

    logger.debug(f"Running LDSC with {m} SNPs, N={n}, intercept={intercept}")

    # Compute summary statistics
    mean_chi2 = sum(chi2_stats) / m
    sorted_chi2 = sorted(chi2_stats)
    median_chi2 = sorted_chi2[m // 2] if m % 2 == 1 else (sorted_chi2[m // 2 - 1] + sorted_chi2[m // 2]) / 2.0
    lambda_gc = median_chi2 / 0.4549  # Expected median of chi2(1)

    # Weighted least squares regression: chi2 = intercept + slope * ld_score
    # Weights: 1 / (ld_score^2) to account for heteroscedasticity
    if HAS_NUMPY:
        result = _ldsc_regression_numpy(chi2_stats, ld_scores, n, m, intercept)
    else:
        result = _ldsc_regression_pure(chi2_stats, ld_scores, n, m, intercept)

    if result is None:
        return {"status": "error", "message": "LDSC regression failed (singular matrix)"}

    slope, slope_se, intercept_val, intercept_se = result

    # Convert slope to h2: h2 = slope * M / N
    h2 = slope * m / n
    h2_se = slope_se * m / n

    # Clamp h2 to [0, 1]
    h2 = max(0.0, min(1.0, h2))

    return {
        "status": "success",
        "h2": float(h2),
        "h2_se": float(h2_se),
        "intercept": float(intercept_val),
        "intercept_se": float(intercept_se),
        "lambda_gc": float(lambda_gc),
        "mean_chi2": float(mean_chi2),
    }


def partitioned_h2(
    chi2_stats: list[float],
    ld_scores_by_category: dict,
    n_samples: int,
) -> dict:
    """Estimate partitioned heritability by functional annotation category.

    Extends LDSC to multiple annotation categories by regressing chi-squared
    statistics on category-specific LD scores. Each category's contribution
    to heritability is estimated simultaneously.

    Args:
        chi2_stats: Chi-squared statistics for each SNP.
        ld_scores_by_category: Dictionary mapping category names (e.g.,
            "promoter", "enhancer", "coding", "intergenic") to lists of
            LD scores for that category. All lists must have the same
            length as chi2_stats.
        n_samples: GWAS sample size.

    Returns:
        Dictionary with keys:
            - status: "success" or "error"
            - per_category: Dict mapping category name to h2, h2_se,
              enrichment, and proportion of h2
            - total_h2: Sum of category h2 values
            - n_categories: Number of categories
            - n_snps: Number of SNPs
    """
    if not chi2_stats:
        return {"status": "error", "message": "Chi2 stats are required"}

    if not ld_scores_by_category:
        return {"status": "error", "message": "LD scores by category are required"}

    m = len(chi2_stats)
    n = n_samples

    # Validate all categories have correct length
    for cat_name, cat_scores in ld_scores_by_category.items():
        if len(cat_scores) != m:
            return {
                "status": "error",
                "message": (f"Category '{cat_name}' has {len(cat_scores)} LD scores, " f"expected {m}"),
            }

    categories = list(ld_scores_by_category.keys())
    n_categories = len(categories)

    logger.info(f"Running partitioned LDSC: {m} SNPs, {n_categories} categories, N={n}")

    if not HAS_NUMPY:
        # Fallback: run separate LDSC per category
        return _partitioned_h2_separate(chi2_stats, ld_scores_by_category, n_samples)

    # Build design matrix: [intercept, cat1_ld, cat2_ld, ...]
    X = np.ones((m, n_categories + 1))
    for c_idx, cat_name in enumerate(categories):
        X[:, c_idx + 1] = np.array(ld_scores_by_category[cat_name])

    y = np.array(chi2_stats, dtype=float)

    # Weights: inverse of predicted variance
    # Initial weights: 1 / (1 + ld_total)^2
    total_ld = np.sum(X[:, 1:], axis=1)
    weights = 1.0 / np.maximum(1.0 + total_ld, 1e-6) ** 2

    # Weighted least squares
    W = np.diag(weights)
    try:
        XtWX = X.T @ W @ X
        XtWy = X.T @ W @ y
        beta = np.linalg.solve(XtWX, XtWy)
    except np.linalg.LinAlgError:
        return {"status": "error", "message": "Partitioned LDSC regression failed"}

    # Standard errors via sandwich estimator
    residuals = y - X @ beta
    sigma2 = float(np.sum(weights * residuals**2) / max(m - n_categories - 1, 1))

    try:
        XtWX_inv = np.linalg.inv(XtWX)
    except np.linalg.LinAlgError:
        XtWX_inv = np.eye(n_categories + 1) * 1e6

    beta_se = np.sqrt(np.maximum(sigma2 * np.diag(XtWX_inv), 0.0))

    # Convert slopes to h2 per category
    per_category: dict[str, dict] = {}
    total_h2 = 0.0

    for c_idx, cat_name in enumerate(categories):
        slope = float(beta[c_idx + 1])
        slope_se_val = float(beta_se[c_idx + 1])

        cat_h2 = slope * m / n
        cat_h2_se = slope_se_val * m / n

        per_category[cat_name] = {
            "h2": float(cat_h2),
            "h2_se": float(cat_h2_se),
            "coefficient": float(slope),
        }
        total_h2 += cat_h2

    # Compute enrichment and proportions
    for cat_name in categories:
        cat_h2 = per_category[cat_name]["h2"]

        # Proportion of total h2
        proportion = cat_h2 / max(abs(total_h2), 1e-12)
        per_category[cat_name]["proportion"] = float(proportion)

        # Enrichment: (proportion of h2) / (proportion of SNPs in category)
        cat_ld_scores = ld_scores_by_category[cat_name]
        n_in_cat = sum(1 for s in cat_ld_scores if s > 0)
        prop_snps = n_in_cat / max(m, 1)
        enrichment = proportion / max(prop_snps, 1e-12) if prop_snps > 0 else 0.0
        per_category[cat_name]["enrichment"] = float(enrichment)

    total_h2 = max(0.0, min(1.0, total_h2))

    return {
        "status": "success",
        "per_category": per_category,
        "total_h2": float(total_h2),
        "n_categories": n_categories,
        "n_snps": m,
    }


def genetic_correlation(
    z1: list[float],
    z2: list[float],
    ld_scores: list[float],
    n1: int,
    n2: int,
    n_overlap: int = 0,
) -> dict:
    """Estimate genetic correlation between two traits using cross-trait LDSC.

    Regresses the product of Z-scores (z1 * z2) on LD scores to estimate
    the genetic covariance, then normalizes by the geometric mean of the
    individual heritabilities to obtain the genetic correlation (rg).

    Args:
        z1: Z-scores for trait 1.
        z2: Z-scores for trait 2.
        ld_scores: LD scores for each SNP.
        n1: Sample size for trait 1.
        n2: Sample size for trait 2.
        n_overlap: Number of overlapping samples between the two studies
            (default 0). Used to correct the intercept for sample overlap.

    Returns:
        Dictionary with keys:
            - status: "success" or "error"
            - rg: Estimated genetic correlation
            - rg_se: Standard error of rg
            - p_value: P-value for rg != 0
            - h2_1: Estimated h2 for trait 1
            - h2_2: Estimated h2 for trait 2
            - intercept: Cross-trait LDSC intercept
    """
    if not z1 or not z2 or not ld_scores:
        return {"status": "error", "message": "Z-scores and LD scores are required"}

    if len(z1) != len(z2) or len(z1) != len(ld_scores):
        return {
            "status": "error",
            "message": (f"Length mismatch: z1={len(z1)}, z2={len(z2)}, " f"ld_scores={len(ld_scores)}"),
        }

    m = len(z1)
    if m < 3:
        return {"status": "error", "message": f"Need at least 3 SNPs, got {m}"}

    logger.debug(f"Estimating genetic correlation: {m} SNPs, N1={n1}, N2={n2}, " f"overlap={n_overlap}")

    # Step 1: Estimate h2 for each trait individually
    chi2_1 = [z * z for z in z1]
    chi2_2 = [z * z for z in z2]

    h2_1_result = estimate_h2_ldsc(chi2_1, ld_scores, n1)
    h2_2_result = estimate_h2_ldsc(chi2_2, ld_scores, n2)

    if h2_1_result["status"] != "success" or h2_2_result["status"] != "success":
        return {
            "status": "error",
            "message": "Failed to estimate individual heritabilities",
        }

    h2_1 = h2_1_result["h2"]
    h2_2 = h2_2_result["h2"]

    # Step 2: Cross-trait LDSC
    # Dependent variable: product of Z-scores
    z_product = [a * b for a, b in zip(z1, z2)]

    # Geometric mean of sample sizes
    n_eff = math.sqrt(n1 * n2)

    # Regression of z_product on ld_scores
    if HAS_NUMPY:
        cross_result = _ldsc_regression_numpy(z_product, ld_scores, int(n_eff), m, True)
    else:
        cross_result = _ldsc_regression_pure(z_product, ld_scores, int(n_eff), m, True)

    if cross_result is None:
        return {"status": "error", "message": "Cross-trait LDSC regression failed"}

    slope, slope_se, intercept_val, intercept_se = cross_result

    # Genetic covariance: cov_g = slope * M / sqrt(N1 * N2)
    cov_g = slope * m / n_eff

    # Genetic correlation: rg = cov_g / sqrt(h2_1 * h2_2)
    denom = math.sqrt(max(h2_1, 1e-12) * max(h2_2, 1e-12))
    rg = cov_g / denom if denom > 1e-12 else 0.0

    # Clamp rg to [-1, 1]
    rg = max(-1.0, min(1.0, rg))

    # Standard error of rg (delta method approximation)
    cov_g_se = slope_se * m / n_eff
    rg_se = cov_g_se / denom if denom > 1e-12 else 1.0

    # P-value (Wald test)
    z_stat = rg / max(rg_se, 1e-12)
    p_value = 2.0 * (1.0 - _normal_cdf(abs(z_stat)))

    return {
        "status": "success",
        "rg": float(rg),
        "rg_se": float(rg_se),
        "p_value": float(p_value),
        "h2_1": float(h2_1),
        "h2_2": float(h2_2),
        "intercept": float(intercept_val),
    }


def haseman_elston_regression(
    grm_offdiag: list[float],
    phenotype_products: list[float],
) -> dict:
    """Estimate heritability using Haseman-Elston regression.

    Regresses pairwise phenotype products (or squared differences) on
    pairwise genetic relatedness values from the GRM off-diagonal.
    The slope estimates the genetic variance component.

    Args:
        grm_offdiag: Off-diagonal elements of the genomic relationship
            matrix, flattened. For n samples, this has n*(n-1)/2 elements.
        phenotype_products: Pairwise products of centered phenotype values,
            in the same order as grm_offdiag. Specifically,
            (y_i - mean(y)) * (y_j - mean(y)) for all pairs i < j.

    Returns:
        Dictionary with keys:
            - status: "success" or "error"
            - h2: Estimated heritability
            - h2_se: Standard error of h2
            - slope: Regression slope (sigma_g^2)
            - intercept: Regression intercept (sigma_e^2)
            - n_pairs: Number of sample pairs used
    """
    if not grm_offdiag or not phenotype_products:
        return {"status": "error", "message": "GRM off-diagonal and phenotype products are required"}

    if len(grm_offdiag) != len(phenotype_products):
        return {
            "status": "error",
            "message": (
                f"Length mismatch: {len(grm_offdiag)} GRM values vs " f"{len(phenotype_products)} phenotype products"
            ),
        }

    n_pairs = len(grm_offdiag)
    if n_pairs < 3:
        return {"status": "error", "message": f"Need at least 3 pairs, got {n_pairs}"}

    logger.debug(f"Running Haseman-Elston regression with {n_pairs} pairs")

    # Simple linear regression: phenotype_product = intercept + slope * grm
    x = grm_offdiag
    y = phenotype_products

    n = n_pairs
    x_mean = sum(x) / n
    y_mean = sum(y) / n

    # Compute slope and intercept
    ss_xy = sum((xi - x_mean) * (yi - y_mean) for xi, yi in zip(x, y))
    ss_xx = sum((xi - x_mean) ** 2 for xi in x)

    if abs(ss_xx) < 1e-20:
        return {"status": "error", "message": "Zero variance in GRM off-diagonal values"}

    slope = ss_xy / ss_xx
    intercept_val = y_mean - slope * x_mean

    # Residual variance and standard error
    y_pred = [intercept_val + slope * xi for xi in x]
    residuals = [yi - yp for yi, yp in zip(y, y_pred)]
    mse = sum(r**2 for r in residuals) / max(n - 2, 1)
    slope_se = math.sqrt(mse / max(ss_xx, 1e-20))

    # h2 = slope / (slope + intercept) if both positive
    # In HE regression, slope = sigma_g^2, intercept approx sigma_e^2
    total_var = abs(slope) + abs(intercept_val)
    if total_var > 1e-12:
        h2 = max(0.0, min(1.0, slope / total_var))
    else:
        h2 = 0.0

    # SE of h2 (delta method)
    if total_var > 1e-12:
        h2_se = slope_se * abs(intercept_val) / (total_var**2)
    else:
        h2_se = 0.0

    return {
        "status": "success",
        "h2": float(h2),
        "h2_se": float(h2_se),
        "slope": float(slope),
        "intercept": float(intercept_val),
        "n_pairs": n_pairs,
    }


def greml_simple(
    grm: list[list[float]] | Any,
    phenotypes: list[float],
) -> dict:
    """Simplified GREML (Genomic REML) variance component estimation.

    Performs iterative REML to estimate the genetic and residual variance
    components from a genomic relationship matrix and phenotype vector.
    Uses eigendecomposition of the GRM for efficient likelihood evaluation,
    and a grid search followed by golden-section refinement.

    This is equivalent to the estimation in the existing
    metainformant.gwas.analysis.heritability module but exposed here as
    part of the heritability submodule for consistency.

    Args:
        grm: Genomic relationship matrix (n x n). Accepts numpy array or
            list of lists.
        phenotypes: Phenotype values for each sample.

    Returns:
        Dictionary with keys:
            - status: "success" or "error"
            - h2: Estimated heritability
            - h2_se: Standard error
            - sigma_g: Genetic variance component
            - sigma_e: Residual variance component
            - log_likelihood: Maximized REML log-likelihood
            - n_samples: Number of samples
            - method: "greml"
    """
    if not HAS_NUMPY:
        return {"status": "error", "message": "numpy is required for GREML estimation"}

    n = len(phenotypes)
    if n < 3:
        return {"status": "error", "message": f"Need at least 3 samples, got {n}"}

    try:
        K = np.asarray(grm, dtype=float)
        y = np.array(phenotypes, dtype=float)
    except (ValueError, TypeError) as exc:
        return {"status": "error", "message": f"Failed to convert inputs: {exc}"}

    if K.shape != (n, n):
        return {
            "status": "error",
            "message": f"GRM shape {K.shape} does not match {n} samples",
        }

    if np.var(y) < 1e-12:
        return {"status": "error", "message": "Phenotype has zero or near-zero variance"}

    logger.debug(f"Running GREML with {n} samples")

    try:
        # Eigendecompose K
        eigenvalues, eigenvectors = np.linalg.eigh(K)
        idx = np.argsort(eigenvalues)[::-1]
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]
        eigenvalues = np.maximum(eigenvalues, 0.0)

        # Rotate phenotypes
        Ut = eigenvectors.T
        y_rot = Ut @ y
        ones_rot = Ut @ np.ones(n)

        # Grid search over h2
        n_grid = 100
        h2_grid = np.linspace(0.001, 0.999, n_grid)
        best_ll = -np.inf
        best_h2 = 0.5

        for h2_candidate in h2_grid:
            ll = _reml_log_likelihood(y_rot, ones_rot, eigenvalues, float(h2_candidate), n)
            if ll > best_ll:
                best_ll = ll
                best_h2 = float(h2_candidate)

        # Refine with golden section search
        step = 1.0 / n_grid
        lo = max(0.001, best_h2 - step)
        hi = min(0.999, best_h2 + step)
        best_h2 = _golden_section_max(
            lambda h2: _reml_log_likelihood(y_rot, ones_rot, eigenvalues, h2, n),
            lo,
            hi,
        )
        best_ll = _reml_log_likelihood(y_rot, ones_rot, eigenvalues, best_h2, n)

        # Variance components
        total_var = float(np.var(y))
        sigma_g = best_h2 * total_var
        sigma_e = (1.0 - best_h2) * total_var

        # Approximate SE
        h2_se = _approximate_h2_se(y_rot, ones_rot, eigenvalues, best_h2, n)

        return {
            "status": "success",
            "h2": float(best_h2),
            "h2_se": float(h2_se),
            "sigma_g": float(sigma_g),
            "sigma_e": float(sigma_e),
            "log_likelihood": float(best_ll),
            "n_samples": n,
            "method": "greml",
        }

    except Exception as exc:
        logger.warning(f"GREML estimation failed: {exc}")
        return {"status": "error", "message": str(exc)}


def compute_liability_h2(
    h2_observed: float,
    prevalence: float,
    sample_prevalence: float,
) -> dict:
    """Convert observed-scale h2 to liability-scale for case-control studies.

    For case-control GWAS, the observed-scale heritability underestimates
    the true liability-scale heritability due to ascertainment. Uses the
    Lee et al. (2012) transformation to correct for population prevalence
    and sample case/control ratio.

    The transformation is:
        h2_liability = h2_observed * K^2 * (1-K)^2 / (P * (1-P) * z^2)

    where K is population prevalence, P is sample prevalence, and z is
    the standard normal density at the liability threshold.

    Args:
        h2_observed: Heritability estimated on the observed (0/1) scale.
        prevalence: Population prevalence of the disease (K).
        sample_prevalence: Proportion of cases in the study sample (P).

    Returns:
        Dictionary with keys:
            - status: "success" or "error"
            - h2_liability: Liability-scale heritability
            - h2_observed: Input observed-scale h2
            - prevalence: Population prevalence used
            - sample_prevalence: Sample prevalence used
            - threshold: Liability threshold
            - conversion_factor: Multiplicative factor applied
    """
    if h2_observed < 0 or h2_observed > 1:
        return {
            "status": "error",
            "message": f"h2_observed must be in [0, 1], got {h2_observed}",
        }

    if prevalence <= 0 or prevalence >= 1:
        return {
            "status": "error",
            "message": f"Prevalence must be in (0, 1), got {prevalence}",
        }

    if sample_prevalence <= 0 or sample_prevalence >= 1:
        return {
            "status": "error",
            "message": f"Sample prevalence must be in (0, 1), got {sample_prevalence}",
        }

    logger.debug(
        f"Converting h2_observed={h2_observed:.4f} to liability scale " f"(K={prevalence}, P={sample_prevalence})"
    )

    K = prevalence
    P = sample_prevalence

    # Liability threshold: Phi^{-1}(1 - K)
    threshold = _normal_quantile(1.0 - K)

    # Standard normal density at the threshold
    z_density = _normal_pdf(threshold)

    if z_density < 1e-20:
        return {"status": "error", "message": "Normal density at threshold is too small"}

    # Lee et al. conversion factor
    numerator = K**2 * (1.0 - K) ** 2
    denominator = P * (1.0 - P) * z_density**2

    if denominator < 1e-20:
        return {"status": "error", "message": "Conversion factor denominator is too small"}

    conversion_factor = numerator / denominator
    h2_liability = h2_observed * conversion_factor

    # Clamp to [0, 1]
    h2_liability = max(0.0, min(1.0, h2_liability))

    return {
        "status": "success",
        "h2_liability": float(h2_liability),
        "h2_observed": float(h2_observed),
        "prevalence": float(prevalence),
        "sample_prevalence": float(sample_prevalence),
        "threshold": float(threshold),
        "conversion_factor": float(conversion_factor),
    }


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _ldsc_regression_numpy(
    y_vals: list[float],
    ld_scores: list[float],
    n_samples: int,
    m_snps: int,
    fit_intercept: bool,
) -> Optional[Tuple[float, float, float, float]]:
    """LDSC regression using numpy weighted least squares.

    Args:
        y_vals: Dependent variable (chi2 or z-product).
        ld_scores: LD scores (independent variable).
        n_samples: Sample size (for weight computation).
        m_snps: Number of SNPs.
        fit_intercept: Whether to fit an intercept.

    Returns:
        Tuple of (slope, slope_se, intercept, intercept_se) or None.
    """
    y_arr = np.array(y_vals, dtype=float)
    ld_arr = np.array(ld_scores, dtype=float)

    # Weights: inverse of approximate variance of chi2
    # Var(chi2) ~ 2 * (1 + N * h2 * l_j / M)^2
    # Simplified: w_j = 1 / max(l_j, 1)^2
    weights = 1.0 / np.maximum(ld_arr, 1.0) ** 2

    if fit_intercept:
        X = np.column_stack([np.ones(m_snps), ld_arr])
    else:
        X = ld_arr.reshape(-1, 1)

    W = np.diag(weights)

    try:
        XtWX = X.T @ W @ X
        XtWy = X.T @ W @ y_arr
        beta = np.linalg.solve(XtWX, XtWy)
    except np.linalg.LinAlgError:
        return None

    # Standard errors
    residuals = y_arr - X @ beta
    dof = m_snps - X.shape[1]
    if dof <= 0:
        dof = 1
    sigma2 = float(np.sum(weights * residuals**2) / dof)

    try:
        XtWX_inv = np.linalg.inv(XtWX)
    except np.linalg.LinAlgError:
        return None

    beta_se = np.sqrt(np.maximum(sigma2 * np.diag(XtWX_inv), 0.0))

    if fit_intercept:
        intercept_val = float(beta[0])
        intercept_se = float(beta_se[0])
        slope = float(beta[1])
        slope_se = float(beta_se[1])
    else:
        intercept_val = 1.0
        intercept_se = 0.0
        slope = float(beta[0])
        slope_se = float(beta_se[0])

    return slope, slope_se, intercept_val, intercept_se


def _ldsc_regression_pure(
    y_vals: list[float],
    ld_scores: list[float],
    n_samples: int,
    m_snps: int,
    fit_intercept: bool,
) -> Optional[Tuple[float, float, float, float]]:
    """Pure Python LDSC regression fallback.

    Args:
        y_vals: Dependent variable.
        ld_scores: LD scores.
        n_samples: Sample size.
        m_snps: Number of SNPs.
        fit_intercept: Whether to fit an intercept.

    Returns:
        Tuple of (slope, slope_se, intercept, intercept_se) or None.
    """
    n = m_snps

    # Compute weights
    weights = [1.0 / max(l, 1.0) ** 2 for l in ld_scores]
    total_w = sum(weights)

    if fit_intercept:
        # Weighted means
        w_sum_x = sum(w * x for w, x in zip(weights, ld_scores))
        w_sum_y = sum(w * y for w, y in zip(weights, y_vals))
        x_mean = w_sum_x / total_w
        y_mean = w_sum_y / total_w

        # Weighted covariance and variance
        ss_xy = sum(w * (x - x_mean) * (y - y_mean) for w, x, y in zip(weights, ld_scores, y_vals))
        ss_xx = sum(w * (x - x_mean) ** 2 for w, x in zip(weights, ld_scores))

        if abs(ss_xx) < 1e-20:
            return None

        slope = ss_xy / ss_xx
        intercept_val = y_mean - slope * x_mean

        # Residuals and SE
        y_pred = [intercept_val + slope * x for x in ld_scores]
        residuals = [y - yp for y, yp in zip(y_vals, y_pred)]
        w_rss = sum(w * r**2 for w, r in zip(weights, residuals))
        mse = w_rss / max(n - 2, 1)

        slope_se = math.sqrt(mse / max(ss_xx, 1e-20))

        # Intercept SE
        sum_w_x2 = sum(w * x**2 for w, x in zip(weights, ld_scores))
        intercept_se = math.sqrt(mse * sum_w_x2 / max(total_w * ss_xx, 1e-20))
    else:
        # No intercept: constrained model
        intercept_val = 1.0
        y_adj = [y - 1.0 for y in y_vals]

        w_sum_xy = sum(w * x * y for w, x, y in zip(weights, ld_scores, y_adj))
        w_sum_xx = sum(w * x**2 for w, x in zip(weights, ld_scores))

        if abs(w_sum_xx) < 1e-20:
            return None

        slope = w_sum_xy / w_sum_xx

        y_pred = [1.0 + slope * x for x in ld_scores]
        residuals = [y - yp for y, yp in zip(y_vals, y_pred)]
        w_rss = sum(w * r**2 for w, r in zip(weights, residuals))
        mse = w_rss / max(n - 1, 1)

        slope_se = math.sqrt(mse / max(w_sum_xx, 1e-20))
        intercept_se = 0.0

    return slope, slope_se, intercept_val, intercept_se


def _partitioned_h2_separate(
    chi2_stats: list[float],
    ld_scores_by_category: dict,
    n_samples: int,
) -> dict:
    """Fallback: estimate partitioned h2 by running separate LDSC per category.

    This is a simplified approach when numpy is not available.

    Args:
        chi2_stats: Chi-squared statistics.
        ld_scores_by_category: LD scores per category.
        n_samples: Sample size.

    Returns:
        Partitioned heritability results dict.
    """
    per_category: dict[str, dict] = {}
    total_h2 = 0.0

    for cat_name, cat_ld in ld_scores_by_category.items():
        result = estimate_h2_ldsc(chi2_stats, cat_ld, n_samples)
        if result["status"] == "success":
            cat_h2 = result["h2"]
            per_category[cat_name] = {
                "h2": cat_h2,
                "h2_se": result["h2_se"],
                "coefficient": 0.0,
                "enrichment": 0.0,
                "proportion": 0.0,
            }
            total_h2 += cat_h2
        else:
            per_category[cat_name] = {
                "h2": 0.0,
                "h2_se": 0.0,
                "coefficient": 0.0,
                "enrichment": 0.0,
                "proportion": 0.0,
            }

    # Compute proportions
    for cat_name in per_category:
        if total_h2 > 1e-12:
            per_category[cat_name]["proportion"] = per_category[cat_name]["h2"] / total_h2

    total_h2 = max(0.0, min(1.0, total_h2))

    return {
        "status": "success",
        "per_category": per_category,
        "total_h2": float(total_h2),
        "n_categories": len(ld_scores_by_category),
        "n_snps": len(chi2_stats),
    }


def _reml_log_likelihood(
    y_rot: Any,
    ones_rot: Any,
    eigenvalues: Any,
    h2: float,
    n: int,
) -> float:
    """Compute the REML log-likelihood for a given h2.

    In the rotated eigenspace of the GRM, the variance matrix is diagonal:
    V_ii = h2 * lambda_i + (1 - h2).

    Args:
        y_rot: Rotated phenotypes (U' y).
        ones_rot: Rotated intercept (U' 1).
        eigenvalues: Eigenvalues of the GRM.
        h2: Candidate heritability in (0, 1).
        n: Number of samples.

    Returns:
        REML log-likelihood value.
    """
    d = h2 * eigenvalues + (1.0 - h2)
    d = np.maximum(d, 1e-10)
    d_inv = 1.0 / d

    # WLS intercept
    ot_dinv_o = float(np.sum(d_inv * ones_rot**2))
    if ot_dinv_o < 1e-30:
        return -1e30
    ot_dinv_y = float(np.sum(d_inv * ones_rot * y_rot))
    beta_hat = ot_dinv_y / ot_dinv_o

    residuals = y_rot - beta_hat * ones_rot
    p = 1

    sigma2 = float(np.sum(d_inv * residuals**2)) / (n - p)
    if sigma2 <= 0:
        return -1e30

    ll = -0.5 * (n - p) * math.log(2 * math.pi * sigma2)
    ll -= 0.5 * float(np.sum(np.log(d)))
    ll -= 0.5 * (n - p)
    ll -= 0.5 * math.log(ot_dinv_o)

    return float(ll)


def _approximate_h2_se(
    y_rot: Any,
    ones_rot: Any,
    eigenvalues: Any,
    h2: float,
    n: int,
    delta: float = 0.005,
) -> float:
    """Approximate SE of h2 from the curvature of the log-likelihood.

    Args:
        y_rot: Rotated phenotypes.
        ones_rot: Rotated intercept.
        eigenvalues: GRM eigenvalues.
        h2: Estimated h2.
        n: Number of samples.
        delta: Finite difference step size.

    Returns:
        Approximate standard error.
    """
    h2_lo = max(0.001, h2 - delta)
    h2_hi = min(0.999, h2 + delta)
    actual_delta = (h2_hi - h2_lo) / 2.0

    ll_lo = _reml_log_likelihood(y_rot, ones_rot, eigenvalues, h2_lo, n)
    ll_mid = _reml_log_likelihood(y_rot, ones_rot, eigenvalues, h2, n)
    ll_hi = _reml_log_likelihood(y_rot, ones_rot, eigenvalues, h2_hi, n)

    d2ll = (ll_hi - 2.0 * ll_mid + ll_lo) / (actual_delta**2)

    if d2ll >= 0:
        return 0.5

    fisher_info = -d2ll
    se = 1.0 / math.sqrt(fisher_info)

    return min(se, 0.5)


def _golden_section_max(
    func: Any,
    a: float,
    b: float,
    tol: float = 1e-6,
    max_iter: int = 100,
) -> float:
    """Golden section search for the maximum of a unimodal function on [a, b].

    Args:
        func: Function to maximize.
        a: Lower bound.
        b: Upper bound.
        tol: Convergence tolerance.
        max_iter: Maximum iterations.

    Returns:
        x at maximum.
    """
    gr = (math.sqrt(5) + 1) / 2

    c = b - (b - a) / gr
    d = a + (b - a) / gr

    for _ in range(max_iter):
        if abs(b - a) < tol:
            break
        if func(c) < func(d):
            a = c
        else:
            b = d
        c = b - (b - a) / gr
        d = a + (b - a) / gr

    return (a + b) / 2.0


def _normal_cdf(x: float) -> float:
    """Approximate standard normal CDF (Abramowitz & Stegun).

    Args:
        x: Input value.

    Returns:
        CDF value.
    """
    a1 = 0.254829592
    a2 = -0.284496736
    a3 = 1.421413741
    a4 = -1.453152027
    a5 = 1.061405429
    p = 0.3275911

    sign = 1 if x >= 0 else -1
    x = abs(x) / math.sqrt(2.0)

    t = 1.0 / (1.0 + p * x)
    y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * math.exp(-x * x)

    return 0.5 * (1.0 + sign * y)


def _normal_pdf(x: float) -> float:
    """Standard normal probability density function.

    Args:
        x: Input value.

    Returns:
        PDF value.
    """
    return math.exp(-0.5 * x * x) / math.sqrt(2.0 * math.pi)


def _normal_quantile(p: float) -> float:
    """Approximate quantile of the standard normal distribution.

    Uses rational approximation from Abramowitz & Stegun 26.2.23.

    Args:
        p: Probability in (0, 1).

    Returns:
        Z-score such that P(Z <= z) = p.
    """
    if p <= 0.0:
        return -10.0
    if p >= 1.0:
        return 10.0
    if p == 0.5:
        return 0.0

    if p < 0.5:
        return -_normal_quantile(1.0 - p)

    t = math.sqrt(-2.0 * math.log(1.0 - p))
    c0 = 2.515517
    c1 = 0.802853
    c2 = 0.010328
    d1 = 1.432788
    d2 = 0.189269
    d3 = 0.001308

    z = t - (c0 + c1 * t + c2 * t * t) / (1.0 + d1 * t + d2 * t * t + d3 * t * t * t)
    return z
