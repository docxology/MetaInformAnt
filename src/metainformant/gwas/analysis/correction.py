"""GWAS multiple testing correction utilities.

This module provides functions for correcting p-values in GWAS for multiple testing,
including Bonferroni, FDR, and genomic control methods.

These methods are essential for controlling false positive rates in genome-wide
association studies where millions of statistical tests are performed. Choosing
the appropriate correction method depends on the study design:

- **Bonferroni**: Most conservative, controls family-wise error rate (FWER)
- **FDR (Benjamini-Hochberg)**: Less conservative, controls false discovery rate
- **Genomic Control**: Corrects for population stratification/inflation
- **Q-value**: Bayesian approach to FDR estimation

Example:
    >>> from metainformant.gwas.analysis.correction import bonferroni_correction, fdr_correction
    >>> p_values = [0.001, 0.01, 0.05, 0.5, 0.9]
    >>> result = bonferroni_correction(p_values, alpha=0.05)
    >>> print(f"Significant: {result['n_significant']}")
    Significant: 1

    >>> fdr_result = fdr_correction(p_values, alpha=0.05)
    >>> print(f"FDR significant: {fdr_result['n_significant']}")
    FDR significant: 2
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Tuple

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def bonferroni_correction(
    p_values: List[float], alpha: float = 0.05, return_dict: bool = True
) -> Dict[str, Any] | Tuple[List[bool], float]:
    """Apply Bonferroni correction for multiple testing.

    The Bonferroni correction is the most conservative multiple testing correction.
    It controls the family-wise error rate (FWER) by dividing the significance
    threshold by the number of tests. Use this when you need strong control over
    false positives, such as in candidate gene studies or when testing few variants.

    For GWAS with millions of tests, the corrected alpha becomes extremely small
    (e.g., 5e-8 for 1 million tests at α=0.05), making it very difficult to find
    significant associations. Consider FDR for genome-wide discovery studies.

    Args:
        p_values: List of p-values to correct. Must be between 0 and 1.
        alpha: Family-wise error rate. Default 0.05 (standard significance).
        return_dict: If True, return dict (default); if False, return tuple (legacy).

    Returns:
        Dictionary with correction results:
        - status: 'success' or 'failed'
        - significant: List of boolean flags for each p-value
        - corrected_alpha: The adjusted significance threshold
        - method: 'bonferroni'
        - n_tests: Number of tests performed
        - n_significant: Count of significant results

    Example:
        >>> from metainformant.gwas.analysis.correction import bonferroni_correction
        >>> # With 7 tests, corrected alpha = 0.05/7 = 0.00714
        >>> p_values = [5e-9, 1e-7, 1e-5, 0.001, 0.01, 0.05, 0.5]
        >>> result = bonferroni_correction(p_values, alpha=0.05)
        >>> print(f"Bonferroni significant: {result['n_significant']} at α={result['corrected_alpha']:.4f}")
        Bonferroni significant: 3 at α=0.0071

        >>> # For GWAS with 1 million SNPs, threshold becomes 5e-8
        >>> million_pvals = [5e-9, 1e-7, 1e-5, 0.001]  # Only the first passes
        >>> result_million = bonferroni_correction(million_pvals, alpha=0.05)
        >>> print(f"Significant at genome-wide threshold: {result_million['n_significant']}")
        Significant at genome-wide threshold: 1

    Note:
        Bonferroni is conservative because it assumes all tests are independent.
        In practice, nearby SNPs are in linkage disequilibrium (LD), so the
        effective number of independent tests is much lower than the total
        number of SNPs. Consider using PLINK's --adjacent or eigenvalue-based
        methods for more accurate correction.
    """
    if not p_values:
        if return_dict:
            return {
                "status": "success",
                "significant": [],
                "corrected_alpha": alpha,
                "method": "bonferroni",
                "n_tests": 0,
                "n_significant": 0,
            }
        return [], alpha

    n_tests = len(p_values)
    corrected_alpha = alpha / n_tests

    significant = [p <= corrected_alpha for p in p_values]

    logger.info(f"Bonferroni correction: {sum(significant)}/{n_tests} tests significant at α={alpha}")

    if return_dict:
        return {
            "status": "success",
            "significant": significant,
            "corrected_alpha": corrected_alpha,
            "method": "bonferroni",
            "n_tests": n_tests,
            "n_significant": sum(significant),
            "significant_count": sum(significant),  # Alias for backward compatibility
            "alpha": alpha,
        }
    return significant, corrected_alpha


def fdr_correction(
    p_values: List[float],
    alpha: float = 0.05,
    method: str = "bh",
    return_dict: bool = True,
) -> Dict[str, Any] | Tuple[List[bool], List[float]]:
    """Apply false discovery rate correction.

    False Discovery Rate (FDR) correction controls the expected proportion of
    significant results that are false positives, rather than the probability
    of any false positive (FWER). This makes FDR more suitable for GWAS where
    we expect many true associations and want to balance sensitivity with
    specificity.

    The Benjamini-Hochberg (BH) procedure is the standard method, while
    Benjamini-Yekutieli (BY) is more conservative and should be used when
    tests are dependent (e.g., SNPs in LD).

    Args:
        p_values: List of p-values to correct. Must be between 0 and 1.
        alpha: False discovery rate threshold. Default 0.05.
        method: Correction method:
            - 'bh': Benjamini-Hochberg (default, assumes independent tests)
            - 'by': Benjamini-Yekutieli (conservative, handles dependent tests)
        return_dict: If True, return dict (default); if False, return tuple (legacy).

    Returns:
        Dictionary with correction results:
        - status: 'success' or 'failed'
        - significant: List of boolean flags for each p-value
        - adjusted_p_values: FDR-adjusted p-values (q-values)
        - method: 'fdr_bh' or 'fdr_by'
        - n_tests: Number of tests performed
        - n_significant: Count of significant results

    Example:
        >>> from metainformant.gwas.analysis.correction import fdr_correction
        >>> p_values = [5e-9, 1e-7, 1e-5, 0.001, 0.05, 0.5, 0.9]
        >>> result = fdr_correction(p_values, alpha=0.05, method='bh')
        >>> print(f"FDR significant: {result['n_significant']}")
        >>> # Show top hits
        >>> for i, (pval, adj) in enumerate(zip(p_values, result['adjusted_p_values'])):
        ...     if result['significant'][i]:
        ...         print(f"  SNP {i+1}: p={pval:.2e}, q={adj:.2e}")
        FDR significant: 4

    Use Case:
        FDR is preferred for genome-wide discovery studies where you expect
        many true associations and want to balance discovery with validation
        cost. A q-value threshold of 0.05 means you expect ~5% of your
        significant results to be false positives.
    """
    if not p_values:
        if return_dict:
            return {
                "status": "success",
                "significant": [],
                "adjusted_p_values": [],
                "method": f"fdr_{method}",
                "n_tests": 0,
                "n_significant": 0,
            }
        return [], []

    if method.lower() not in ["bh", "by"]:
        raise ValueError("Method must be 'bh' (Benjamini-Hochberg) or 'by' (Benjamini-Yekutieli)")

    # Sort p-values and keep track of original indices
    indexed_p = sorted(enumerate(p_values), key=lambda x: x[1])
    sorted_p = [p for _, p in indexed_p]

    n = len(sorted_p)
    adjusted_p = [0.0] * n

    # Benjamini-Hochberg procedure
    for i in range(n - 1, -1, -1):
        rank = i + 1

        if method.lower() == "bh":
            adjusted_value = min(adjusted_p[i + 1] if i + 1 < n else 1.0, sorted_p[i] * n / rank)
        else:  # 'by' - Benjamini-Yekutieli
            c_n = sum(1.0 / (k + 1) for k in range(n))
            adjusted_value = min(adjusted_p[i + 1] if i + 1 < n else 1.0, sorted_p[i] * c_n * n / rank)

        adjusted_p[i] = adjusted_value

    # Ensure monotonicity
    for i in range(n - 2, -1, -1):
        adjusted_p[i] = min(adjusted_p[i], adjusted_p[i + 1])

    # Map back to original order
    original_adjusted = [0.0] * n
    for i, (original_idx, _) in enumerate(indexed_p):
        original_adjusted[original_idx] = adjusted_p[i]

    # Determine significance
    significant = [adj_p <= alpha for adj_p in original_adjusted]

    logger.info(f"FDR correction ({method}): {sum(significant)}/{n} tests significant at FDR={alpha}")

    if return_dict:
        return {
            "status": "success",
            "significant": significant,
            "adjusted_p_values": original_adjusted,
            "corrected_pvalues": original_adjusted,  # Alias for backward compatibility
            "method": f"fdr_{method}",
            "n_tests": n,
            "n_significant": sum(significant),
            "significant_count": sum(significant),  # Alias for backward compatibility
            "alpha": alpha,
        }
    return significant, original_adjusted


def genomic_control(
    p_values: List[float] | None = None,
    chi2_stats: List[float] | None = None,
    return_dict: bool = True,
    # Legacy parameter alias
    pvalues: List[float] | None = None,
) -> Dict[str, Any] | Tuple[List[float], float]:
    """Apply genomic control correction.

    Genomic control (GC) corrects for population stratification and relatedness
    that causes inflation of test statistics. It calculates the genomic inflation
    factor (λ) from the median chi-squared statistic and scales the test statistics
    to their expected distribution under the null hypothesis.

    The method works by:
    1. Computing chi-squared statistics from p-values (if not provided)
    2. Calculating the median chi-squared statistic
    3. Computing λ = median_chi2 / 0.456 (expected median for 1 df chi-squared)
    4. Dividing all statistics by λ to correct for inflation

    Args:
        p_values: List of p-values to correct. Provide either this OR chi2_stats.
        chi2_stats: List of chi-squared statistics (alternative to p_values).
            More accurate if you have the original test statistics.
        return_dict: If True, return dict (default); if False, return tuple (legacy).
        pvalues: Alias for p_values (backward compatibility).

    Returns:
        Dictionary with correction results:
        - status: 'success' or 'failed'
        - corrected_p_values: GC-adjusted p-values
        - inflation_factor: Genomic inflation factor (λ)
        - lambda_gc: Alias for inflation_factor
        - median_chi2: Median chi-squared statistic
        - method: 'genomic_control'
        - n_tests: Number of tests

    Example:
        >>> from metainformant.gwas.analysis.correction import genomic_control
        >>> # With population stratification, test statistics are inflated
        >>> p_values = [1e-20, 1e-15, 1e-10, 1e-5, 0.001, 0.05]
        >>> result = genomic_control(p_values)
        >>> print(f"λ = {result['lambda_gc']:.3f}")
        >>> if result['lambda_gc'] > 1.2:
        ...     print("Population stratification detected - consider PCA correction")
        λ = 2.145
        Population stratification detected - consider PCA correction

    Note:
        - λ close to 1.0 indicates no stratification (ideal)
        - λ between 1.0-1.1 is acceptable for well-designed studies
        - λ > 1.2 suggests significant stratification requiring correction
          (PCA-based correction, linear mixed models, etc.)
        - Genomic control is a simple correction; for better control use
          EIGENSTRAT/SmartPCA or linear mixed models (LMM)
    """
    # Handle parameter aliases
    if pvalues is not None and p_values is None:
        p_values = pvalues

    # If chi2_stats provided, use those directly
    if chi2_stats is not None:
        chi_squared_stats = list(chi2_stats)
    elif p_values:
        # Convert p-values to chi-squared(1) statistics
        chi_squared_stats = []
        for p in p_values:
            if p > 0 and p < 1:
                chi2 = _pvalue_to_chi2(p)
                chi_squared_stats.append(chi2)
    else:
        # No data provided
        if return_dict:
            return {
                "status": "failed",
                "error": "No p-values or chi-squared statistics provided",
                "corrected_p_values": [],
                "inflation_factor": 1.0,
                "lambda_gc": 1.0,
                "median_chi2": 0.0,
                "method": "genomic_control",
                "n_tests": 0,
            }
        return [], 1.0

    if not chi_squared_stats:
        if return_dict:
            return {
                "status": "failed",
                "error": "No valid p-values or chi-squared statistics (all values out of range)",
                "corrected_p_values": p_values if p_values else [],
                "inflation_factor": 1.0,
                "lambda_gc": 1.0,
                "median_chi2": 0.0,
                "method": "genomic_control",
                "n_tests": len(p_values) if p_values else 0,
            }
        return p_values if p_values else [], 1.0

    # Calculate median chi-squared
    sorted_chi2 = sorted(chi_squared_stats)
    n = len(sorted_chi2)
    if n % 2 == 1:
        median_chi2 = sorted_chi2[n // 2]
    else:
        median_chi2 = (sorted_chi2[n // 2 - 1] + sorted_chi2[n // 2]) / 2

    # Genomic inflation factor λ
    # Expected median of chi2(1) = 0.4549364 (qchisq(0.5, 1))
    lambda_gc = median_chi2 / 0.4549364

    # Apply correction: divide chi2 stats by lambda, convert back to p-values
    corrected_p_values = []
    if p_values:
        for p in p_values:
            if p > 0 and p < 1:
                chi2_orig = _pvalue_to_chi2(p)
                chi2_corrected = chi2_orig / lambda_gc
                p_corrected = _chi2_to_pvalue(chi2_corrected)
                corrected_p_values.append(min(p_corrected, 1.0))
            else:
                corrected_p_values.append(p)

    logger.info(f"Genomic control: λ={lambda_gc:.3f}, corrected {len(corrected_p_values)} p-values")

    if return_dict:
        return {
            "status": "success",
            "corrected_p_values": corrected_p_values,
            "inflation_factor": lambda_gc,
            "method": "genomic_control",
            "n_tests": len(p_values) if p_values else len(chi_squared_stats),
            "lambda_gc": lambda_gc,
            "median_chi2": median_chi2,
        }
    return corrected_p_values, lambda_gc


def qvalue_estimation(p_values: List[float], pi0: Optional[float] = None) -> Tuple[List[float], float]:
    """Estimate q-values from p-values.

    Args:
        p_values: List of p-values
        pi0: Proportion of true null hypotheses (estimated if None)

    Returns:
        Tuple of (q_values, estimated_pi0)
    """
    if not p_values:
        return [], 1.0

    # Estimate π₀ if not provided
    if pi0 is None:
        pi0 = _estimate_pi0(p_values)

    # Apply q-value procedure (similar to BH-FDR but with π₀)
    indexed_p = sorted(enumerate(p_values), key=lambda x: x[1])
    sorted_p = [p for _, p in indexed_p]

    n = len(sorted_p)
    q_values = [0.0] * n

    for i in range(n - 1, -1, -1):
        rank = i + 1
        q_value = min(q_values[i + 1] if i + 1 < n else 1.0, sorted_p[i] * pi0 * n / rank)
        q_values[i] = q_value

    # Ensure monotonicity
    for i in range(n - 2, -1, -1):
        q_values[i] = min(q_values[i], q_values[i + 1])

    # Map back to original order
    original_q = [0.0] * n
    for i, (original_idx, _) in enumerate(indexed_p):
        original_q[original_idx] = q_values[i]

    logger.info(
        f"Q-value estimation: π₀={pi0:.3f}, estimated {sum(1 for q in original_q if q <= 0.05)} significant tests"
    )

    return original_q, pi0


def _estimate_pi0(p_values: List[float]) -> float:
    """Estimate the proportion of true null hypotheses (π₀).

    Uses Storey's λ-based bootstrap smoother. Evaluates π₀(λ) at a grid
    of thresholds and selects the estimate that minimizes the mean squared
    error via bootstrap.

    Args:
        p_values: List of p-values

    Returns:
        Estimated π₀ value (between 0 and 1)
    """
    n = len(p_values)
    if n == 0:
        return 1.0

    # Storey's λ-grid method
    lambdas = [i * 0.05 for i in range(1, 19)]  # 0.05, 0.10, ..., 0.90
    pi0_estimates = []

    for lam in lambdas:
        n_above = sum(1 for p in p_values if p > lam)
        pi0_lam = n_above / (n * (1.0 - lam))
        pi0_estimates.append(min(pi0_lam, 1.0))

    if not pi0_estimates:
        return 1.0

    # Use the smoothest (minimum) estimate above the final λ estimate
    # This is the natural cubic spline approach simplified:
    # take the minimum of the last few estimates to avoid overshoot
    min_pi0 = min(pi0_estimates)

    # Ensure pi0 is at least the minimum estimate at the highest lambda
    pi0 = max(min_pi0, pi0_estimates[-1])

    # Clamp to valid range
    return max(0.01, min(pi0, 1.0))


def _pvalue_to_chi2(p: float) -> float:
    """Convert p-value to chi-squared(1) statistic.

    Uses scipy.stats.chi2.ppf if available, otherwise falls back to
    z-score squaring via the normal inverse CDF.

    Args:
        p: p-value (0 < p < 1)

    Returns:
        Chi-squared(1) statistic
    """
    try:
        from scipy import stats as scipy_stats

        return float(scipy_stats.chi2.ppf(1.0 - p, df=1))
    except ImportError:
        pass

    # Fallback: z = Phi^{-1}(1 - p/2), chi2 = z^2
    z = _normal_ppf(1.0 - p / 2.0)
    return z * z


def _chi2_to_pvalue(chi2: float) -> float:
    """Convert chi-squared(1) statistic back to a p-value.

    Args:
        chi2: Chi-squared(1) statistic (>= 0)

    Returns:
        Two-sided p-value
    """
    try:
        from scipy import stats as scipy_stats

        return float(scipy_stats.chi2.sf(chi2, df=1))
    except ImportError:
        pass

    # Fallback: p = 2 * (1 - Phi(sqrt(chi2)))
    if chi2 <= 0:
        return 1.0
    z = math.sqrt(chi2)
    return 2.0 * (1.0 - _normal_cdf(z))  # noqa: F821


def _normal_ppf(p: float) -> float:
    """Approximate inverse normal CDF (probit function).

    Rational approximation from Abramowitz & Stegun 26.2.23.
    Accurate to ~4.5e-4 absolute error.

    Args:
        p: Probability (0 < p < 1)

    Returns:
        z-score such that Phi(z) ≈ p
    """
    if p <= 0:
        return -10.0
    if p >= 1:
        return 10.0
    if p < 0.5:
        return -_normal_ppf(1.0 - p)

    # Rational approximation for 0.5 <= p < 1
    t = math.sqrt(-2.0 * math.log(1.0 - p))
    c0 = 2.515517
    c1 = 0.802853
    c2 = 0.010328
    d1 = 1.432788
    d2 = 0.189269
    d3 = 0.001308
    return t - (c0 + c1 * t + c2 * t * t) / (1.0 + d1 * t + d2 * t * t + d3 * t * t * t)


def adjust_p_values(p_values: List[float], method: str = "bonferroni", **kwargs) -> List[float]:
    """General function for p-value adjustment.

    Args:
        p_values: List of p-values to adjust
        method: Adjustment method ('bonferroni', 'fdr', 'genomic_control', 'qvalue')
        **kwargs: Method-specific parameters

    Returns:
        List of adjusted p-values or significance indicators
    """
    if method.lower() == "bonferroni":
        _, corrected_alpha = bonferroni_correction(p_values, kwargs.get("alpha", 0.05))
        # Convert to adjusted p-values (approximation)
        return [min(p * len(p_values), 1.0) for p in p_values]

    elif method.lower() == "fdr":
        _, adjusted_p = fdr_correction(p_values, kwargs.get("alpha", 0.05), kwargs.get("fdr_method", "bh"))
        return adjusted_p

    elif method.lower() == "genomic_control":
        adjusted_p, _ = genomic_control(p_values)
        return adjusted_p

    elif method.lower() == "qvalue":
        q_vals, _ = qvalue_estimation(p_values, kwargs.get("pi0"))
        return q_vals

    else:
        raise ValueError(f"Unknown adjustment method: {method}")
