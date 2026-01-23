"""GWAS multiple testing correction utilities.

This module provides functions for correcting p-values in GWAS for multiple testing,
including Bonferroni, FDR, and genomic control methods.
"""

from __future__ import annotations

import math
from typing import List, Tuple

from metainformant.core import logging

logger = logging.get_logger(__name__)


def bonferroni_correction(p_values: List[float], alpha: float = 0.05) -> Tuple[List[bool], float]:
    """Apply Bonferroni correction for multiple testing.

    Args:
        p_values: List of p-values to correct
        alpha: Family-wise error rate

    Returns:
        Tuple of (significant_flags, corrected_alpha)
    """
    if not p_values:
        return [], alpha

    n_tests = len(p_values)
    corrected_alpha = alpha / n_tests

    significant = [p <= corrected_alpha for p in p_values]

    logger.info(f"Bonferroni correction: {sum(significant)}/{n_tests} tests significant at α={alpha}")

    return significant, corrected_alpha


def fdr_correction(p_values: List[float], alpha: float = 0.05, method: str = "bh") -> Tuple[List[bool], List[float]]:
    """Apply false discovery rate correction.

    Args:
        p_values: List of p-values to correct
        alpha: False discovery rate threshold
        method: Correction method ('bh' for Benjamini-Hochberg, 'by' for Benjamini-Yekutieli)

    Returns:
        Tuple of (significant_flags, adjusted_p_values)
    """
    if not p_values:
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

    return significant, original_adjusted


def genomic_control(p_values: List[float]) -> Tuple[List[float], float]:
    """Apply genomic control correction.

    Args:
        p_values: List of p-values to correct

    Returns:
        Tuple of (corrected_p_values, inflation_factor)
    """
    if not p_values:
        return [], 1.0

    # Calculate median chi-squared statistic
    chi_squared_stats = []
    for p in p_values:
        if p > 0 and p < 1:
            # Convert p-value to chi-squared statistic (1 df)
            chi2 = -2 * math.log(p)
            chi_squared_stats.append(chi2)

    if not chi_squared_stats:
        return p_values, 1.0

    # Calculate median chi-squared
    sorted_chi2 = sorted(chi_squared_stats)
    n = len(sorted_chi2)
    if n % 2 == 1:
        median_chi2 = sorted_chi2[n // 2]
    else:
        median_chi2 = (sorted_chi2[n // 2 - 1] + sorted_chi2[n // 2]) / 2

    # Genomic inflation factor λ
    lambda_gc = median_chi2 / 0.456  # Expected median for 1 df chi-squared

    # Apply correction
    corrected_p_values = []
    for p in p_values:
        if p > 0 and p < 1:
            # Correct chi-squared and convert back to p-value
            chi2_corrected = (-2 * math.log(p)) / lambda_gc
            p_corrected = math.exp(-chi2_corrected / 2)
            corrected_p_values.append(min(p_corrected, 1.0))  # Cap at 1.0
        else:
            corrected_p_values.append(p)

    logger.info(f"Genomic control: λ={lambda_gc:.3f}, corrected {len(corrected_p_values)} p-values")

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

    Args:
        p_values: List of p-values

    Returns:
        Estimated π₀ value
    """
    # Conservative approach: assume π₀ = 1 (all nulls true)
    # More sophisticated methods would use lambda tuning
    return 1.0


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
