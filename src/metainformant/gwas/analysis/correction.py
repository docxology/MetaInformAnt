"""GWAS multiple testing correction utilities.

This module provides functions for correcting p-values in GWAS for multiple testing,
including Bonferroni, FDR, and genomic control methods.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Tuple

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

try:
    from scipy import stats as _scipy_stats

    HAS_SCIPY = True
except ImportError:  # pragma: no cover - exercised only in lean environments
    _scipy_stats = None
    HAS_SCIPY = False

EXPECTED_MEDIAN_CHI2_1DF = 0.454936423119572


def _valid_p_value(p_value: Any) -> bool:
    """Return True for finite p-values in the closed interval (0, 1]."""
    try:
        p = float(p_value)
    except (TypeError, ValueError):
        return False
    return math.isfinite(p) and 0 < p <= 1


def _chi2_from_p_value(p_value: float) -> float:
    """Convert a two-sided 1-df association p-value to a chi-square statistic."""
    p = min(max(float(p_value), 1e-300), 1.0)
    if HAS_SCIPY and _scipy_stats is not None:
        chi2 = float(_scipy_stats.chi2.isf(p, 1))
    else:
        # Wilson-Hilferty approximation to the 1-df chi-square inverse survival.
        # This is a fallback only; scipy is preferred for production GWAS runs.
        z = math.sqrt(2.0) * _erfcinv_approx(p)
        chi2 = z * z
    if not math.isfinite(chi2):
        return 0.0 if p >= 1.0 else 1e6
    return max(chi2, 0.0)


def _p_value_from_chi2(chi2_stat: float) -> float:
    """Convert a 1-df chi-square statistic back to an upper-tail p-value."""
    chi2 = max(float(chi2_stat), 0.0)
    if HAS_SCIPY and _scipy_stats is not None:
        p_value = float(_scipy_stats.chi2.sf(chi2, 1))
    else:
        # Exact 1-df survival expressed through erfc, using the stdlib fallback.
        p_value = math.erfc(math.sqrt(chi2 / 2.0))
    if not math.isfinite(p_value):
        return 0.0
    return min(max(p_value, 0.0), 1.0)


def _erfcinv_approx(y: float) -> float:
    """Approximate erfc inverse for scipy-free environments."""
    # Mike Giles-style approximation, accurate enough for fallback plotting/QC.
    if y <= 0:
        return float("inf")
    if y >= 2:
        return float("-inf")
    z = y if y < 1 else 2 - y
    t = math.sqrt(-2.0 * math.log(z / 2.0))
    x = -0.70711 * (
        (2.30753 + t * 0.27061) / (1.0 + t * (0.99229 + t * 0.04481)) - t
    )
    for _ in range(2):
        err = math.erfc(x) - z
        x += err / (1.1283791670955126 * math.exp(-(x * x)) - x * err)
    return x if y < 1 else -x


def bonferroni_correction(
    p_values: List[float], alpha: float = 0.05, return_dict: bool = True
) -> Dict[str, Any] | Tuple[List[bool], float]:
    """Apply Bonferroni correction for multiple testing.

    Args:
        p_values: List of p-values to correct
        alpha: Family-wise error rate
        return_dict: If True, return dict (default); if False, return tuple (legacy)

    Returns:
        Dictionary with correction results or tuple of (significant_flags, corrected_alpha)
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
    p_values: List[float], alpha: float = 0.05, method: str = "bh", return_dict: bool = True
) -> Dict[str, Any] | Tuple[List[bool], List[float]]:
    """Apply false discovery rate correction.

    Args:
        p_values: List of p-values to correct
        alpha: False discovery rate threshold
        method: Correction method ('bh' for Benjamini-Hochberg, 'by' for Benjamini-Yekutieli)
        return_dict: If True, return dict (default); if False, return tuple (legacy)

    Returns:
        Dictionary with correction results or tuple of (significant_flags, adjusted_p_values)
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

    Args:
        p_values: List of p-values to correct
        chi2_stats: List of chi-squared statistics (alternative to p_values)
        return_dict: If True, return dict (default); if False, return tuple (legacy)
        pvalues: Alias for p_values (backward compatibility)

    Returns:
        Dictionary with correction results or tuple of (corrected_p_values, inflation_factor)
    """
    # Handle parameter aliases
    if pvalues is not None and p_values is None:
        p_values = pvalues

    # If chi2_stats provided, use those directly
    if chi2_stats is not None:
        chi_squared_stats = []
        for stat in chi2_stats:
            try:
                chi2 = float(stat)
            except (TypeError, ValueError):
                continue
            if math.isfinite(chi2) and chi2 >= 0:
                chi_squared_stats.append(chi2)
        p_values_provided = False
    elif p_values:
        p_values_provided = True
        # Convert association p-values to 1-df chi-square statistics. GWAS
        # lambda GC is based on the chi-square survival distribution, not the
        # -2 log(p) transform used by Fisher's method.
        chi_squared_stats = [_chi2_from_p_value(float(p)) for p in p_values if _valid_p_value(p)]
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

    # Genomic inflation factor λ. The exact 1-df null median is
    # scipy.stats.chi2.ppf(0.5, 1) = 0.4549364231...
    lambda_gc = median_chi2 / EXPECTED_MEDIAN_CHI2_1DF
    correction_lambda = lambda_gc if math.isfinite(lambda_gc) and lambda_gc > 0 else 1.0

    # Apply correction
    corrected_p_values = []
    if p_values:
        for p in p_values:
            if _valid_p_value(p):
                # Correct the 1-df chi-square statistic and convert back to a
                # 1-df upper-tail p-value. Keep p=1 as 1.
                chi2_corrected = _chi2_from_p_value(float(p)) / correction_lambda
                corrected_p_values.append(_p_value_from_chi2(chi2_corrected))
            else:
                corrected_p_values.append(p)

    logger.info(f"Genomic control: λ={lambda_gc:.3f}, corrected {len(corrected_p_values)} p-values")

    if return_dict:
        return {
            "status": "success",
            "corrected_p_values": corrected_p_values,
            "corrected_pvalues": corrected_p_values,
            "inflation_factor": lambda_gc,
            "method": "genomic_control",
            "n_tests": len(p_values) if p_values else len(chi_squared_stats),
            "lambda_gc": lambda_gc,
            "median_chi2": median_chi2,
            "expected_median_chi2": EXPECTED_MEDIAN_CHI2_1DF,
            "p_values_provided": p_values_provided,
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
        bonferroni_correction(p_values, kwargs.get("alpha", 0.05), return_dict=False)
        # Convert to adjusted p-values (approximation)
        return [min(p * len(p_values), 1.0) for p in p_values]

    elif method.lower() == "fdr":
        _, adjusted_p = fdr_correction(
            p_values,
            kwargs.get("alpha", 0.05),
            kwargs.get("fdr_method", "bh"),
            return_dict=False,
        )
        return adjusted_p

    elif method.lower() == "genomic_control":
        adjusted_p, _ = genomic_control(p_values, return_dict=False)
        return adjusted_p

    elif method.lower() == "qvalue":
        q_vals, _ = qvalue_estimation(p_values, kwargs.get("pi0"))
        return q_vals

    else:
        raise ValueError(f"Unknown adjustment method: {method}")
