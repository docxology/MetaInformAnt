"""GWAS multiple testing correction utilities.

This module provides functions for correcting p-values in GWAS for multiple testing,
including Bonferroni, FDR, and genomic control methods.
"""

from __future__ import annotations

import math
from typing import Any, Dict, Iterable, List, Optional, Tuple, cast

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
    x = -0.70711 * ((2.30753 + t * 0.27061) / (1.0 + t * (0.99229 + t * 0.04481)) - t)
    for _ in range(2):
        err = math.erfc(x) - z
        x += err / (1.1283791670955126 * math.exp(-(x * x)) - x * err)
    return x if y < 1 else -x


def lambda_gc_from_p_values(p_values: Iterable[float]) -> Optional[float]:
    """Compute the genomic inflation factor λ_GC from association p-values (1 df).

    λ_GC = median(χ²(1) statistics derived from p) / median(χ²(1) under null),
    where the exact 1-df null median is ``EXPECTED_MEDIAN_CHI2_1DF``. This is the
    canonical conversion used across the GWAS domain; callers must not divide
    median p-values (or -log10 p-values) by a chi-square quantile directly.

    Args:
        p_values: Iterable of association p-values (values outside (0, 1] are ignored).

    Returns:
        λ_GC as a float, or None when no valid p-values are supplied.
    """
    chi2_stats = [_chi2_from_p_value(float(p)) for p in p_values if _valid_p_value(p)]
    if not chi2_stats:
        return None
    chi2_stats.sort()
    n = len(chi2_stats)
    median_chi2 = (
        chi2_stats[n // 2]
        if n % 2 == 1
        else (chi2_stats[n // 2 - 1] + chi2_stats[n // 2]) / 2.0
    )
    lambda_gc = median_chi2 / EXPECTED_MEDIAN_CHI2_1DF
    if not math.isfinite(lambda_gc) or lambda_gc <= 0:
        return None
    return lambda_gc


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

    logger.info(
        f"Bonferroni correction: {sum(significant)}/{n_tests} tests significant at α={alpha}"
    )

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
        raise ValueError(
            "Method must be 'bh' (Benjamini-Hochberg) or 'by' (Benjamini-Yekutieli)"
        )

    # Sort p-values and keep track of original indices
    indexed_p = sorted(enumerate(p_values), key=lambda x: x[1])
    sorted_p = [p for _, p in indexed_p]

    n = len(sorted_p)
    adjusted_p = [0.0] * n

    # Benjamini-Hochberg procedure
    for i in range(n - 1, -1, -1):
        rank = i + 1

        if method.lower() == "bh":
            adjusted_value = min(
                adjusted_p[i + 1] if i + 1 < n else 1.0, sorted_p[i] * n / rank
            )
        else:  # 'by' - Benjamini-Yekutieli
            c_n = sum(1.0 / (k + 1) for k in range(n))
            adjusted_value = min(
                adjusted_p[i + 1] if i + 1 < n else 1.0, sorted_p[i] * c_n * n / rank
            )

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

    logger.info(
        f"FDR correction ({method}): {sum(significant)}/{n} tests significant at FDR={alpha}"
    )

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
        Dictionary with correction results or tuple of (corrected_p_values, inflation_factor).
        ``corrected_p_values`` is aligned to the input list for both p-value
        and chi-square inputs (chi-square statistics are converted to 1-df
        upper-tail p-values, lambda-corrected, and converted back).
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
        chi_squared_stats = [
            _chi2_from_p_value(float(p)) for p in p_values if _valid_p_value(p)
        ]
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
    corrected_p_values: List[float] = []
    if p_values:
        for p in p_values:
            if _valid_p_value(p):
                # Correct the 1-df chi-square statistic and convert back to a
                # 1-df upper-tail p-value. Keep p=1 as 1.
                chi2_corrected = _chi2_from_p_value(float(p)) / correction_lambda
                corrected_p_values.append(_p_value_from_chi2(chi2_corrected))
            else:
                corrected_p_values.append(p)
    elif chi2_stats is not None:
        # Chi-square input: convert each statistic to a 1-df p-value, apply
        # the genomic-control lambda, and convert back. Output stays aligned
        # with the input list; uninterpretable entries become NaN.
        for stat in chi2_stats:
            try:
                chi2 = float(stat)
            except (TypeError, ValueError):
                corrected_p_values.append(float("nan"))
                continue
            if math.isfinite(chi2) and chi2 >= 0:
                corrected_p_values.append(_p_value_from_chi2(chi2 / correction_lambda))
            else:
                corrected_p_values.append(float("nan"))

    logger.info(
        f"Genomic control: λ={lambda_gc:.3f}, corrected {len(corrected_p_values)} p-values"
    )

    if return_dict:
        return {
            "status": "success",
            "corrected_p_values": corrected_p_values,
            "corrected_pvalues": corrected_p_values,
            "inflation_factor": lambda_gc,
            "method": "genomic_control",
            "n_tests": len(chi2_stats)
            if chi2_stats is not None
            else len(p_values or []),
            "lambda_gc": lambda_gc,
            "median_chi2": median_chi2,
            "expected_median_chi2": EXPECTED_MEDIAN_CHI2_1DF,
            "p_values_provided": p_values_provided,
        }
    return corrected_p_values, lambda_gc


def qvalue_estimation(
    p_values: List[float], pi0: Optional[float] = None
) -> Tuple[List[float], float]:
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
        q_value = min(
            q_values[i + 1] if i + 1 < n else 1.0, sorted_p[i] * pi0 * n / rank
        )
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


def _linear_fit_at_last_lambda(curve: List[Tuple[float, float]]) -> float:
    """Evaluate a weighted least-squares linear fit at the largest lambda.

    Fits y = a + b*x through the (lambda, pi0(lambda)) points via the
    weighted normal equations and evaluates the fit at the largest tuning
    parameter in the curve (interpolation at the final knot, following the
    ``qvalue`` smoother convention rather than extrapolating beyond the
    data). Points are weighted by ``w = (1 - lambda)^2``: raw tail-based
    pi0(lambda) estimates grow noisier as lambda approaches 1
    (variance ~ lambda / ((1 - lambda) * n)), so down-weighting them keeps
    the smoothed estimate stable for both flat-null and signal-mixture
    histograms. Returns NaN when the system is singular.
    """
    s = [0.0] * 3  # weighted sums of lambda^0..lambda^2
    t = [0.0] * 2  # weighted sums of pi0 * lambda^0..pi0 * lambda^1
    for x, y in curve:
        weight = (1.0 - x) ** 2
        power = 1.0
        for k in range(3):
            s[k] += weight * power
            if k < 2:
                t[k] += weight * y * power
            power *= x

    # Solve the 2x2 weighted normal system [[s0, s1], [s1, s2]] y = t by
    # elimination with partial pivoting.
    a00, a01, a11, b0, b1 = s[0], s[1], s[2], t[0], t[1]
    if abs(a01) > abs(a00):
        a00, a01, a11, b0, b1 = a01, a11, a00, b1, b0
    if abs(a00) < 1e-12:
        return float("nan")
    factor = a01 / a00
    a11p = a11 - factor * a01
    b1p = b1 - factor * b0
    if abs(a11p) < 1e-12:
        return float("nan")
    b = b1p / a11p
    a = (b0 - a01 * b) / a00
    x_last = curve[-1][0]
    return a + b * x_last


def _estimate_pi0(p_values: List[float]) -> float:
    """Estimate the proportion of true null hypotheses (pi0).

    Storey's lambda-tuning estimator over the p-value histogram: for tuning
    parameters lambda in [0.05, 0.95], pi0(lambda) = #(p > lambda) / ((1 -
    lambda) * n). A weighted linear least-squares smoother over the
    pi0(lambda) curve, evaluated at the largest lambda knot, damps the noise of
    the raw tail estimates. The result is bounded to (0, 1]; a flat null histogram
    (uniform p-values) yields pi0 ~ 1, under which the q-value procedure is
    exactly Benjamini-Hochberg.

    Args:
        p_values: List of p-values

    Returns:
        Estimated pi0 in (0, 1]
    """
    valid = [float(p) for p in p_values if _valid_p_value(p)]
    if not valid:
        return 1.0

    n = len(valid)
    lambdas = [0.05 * k for k in range(1, 20)]  # 0.05, 0.10, ..., 0.95
    curve: List[Tuple[float, float]] = []
    for lam in lambdas:
        tail = sum(1 for p in valid if p > lam)
        curve.append((lam, tail / ((1.0 - lam) * n)))

    pi0_fit = _linear_fit_at_last_lambda(curve)
    curve_min = min(value for _, value in curve)
    if math.isfinite(pi0_fit):
        # Keep the smoothed estimate within the observed curve's range and 1.
        pi0 = min(max(pi0_fit, curve_min), 1.0)
    else:
        pi0 = min(curve_min, 1.0)
    return max(pi0, 1.0 / (n + 1.0))


def adjust_p_values(
    p_values: List[float], method: str = "bonferroni", **kwargs: Any
) -> List[float]:
    """General function for p-value adjustment.

    Args:
        p_values: List of p-values to adjust
        method: Adjustment method ('bonferroni', 'fdr', 'genomic_control', 'qvalue')
        **kwargs: Method-specific parameters

    Returns:
        List of adjusted p-values or significance indicators
    """
    if method.lower() == "bonferroni":
        # Bonferroni-adjusted p-value (approximation): p * n, capped at 1.0
        return [min(p * len(p_values), 1.0) for p in p_values]

    elif method.lower() == "fdr":
        _, adjusted_p = cast(
            Tuple[List[bool], List[float]],
            fdr_correction(
                p_values,
                kwargs.get("alpha", 0.05),
                kwargs.get("fdr_method", "bh"),
                return_dict=False,
            ),
        )
        return adjusted_p

    elif method.lower() == "genomic_control":
        adjusted_p, _ = cast(
            Tuple[List[float], float], genomic_control(p_values, return_dict=False)
        )
        return adjusted_p

    elif method.lower() == "qvalue":
        q_vals, _ = qvalue_estimation(p_values, kwargs.get("pi0"))
        return q_vals

    else:
        raise ValueError(f"Unknown adjustment method: {method}")
