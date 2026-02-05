"""Information-theoretic hypothesis testing for biological data.

This module implements statistical hypothesis testing methods based on
information theory, including permutation tests for mutual information,
independence tests, bootstrap confidence intervals for entropy, and
significance filtering with multiple testing correction.
"""

from __future__ import annotations

import math
from collections import Counter
from typing import Any, Dict, List, Sequence

import numpy as np

from metainformant.core import logging, validation

logger = logging.get_logger(__name__)

# Graceful import of syntactic module functions
try:
    from .syntactic import (
        mutual_information,
        shannon_entropy,
        shannon_entropy_from_counts,
        transfer_entropy,
    )

    HAS_SYNTACTIC = True
except ImportError:
    HAS_SYNTACTIC = False
    logger.warning("Could not import syntactic module; hypothesis tests will be unavailable")


def mi_permutation_test(
    x: Sequence[Any],
    y: Sequence[Any],
    n_permutations: int = 1000,
    seed: int = 42,
) -> Dict[str, Any]:
    """Permutation test for mutual information significance.

    Assesses whether the observed mutual information between two sequences
    is significantly greater than expected under the null hypothesis of
    independence. The null distribution is constructed by randomly shuffling
    the second sequence and recomputing MI for each permutation.

    Args:
        x: First sequence of observations
        y: Second sequence of observations (same length as x)
        n_permutations: Number of random permutations (default 1000)
        seed: Random seed for reproducibility

    Returns:
        Dictionary with keys:
            - observed_mi: The observed mutual information value
            - p_value: Proportion of permuted MI values >= observed MI
            - null_distribution: List of MI values from permutations
            - significant: Whether p_value < 0.05

    Raises:
        ValueError: If sequences have different lengths or n_permutations < 1
        RuntimeError: If syntactic module is not available
    """
    if not HAS_SYNTACTIC:
        raise RuntimeError("Syntactic module required for mi_permutation_test")

    validation.validate_type(x, (list, tuple), "x")
    validation.validate_type(y, (list, tuple), "y")

    if len(x) != len(y):
        raise ValueError("Sequences must have the same length")
    if n_permutations < 1:
        raise ValueError("n_permutations must be >= 1")

    rng = np.random.default_rng(seed)

    # Compute observed mutual information
    observed_mi = mutual_information(x, y)

    # Build null distribution by shuffling y
    y_array = np.array(y)
    null_distribution: List[float] = []

    for _ in range(n_permutations):
        y_permuted = y_array.copy()
        rng.shuffle(y_permuted)
        mi_perm = mutual_information(x, list(y_permuted))
        null_distribution.append(float(mi_perm))

    # p-value: proportion of null MI values >= observed MI
    n_extreme = sum(1 for mi_null in null_distribution if mi_null >= observed_mi)
    p_value = (n_extreme + 1) / (n_permutations + 1)

    logger.debug(
        "MI permutation test: observed_mi=%.4f, p_value=%.4f, n_permutations=%d",
        observed_mi,
        p_value,
        n_permutations,
    )

    return {
        "observed_mi": float(observed_mi),
        "p_value": float(p_value),
        "null_distribution": null_distribution,
        "significant": bool(p_value < 0.05),
    }


def independence_test(
    x: Sequence[Any],
    y: Sequence[Any],
    method: str = "chi_squared",
) -> Dict[str, Any]:
    """Test independence between two categorical sequences.

    Constructs a contingency table from the observed joint and marginal
    frequencies, then applies the specified statistical test to determine
    whether x and y are independent.

    Args:
        x: First sequence of categorical observations
        y: Second sequence of categorical observations (same length as x)
        method: Test method. Currently supports "chi_squared"

    Returns:
        Dictionary with keys:
            - statistic: The test statistic value
            - p_value: Approximate p-value from chi-squared distribution
            - method: The method used
            - dof: Degrees of freedom
            - significant: Whether p_value < 0.05

    Raises:
        ValueError: If sequences have different lengths or method is unknown
    """
    validation.validate_type(x, (list, tuple), "x")
    validation.validate_type(y, (list, tuple), "y")

    if len(x) != len(y):
        raise ValueError("Sequences must have the same length")

    if method != "chi_squared":
        raise ValueError(f"Unknown method: {method}. Supported: 'chi_squared'")

    n = len(x)

    # Build contingency table
    x_labels = sorted(set(x), key=str)
    y_labels = sorted(set(y), key=str)

    x_index = {label: i for i, label in enumerate(x_labels)}
    y_index = {label: i for i, label in enumerate(y_labels)}

    n_rows = len(x_labels)
    n_cols = len(y_labels)

    # Observed counts
    observed = np.zeros((n_rows, n_cols), dtype=float)
    for xi, yi in zip(x, y):
        observed[x_index[xi], y_index[yi]] += 1.0

    # Marginal totals
    row_totals = observed.sum(axis=1)
    col_totals = observed.sum(axis=0)

    # Expected counts under independence: E_ij = (row_i * col_j) / n
    expected = np.outer(row_totals, col_totals) / n

    # Chi-squared statistic: sum((O - E)^2 / E) for cells where E > 0
    chi2 = 0.0
    for i in range(n_rows):
        for j in range(n_cols):
            if expected[i, j] > 0:
                chi2 += (observed[i, j] - expected[i, j]) ** 2 / expected[i, j]

    # Degrees of freedom
    dof = (n_rows - 1) * (n_cols - 1)

    # Approximate p-value using the regularized incomplete gamma function
    # P(X >= chi2) = 1 - gamma_cdf(chi2, k/2, 2) = 1 - regularized_gamma(k/2, chi2/2)
    p_value = _chi2_survival(chi2, dof)

    logger.debug(
        "Independence test (%s): statistic=%.4f, dof=%d, p_value=%.4f",
        method,
        chi2,
        dof,
        p_value,
    )

    return {
        "statistic": float(chi2),
        "p_value": float(p_value),
        "method": method,
        "dof": int(dof),
        "significant": bool(p_value < 0.05),
    }


def entropy_confidence_interval(
    data: Sequence[Any],
    confidence: float = 0.95,
    n_bootstrap: int = 1000,
    seed: int = 42,
) -> Dict[str, Any]:
    """Bootstrap confidence interval for Shannon entropy estimate.

    Resamples the input data with replacement to construct an empirical
    distribution of entropy estimates, then derives the confidence interval
    from the quantiles of that distribution.

    Args:
        data: Sequence of observations (categorical values)
        confidence: Confidence level (default 0.95 for 95% CI)
        n_bootstrap: Number of bootstrap resamples
        seed: Random seed for reproducibility

    Returns:
        Dictionary with keys:
            - entropy: Point estimate of Shannon entropy (bits)
            - ci_lower: Lower bound of the confidence interval
            - ci_upper: Upper bound of the confidence interval
            - confidence: The confidence level used
            - std: Standard deviation of the bootstrap entropy distribution

    Raises:
        ValueError: If data is empty, confidence not in (0,1), or n_bootstrap < 1
        RuntimeError: If syntactic module is not available
    """
    if not HAS_SYNTACTIC:
        raise RuntimeError("Syntactic module required for entropy_confidence_interval")

    validation.validate_type(data, (list, tuple), "data")

    if len(data) == 0:
        raise ValueError("Data sequence must not be empty")
    if not (0.0 < confidence < 1.0):
        raise ValueError("Confidence must be between 0 and 1 (exclusive)")
    if n_bootstrap < 1:
        raise ValueError("n_bootstrap must be >= 1")

    rng = np.random.default_rng(seed)

    # Point estimate
    point_entropy = shannon_entropy_from_counts(Counter(data))

    # Bootstrap resampling
    data_array = np.array(data)
    n = len(data_array)
    bootstrap_entropies: List[float] = []

    for _ in range(n_bootstrap):
        indices = rng.integers(0, n, size=n)
        resample = data_array[indices]
        h = shannon_entropy_from_counts(Counter(resample.tolist()))
        bootstrap_entropies.append(float(h))

    # Compute confidence interval from quantiles
    alpha = 1.0 - confidence
    ci_lower = float(np.percentile(bootstrap_entropies, 100 * alpha / 2))
    ci_upper = float(np.percentile(bootstrap_entropies, 100 * (1 - alpha / 2)))
    std = float(np.std(bootstrap_entropies, ddof=1))

    logger.debug(
        "Entropy CI: H=%.4f, [%.4f, %.4f] (%.0f%% CI, std=%.4f)",
        point_entropy,
        ci_lower,
        ci_upper,
        confidence * 100,
        std,
    )

    return {
        "entropy": float(point_entropy),
        "ci_lower": ci_lower,
        "ci_upper": ci_upper,
        "confidence": confidence,
        "std": std,
    }


def information_significance_filter(
    variables: List[Sequence[Any]],
    target: Sequence[Any],
    alpha: float = 0.05,
    correction: str = "bonferroni",
) -> Dict[str, Any]:
    """Filter variables by mutual information significance with multiple testing correction.

    For each variable, performs a permutation test of MI against the target,
    then applies the specified multiple testing correction to control the
    family-wise error rate or false discovery rate.

    Args:
        variables: List of variable sequences to test against target
        target: Target sequence
        alpha: Significance threshold (default 0.05)
        correction: Multiple testing correction method.
            "bonferroni": Bonferroni correction (controls FWER)
            "fdr_bh": Benjamini-Hochberg procedure (controls FDR)
            "none": No correction

    Returns:
        Dictionary with keys:
            - significant_indices: List of indices of significant variables
            - p_values: List of raw p-values for each variable
            - adjusted_p_values: List of adjusted p-values
            - mi_values: List of MI values for each variable
            - alpha: The significance threshold used
            - correction: The correction method used

    Raises:
        ValueError: If variables list is empty, sequences have different lengths,
            or correction method is unknown
        RuntimeError: If syntactic module is not available
    """
    if not HAS_SYNTACTIC:
        raise RuntimeError("Syntactic module required for information_significance_filter")

    validation.validate_type(variables, list, "variables")
    validation.validate_type(target, (list, tuple), "target")

    if len(variables) == 0:
        raise ValueError("Variables list must not be empty")

    if correction not in ("bonferroni", "fdr_bh", "none"):
        raise ValueError(f"Unknown correction method: {correction}. Supported: 'bonferroni', 'fdr_bh', 'none'")

    n_vars = len(variables)
    p_values: List[float] = []
    mi_values: List[float] = []

    # Compute MI and p-value for each variable against target
    for i, var in enumerate(variables):
        result = mi_permutation_test(list(var), list(target), n_permutations=500, seed=42 + i)
        p_values.append(result["p_value"])
        mi_values.append(result["observed_mi"])

    # Apply multiple testing correction
    if correction == "bonferroni":
        adjusted_p_values = [min(p * n_vars, 1.0) for p in p_values]
    elif correction == "fdr_bh":
        adjusted_p_values = _benjamini_hochberg(p_values)
    else:
        adjusted_p_values = list(p_values)

    # Determine significant variables
    significant_indices = [i for i, adj_p in enumerate(adjusted_p_values) if adj_p < alpha]

    logger.debug(
        "Significance filter: %d/%d variables significant (alpha=%.2f, correction=%s)",
        len(significant_indices),
        n_vars,
        alpha,
        correction,
    )

    return {
        "significant_indices": significant_indices,
        "p_values": p_values,
        "adjusted_p_values": adjusted_p_values,
        "mi_values": mi_values,
        "alpha": alpha,
        "correction": correction,
    }


def entropy_rate_test(
    x: Sequence[Any],
    y: Sequence[Any],
    lag: int = 1,
    n_permutations: int = 500,
    seed: int = 42,
) -> Dict[str, Any]:
    """Test whether transfer entropy from x to y is significant.

    Uses a permutation test to assess whether the observed transfer entropy
    from sequence x to sequence y exceeds what would be expected by chance.
    The null distribution is constructed by shuffling x and recomputing
    transfer entropy for each permutation.

    Args:
        x: Source sequence
        y: Target sequence (same length as x)
        lag: Time lag for transfer entropy computation (default 1)
        n_permutations: Number of permutations for the null distribution
        seed: Random seed for reproducibility

    Returns:
        Dictionary with keys:
            - observed_te: The observed transfer entropy value
            - p_value: Proportion of permuted TE values >= observed TE
            - significant: Whether p_value < 0.05
            - lag: The lag value used

    Raises:
        ValueError: If sequences have different lengths, lag < 1, or sequences
            are too short for the given lag
        RuntimeError: If syntactic module is not available
    """
    if not HAS_SYNTACTIC:
        raise RuntimeError("Syntactic module required for entropy_rate_test")

    validation.validate_type(x, (list, tuple), "x")
    validation.validate_type(y, (list, tuple), "y")

    if len(x) != len(y):
        raise ValueError("Sequences must have the same length")
    if lag < 1:
        raise ValueError("Lag must be >= 1")
    if len(x) <= lag + 1:
        raise ValueError("Sequences too short for given lag")
    if n_permutations < 1:
        raise ValueError("n_permutations must be >= 1")

    rng = np.random.default_rng(seed)

    # Compute observed transfer entropy
    observed_te = transfer_entropy(x, y, lag=lag)

    # Build null distribution by shuffling x
    x_array = np.array(x)
    null_distribution: List[float] = []

    for _ in range(n_permutations):
        x_permuted = x_array.copy()
        rng.shuffle(x_permuted)
        te_perm = transfer_entropy(list(x_permuted), y, lag=lag)
        null_distribution.append(float(te_perm))

    # p-value: proportion of null TE values >= observed TE
    n_extreme = sum(1 for te_null in null_distribution if te_null >= observed_te)
    p_value = (n_extreme + 1) / (n_permutations + 1)

    logger.debug(
        "Transfer entropy test: observed_te=%.4f, p_value=%.4f, lag=%d",
        observed_te,
        p_value,
        lag,
    )

    return {
        "observed_te": float(observed_te),
        "p_value": float(p_value),
        "significant": bool(p_value < 0.05),
        "lag": int(lag),
    }


# ---------------------------------------------------------------------------
# Internal helper functions
# ---------------------------------------------------------------------------


def _chi2_survival(x: float, k: int) -> float:
    """Compute the survival function (1 - CDF) of chi-squared distribution.

    Uses the regularized incomplete gamma function to compute
    P(X >= x) for X ~ chi-squared(k).

    Args:
        x: The chi-squared statistic value
        k: Degrees of freedom

    Returns:
        p-value (survival probability)
    """
    if k <= 0:
        return 1.0
    if x <= 0:
        return 1.0

    # P(X >= x) = 1 - gamma_cdf(x, k) = 1 - regularized_lower_gamma(k/2, x/2)
    # which equals the regularized upper incomplete gamma function Q(k/2, x/2)
    a = k / 2.0
    z = x / 2.0

    return _regularized_upper_gamma(a, z)


def _regularized_upper_gamma(a: float, x: float) -> float:
    """Compute the regularized upper incomplete gamma function Q(a, x).

    Q(a, x) = 1 - P(a, x) where P is the regularized lower incomplete gamma.
    Uses series expansion for small x and continued fraction for large x.

    Args:
        a: Shape parameter (must be > 0)
        x: Integration lower bound (must be >= 0)

    Returns:
        Q(a, x) value in [0, 1]
    """
    if x < 0:
        return 1.0
    if x == 0:
        return 1.0

    if x < a + 1.0:
        # Use series expansion for P(a, x), then Q = 1 - P
        return 1.0 - _lower_gamma_series(a, x)
    else:
        # Use continued fraction for Q(a, x)
        return _upper_gamma_cf(a, x)


def _lower_gamma_series(a: float, x: float) -> float:
    """Regularized lower incomplete gamma via series expansion.

    P(a, x) = e^{-x} * x^a * sum_{n=0}^{inf} x^n / Gamma(a+n+1)

    Args:
        a: Shape parameter
        x: Upper integration limit

    Returns:
        P(a, x)
    """
    max_iter = 200
    eps = 1e-12

    ap = a
    term = 1.0 / a
    total = term

    for _ in range(1, max_iter):
        ap += 1.0
        term *= x / ap
        total += term
        if abs(term) < abs(total) * eps:
            break

    # P(a,x) = total * exp(-x + a*ln(x) - ln(Gamma(a)))
    log_val = -x + a * math.log(x) - math.lgamma(a)
    try:
        return total * math.exp(log_val)
    except OverflowError:
        return 1.0


def _upper_gamma_cf(a: float, x: float) -> float:
    """Regularized upper incomplete gamma via Lentz continued fraction.

    Args:
        a: Shape parameter
        x: Lower integration limit

    Returns:
        Q(a, x)
    """
    max_iter = 200
    eps = 1e-12
    tiny = 1e-30

    # Modified Lentz method for the continued fraction
    # Q(a,x) = exp(-x + a*ln(x) - lgamma(a)) * cf
    # where cf = 1/(x + 1-a - 1*( 1-a)/(x + 3-a - 2*(2-a)/(x + 5-a - ...)))

    b = x + 1.0 - a
    c = 1.0 / tiny
    d = 1.0 / b if abs(b) > tiny else 1.0 / tiny
    f = d

    for i in range(1, max_iter):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < tiny:
            d = tiny
        c = b + an / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        f *= delta
        if abs(delta - 1.0) < eps:
            break

    log_val = -x + a * math.log(x) - math.lgamma(a)
    try:
        return f * math.exp(log_val)
    except OverflowError:
        return 0.0


def _benjamini_hochberg(p_values: List[float]) -> List[float]:
    """Apply Benjamini-Hochberg FDR correction to a list of p-values.

    Args:
        p_values: List of raw p-values

    Returns:
        List of adjusted p-values (same order as input)
    """
    n = len(p_values)
    if n == 0:
        return []

    # Sort p-values and track original indices
    indexed = sorted(enumerate(p_values), key=lambda pair: pair[1])

    adjusted = [0.0] * n
    cumulative_min = 1.0

    # Traverse from largest to smallest p-value
    for rank_idx in range(n - 1, -1, -1):
        original_idx, p = indexed[rank_idx]
        rank = rank_idx + 1  # 1-based rank
        adjusted_p = p * n / rank
        cumulative_min = min(cumulative_min, adjusted_p)
        adjusted[original_idx] = min(cumulative_min, 1.0)

    return adjusted
