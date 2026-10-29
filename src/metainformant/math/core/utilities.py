"""Mathematical utility functions.

This module provides general mathematical utility functions for statistical analysis.
"""

from __future__ import annotations

import math
from typing import List, Sequence

import numpy as np

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def correlation(x: Sequence[float] | np.ndarray, y: Sequence[float] | np.ndarray) -> float:
    """Calculate Pearson correlation coefficient.

    Args:
        x: First variable values
        y: Second variable values

    Returns:
        Pearson correlation coefficient
    """
    if len(x) != len(y):
        raise ValueError("Input lists must have the same length")

    if len(x) < 2:
        if len(x) == 1 and x[0] != y[0]:
            return 0.0
        raise ValueError("Input lists must contain at least 2 elements")

    # Convert to numpy arrays
    x_arr = np.asarray(x)
    y_arr = np.asarray(y)

    if np.var(x_arr) == 0 or np.var(y_arr) == 0:
        return 0.0

    return float(np.corrcoef(x_arr, y_arr)[0, 1])


def correlation_coefficient(x: List[float], y: List[float]) -> float:
    """Calculate Pearson correlation coefficient (alias for correlation).

    Args:
        x: First variable values
        y: Second variable values

    Returns:
        Pearson correlation coefficient
    """
    return correlation(x, y)


def linear_regression(x: Sequence[float] | np.ndarray, y: Sequence[float] | np.ndarray) -> tuple[float, float, float]:
    """Perform linear regression and return slope, intercept, and r-squared.

    Args:
        x: Independent variable values
        y: Dependent variable values

    Returns:
        Tuple of (slope, intercept, r_squared)
    """
    if len(x) != len(y):
        raise ValueError("Input lists must have the same length")

    if len(x) < 2:
        raise ValueError("Input lists must contain at least 2 elements")

    # Convert to numpy arrays
    x_arr = np.asarray(x)
    y_arr = np.asarray(y)

    if np.var(x_arr) == 0:
        return 0.0, 0.0, 0.0

    # Perform linear regression
    slope, intercept = np.polyfit(x_arr, y_arr, 1)

    # Calculate R-squared
    y_pred = slope * x_arr + intercept
    ss_res = np.sum((y_arr - y_pred) ** 2)
    ss_tot = np.sum((y_arr - np.mean(y_arr)) ** 2)

    if ss_tot == 0:
        r_squared = 1.0 if ss_res == 0 else 0.0
    else:
        r_squared = 1 - (ss_res / ss_tot)

    return slope, intercept, r_squared


def r_squared(
    x_or_pA: Sequence[float] | np.ndarray | float,
    y_or_pa: Sequence[float] | np.ndarray | float | None = None,
    pB: float | None = None,
    pb: float | None = None,
    pAB: float | None = None,
) -> float:
    """Calculate R-squared for linear regression or LD.

    Can be called with:
    - Two lists: r_squared(x, y) for linear regression R-squared
    - Allele frequencies: r_squared(pA, pa, pB, pb, pAB) for LD r²

    Args:
        x_or_pA: Either x values for regression or pA allele frequency
        y_or_pa: Either y values for regression or pa allele frequency
        pB: Allele frequency for B (if using LD mode)
        pb: Allele frequency for b (if using LD mode)
        pAB: Haplotype frequency for AB (if using LD mode)

    Returns:
        R-squared value

    Raises:
        TypeError: If the first two arguments are not numbers in LD mode, or not lists (tuples or
            arrays) in regression mode
        ValueError: If ``y`` is missing in regression mode, or the regression inputs are unusable
    """
    # Check if we're in LD mode (allele frequencies)
    if pB is not None and pb is not None and pAB is not None:
        for name, value in (("pA", x_or_pA), ("pa", y_or_pa)):
            if isinstance(value, bool) or not isinstance(value, (int, float, np.integer, np.floating)):
                raise TypeError(f"{name} must be an allele frequency (a number) in LD mode, got {type(value).__name__}")
        pA = float(x_or_pA)  # type: ignore[arg-type]
        pa = float(y_or_pa)  # type: ignore[arg-type]
        # Calculate D = P(AB) - P(A)*P(B)
        D = pAB - pA * pB
        # r² = D² / (pA * pa * pB * pb)
        denom = pA * pa * pB * pb
        if denom == 0:
            return 0.0
        return (D * D) / denom

    # Linear regression mode
    x = x_or_pA
    y = y_or_pa
    if y is None:
        raise ValueError("y is required for linear regression mode")
    for name, value in (("x", x), ("y", y)):
        if not isinstance(value, (list, tuple, np.ndarray)):
            raise TypeError(f"{name} must be a list, tuple or array in regression mode, got {type(value).__name__}")
    _, _, r2 = linear_regression(x, y)  # type: ignore[arg-type]
    return r2


def fisher_exact_test(a: int, b: int, c: int, d: int) -> tuple[float, float]:
    """Perform Fisher's exact test on a 2x2 contingency table.

    Fisher's exact test is used to determine if there is a significant association
    between two categorical variables in a 2x2 contingency table.

    Args:
        a: Count in cell (1,1)
        b: Count in cell (1,2)
        c: Count in cell (2,1)
        d: Count in cell (2,2)

    Returns:
        Tuple of (odds_ratio, p_value). The odds ratio is the sample odds ratio ``a*d / (b*c)``
        (infinite when ``b*c`` is 0) and the p-value is the exact two-sided conditional
        (hypergeometric) probability from ``scipy.stats.fisher_exact``. A table with an empty
        row or column gives ``(nan, 1.0)``.

    Raises:
        ImportError: If scipy, a required dependency, is not installed. There is deliberately no
            approximate fallback: a p-value is either the exact one or not returned at all.

    Examples:
        >>> # Test association between treatment and outcome
        >>> odds_ratio, p_val = fisher_exact_test(8, 2, 1, 9)
        >>> print(f"p-value: {p_val:.4f}, odds ratio: {odds_ratio:.2f}")
    """
    try:
        from scipy.stats import fisher_exact
    except ImportError as exc:
        raise ImportError(
            "scipy is required for fisher_exact_test() and is a core metainformant dependency; "
            "install it with: uv pip install scipy"
        ) from exc

    odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative="two-sided")
    return float(odds_ratio), float(p_value)


def covariance(x: List[float], y: List[float]) -> float:
    """Calculate the population covariance between two variables.

    The denominator is the number of pairs n, not n - 1 (``numpy.cov(..., ddof=0)``), which matches
    the population variance of :func:`metainformant.math.population_genetics.statistics.variance`.

    Args:
        x: First variable values
        y: Second variable values

    Returns:
        Population covariance; 0.0 for fewer than two pairs

    Examples:
        >>> x = [1, 2, 3, 4, 5]
        >>> y = [2, 4, 6, 8, 10]
        >>> cov = covariance(x, y)
        >>> print(f"Covariance: {cov}")
        Covariance: 4.0
    """
    if len(x) != len(y):
        raise ValueError("Input lists must have equal length")

    if len(x) < 2:
        return 0.0

    x_arr = np.array(x)
    y_arr = np.array(y)

    return float(np.cov(x_arr, y_arr, ddof=0)[0, 1])


def shannon_entropy(values: List[float], base: float = 2.0) -> float:
    """Calculate Shannon entropy from a list of probability values or counts.

    Args:
        values: List of probability values (will be normalized if they don't sum to 1)
               or frequency counts
        base: Logarithm base (2.0 for bits, math.e for nats)

    Returns:
        Shannon entropy value

    Examples:
        >>> shannon_entropy([0.25, 0.25, 0.25, 0.25])  # Uniform distribution
        2.0
        >>> shannon_entropy([1.0, 0.0, 0.0])  # Deterministic
        0.0
    """
    if not values:
        return 0.0

    # Convert to numpy array
    values_arr = np.array(values, dtype=float)

    # Handle zero or negative values
    if np.any(values_arr < 0):
        raise ValueError("All values must be non-negative")

    # Normalize if not already a probability distribution
    total = np.sum(values_arr)
    if total > 0:
        probs = values_arr / total
    else:
        return 0.0

    # Remove zero probabilities to avoid log(0)
    probs = probs[probs > 0]

    if len(probs) == 0:
        return 0.0

    # Calculate entropy
    entropy = -np.sum(probs * np.log(probs))

    # Convert to specified base
    if base != math.e:
        entropy = entropy / math.log(base)

    return float(entropy)


def jensen_shannon_divergence(p: List[float], q: List[float], base: float = 2.0) -> float:
    """Calculate Jensen-Shannon divergence between two probability distributions.

    The Jensen-Shannon divergence is a symmetrized version of the Kullback-Leibler
    divergence, bounded between 0 and 1 (in bits).

    Args:
        p: First probability distribution
        q: Second probability distribution
        base: Logarithm base (2.0 for bits, math.e for nats)

    Returns:
        Jensen-Shannon divergence value

    Examples:
        >>> p = [0.5, 0.5]
        >>> q = [0.1, 0.9]
        >>> jsd = jensen_shannon_divergence(p, q)
        >>> print(f"JSD: {jsd:.3f}")
    """
    if len(p) != len(q):
        raise ValueError("Probability distributions must have the same length")

    # Convert to numpy arrays
    p_arr = np.array(p, dtype=float)
    q_arr = np.array(q, dtype=float)

    # Normalize distributions
    p_arr = p_arr / np.sum(p_arr) if np.sum(p_arr) > 0 else p_arr
    q_arr = q_arr / np.sum(q_arr) if np.sum(q_arr) > 0 else q_arr

    # Calculate mixture distribution
    m = 0.5 * (p_arr + q_arr)

    # Calculate KL divergences
    kl_pm = _kl_divergence(p_arr, m, base)
    kl_qm = _kl_divergence(q_arr, m, base)

    # Jensen-Shannon divergence
    jsd = 0.5 * (kl_pm + kl_qm)

    return float(jsd)


def _kl_divergence(p: np.ndarray, q: np.ndarray, base: float = 2.0) -> float:
    """Calculate Kullback-Leibler divergence between two distributions.

    Args:
        p: First distribution
        q: Second distribution
        base: Logarithm base

    Returns:
        KL divergence value
    """
    # Avoid division by zero and log(0)
    mask = (p > 0) & (q > 0)
    if not np.any(mask):
        return 0.0

    p_masked = p[mask]
    q_masked = q[mask]

    kl = np.sum(p_masked * np.log(p_masked / q_masked))

    # Convert to specified base
    if base != math.e:
        kl = kl / math.log(base)

    return float(kl)
