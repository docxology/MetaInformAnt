"""Bias-corrected entropy and information estimation methods.

This module implements various estimation methods for information-theoretic
quantities with bias correction to improve accuracy on finite samples.
"""

from __future__ import annotations

import math
from collections import Counter
from typing import Any, Dict, List, Optional, Union

import numpy as np

from metainformant.core.data import validation
from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def entropy_estimator(
    counts: Union[Dict[Any, int], List[int]], method: str = "plugin", bias_correction: bool = True
) -> float:
    """Estimate Shannon entropy with various methods and bias correction.

    Args:
        counts: Either dict mapping items to counts or list of counts
        method: Estimation method ('plugin', 'miller_madow', 'chao_shen', 'jackknife')
        bias_correction: Whether to apply bias correction

    Returns:
        Entropy estimate in bits

    Raises:
        ValueError: If invalid method or insufficient data
    """
    if isinstance(counts, dict):
        count_values = list(counts.values())
    else:
        count_values = list(counts)

    validation.validate_type(count_values, list, "counts")

    count_array = np.array(count_values, dtype=int)
    if np.any(count_array < 0):
        raise ValueError("Counts cannot be negative")

    total = np.sum(count_array)
    if total == 0:
        return 0.0

    if method == "plugin":
        return _plugin_entropy_estimator(count_array, total, bias_correction)
    elif method == "miller_madow":
        return _miller_madow_entropy_estimator(count_array, total)
    elif method == "chao_shen":
        return _chao_shen_entropy_estimator(count_array, total)
    elif method == "jackknife":
        return _jackknife_entropy_estimator(count_array, total)
    else:
        raise ValueError(f"Unknown entropy estimation method: {method}")


def _plugin_entropy_estimator(counts: np.ndarray, total: int, bias_correction: bool) -> float:
    """Plugin (maximum likelihood) entropy estimator."""
    # Convert to probabilities
    probs = counts / total
    probs = probs[probs > 0]  # Remove zeros

    if len(probs) == 0:
        return 0.0

    # Plugin entropy: -sum(p * log2(p))
    entropy = -np.sum(probs * np.log2(probs))

    # Bias correction (Miller, 1955): subtract (k-1)/(2n)
    if bias_correction and total > 1:
        k = len(probs)  # Number of non-zero categories
        correction = (k - 1) / (2 * total)
        entropy -= correction

    return max(0.0, entropy)


def _miller_madow_entropy_estimator(counts: np.ndarray, total: int) -> float:
    """Miller-Madow entropy estimator (bias-corrected plugin)."""
    # Convert to probabilities
    probs = counts / total
    probs = probs[probs > 0]  # Remove zeros

    if len(probs) == 0:
        return 0.0

    # Miller-Madow correction: subtract (k-1)/(2n) where k is number of non-zero categories
    k = len(probs)
    correction = (k - 1) / (2 * total) if total > 1 else 0

    entropy = -np.sum(probs * np.log2(probs)) - correction

    return max(0.0, entropy)


def _chao_shen_entropy_estimator(counts: np.ndarray, total: int) -> float:
    """Chao-Shen entropy estimator for sparse data."""
    # Sort counts in descending order
    sorted_counts = np.sort(counts)[::-1]
    k = len(sorted_counts)

    if k == 0:
        return 0.0

    # Chao-Shen estimator components
    C = 1 - (np.sum(sorted_counts == 1) / total) if total > 0 else 0

    if C == 0:
        # No singleton counts, use plugin estimator
        probs = sorted_counts / total
        return -np.sum(probs * np.log2(probs))

    # Weighted sum for Chao-Shen
    weights = np.array([1 / (i * (i - 1)) for i in range(2, k + 2)])

    # Calculate lambda terms
    lambda_terms = []
    for i in range(k):
        if sorted_counts[i] > 1:
            term = (sorted_counts[i] / total) * np.log2(sorted_counts[i] / total)
            lambda_terms.append(term)

    if lambda_terms:
        lambda_sum = np.sum(lambda_terms)
    else:
        lambda_sum = 0.0

    # Chao-Shen formula
    entropy = -lambda_sum + (np.sum(sorted_counts == 1) / total) * np.log2(total) * C

    return max(0.0, entropy)


def _jackknife_entropy_estimator(counts: np.ndarray, total: int) -> float:
    """Jackknife entropy estimator for bias reduction."""
    counts = counts[counts > 0]  # Remove zeros
    k = len(counts)

    if k <= 1:
        return 0.0

    # Plugin entropy
    h_plugin = entropy_estimator(counts, method="plugin", bias_correction=False)

    # Jackknife: remove one category at a time
    h_jackknife_terms = []

    for i in range(k):
        # Remove i-th category
        reduced_counts = np.delete(counts, i)
        reduced_total = total - counts[i]

        if reduced_total > 0 and len(reduced_counts) > 0:
            reduced_probs = reduced_counts / reduced_total
            h_reduced = -np.sum(reduced_probs * np.log2(reduced_probs))
            h_jackknife_terms.append(h_reduced)

    if not h_jackknife_terms:
        return h_plugin

    # Jackknife estimate: k * h_plugin - ((k-1)/k) * sum(h_{-i})
    h_avg_reduced = np.mean(h_jackknife_terms)
    h_jackknife = k * h_plugin - (k - 1) * h_avg_reduced

    return max(0.0, h_jackknife)


def mutual_information_estimator(
    x: List[Any], y: List[Any], method: str = "plugin", bias_correction: bool = True
) -> float:
    """Estimate mutual information I(X;Y) with bias correction.

    Args:
        x: Samples from first variable
        y: Samples from second variable
        method: Estimation method ('plugin', 'miller_madow')
        bias_correction: Whether to apply bias correction

    Returns:
        Mutual information estimate

    Raises:
        ValueError: If sequences have different lengths
    """
    validation.validate_type(x, list, "x")
    validation.validate_type(y, list, "y")

    if len(x) != len(y):
        raise ValueError("Sequences must have the same length")

    # Convert to count dictionaries
    # Joint counts
    joint_counts = Counter(zip(x, y))

    # Marginal counts
    x_counts = Counter(x)
    y_counts = Counter(y)

    total = len(x)

    # I(X;Y) = H(X) + H(Y) - H(X,Y)
    h_x = entropy_estimator(x_counts, method=method, bias_correction=bias_correction)
    h_y = entropy_estimator(y_counts, method=method, bias_correction=bias_correction)
    h_xy = entropy_estimator(joint_counts, method=method, bias_correction=bias_correction)

    mi = h_x + h_y - h_xy
    return max(0.0, mi)  # Ensure non-negative


def kl_divergence_estimator(p: List[Any], q: List[Any], method: str = "plugin", bias_correction: bool = True) -> float:
    """Estimate KL divergence D_KL(P||Q) with bias correction.

    Args:
        p: Samples from distribution P
        q: Samples from distribution Q
        method: Estimation method
        bias_correction: Whether to apply bias correction

    Returns:
        KL divergence estimate

    Raises:
        ValueError: If sample lists have different lengths
    """
    if len(p) != len(q):
        raise ValueError("Sample lists must have the same length")

    # Convert to probability distributions using counts
    p_counts = Counter(p)
    q_counts = Counter(q)

    total = len(p)

    # Convert to probabilities
    p_probs = {item: count / total for item, count in p_counts.items()}
    q_probs = {item: count / total for item, count in q_counts.items()}

    # All possible items
    all_items = set(p_probs.keys()) | set(q_probs.keys())

    # KL divergence: sum(p * log(p/q))
    kl_div = 0.0
    for item in all_items:
        p_prob = p_probs.get(item, 0.0)
        q_prob = q_probs.get(item, 0.0)

        if p_prob > 0:
            if q_prob > 0:
                kl_div += p_prob * math.log2(p_prob / q_prob)
            else:
                return float("inf")  # Infinite divergence

    return max(0.0, kl_div)


def bias_correction(entropy: float, sample_size: int, alphabet_size: int) -> float:
    """Apply general bias correction to entropy estimate.

    Args:
        entropy: Raw entropy estimate
        sample_size: Number of samples (n)
        alphabet_size: Size of alphabet (d)

    Returns:
        Bias-corrected entropy

    Raises:
        ValueError: If parameters are invalid
    """
    if sample_size <= 0:
        raise ValueError("Sample size must be positive")
    if alphabet_size <= 0:
        raise ValueError("Alphabet size must be positive")

    # General bias correction: entropy - (d-1)/(2n)
    correction = (alphabet_size - 1) / (2 * sample_size)
    corrected_entropy = entropy - correction

    return max(0.0, corrected_entropy)


def entropy_bootstrap_confidence(
    counts: Union[Dict[Any, int], List[int]],
    method: str = "plugin",
    n_bootstraps: int = 1000,
    confidence_level: float = 0.95,
    random_state: Optional[int] = None,
) -> Dict[str, float]:
    """Calculate bootstrap confidence interval for entropy estimate.

    Args:
        counts: Count data
        method: Estimation method
        n_bootstraps: Number of bootstrap samples
        confidence_level: Confidence level (0-1)
        random_state: Random state for reproducibility

    Returns:
        Dictionary with entropy estimate and confidence interval

    Raises:
        ValueError: If invalid parameters
    """
    if not (0 < confidence_level < 1):
        raise ValueError("Confidence level must be between 0 and 1")

    if isinstance(counts, dict):
        items = []
        weights = []
        for item, count in counts.items():
            items.extend([item] * count)
            weights.append(count)
    else:
        items = []
        for i, count in enumerate(counts):
            items.extend([i] * count)

    if not items:
        return {
            "entropy": 0.0,
            "ci_lower": 0.0,
            "ci_upper": 0.0,
            "confidence_level": confidence_level,
        }

    np.random.seed(random_state)
    bootstrap_entropies = []

    # Generate bootstrap samples
    for _ in range(n_bootstraps):
        # Bootstrap resampling
        bootstrap_sample = np.random.choice(items, size=len(items), replace=True)

        # Convert back to counts
        bootstrap_counts = Counter(bootstrap_sample)

        # Estimate entropy
        entropy_est = entropy_estimator(bootstrap_counts, method=method, bias_correction=True)
        bootstrap_entropies.append(entropy_est)

    bootstrap_entropies = np.array(bootstrap_entropies)

    # Calculate confidence interval
    alpha = 1 - confidence_level
    ci_lower = np.percentile(bootstrap_entropies, alpha / 2 * 100)
    ci_upper = np.percentile(bootstrap_entropies, (1 - alpha / 2) * 100)

    # Main estimate (using original data)
    if isinstance(counts, dict):
        main_counts = counts
    else:
        main_counts = Counter(items)

    main_entropy = entropy_estimator(main_counts, method=method, bias_correction=True)

    return {
        "entropy": main_entropy,
        "ci_lower": ci_lower,
        "ci_upper": ci_upper,
        "confidence_level": confidence_level,
        "n_bootstraps": n_bootstraps,
    }


def effective_sample_size_correction(entropy: float, sample_size: int, alphabet_size: int) -> float:
    """Apply effective sample size correction for entropy estimation.

    This correction accounts for the fact that the effective sample size
    may be smaller than the nominal sample size due to dependencies.

    Args:
        entropy: Raw entropy estimate
        sample_size: Nominal sample size
        alphabet_size: Alphabet size

    Returns:
        Corrected entropy estimate
    """
    if sample_size <= 1:
        return entropy

    # Correction factor based on alphabet size and sample size
    # This is a heuristic correction
    correction_factor = 1.0 - (alphabet_size - 1) / (2 * sample_size * math.log(2))

    corrected_entropy = entropy * correction_factor

    return max(0.0, corrected_entropy)


def panzeri_treves_bias_correction(
    entropy: float, sample_size: int, alphabet_size: int, response_frequencies: Optional[np.ndarray] = None
) -> float:
    """Apply Panzeri-Treves bias correction for entropy estimation.

    This is a more sophisticated bias correction that accounts for
    response frequencies and provides better small-sample performance.

    Args:
        entropy: Raw entropy estimate
        sample_size: Sample size (n)
        alphabet_size: Number of possible responses (d)
        response_frequencies: Array of response frequencies (optional)

    Returns:
        Bias-corrected entropy

    References:
        Panzeri & Treves (1996). Analytical estimates of limited sampling biases in different information measures.
    """
    if sample_size <= 1:
        return entropy

    if response_frequencies is not None:
        # Use provided frequencies
        freq_array = np.array(response_frequencies)
    else:
        # Assume uniform distribution as fallback
        freq_array = np.ones(alphabet_size) / alphabet_size

    # Panzeri-Treves correction
    # Sum over responses: (f_r - 1/n) / (n * log(2))
    correction_terms = []
    for freq in freq_array:
        if freq > 0:
            prob = freq / sample_size
            if prob < 1.0:  # Avoid log(0)
                term = (freq * (1 - prob)) / (sample_size * math.log(2))
                correction_terms.append(term)

    correction = sum(correction_terms) if correction_terms else 0.0

    corrected_entropy = entropy - correction

    return max(0.0, corrected_entropy)


def entropy_rate_estimator(sequence: List[Any], order: int = 1, method: str = "plugin") -> float:
    """Estimate entropy rate of a sequence.

    The entropy rate is the limit of n-block entropy divided by n as n→∞.
    For Markov chains, this equals the conditional entropy H(X_{n+1}|X_n).

    Args:
        sequence: Input sequence
        order: Markov order (1 for first-order Markov)
        method: Entropy estimation method

    Returns:
        Entropy rate estimate

    Raises:
        ValueError: If sequence is too short or order is invalid
    """
    validation.validate_type(sequence, list, "sequence")

    n = len(sequence)
    if n <= order + 1:
        raise ValueError(f"Sequence too short for order {order}")

    if order < 1:
        raise ValueError("Order must be >= 1")

    # For entropy rate, we estimate H(X_{n+1} | X_1^n)
    # Using the chain rule: H(X_1^{n+1}) - H(X_1^n)

    # Build n+1 blocks
    blocks_n1 = [tuple(sequence[i : i + order + 1]) for i in range(n - order)]
    # Build n blocks
    blocks_n = [tuple(sequence[i : i + order]) for i in range(n - order + 1)]

    # Estimate entropies
    h_n1 = entropy_estimator(Counter(blocks_n1), method=method)
    h_n = entropy_estimator(Counter(blocks_n), method=method)

    # Entropy rate: H(X_{n+1}|X_1^n) = H(X_1^{n+1}) - H(X_1^n)
    entropy_rate = h_n1 - h_n

    return max(0.0, entropy_rate)
