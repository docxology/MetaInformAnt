"""Syntactic information theory measures for biological sequences.

This module implements core information-theoretic measures including Shannon entropy,
mutual information, KL divergence, and higher-order measures for analyzing
biological sequences and data patterns.
"""

from __future__ import annotations

import math
from collections import Counter
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from metainformant.core import logging, validation

logger = logging.get_logger(__name__)


def shannon_entropy(probs: Sequence[float], base: float = 2.0) -> float:
    """Calculate Shannon entropy from probability distribution.

    Args:
        probs: Sequence of probabilities (must sum to 1.0)
        base: Logarithm base (2.0 for bits, math.e for nats)

    Returns:
        Shannon entropy value

    Raises:
        ValueError: If probabilities don't sum to 1.0 or contain negative values
    """
    validation.validate_type(probs, (list, tuple, np.ndarray), "probs")

    # Convert to numpy array for easier handling
    probs_array = np.array(probs, dtype=float)

    # Validate probabilities
    if np.any(probs_array < 0):
        raise ValueError("Probabilities cannot be negative")
    if not np.isclose(np.sum(probs_array), 1.0, atol=1e-6):
        raise ValueError("Probabilities must sum to 1.0")

    # Remove zero probabilities to avoid log(0)
    probs_array = probs_array[probs_array > 0]

    # Calculate entropy
    entropy = -np.sum(probs_array * np.log(probs_array) / np.log(base))

    return float(entropy)


def shannon_entropy_from_counts(
    counts: Union[Sequence[int], Dict[Any, int]]
) -> float:
    """Calculate Shannon entropy from frequency counts.

    Args:
        counts: Either a sequence of counts or a dictionary mapping items to counts

    Returns:
        Shannon entropy in bits

    Raises:
        ValueError: If counts contain negative values or all counts are zero
    """
    if isinstance(counts, dict):
        count_values = list(counts.values())
    else:
        count_values = list(counts)

    validation.validate_type(count_values, list, "counts")

    # Convert to numpy array
    count_array = np.array(count_values, dtype=int)

    if np.any(count_array < 0):
        raise ValueError("Counts cannot be negative")

    total = np.sum(count_array)
    if total == 0:
        return 0.0

    # Convert to probabilities
    probs = count_array / total

    return shannon_entropy(probs, base=2.0)


def joint_entropy(x: Sequence[Any], y: Sequence[Any], base: float = 2.0) -> float:
    """Calculate joint entropy H(X,Y).

    Args:
        x: First sequence
        y: Second sequence (same length as x)
        base: Logarithm base

    Returns:
        Joint entropy value

    Raises:
        ValueError: If sequences have different lengths
    """
    validation.validate_type(x, (list, tuple), "x")
    validation.validate_type(y, (list, tuple), "y")

    if len(x) != len(y):
        raise ValueError("Sequences must have the same length")

    # Count joint occurrences
    joint_counts = Counter(zip(x, y))

    # Calculate total
    total = len(x)

    # Convert to probabilities
    joint_probs = [count / total for count in joint_counts.values()]

    return shannon_entropy(joint_probs, base=base)


def conditional_entropy(
    x: Sequence[Any],
    y: Sequence[Any],
    base: float = 2.0
) -> float:
    """Calculate conditional entropy H(X|Y).

    Args:
        x: First sequence
        y: Second sequence (condition)
        base: Logarithm base

    Returns:
        Conditional entropy H(X|Y)

    Raises:
        ValueError: If sequences have different lengths
    """
    validation.validate_type(x, (list, tuple), "x")
    validation.validate_type(y, (list, tuple), "y")

    if len(x) != len(y):
        raise ValueError("Sequences must have the same length")

    # Calculate joint entropy H(X,Y)
    h_xy = joint_entropy(x, y, base=base)

    # Calculate entropy of Y H(Y)
    h_y = shannon_entropy_from_counts(Counter(y))

    # H(X|Y) = H(X,Y) - H(Y)
    return h_xy - h_y


def mutual_information(
    x: Sequence[Any],
    y: Sequence[Any],
    base: float = 2.0
) -> float:
    """Calculate mutual information I(X;Y).

    Args:
        x: First sequence
        y: Second sequence
        base: Logarithm base

    Returns:
        Mutual information I(X;Y)

    Raises:
        ValueError: If sequences have different lengths
    """
    validation.validate_type(x, (list, tuple), "x")
    validation.validate_type(y, (list, tuple), "y")

    if len(x) != len(y):
        raise ValueError("Sequences must have the same length")

    # I(X;Y) = H(X) + H(Y) - H(X,Y)
    h_x = shannon_entropy_from_counts(Counter(x))
    h_y = shannon_entropy_from_counts(Counter(y))
    h_xy = joint_entropy(x, y, base=base)

    return h_x + h_y - h_xy


def conditional_mutual_information(
    x: Sequence[Any],
    y: Sequence[Any],
    z: Sequence[Any],
    base: float = 2.0
) -> float:
    """Calculate conditional mutual information I(X;Y|Z).

    Args:
        x: First sequence
        y: Second sequence
        z: Conditioning sequence
        base: Logarithm base

    Returns:
        Conditional mutual information I(X;Y|Z)

    Raises:
        ValueError: If sequences have different lengths
    """
    validation.validate_type(x, (list, tuple), "x")
    validation.validate_type(y, (list, tuple), "y")
    validation.validate_type(z, (list, tuple), "z")

    if not (len(x) == len(y) == len(z)):
        raise ValueError("All sequences must have the same length")

    # I(X;Y|Z) = H(X,Z) + H(Y,Z) - H(Z) - H(X,Y,Z)
    h_xz = joint_entropy(x, z, base=base)
    h_yz = joint_entropy(y, z, base=base)
    h_z = shannon_entropy_from_counts(Counter(z))

    # Calculate joint entropy H(X,Y,Z)
    xyz_counts = Counter(zip(x, y, z))
    total = len(x)
    xyz_probs = [count / total for count in xyz_counts.values()]
    h_xyz = shannon_entropy(xyz_probs, base=base)

    return h_xz + h_yz - h_z - h_xyz


def kl_divergence(
    p: Sequence[float],
    q: Sequence[float],
    base: float = 2.0
) -> float:
    """Calculate Kullback-Leibler divergence D_KL(p||q).

    Args:
        p: First probability distribution
        q: Second probability distribution
        base: Logarithm base

    Returns:
        KL divergence value

    Raises:
        ValueError: If distributions have different lengths or don't sum to 1
    """
    validation.validate_type(p, (list, tuple, np.ndarray), "p")
    validation.validate_type(q, (list, tuple, np.ndarray), "q")

    p_array = np.array(p, dtype=float)
    q_array = np.array(q, dtype=float)

    if len(p_array) != len(q_array):
        raise ValueError("Probability distributions must have the same length")

    if not np.isclose(np.sum(p_array), 1.0, atol=1e-6):
        raise ValueError("First distribution must sum to 1.0")
    if not np.isclose(np.sum(q_array), 1.0, atol=1e-6):
        raise ValueError("Second distribution must sum to 1.0")

    # KL divergence: sum(p * log(p/q))
    # Handle zero probabilities carefully
    kl_sum = 0.0
    for pi, qi in zip(p_array, q_array):
        if pi > 0:
            if qi > 0:
                kl_sum += pi * (math.log(pi / qi) / math.log(base))
            else:
                return float('inf')  # Infinite divergence if q has zero where p doesn't

    return kl_sum


def cross_entropy(
    p: Sequence[float],
    q: Sequence[float],
    base: float = 2.0
) -> float:
    """Calculate cross-entropy H(p, q).

    Args:
        p: True probability distribution
        q: Predicted probability distribution
        base: Logarithm base

    Returns:
        Cross-entropy value

    Raises:
        ValueError: If distributions have different lengths or don't sum to 1
    """
    validation.validate_type(p, (list, tuple, np.ndarray), "p")
    validation.validate_type(q, (list, tuple, np.ndarray), "q")

    p_array = np.array(p, dtype=float)
    q_array = np.array(q, dtype=float)

    if len(p_array) != len(q_array):
        raise ValueError("Probability distributions must have the same length")

    if not np.isclose(np.sum(p_array), 1.0, atol=1e-6):
        raise ValueError("True distribution must sum to 1.0")

    # Cross-entropy: -sum(p * log(q))
    ce_sum = 0.0
    for pi, qi in zip(p_array, q_array):
        if pi > 0:
            if qi > 0:
                ce_sum -= pi * (math.log(qi) / math.log(base))
            else:
                return float('inf')  # Infinite if q has zero probability

    return ce_sum


def jensen_shannon_divergence(
    p: Sequence[float],
    q: Sequence[float],
    base: float = 2.0
) -> float:
    """Calculate Jensen-Shannon divergence JSD(p||q).

    Args:
        p: First probability distribution
        q: Second probability distribution
        base: Logarithm base

    Returns:
        Jensen-Shannon divergence value
    """
    p_array = np.array(p, dtype=float)
    q_array = np.array(q, dtype=float)

    # Jensen-Shannon divergence: (1/2) * [KL(p||m) + KL(q||m)]
    # where m = (p + q) / 2
    m = (p_array + q_array) / 2

    kl_pm = kl_divergence(p_array, m, base=base)
    kl_qm = kl_divergence(q_array, m, base=base)

    return 0.5 * (kl_pm + kl_qm)


def total_correlation(
    variables: List[Sequence[Any]],
    base: float = 2.0
) -> float:
    """Calculate total correlation (multi-information) for multiple variables.

    Args:
        variables: List of sequences (all same length)
        base: Logarithm base

    Returns:
        Total correlation value

    Raises:
        ValueError: If variables have different lengths
    """
    if not variables:
        return 0.0

    validation.validate_type(variables, list, "variables")

    # Check all variables have same length
    lengths = [len(var) for var in variables]
    if len(set(lengths)) != 1:
        raise ValueError("All variable sequences must have the same length")

    n_vars = len(variables)

    if n_vars == 1:
        return 0.0

    # Total correlation: sum(H(X_i)) - H(X_1, X_2, ..., X_n)
    individual_entropies = sum(
        shannon_entropy_from_counts(Counter(var))
        for var in variables
    )

    # Joint entropy of all variables
    joint_counts = Counter(zip(*variables))
    joint_entropy_val = shannon_entropy_from_counts(joint_counts)

    return individual_entropies - joint_entropy_val


def transfer_entropy(
    x: Sequence[Any],
    y: Sequence[Any],
    lag: int = 1,
    base: float = 2.0
) -> float:
    """Calculate transfer entropy from X to Y.

    Transfer entropy measures the amount of information that X provides
    about future states of Y that is not available from past states of Y.

    Args:
        x: Source sequence
        y: Target sequence
        lag: Time lag for transfer
        base: Logarithm base

    Returns:
        Transfer entropy T(X→Y)

    Raises:
        ValueError: If sequences have different lengths or lag is invalid
    """
    validation.validate_type(x, (list, tuple), "x")
    validation.validate_type(y, (list, tuple), "y")

    if len(x) != len(y):
        raise ValueError("Sequences must have the same length")

    if lag < 1:
        raise ValueError("Lag must be >= 1")

    n = len(x)
    if n <= lag + 1:
        raise ValueError("Sequences too short for given lag")

    # Transfer entropy: H(Y_{t+1} | Y_t) - H(Y_{t+1} | Y_t, X_t)
    # This is equivalent to: I(Y_{t+1}; X_t | Y_t)

    # Build sequences for conditional mutual information
    y_future = y[lag:]      # Y_{t+1}
    y_past = y[:-lag]       # Y_t
    x_past = x[:-lag]       # X_t

    # Ensure all sequences have same length
    min_len = min(len(y_future), len(y_past), len(x_past))
    y_future = y_future[:min_len]
    y_past = y_past[:min_len]
    x_past = x_past[:min_len]

    # Calculate I(Y_{t+1}; X_t | Y_t)
    return conditional_mutual_information(y_future, x_past, y_past, base=base)


def renyi_entropy(
    probs: Sequence[float],
    alpha: float = 2.0,
    base: float = 2.0
) -> float:
    """Calculate Rényi entropy of order alpha.

    Args:
        probs: Probability distribution
        alpha: Order parameter (alpha != 1)
        base: Logarithm base

    Returns:
        Rényi entropy value

    Raises:
        ValueError: If alpha = 1 or invalid probabilities
    """
    if alpha == 1:
        raise ValueError("Rényi entropy is undefined for alpha = 1")

    validation.validate_type(probs, (list, tuple, np.ndarray), "probs")

    probs_array = np.array(probs, dtype=float)

    if np.any(probs_array < 0):
        raise ValueError("Probabilities cannot be negative")
    if not np.isclose(np.sum(probs_array), 1.0, atol=1e-6):
        raise ValueError("Probabilities must sum to 1.0")

    # Remove zero probabilities
    probs_array = probs_array[probs_array > 0]

    if alpha == 0:
        # Hartley entropy (log of support size)
        return math.log(len(probs_array), base)
    elif alpha == float('inf'):
        # Min entropy
        return -math.log(np.max(probs_array), base)
    else:
        # General Rényi entropy: (1/(1-alpha)) * log(sum(p^alpha))
        sum_p_alpha = np.sum(probs_array ** alpha)
        if sum_p_alpha > 0:
            return (1.0 / (1.0 - alpha)) * math.log(sum_p_alpha, base)
        else:
            return 0.0


def tsallis_entropy(
    probs: Sequence[float],
    q: float = 2.0,
    base: float = 2.0
) -> float:
    """Calculate Tsallis entropy of order q.

    Args:
        probs: Probability distribution
        q: Entropic index (q != 1)
        base: Logarithm base

    Returns:
        Tsallis entropy value

    Raises:
        ValueError: If q = 1 or invalid probabilities
    """
    if q == 1:
        raise ValueError("Tsallis entropy is undefined for q = 1")

    validation.validate_type(probs, (list, tuple, np.ndarray), "probs")

    probs_array = np.array(probs, dtype=float)

    if np.any(probs_array < 0):
        raise ValueError("Probabilities cannot be negative")
    if not np.isclose(np.sum(probs_array), 1.0, atol=1e-6):
        raise ValueError("Probabilities must sum to 1.0")

    # Remove zero probabilities
    probs_array = probs_array[probs_array > 0]

    if q == 0:
        # Hartley entropy
        return math.log(len(probs_array), base)
    else:
        # Tsallis entropy: (1/(q-1)) * (sum(p^q) - 1)
        sum_p_q = np.sum(probs_array ** q)
        return (1.0 / (q - 1.0)) * (sum_p_q - 1.0)


def normalized_mutual_information(
    x: Sequence[Any],
    y: Sequence[Any],
    method: str = "arithmetic",
    base: float = 2.0
) -> float:
    """Calculate normalized mutual information.

    Args:
        x: First sequence
        y: Second sequence
        method: Normalization method ('arithmetic', 'geometric', 'max', 'min')
        base: Logarithm base

    Returns:
        Normalized mutual information (0-1 range)

    Raises:
        ValueError: If sequences have different lengths or invalid method
    """
    validation.validate_type(x, (list, tuple), "x")
    validation.validate_type(y, (list, tuple), "y")

    if len(x) != len(y):
        raise ValueError("Sequences must have the same length")

    # Calculate mutual information
    mi = mutual_information(x, y, base=base)

    # Calculate entropies for normalization
    h_x = shannon_entropy_from_counts(Counter(x))
    h_y = shannon_entropy_from_counts(Counter(y))

    # Choose normalization method
    if method == "arithmetic":
        # NMI = 2 * I(X;Y) / (H(X) + H(Y))
        denominator = h_x + h_y
    elif method == "geometric":
        # NMI = 2 * I(X;Y) / (H(X) + H(Y)) * 2 / (H(X) + H(Y)) wait, that's wrong
        # Actually geometric mean normalization
        denominator = math.sqrt(h_x * h_y) * 2
    elif method == "max":
        # NMI = I(X;Y) / max(H(X), H(Y))
        denominator = max(h_x, h_y)
    elif method == "min":
        # NMI = I(X;Y) / min(H(X), H(Y))
        denominator = min(h_x, h_y)
    else:
        raise ValueError(f"Unknown normalization method: {method}")

    if denominator == 0:
        return 0.0

    return 2.0 * mi / denominator if method in ["arithmetic", "geometric"] else mi / denominator


def information_coefficient(
    x: Sequence[Any],
    y: Sequence[Any],
    base: float = 2.0
) -> float:
    """Calculate information coefficient (IC).

    This is a measure of the total information shared between two variables,
    similar to mutual information but normalized.

    Args:
        x: First sequence
        y: Second sequence
        base: Logarithm base

    Returns:
        Information coefficient value

    Raises:
        ValueError: If sequences have different lengths
    """
    validation.validate_type(x, (list, tuple), "x")
    validation.validate_type(y, (list, tuple), "y")

    if len(x) != len(y):
        raise ValueError("Sequences must have the same length")

    # Information coefficient: I(X;Y) / [H(X) + H(Y) - I(X;Y)]
    # This gives a value between 0 and 1
    mi = mutual_information(x, y, base=base)
    h_x = shannon_entropy_from_counts(Counter(x))
    h_y = shannon_entropy_from_counts(Counter(y))

    denominator = h_x + h_y - mi

    if denominator == 0:
        return 0.0

    return mi / denominator






