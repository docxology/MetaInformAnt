"""Syntactic information theory methods.

This module implements fundamental information-theoretic measures including
Shannon entropy, mutual information, Kullback-Leibler divergence, and related
quantities for analyzing biological data.
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from typing import Any, Sequence

import numpy as np


def shannon_entropy(probs: Sequence[float], base: float = 2.0) -> float:
    """Calculate Shannon entropy of a probability distribution.
    
    Measures the information content or uncertainty in a probability
    distribution. Higher entropy indicates greater uncertainty/randomness.
    
    Args:
        probs: Probability distribution (should sum to 1, but function
            handles non-normalized inputs by normalizing)
        base: Logarithm base (2 for bits, e for nats, 10 for dits)
        
    Returns:
        Shannon entropy. Formula: H = -Σ p_i × log_base(p_i) for p_i > 0
        
    Examples:
        >>> shannon_entropy([0.5, 0.5])  # Maximum entropy for 2 outcomes
        1.0
        >>> shannon_entropy([1.0, 0.0])  # Certainty (zero entropy)
        0.0
        >>> shannon_entropy([0.25, 0.25, 0.25, 0.25])  # Maximum for 4 outcomes
        2.0
        
    References:
        Shannon, C. E. (1948). A mathematical theory of communication.
        Bell System Technical Journal, 27(3), 379-423.
    """
    if not probs:
        return 0.0
    
    # Normalize probabilities
    total = sum(p for p in probs if p > 0)
    if total == 0:
        return 0.0
    
    entropy = 0.0
    for p in probs:
        if p > 0:
            p_norm = p / total
            entropy -= p_norm * math.log(p_norm, base)
    
    return entropy


def shannon_entropy_from_counts(counts: Sequence[int] | dict[Any, int]) -> float:
    """Calculate Shannon entropy from counts or frequency dictionary.
    
    Args:
        counts: Counts or frequency dictionary
        
    Returns:
        Shannon entropy in bits
    """
    if isinstance(counts, dict):
        counts_list = list(counts.values())
    else:
        counts_list = list(counts)
    
    total = sum(counts_list)
    if total == 0:
        return 0.0
    
    probs = [c / total for c in counts_list if c > 0]
    return shannon_entropy(probs, base=2.0)


def joint_entropy(
    joint_probs: dict[tuple[Any, Any], float] | np.ndarray,
    base: float = 2.0
) -> float:
    """Calculate joint entropy of two random variables.
    
    Measures the total uncertainty in the joint distribution of X and Y.
    
    Args:
        joint_probs: Joint probability distribution as dictionary of
            (x, y) -> probability pairs, or 2D array
        base: Logarithm base
        
    Returns:
        Joint entropy H(X, Y)
        
    Examples:
        >>> joint_entropy({(0, 0): 0.25, (0, 1): 0.25, (1, 0): 0.25, (1, 1): 0.25})
        2.0  # Independent variables, maximum entropy
    """
    if isinstance(joint_probs, np.ndarray):
        entropy = 0.0
        total = np.sum(joint_probs)
        if total == 0:
            return 0.0
        for p in joint_probs.flatten():
            if p > 0:
                p_norm = p / total
                entropy -= p_norm * math.log(p_norm, base)
        return entropy
    
    total = sum(joint_probs.values())
    if total == 0:
        return 0.0
    
    entropy = 0.0
    for p in joint_probs.values():
        if p > 0:
            p_norm = p / total
            entropy -= p_norm * math.log(p_norm, base)
    
    return entropy


def conditional_entropy(
    x: Sequence[Any],
    y: Sequence[Any],
    base: float = 2.0
) -> float:
    """Calculate conditional entropy H(X|Y).
    
    Measures the uncertainty in X given knowledge of Y.
    
    Args:
        x: Sequence of X values
        y: Sequence of Y values (must match length of x)
        base: Logarithm base
        
    Returns:
        Conditional entropy H(X|Y)
        
    Examples:
        >>> x = [0, 1, 0, 1]
        >>> y = [0, 0, 1, 1]  # X = Y
        >>> conditional_entropy(x, y)  # Should be 0 (no uncertainty given Y)
        0.0
    """
    if len(x) != len(y):
        raise ValueError("X and Y must have the same length")
    
    if not x:
        return 0.0
    
    # Count joint and marginal distributions
    joint_counts: dict[tuple[Any, Any], int] = defaultdict(int)
    y_counts: dict[Any, int] = defaultdict(int)
    
    for xi, yi in zip(x, y):
        joint_counts[(xi, yi)] += 1
        y_counts[yi] += 1
    
    n = len(x)
    entropy = 0.0
    
    for y_val, y_count in y_counts.items():
        p_y = y_count / n
        if p_y == 0:
            continue
        
        # Calculate H(X|Y=y)
        conditional_entropy_y = 0.0
        for (xi, yi), joint_count in joint_counts.items():
            if yi == y_val:
                p_x_given_y = joint_count / y_count
                if p_x_given_y > 0:
                    conditional_entropy_y -= p_x_given_y * math.log(p_x_given_y, base)
        
        entropy += p_y * conditional_entropy_y
    
    return entropy


def mutual_information(
    x: Sequence[Any],
    y: Sequence[Any],
    base: float = 2.0
) -> float:
    """Calculate mutual information I(X; Y).
    
    Measures the amount of information shared between X and Y.
    I(X; Y) = H(X) + H(Y) - H(X, Y) = H(X) - H(X|Y)
    
    Args:
        x: Sequence of X values
        y: Sequence of Y values (must match length of x)
        base: Logarithm base
        
    Returns:
        Mutual information in bits (or nats/dits depending on base)
        
    Examples:
        >>> x = [0, 1, 0, 1]
        >>> y = [0, 1, 0, 1]  # Perfect correlation
        >>> mutual_information(x, y)  # Should equal H(X) = 1.0
        1.0
        
    References:
        Cover, T. M., & Thomas, J. A. (2006). Elements of Information Theory.
        John Wiley & Sons.
    """
    if len(x) != len(y):
        raise ValueError("X and Y must have the same length")
    
    if not x:
        return 0.0
    
    # Count distributions
    x_counts = Counter(x)
    y_counts = Counter(y)
    joint_counts: dict[tuple[Any, Any], int] = defaultdict(int)
    
    for xi, yi in zip(x, y):
        joint_counts[(xi, yi)] += 1
    
    n = len(x)
    
    # Calculate H(X)
    h_x = shannon_entropy_from_counts(x_counts)
    
    # Calculate H(Y)
    h_y = shannon_entropy_from_counts(y_counts)
    
    # Calculate H(X, Y)
    joint_probs = {k: v / n for k, v in joint_counts.items()}
    h_xy = joint_entropy(joint_probs, base=base)
    
    # I(X; Y) = H(X) + H(Y) - H(X, Y)
    mi = h_x + h_y - h_xy
    
    return max(0.0, mi)  # Mutual information is non-negative


def conditional_mutual_information(
    x: Sequence[Any],
    y: Sequence[Any],
    z: Sequence[Any],
    base: float = 2.0
) -> float:
    """Calculate conditional mutual information I(X; Y|Z).
    
    Measures the information shared between X and Y given Z.
    I(X; Y|Z) = H(X|Z) + H(Y|Z) - H(X, Y|Z)
    
    Args:
        x: Sequence of X values
        y: Sequence of Y values
        z: Sequence of Z values (all must have same length)
        base: Logarithm base
        
    Returns:
        Conditional mutual information
    """
    if not (len(x) == len(y) == len(z)):
        raise ValueError("X, Y, and Z must have the same length")
    
    if not x:
        return 0.0
    
    # Calculate H(X|Z)
    h_x_given_z = conditional_entropy(x, z, base=base)
    
    # Calculate H(Y|Z)
    h_y_given_z = conditional_entropy(y, z, base=base)
    
    # Calculate H(X, Y|Z) by creating combined XY sequence
    xy = [(xi, yi) for xi, yi in zip(x, y)]
    h_xy_given_z = conditional_entropy(xy, z, base=base)
    
    # I(X; Y|Z) = H(X|Z) + H(Y|Z) - H(X, Y|Z)
    cmi = h_x_given_z + h_y_given_z - h_xy_given_z
    
    return max(0.0, cmi)


def kl_divergence(
    p: Sequence[float],
    q: Sequence[float],
    base: float = 2.0
) -> float:
    """Calculate Kullback-Leibler divergence D_KL(P||Q).
    
    Measures how different distribution Q is from distribution P.
    Also known as relative entropy.
    
    Args:
        p: Probability distribution P
        q: Probability distribution Q (must have same length as p)
        base: Logarithm base
        
    Returns:
        KL divergence. Returns inf if Q has zeros where P has positives.
        
    Examples:
        >>> p = [0.5, 0.5]
        >>> q = [0.5, 0.5]
        >>> kl_divergence(p, q)  # Identical distributions
        0.0
        
    References:
        Kullback, S., & Leibler, R. A. (1951). On information and sufficiency.
        The Annals of Mathematical Statistics, 22(1), 79-86.
    """
    if len(p) != len(q):
        raise ValueError("P and Q must have the same length")
    
    # Normalize distributions
    p_sum = sum(p)
    q_sum = sum(q)
    
    if p_sum == 0:
        return 0.0
    if q_sum == 0:
        return float('inf')
    
    p_norm = [pi / p_sum for pi in p]
    q_norm = [qi / q_sum for qi in q]
    
    kl = 0.0
    for pi, qi in zip(p_norm, q_norm):
        if pi > 0:
            if qi == 0:
                return float('inf')
            kl += pi * math.log(pi / qi, base)
    
    return kl


def cross_entropy(
    p: Sequence[float],
    q: Sequence[float],
    base: float = 2.0
) -> float:
    """Calculate cross-entropy H(P, Q).
    
    Measures the average number of bits needed to encode events from
    distribution P using a code optimized for distribution Q.
    
    Args:
        p: True probability distribution P
        q: Approximate probability distribution Q
        base: Logarithm base
        
    Returns:
        Cross-entropy H(P, Q) = -Σ p_i × log(q_i)
    """
    if len(p) != len(q):
        raise ValueError("P and Q must have the same length")
    
    # Normalize distributions
    p_sum = sum(p)
    q_sum = sum(q)
    
    if p_sum == 0:
        return 0.0
    if q_sum == 0:
        return float('inf')
    
    p_norm = [pi / p_sum for pi in p]
    q_norm = [qi / q_sum for qi in q]
    
    entropy = 0.0
    for pi, qi in zip(p_norm, q_norm):
        if pi > 0:
            if qi == 0:
                return float('inf')
            entropy -= pi * math.log(qi, base)
    
    return entropy


def total_correlation(
    variables: list[Sequence[Any]],
    base: float = 2.0
) -> float:
    """Calculate total correlation (multivariate mutual information).
    
    Measures the total amount of dependence among multiple variables.
    Also known as multi-information.
    
    Args:
        variables: List of sequences, each representing a variable
        base: Logarithm base
        
    Returns:
        Total correlation. For n variables: TC = Σ H(X_i) - H(X_1, ..., X_n)
    """
    if not variables:
        return 0.0
    
    n_vars = len(variables)
    if n_vars < 2:
        return 0.0
    
    # Check all have same length
    lengths = [len(v) for v in variables]
    if len(set(lengths)) != 1:
        raise ValueError("All variables must have the same length")
    
    n = lengths[0]
    if n == 0:
        return 0.0
    
    # Calculate sum of individual entropies
    sum_individual_entropies = 0.0
    for var in variables:
        counts = Counter(var)
        sum_individual_entropies += shannon_entropy_from_counts(counts)
    
    # Calculate joint entropy
    # Create tuples of all variables
    joint_tuples = [tuple(v[i] for v in variables) for i in range(n)]
    joint_counts = Counter(joint_tuples)
    joint_probs = {k: v / n for k, v in joint_counts.items()}
    h_joint = shannon_entropy(list(joint_probs.values()), base=base)
    
    # Total correlation
    tc = sum_individual_entropies - h_joint
    
    return max(0.0, tc)


def transfer_entropy(
    x: Sequence[Any],
    y: Sequence[Any],
    lag: int = 1,
    base: float = 2.0
) -> float:
    """Calculate transfer entropy T(Y -> X).
    
    Measures the information transferred from Y to X, accounting for
    X's own past. T(Y -> X) = H(X_t | X_{t-1}) - H(X_t | X_{t-1}, Y_{t-1})
    
    Args:
        x: Time series X
        y: Time series Y (must have same length as x)
        lag: Time lag for past values
        base: Logarithm base
        
    Returns:
        Transfer entropy in bits
        
    Examples:
        >>> x = [0, 1, 0, 1, 0, 1]
        >>> y = [0, 0, 1, 1, 0, 0]  # Y doesn't predict X
        >>> te = transfer_entropy(x, y, lag=1)
        >>> te < 0.5  # Should be low
        True
        
    References:
        Schreiber, T. (2000). Measuring information transfer.
        Physical Review Letters, 85(2), 461.
    """
    if len(x) != len(y):
        raise ValueError("X and Y must have the same length")
    
    if len(x) <= lag:
        return 0.0
    
    # Create sequences for conditional entropy calculation
    x_present = x[lag:]
    x_past = x[:-lag]
    y_past = y[:-lag]
    
    # H(X_t | X_{t-1})
    h_x_given_x_past = conditional_entropy(x_present, x_past, base=base)
    
    # H(X_t | X_{t-1}, Y_{t-1})
    # Combine X_past and Y_past into tuples
    xy_past = [(xp, yp) for xp, yp in zip(x_past, y_past)]
    h_x_given_xy_past = conditional_entropy(x_present, xy_past, base=base)
    
    # Transfer entropy
    te = h_x_given_x_past - h_x_given_xy_past
    
    return max(0.0, te)


def jensen_shannon_divergence(
    p: Sequence[float],
    q: Sequence[float],
    base: float = 2.0
) -> float:
    """Calculate Jensen-Shannon divergence between two probability distributions.
    
    JS divergence is a symmetrized version of Kullback-Leibler divergence,
    measuring the distance between two probability distributions. Always
    bounded in [0, 1] when using log base 2.
    
    Args:
        p: First probability distribution (list of probabilities)
        q: Second probability distribution (must have same length as p)
        base: Logarithm base
        
    Returns:
        Jensen-Shannon divergence. Returns 0.0 if distributions are identical
        or if either distribution sums to zero.
        Formula: JS(P||Q) = H(M) - [H(P) + H(Q)]/2 where M = (P+Q)/2
        
    Examples:
        >>> p = [0.5, 0.5]
        >>> q = [0.5, 0.5]
        >>> jensen_shannon_divergence(p, q)
        0.0  # Identical distributions
        
        >>> p = [1.0, 0.0]
        >>> q = [0.0, 1.0]
        >>> js = jensen_shannon_divergence(p, q)
        >>> js > 0.0  # Divergent distributions
        True
        
    Raises:
        ValueError: If distributions have different lengths
        
    References:
        Lin, J. (1991). Divergence measures based on the Shannon entropy.
        IEEE Transactions on Information Theory, 37(1), 145-151.
    """
    if len(p) != len(q):
        raise ValueError("Distributions must have same length")
    
    # Normalize distributions
    p_sum = sum(p)
    q_sum = sum(q)
    
    if p_sum == 0 or q_sum == 0:
        return 0.0
    
    p_norm = [pi / p_sum for pi in p]
    q_norm = [qi / q_sum for qi in q]
    
    # Calculate midpoint distribution
    m = [(pi + qi) / 2 for pi, qi in zip(p_norm, q_norm)]
    
    # Calculate entropies
    h_p = shannon_entropy(p_norm, base=base)
    h_q = shannon_entropy(q_norm, base=base)
    h_m = shannon_entropy(m, base=base)
    
    return h_m - (h_p + h_q) / 2


def renyi_entropy(
    probs: Sequence[float],
    alpha: float = 2.0,
    base: float = 2.0
) -> float:
    """Calculate Rényi entropy of order α.
    
    Rényi entropy is a generalization of Shannon entropy. When α=1, it equals
    Shannon entropy. When α=2, it's called collision entropy.
    
    Args:
        probs: Probability distribution
        alpha: Order parameter (α > 0, α ≠ 1). α=1 uses Shannon entropy.
        base: Logarithm base
        
    Returns:
        Rényi entropy. Formula: H_α(X) = (1/(1-α)) × log(Σ p_i^α)
        
    Examples:
        >>> probs = [0.5, 0.5]
        >>> renyi_entropy(probs, alpha=2.0)  # Collision entropy
        1.0
        >>> abs(renyi_entropy(probs, alpha=1.0) - shannon_entropy(probs)) < 1e-10
        True  # α=1 equals Shannon entropy
        
    References:
        Rényi, A. (1961). On measures of entropy and information.
        Proceedings of the Fourth Berkeley Symposium on Mathematical Statistics
        and Probability, 1, 547-561.
    """
    if not probs:
        return 0.0
    
    # Normalize probabilities
    total = sum(p for p in probs if p > 0)
    if total == 0:
        return 0.0
    
    p_norm = [p / total for p in probs if p > 0]
    
    # Special case: α = 1 (Shannon entropy)
    if abs(alpha - 1.0) < 1e-10:
        return shannon_entropy(p_norm, base=base)
    
    # General case: H_α = (1/(1-α)) × log(Σ p_i^α)
    if alpha < 0:
        raise ValueError("Alpha must be non-negative")
    
    sum_power = sum(p ** alpha for p in p_norm)
    if sum_power == 0:
        return 0.0
    
    return (1.0 / (1.0 - alpha)) * math.log(sum_power, base)


def tsallis_entropy(
    probs: Sequence[float],
    q: float = 2.0,
    base: float = 2.0
) -> float:
    """Calculate Tsallis entropy (non-extensive entropy).
    
    Tsallis entropy is a generalization of Shannon entropy used in
    non-extensive statistical mechanics. When q=1, it equals Shannon entropy.
    
    Args:
        probs: Probability distribution
        q: Entropic index (q > 0, q ≠ 1). q=1 uses Shannon entropy.
        base: Logarithm base (for normalization)
        
    Returns:
        Tsallis entropy. Formula: S_q = (1/(q-1)) × (1 - Σ p_i^q)
        
    Examples:
        >>> probs = [0.5, 0.5]
        >>> tsallis_entropy(probs, q=2.0)
        0.5  # For q=2, S_2 = 1 - Σ p_i^2
        >>> abs(tsallis_entropy(probs, q=1.0) - shannon_entropy(probs)) < 1e-10
        True  # q=1 equals Shannon entropy
        
    References:
        Tsallis, C. (1988). Possible generalization of Boltzmann-Gibbs statistics.
        Journal of Statistical Physics, 52(1-2), 479-487.
    """
    if not probs:
        return 0.0
    
    # Normalize probabilities
    total = sum(p for p in probs if p > 0)
    if total == 0:
        return 0.0
    
    p_norm = [p / total for p in probs if p > 0]
    
    # Special case: q = 1 (Shannon entropy)
    if abs(q - 1.0) < 1e-10:
        return shannon_entropy(p_norm, base=base)
    
    if q < 0:
        raise ValueError("q must be non-negative")
    
    # Tsallis entropy: S_q = (1/(q-1)) × (1 - Σ p_i^q)
    sum_power = sum(p ** q for p in p_norm)
    return (1.0 / (q - 1.0)) * (1.0 - sum_power)


def normalized_mutual_information(
    x: Sequence[Any],
    y: Sequence[Any],
    method: str = "arithmetic",
    base: float = 2.0
) -> float:
    """Calculate normalized mutual information (NMI).
    
    NMI normalizes mutual information to [0, 1] range for easier comparison.
    Common normalization methods: arithmetic mean, geometric mean, min, max.
    
    Args:
        x: Sequence of X values
        y: Sequence of Y values (must match length of x)
        method: Normalization method ("arithmetic", "geometric", "min", "max")
        base: Logarithm base
        
    Returns:
        Normalized mutual information in [0, 1]
        
    Examples:
        >>> x = [0, 1, 0, 1]
        >>> y = [0, 1, 0, 1]  # Perfect correlation
        >>> nmi = normalized_mutual_information(x, y)
        >>> abs(nmi - 1.0) < 1e-10  # Should be 1.0
        True
    """
    if len(x) != len(y):
        raise ValueError("X and Y must have the same length")
    
    if not x:
        return 0.0
    
    # Calculate MI and individual entropies
    mi = mutual_information(x, y, base=base)
    
    x_counts = Counter(x)
    y_counts = Counter(y)
    h_x = shannon_entropy_from_counts(x_counts)
    h_y = shannon_entropy_from_counts(y_counts)
    
    if h_x == 0.0 and h_y == 0.0:
        return 1.0 if mi > 0 else 0.0
    
    # Normalize by chosen method
    if method == "arithmetic":
        denominator = (h_x + h_y) / 2.0
    elif method == "geometric":
        denominator = math.sqrt(h_x * h_y) if h_x > 0 and h_y > 0 else 0.0
    elif method == "min":
        denominator = min(h_x, h_y)
    elif method == "max":
        denominator = max(h_x, h_y)
    else:
        raise ValueError(f"Unknown normalization method: {method}")
    
    if denominator == 0:
        return 0.0
    
    return mi / denominator


def information_coefficient(
    x: Sequence[float],
    y: Sequence[float],
    grid_resolution: int = 10
) -> float:
    """Calculate information coefficient (MIC-like measure).
    
    Estimates mutual information for continuous variables by discretization.
    Similar to Maximal Information Coefficient (MIC) but simpler.
    
    Args:
        x: Sequence of continuous X values
        y: Sequence of continuous Y values (must match length of x)
        grid_resolution: Number of bins for discretization
        
    Returns:
        Information coefficient in [0, 1] (normalized MI)
        
    Examples:
        >>> import numpy as np
        >>> x = np.random.randn(100)
        >>> y = x + np.random.randn(100) * 0.1  # Strong correlation
        >>> ic = information_coefficient(x, y)
        >>> ic > 0.5  # Should have high IC
        True
        
    References:
        Reshef, D. N., et al. (2011). Detecting novel associations in large datasets.
        Science, 334(6062), 1518-1524.
    """
    if len(x) != len(y):
        raise ValueError("X and Y must have the same length")
    
    if len(x) < 2:
        return 0.0
    
    # Discretize continuous variables
    x_min, x_max = min(x), max(x)
    y_min, y_max = min(y), max(y)
    
    if x_max == x_min or y_max == y_min:
        return 0.0
    
    x_bins = [int((xi - x_min) / (x_max - x_min) * grid_resolution) for xi in x]
    y_bins = [int((yi - y_min) / (y_max - y_min) * grid_resolution) for yi in y]
    
    # Calculate normalized MI
    return normalized_mutual_information(x_bins, y_bins, method="geometric")

