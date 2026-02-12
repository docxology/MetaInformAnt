"""Advanced information-theoretic measures for biological data.

This module implements advanced measures including Fisher information,
variation of information, interaction information, and other quantities
useful for biological sequence and systems analysis.
"""

from __future__ import annotations

import math
from collections import Counter
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from metainformant.core.data import validation
from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def fisher_information(
    samples: np.ndarray,
    log_likelihood_grad: Optional[Callable] = None,
    method: str = "empirical",
) -> float:
    """Calculate scalar Fisher information for a 1D parameter.

    Fisher information quantifies the amount of information that an
    observable random variable carries about an unknown parameter.

    Args:
        samples: 1D array of observed samples
        log_likelihood_grad: Optional gradient function of log-likelihood.
            If None, uses empirical estimation from score function.
        method: Estimation method ('empirical', 'parametric_normal')

    Returns:
        Fisher information value

    Raises:
        ValueError: If samples array is too small
    """
    samples = np.asarray(samples).flatten()
    if len(samples) < 5:
        raise ValueError("Need at least 5 samples for Fisher information estimation")

    if method == "parametric_normal":
        # For normal distribution: I(mu) = n/sigma^2
        variance = np.var(samples, ddof=1)
        if variance <= 0:
            return float("inf")
        return 1.0 / variance

    elif method == "empirical":
        if log_likelihood_grad is not None:
            # Use provided gradient function
            grads = np.array([log_likelihood_grad(s) for s in samples])
            return float(np.mean(grads**2))
        else:
            # Empirical estimation using KDE score function
            from scipy.stats import gaussian_kde

            kde = gaussian_kde(samples)
            # Score function: d/dx log f(x) = f'(x)/f(x)
            h = 1e-5 * np.std(samples)
            if h == 0:
                return 0.0
            scores = []
            for s in samples:
                f_plus = kde.evaluate(s + h)[0]
                f_minus = kde.evaluate(s - h)[0]
                f_val = kde.evaluate(s)[0]
                if f_val > 0:
                    score = (f_plus - f_minus) / (2 * h * f_val) * f_val
                    scores.append(score**2 / f_val)
            if scores:
                return float(np.mean(scores))
            return 0.0
    else:
        raise ValueError(f"Unknown method: {method}")


def fisher_information_matrix(
    samples: np.ndarray,
    method: str = "empirical",
) -> np.ndarray:
    """Calculate the Fisher information matrix for multivariate data.

    For d-dimensional data, returns a d x d matrix where each entry
    I_{ij} measures the joint information about parameters i and j.

    Args:
        samples: 2D array (n_samples, n_dimensions)
        method: Estimation method ('empirical', 'parametric_normal')

    Returns:
        Fisher information matrix (d x d)

    Raises:
        ValueError: If samples has wrong shape
    """
    samples = np.asarray(samples)
    if samples.ndim == 1:
        samples = samples.reshape(-1, 1)
    if samples.ndim != 2:
        raise ValueError("Samples must be 1D or 2D array")

    n_samples, n_dims = samples.shape
    if n_samples < n_dims + 2:
        raise ValueError(f"Need at least {n_dims + 2} samples for {n_dims}-dimensional FIM")

    if method == "parametric_normal":
        # For multivariate normal: FIM = inverse of covariance matrix
        cov = np.cov(samples, rowvar=False)
        try:
            fim = np.linalg.inv(cov)
        except np.linalg.LinAlgError:
            fim = np.linalg.pinv(cov)
        return fim

    elif method == "empirical":
        # Empirical FIM via numerical score function estimation
        from scipy.stats import gaussian_kde

        try:
            kde = gaussian_kde(samples.T)
        except np.linalg.LinAlgError:
            # Singular covariance — fall back to parametric
            return fisher_information_matrix(samples, method="parametric_normal")

        fim = np.zeros((n_dims, n_dims))
        h = 1e-5 * np.std(samples, axis=0)
        h = np.maximum(h, 1e-10)

        for idx in range(min(n_samples, 500)):  # Subsample for performance
            x = samples[idx]
            f_val = kde.evaluate(x.reshape(-1, 1))[0]
            if f_val <= 0:
                continue

            score = np.zeros(n_dims)
            for d in range(n_dims):
                x_plus = x.copy()
                x_minus = x.copy()
                x_plus[d] += h[d]
                x_minus[d] -= h[d]
                f_plus = kde.evaluate(x_plus.reshape(-1, 1))[0]
                f_minus = kde.evaluate(x_minus.reshape(-1, 1))[0]
                score[d] = (f_plus - f_minus) / (2 * h[d] * f_val)

            fim += np.outer(score, score)

        n_used = min(n_samples, 500)
        fim /= n_used
        return fim

    else:
        raise ValueError(f"Unknown method: {method}")


def relative_information_gain(
    prior_probs: Sequence[float],
    posterior_probs: Sequence[float],
    base: float = 2.0,
) -> float:
    """Calculate relative information gain (KL divergence from prior to posterior).

    Measures how much information was gained by updating from prior
    to posterior distribution. Used in Bayesian analysis of biological data.

    Args:
        prior_probs: Prior probability distribution
        posterior_probs: Posterior probability distribution
        base: Logarithm base

    Returns:
        Information gain in bits (or nats if base=e)

    Raises:
        ValueError: If distributions have different lengths or don't sum to 1
    """
    from metainformant.information.metrics.core.syntactic import kl_divergence

    return kl_divergence(posterior_probs, prior_probs, base=base)


def variation_of_information(
    x: Sequence[Any],
    y: Sequence[Any],
    base: float = 2.0,
) -> float:
    """Calculate variation of information VI(X, Y).

    VI is a true metric on the space of clusterings/partitions.
    VI(X,Y) = H(X|Y) + H(Y|X) = H(X,Y) - I(X;Y)

    Useful for comparing clustering results in single-cell analysis
    or community detection in biological networks.

    Args:
        x: First sequence (e.g., cluster labels)
        y: Second sequence (e.g., cluster labels)
        base: Logarithm base

    Returns:
        Variation of information (non-negative, 0 = identical)

    Raises:
        ValueError: If sequences have different lengths
    """
    from metainformant.information.metrics.core.syntactic import joint_entropy, mutual_information

    h_xy = joint_entropy(x, y, base=base)
    mi = mutual_information(x, y, base=base)

    # VI = H(X,Y) - I(X;Y) = H(X|Y) + H(Y|X)
    # Equivalently: 2*H(X,Y) - H(X) - H(Y)
    vi = h_xy - mi
    return max(0.0, vi)


def interaction_information(
    x: Sequence[Any],
    y: Sequence[Any],
    z: Sequence[Any],
    base: float = 2.0,
) -> float:
    """Calculate interaction information II(X;Y;Z).

    Interaction information generalizes mutual information to three variables.
    It can be positive (synergy) or negative (redundancy).

    II(X;Y;Z) = I(X;Y|Z) - I(X;Y)
              = I(X;Y) + I(X;Z) + I(Y;Z) - I(X;Y;Z_joint) ... (inclusion-exclusion)

    In biology, this reveals whether the combined effect of two variables
    on a third is synergistic or redundant.

    Args:
        x: First variable sequence
        y: Second variable sequence
        z: Third variable sequence
        base: Logarithm base

    Returns:
        Interaction information (positive = synergy, negative = redundancy)

    Raises:
        ValueError: If sequences have different lengths
    """
    from metainformant.information.metrics.core.syntactic import conditional_mutual_information, mutual_information

    # II(X;Y;Z) = I(X;Y|Z) - I(X;Y)
    cmi = conditional_mutual_information(x, y, z, base=base)
    mi = mutual_information(x, y, base=base)

    return cmi - mi


def binding_information(
    x: Sequence[Any],
    y: Sequence[Any],
    base: float = 2.0,
) -> float:
    """Calculate binding information B(X;Y).

    Binding information measures the minimum information that must be
    communicated to reconstruct the joint distribution from the marginals.
    B(X;Y) = I(X;Y) for pairs but generalizes differently for >2 variables.

    Relevant to modeling molecular binding where the joint state of
    two molecules determines function.

    Args:
        x: First variable sequence
        y: Second variable sequence
        base: Logarithm base

    Returns:
        Binding information value

    Raises:
        ValueError: If sequences have different lengths
    """
    from metainformant.information.metrics.core.syntactic import mutual_information

    return mutual_information(x, y, base=base)


def lautum_information(
    x: Sequence[Any],
    y: Sequence[Any],
    base: float = 2.0,
) -> float:
    """Calculate lautum information L(X;Y).

    Lautum information is the 'reverse' of mutual information:
    L(X;Y) = D_KL(P(X)P(Y) || P(X,Y))

    While MI = D_KL(P(X,Y) || P(X)P(Y)), lautum reverses the arguments.
    The name 'lautum' is 'mutual' backwards.

    Args:
        x: First variable sequence
        y: Second variable sequence
        base: Logarithm base

    Returns:
        Lautum information value (non-negative)

    Raises:
        ValueError: If sequences have different lengths
    """
    validation.validate_type(x, (list, tuple), "x")
    validation.validate_type(y, (list, tuple), "y")

    if len(x) != len(y):
        raise ValueError("Sequences must have the same length")

    n = len(x)
    # Compute joint and marginal distributions
    joint_counts = Counter(zip(x, y))
    x_counts = Counter(x)
    y_counts = Counter(y)

    # L(X;Y) = D_KL(P(X)P(Y) || P(X,Y))
    # = sum_{x,y} P(x)P(y) log(P(x)P(y) / P(x,y))
    lautum = 0.0
    for xi in x_counts:
        for yi in y_counts:
            p_x = x_counts[xi] / n
            p_y = y_counts[yi] / n
            p_marginal_product = p_x * p_y
            p_joint = joint_counts.get((xi, yi), 0) / n

            if p_marginal_product > 0 and p_joint > 0:
                lautum += p_marginal_product * math.log(p_marginal_product / p_joint) / math.log(base)
            elif p_marginal_product > 0 and p_joint == 0:
                return float("inf")

    return max(0.0, lautum)
