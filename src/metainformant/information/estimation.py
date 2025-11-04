"""Entropy and information estimation methods with bias correction.

This module provides methods for estimating entropy and mutual information
from samples, including various bias correction techniques for small
sample sizes.
"""

from __future__ import annotations

import math
from collections import Counter
from typing import Any

import numpy as np

from metainformant.information.syntactic import (
    mutual_information,
    shannon_entropy,
    shannon_entropy_from_counts,
)


def entropy_estimator(
    counts: dict[Any, int] | list[int],
    method: str = "plugin",
    bias_correction: bool = True
) -> float:
    """Estimate entropy from counts using various methods.
    
    Args:
        counts: Dictionary of counts or list of counts
        method: Estimation method ("plugin", "miller_madow", "chao_shen", "jackknife")
        bias_correction: Whether to apply bias correction
        
    Returns:
        Entropy estimate in bits
        
    Examples:
        >>> counts = {"A": 50, "T": 30, "G": 20}
        >>> h = entropy_estimator(counts, method="miller_madow")
        >>> h > 0
        True
    """
    if isinstance(counts, dict):
        counts_list = list(counts.values())
        total = sum(counts_list)
        non_zero = sum(1 for c in counts_list if c > 0)
    else:
        counts_list = list(counts)
        total = sum(counts_list)
        non_zero = sum(1 for c in counts_list if c > 0)
    
    if total == 0:
        return 0.0
    
    if method == "plugin":
        # Plug-in estimator (standard Shannon entropy)
        return shannon_entropy_from_counts(counts)
    
    elif method == "miller_madow":
        # Miller-Madow bias correction
        plugin_est = shannon_entropy_from_counts(counts)
        if bias_correction:
            # Bias ≈ (m-1)/(2N) where m is number of symbols
            bias = (non_zero - 1) / (2 * total)
            return plugin_est + bias
        return plugin_est
    
    elif method == "chao_shen":
        # Chao-Shen estimator with Horvitz-Thompson adjustment
        if isinstance(counts, dict):
            counts_dict = counts
        else:
            # Convert to dict with indices
            counts_dict = {i: c for i, c in enumerate(counts) if c > 0}
        
        # Estimate coverage
        n = total
        m = len(counts_dict)
        coverage = 1.0 - (m == 0)
        
        if m > 0:
            # Coverage estimate: C = 1 - f1/n where f1 is singletons
            f1 = sum(1 for c in counts_dict.values() if c == 1)
            coverage = 1.0 - f1 / n if n > 0 else 1.0
        
        # Chao-Shen estimator
        entropy = 0.0
        for count in counts_dict.values():
            if count > 0:
                p = count / n
                if coverage > 0:
                    p_adj = coverage * p
                    if p_adj > 0:
                        entropy -= p_adj * math.log2(p_adj) / coverage
        
        return entropy
    
    elif method == "jackknife":
        # Jackknife bias correction (simplified)
        plugin_est = shannon_entropy_from_counts(counts)
        if not bias_correction or total < 2:
            return plugin_est
        
        # Simplified jackknife: subtract approximate bias
        bias = (non_zero - 1) / (2 * total)
        return plugin_est - bias
    
    else:
        raise ValueError(f"Unknown method: {method}")


def mutual_information_estimator(
    x: list[Any],
    y: list[Any],
    method: str = "plugin",
    bias_correction: bool = True
) -> float:
    """Estimate mutual information with bias correction.
    
    Args:
        x: Sequence of X values
        y: Sequence of Y values
        method: Estimation method ("plugin", "shrinkage", "jackknife")
        bias_correction: Whether to apply bias correction
        
    Returns:
        Mutual information estimate in bits
        
    Examples:
        >>> x = [0, 1, 0, 1, 0, 1]
        >>> y = [0, 1, 0, 1, 0, 1]  # Perfect correlation
        >>> mi = mutual_information_estimator(x, y)
        >>> mi > 0.5  # Should be high
        True
    """
    if len(x) != len(y):
        raise ValueError("X and Y must have the same length")
    
    if method == "plugin":
        # Standard plug-in estimator
        mi = mutual_information(x, y)
        
        if bias_correction:
            # Simple bias correction for MI
            # Bias tends to be positive for small samples
            n = len(x)
            x_unique = len(set(x))
            y_unique = len(set(y))
            
            # Approximate bias: (|X|-1)(|Y|-1)/(2N ln(2))
            if n > 0:
                bias = (x_unique - 1) * (y_unique - 1) / (2 * n * math.log(2))
                mi = max(0.0, mi - bias)
        
        return mi
    
    elif method == "shrinkage":
        # Shrinkage estimator (simplified)
        plugin_mi = mutual_information(x, y)
        
        if not bias_correction:
            return plugin_mi
        
        # Shrinkage towards zero
        n = len(x)
        lambda_shrink = 1.0 / (1.0 + n * 0.1)  # Simple shrinkage factor
        
        return lambda_shrink * plugin_mi
    
    elif method == "jackknife":
        # Jackknife estimator (simplified)
        plugin_mi = mutual_information(x, y)
        
        if not bias_correction or len(x) < 2:
            return plugin_mi
        
        # Simplified jackknife bias correction
        n = len(x)
        x_unique = len(set(x))
        y_unique = len(set(y))
        bias = (x_unique - 1) * (y_unique - 1) / (2 * n * math.log(2))
        
        return max(0.0, plugin_mi - bias)
    
    else:
        raise ValueError(f"Unknown method: {method}")


def kl_divergence_estimator(
    p_samples: list[Any],
    q_samples: list[Any],
    method: str = "plugin",
    bias_correction: bool = True
) -> float:
    """Estimate KL divergence from samples with bias correction.
    
    Args:
        p_samples: Samples from distribution P
        q_samples: Samples from distribution Q
        method: Estimation method ("plugin")
        bias_correction: Whether to apply bias correction
        
    Returns:
        KL divergence estimate in bits
    """
    if len(p_samples) != len(q_samples):
        raise ValueError("P and Q samples must have the same length")
    
    # Convert to probability distributions
    p_counts = Counter(p_samples)
    q_counts = Counter(q_samples)
    
    # Get all unique values
    all_values = set(p_counts.keys()) | set(q_counts.keys())
    
    # Create probability distributions
    p_total = sum(p_counts.values())
    q_total = sum(q_counts.values())
    
    if p_total == 0 or q_total == 0:
        return 0.0
    
    # Calculate KL divergence
    kl = 0.0
    for value in all_values:
        p_val = p_counts.get(value, 0) / p_total
        q_val = q_counts.get(value, 0) / q_total
        
        if p_val > 0:
            if q_val == 0:
                return float('inf')
            kl += p_val * math.log2(p_val / q_val)
    
    if bias_correction and method == "plugin":
        # Simple bias correction
        n_p = len(p_samples)
        n_q = len(q_samples)
        n_p_unique = len(p_counts)
        n_q_unique = len(q_counts)
        
        # Approximate bias correction
        bias = (n_p_unique - 1) / (2 * n_p * math.log(2))
        kl = max(0.0, kl - bias)
    
    return kl


def bias_correction(
    estimate: float,
    n_samples: int,
    n_bins: int,
    measure: str = "entropy"
) -> float:
    """Apply bias correction to entropy/MI estimates.
    
    Args:
        estimate: Uncorrected estimate
        n_samples: Number of samples
        n_bins: Number of bins/symbols
        measure: Type of measure ("entropy", "mutual_information")
        
    Returns:
        Bias-corrected estimate
    """
    if n_samples == 0:
        return 0.0
    
    if measure == "entropy":
        # Miller-Madow bias: (m-1)/(2N)
        bias = (n_bins - 1) / (2 * n_samples)
        return estimate + bias
    
    elif measure == "mutual_information":
        # MI bias correction (simplified)
        # Bias ≈ (|X|-1)(|Y|-1)/(2N ln(2))
        # Here we use n_bins as proxy for |X|×|Y|
        bias = (n_bins - 1) / (2 * n_samples * math.log(2))
        return max(0.0, estimate - bias)
    
    else:
        return estimate

