"""Continuous information theory methods.

This module implements information-theoretic measures for continuous
probability distributions, including differential entropy and methods
for estimating entropy from continuous data samples.
"""

from __future__ import annotations

import math
from collections import Counter
from typing import Any

import numpy as np


def differential_entropy(
    samples: np.ndarray,
    method: str = "histogram",
    bins: int | None = None
) -> float:
    """Calculate differential entropy for continuous data.
    
    Differential entropy extends Shannon entropy to continuous distributions.
    Estimated from samples using histogram or kernel density estimation.
    
    Args:
        samples: Array of continuous samples
        method: Estimation method ("histogram" or "kde")
        bins: Number of bins for histogram method (auto if None)
        
    Returns:
        Differential entropy estimate in nats (natural log base)
        
    Examples:
        >>> import numpy as np
        >>> samples = np.random.normal(0, 1, 1000)
        >>> h = differential_entropy(samples)
        >>> h > 0  # Should be positive
        True
        
    Raises:
        ValueError: If samples array is empty or method is invalid
        
    References:
        Cover, T. M., & Thomas, J. A. (2006). Elements of Information Theory.
        John Wiley & Sons.
    """
    if len(samples) == 0:
        raise ValueError("Samples array cannot be empty")
    
    valid_methods = ["histogram", "kde"]
    if method not in valid_methods:
        raise ValueError(f"Method must be one of {valid_methods}, got {method}")
    
    samples = np.asarray(samples).flatten()
    
    if method == "histogram":
        # Use histogram-based estimation
        if bins is None:
            # Sturges' rule for bin count
            bins = int(1 + math.log2(len(samples)))
        
        # Create histogram
        hist, bin_edges = np.histogram(samples, bins=bins)
        
        # Normalize to get probability density
        bin_width = bin_edges[1] - bin_edges[0]
        probs = hist / (len(samples) * bin_width)
        
        # Calculate differential entropy: h(X) = -∫ p(x) log p(x) dx
        entropy = 0.0
        for p in probs:
            if p > 0:
                entropy -= p * math.log(p) * bin_width
        
        return entropy
    
    elif method == "kde":
        # Kernel density estimation (simplified)
        # For full implementation, would use scipy.stats.gaussian_kde
        # Here we use a simple Gaussian kernel approximation
        std = np.std(samples)
        if std == 0:
            return 0.0
        
        # Approximation for Gaussian: h(X) ≈ 0.5 * log(2πeσ²)
        return 0.5 * math.log(2 * math.pi * math.e * std * std)
    
    else:
        raise ValueError(f"Unknown method: {method}")


def mutual_information_continuous(
    x: np.ndarray,
    y: np.ndarray,
    method: str = "histogram",
    bins: int | None = None
) -> float:
    """Calculate mutual information for continuous variables.
    
    Estimates MI between continuous variables using binning or KDE.
    
    Args:
        x: Array of X samples
        y: Array of Y samples (must match length of x)
        method: Estimation method ("histogram" or "kde")
        bins: Number of bins for histogram method
        
    Returns:
        Mutual information estimate in nats
        
    Examples:
        >>> import numpy as np
        >>> x = np.random.randn(1000)
        >>> y = x + np.random.randn(1000) * 0.1  # Strong correlation
        >>> mi = mutual_information_continuous(x, y)
        >>> mi > 0  # Should be positive
        True
        
    Raises:
        ValueError: If x or y are empty, have different lengths, or method is invalid
    """
    x = np.asarray(x).flatten()
    y = np.asarray(y).flatten()
    
    if len(x) == 0:
        raise ValueError("X array cannot be empty")
    
    if len(y) == 0:
        raise ValueError("Y array cannot be empty")
    
    if len(x) != len(y):
        raise ValueError(f"X and Y must have the same length, got {len(x)} and {len(y)}")
    
    if len(x) < 2:
        raise ValueError("At least 2 samples required for MI calculation")
    
    valid_methods = ["histogram", "kde"]
    if method not in valid_methods:
        raise ValueError(f"Method must be one of {valid_methods}, got {method}")
    
    if method == "histogram":
        if bins is None:
            bins = int(1 + math.log2(len(x)))
        
        # Create 2D histogram
        hist_2d, x_edges, y_edges = np.histogram2d(x, y, bins=bins)
        
        # Normalize
        bin_area = (x_edges[1] - x_edges[0]) * (y_edges[1] - y_edges[0])
        joint_probs = hist_2d / (len(x) * bin_area)
        
        # Marginal distributions
        x_probs = np.sum(joint_probs, axis=1) * (x_edges[1] - x_edges[0])
        y_probs = np.sum(joint_probs, axis=0) * (y_edges[1] - y_edges[0])
        
        # Calculate MI: I(X;Y) = ∫∫ p(x,y) log(p(x,y)/(p(x)p(y))) dx dy
        mi = 0.0
        for i in range(len(x_probs)):
            for j in range(len(y_probs)):
                if joint_probs[i, j] > 0 and x_probs[i] > 0 and y_probs[j] > 0:
                    mi += joint_probs[i, j] * bin_area * math.log(
                        joint_probs[i, j] / (x_probs[i] * y_probs[j] / bin_area)
                    )
        
        return max(0.0, mi)
    
    elif method == "kde":
        # For KDE, use differential entropy approximation
        h_x = differential_entropy(x, method="kde")
        h_y = differential_entropy(y, method="kde")
        
        # Joint entropy approximation (simplified)
        # For correlated variables, joint entropy < sum of individual entropies
        # This is a rough approximation
        h_xy = h_x + h_y - 0.5 * abs(np.corrcoef(x, y)[0, 1]) * min(h_x, h_y)
        
        return max(0.0, h_x + h_y - h_xy)
    
    else:
        raise ValueError(f"Unknown method: {method}")


def kl_divergence_continuous(
    p_samples: np.ndarray,
    q_samples: np.ndarray,
    method: str = "histogram",
    bins: int | None = None
) -> float:
    """Estimate KL divergence between continuous distributions from samples.
    
    Args:
        p_samples: Samples from distribution P
        q_samples: Samples from distribution Q
        method: Estimation method ("histogram")
        bins: Number of bins for histogram method
        
    Returns:
        KL divergence estimate in nats
        
    Examples:
        >>> import numpy as np
        >>> p_samples = np.random.normal(0, 1, 1000)
        >>> q_samples = np.random.normal(0, 1, 1000)  # Same distribution
        >>> kl = kl_divergence_continuous(p_samples, q_samples)
        >>> abs(kl) < 0.5  # Should be close to zero
        True
    """
    p_samples = np.asarray(p_samples).flatten()
    q_samples = np.asarray(q_samples).flatten()
    
    if method == "histogram":
        if bins is None:
            bins = int(1 + math.log2(len(p_samples)))
        
        # Find common range
        all_samples = np.concatenate([p_samples, q_samples])
        min_val = np.min(all_samples)
        max_val = np.max(all_samples)
        
        if max_val == min_val:
            return 0.0
        
        # Create histograms on same bins
        p_hist, _ = np.histogram(p_samples, bins=bins, range=(min_val, max_val))
        q_hist, bin_edges = np.histogram(q_samples, bins=bins, range=(min_val, max_val))
        
        # Normalize to probabilities
        bin_width = bin_edges[1] - bin_edges[0]
        p_probs = p_hist / len(p_samples)
        q_probs = q_hist / len(q_samples)
        
        # Calculate KL divergence
        kl = 0.0
        for i in range(bins):
            if p_probs[i] > 0:
                if q_probs[i] == 0:
                    return float('inf')
                kl += p_probs[i] * math.log(p_probs[i] / q_probs[i])
        
        return kl
    
    else:
        raise ValueError(f"Unknown method: {method}")


def entropy_estimation(
    samples: np.ndarray,
    method: str = "plugin",
    bins: int | None = None
) -> float:
    """Estimate entropy from continuous samples using various methods.
    
    Args:
        samples: Array of continuous samples
        method: Estimation method ("plugin", "miller_madow", "chao_shen")
        bins: Number of bins for discretization (if applicable)
        
    Returns:
        Entropy estimate in bits
        
    Examples:
        >>> import numpy as np
        >>> samples = np.random.normal(0, 1, 1000)
        >>> h = entropy_estimation(samples, method="plugin")
        >>> h > 0
        True
    """
    samples = np.asarray(samples).flatten()
    
    if len(samples) == 0:
        return 0.0
    
    if method == "plugin":
        # Plug-in estimator: discretize and use Shannon entropy
        if bins is None:
            bins = int(1 + math.log2(len(samples)))
        
        hist, _ = np.histogram(samples, bins=bins)
        probs = hist / len(samples)
        
        from metainformant.information.syntactic import shannon_entropy_from_counts
        
        return shannon_entropy_from_counts(hist)
    
    elif method == "miller_madow":
        # Miller-Madow bias correction
        # First get plug-in estimate
        plugin_est = entropy_estimation(samples, method="plugin", bins=bins)
        
        if bins is None:
            bins = int(1 + math.log2(len(samples)))
        
        hist, _ = np.histogram(samples, bins=bins)
        non_zero_bins = np.sum(hist > 0)
        
        # Bias correction: H_MM = H_plugin + (m-1)/(2N)
        # where m is number of non-zero bins
        bias_correction = (non_zero_bins - 1) / (2 * len(samples))
        
        return plugin_est + bias_correction
    
    elif method == "chao_shen":
        # Chao-Shen estimator (for discrete data, but applicable here)
        if bins is None:
            bins = int(1 + math.log2(len(samples)))
        
        hist, _ = np.histogram(samples, bins=bins)
        probs = hist / len(samples)
        
        # Horvitz-Thompson type estimator
        n = len(samples)
        entropy = 0.0
        
        for count in hist:
            if count > 0:
                p = count / n
                # Coverage adjustment
                coverage = 1 - (hist == 0).sum() / bins
                if coverage > 0:
                    p_adj = coverage * p
                    if p_adj > 0:
                        entropy -= p_adj * math.log2(p_adj) / coverage
        
        return entropy
    
    else:
        raise ValueError(f"Unknown method: {method}")

