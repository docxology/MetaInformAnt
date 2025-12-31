"""Continuous information theory measures for numerical data.

This module implements information-theoretic measures for continuous-valued
data including differential entropy, continuous mutual information, and
continuous divergence measures.
"""

from __future__ import annotations

import numpy as np
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from typing import Optional, Tuple, Union

from metainformant.core import logging, validation

logger = logging.get_logger(__name__)


def differential_entropy(
    samples: np.ndarray,
    method: str = "histogram",
    bins: Optional[int] = None
) -> float:
    """Calculate differential entropy of continuous data.

    Args:
        samples: 1D array of samples from the distribution
        method: Estimation method ('histogram', 'kde', 'knn')
        bins: Number of bins for histogram method (auto if None)

    Returns:
        Differential entropy estimate

    Raises:
        ValueError: If samples is not 1D or has insufficient data
    """
    validation.validate_type(samples, np.ndarray, "samples")

    samples = np.asarray(samples).flatten()
    n_samples = len(samples)

    if n_samples < 10:
        raise ValueError("Need at least 10 samples for entropy estimation")

    if method == "histogram":
        return _differential_entropy_histogram(samples, bins)
    elif method == "kde":
        return _differential_entropy_kde(samples)
    elif method == "knn":
        return _differential_entropy_knn(samples)
    else:
        raise ValueError(f"Unknown method: {method}")


def _differential_entropy_histogram(
    samples: np.ndarray,
    bins: Optional[int] = None
) -> float:
    """Estimate differential entropy using histogram method."""
    if bins is None:
        # Use Sturges' rule for number of bins
        bins = int(np.ceil(np.log2(len(samples)) + 1))

    # Create histogram
    hist, bin_edges = np.histogram(samples, bins=bins, density=True)

    # Remove zero bins
    hist = hist[hist > 0]
    bin_width = bin_edges[1] - bin_edges[0]

    if len(hist) == 0:
        return 0.0

    # Differential entropy: -sum(p * log(p * bin_width))
    # where p is the probability density
    probs = hist * bin_width
    probs = probs[probs > 0]  # Remove any remaining zeros

    if len(probs) == 0:
        return 0.0

    return -np.sum(probs * np.log(probs))


def _differential_entropy_kde(samples: np.ndarray) -> float:
    """Estimate differential entropy using kernel density estimation."""
    try:
        from sklearn.neighbors import KernelDensity
    except ImportError:
        raise ImportError("scikit-learn required for KDE entropy estimation")

    # Fit KDE
    kde = KernelDensity(bandwidth='scott', kernel='gaussian')
    kde.fit(samples.reshape(-1, 1))

    # Evaluate on a grid
    n_grid = 1000
    sample_range = np.ptp(samples)
    grid_min = np.min(samples) - 0.1 * sample_range
    grid_max = np.max(samples) + 0.1 * sample_range

    grid = np.linspace(grid_min, grid_max, n_grid).reshape(-1, 1)
    log_density = kde.score_samples(grid)

    # Numerical integration of -∫ p(x) log p(x) dx
    density = np.exp(log_density)
    grid_spacing = (grid_max - grid_min) / (n_grid - 1)

    # Trapezoidal integration
    integrand = -density * log_density
    entropy = np.trapz(integrand, dx=grid_spacing)

    return entropy


def _differential_entropy_knn(samples: np.ndarray, k: int = 3) -> float:
    """Estimate differential entropy using k-nearest neighbors method.

    Based on Kozachenko-Leonenko estimator.
    """
    samples = samples.reshape(-1, 1)
    n_samples = len(samples)

    if n_samples <= k:
        raise ValueError(f"Need more than {k} samples for k-NN entropy estimation")

    # Compute distances to k-th nearest neighbor
    distances = np.zeros(n_samples)

    for i in range(n_samples):
        # Distance to all other points
        diff = samples - samples[i]
        dists = np.sqrt(np.sum(diff**2, axis=1))

        # Sort distances (excluding self, which is 0)
        sorted_dists = np.sort(dists)[1:]  # Remove self-distance

        if len(sorted_dists) >= k:
            distances[i] = sorted_dists[k-1]
        else:
            distances[i] = sorted_dists[-1]  # Use furthest available

    # Kozachenko-Leonenko estimator
    # h = (d/n) * sum(log(r_i)) + log(vol_d) + gamma_constant
    # For 1D: vol_d = 2, gamma_constant ≈ 0.5772 (Euler-Mascheroni)
    d = 1  # Dimensionality
    vol_d = 2  # Volume of d-dimensional unit ball
    gamma_constant = 0.57721566490153286060651209008240243104215933593992

    if np.any(distances <= 0):
        # Handle zero distances (likely due to identical points)
        distances = np.maximum(distances, np.finfo(float).eps)

    log_distances = np.log(distances)
    entropy = (d / n_samples) * np.sum(log_distances) + np.log(vol_d) + gamma_constant

    return entropy


def mutual_information_continuous(
    x: np.ndarray,
    y: np.ndarray,
    method: str = "histogram",
    bins: Optional[int] = None
) -> float:
    """Calculate mutual information between two continuous variables.

    Args:
        x: Samples from first variable
        y: Samples from second variable
        method: Estimation method ('histogram', 'kde', 'knn')
        bins: Number of bins for histogram method

    Returns:
        Mutual information estimate

    Raises:
        ValueError: If input arrays have different lengths
    """
    x = np.asarray(x).flatten()
    y = np.asarray(y).flatten()

    if len(x) != len(y):
        raise ValueError("Input arrays must have the same length")

    if len(x) < 10:
        raise ValueError("Need at least 10 samples for MI estimation")

    # MI = H(X) + H(Y) - H(X,Y)
    h_x = differential_entropy(x, method=method, bins=bins)
    h_y = differential_entropy(y, method=method, bins=bins)

    # Joint entropy
    xy = np.column_stack([x, y])
    h_xy = differential_entropy(xy, method=method, bins=bins)

    mi = h_x + h_y - h_xy
    return max(0.0, mi)  # MI cannot be negative due to estimation errors


def kl_divergence_continuous(
    p_samples: np.ndarray,
    q_samples: np.ndarray,
    method: str = "histogram",
    bins: Optional[int] = None
) -> float:
    """Calculate KL divergence between two continuous distributions.

    Args:
        p_samples: Samples from distribution P
        q_samples: Samples from distribution Q
        method: Estimation method ('histogram', 'kde')
        bins: Number of bins for histogram method

    Returns:
        KL divergence estimate D_KL(P||Q)

    Raises:
        ValueError: If sample arrays are too small
    """
    p_samples = np.asarray(p_samples).flatten()
    q_samples = np.asarray(q_samples).flatten()

    if len(p_samples) < 10 or len(q_samples) < 10:
        raise ValueError("Need at least 10 samples for each distribution")

    if method == "histogram":
        return _kl_divergence_histogram(p_samples, q_samples, bins)
    elif method == "kde":
        return _kl_divergence_kde(p_samples, q_samples)
    else:
        raise ValueError(f"Unknown method: {method}")


def _kl_divergence_histogram(
    p_samples: np.ndarray,
    q_samples: np.ndarray,
    bins: Optional[int] = None
) -> float:
    """Estimate KL divergence using histogram method."""
    if bins is None:
        # Use combined data range for bins
        all_samples = np.concatenate([p_samples, q_samples])
        bins = int(np.ceil(np.log2(len(all_samples)) + 1))

    # Create histograms
    combined_range = (min(p_samples.min(), q_samples.min()),
                     max(p_samples.max(), q_samples.max()))

    hist_p, bin_edges = np.histogram(p_samples, bins=bins, range=combined_range, density=True)
    hist_q, _ = np.histogram(q_samples, bins=bins, range=combined_range, density=True)

    # Avoid division by zero and log of zero
    hist_p = np.maximum(hist_p, np.finfo(float).eps)
    hist_q = np.maximum(hist_q, np.finfo(float).eps)

    # KL divergence: sum(p * log(p/q))
    kl_div = np.sum(hist_p * np.log(hist_p / hist_q))

    return max(0.0, kl_div)  # Ensure non-negative


def _kl_divergence_kde(
    p_samples: np.ndarray,
    q_samples: np.ndarray
) -> float:
    """Estimate KL divergence using KDE method."""
    try:
        from sklearn.neighbors import KernelDensity
    except ImportError:
        raise ImportError("scikit-learn required for KDE KL divergence")

    # Fit KDEs
    kde_p = KernelDensity(bandwidth='scott', kernel='gaussian')
    kde_q = KernelDensity(bandwidth='scott', kernel='gaussian')

    kde_p.fit(p_samples.reshape(-1, 1))
    kde_q.fit(q_samples.reshape(-1, 1))

    # Evaluate on a grid covering both distributions
    all_samples = np.concatenate([p_samples, q_samples])
    sample_range = np.ptp(all_samples)
    grid_min = np.min(all_samples) - 0.1 * sample_range
    grid_max = np.max(all_samples) + 0.1 * sample_range

    n_grid = 1000
    grid = np.linspace(grid_min, grid_max, n_grid).reshape(-1, 1)

    # Get log densities
    log_density_p = kde_p.score_samples(grid)
    log_density_q = kde_q.score_samples(grid)

    # Convert to densities
    density_p = np.exp(log_density_p)
    density_q = np.exp(log_density_q)

    # KL divergence: ∫ p(x) log(p(x)/q(x)) dx
    # Approximate with numerical integration
    ratio = density_p / np.maximum(density_q, np.finfo(float).eps)
    integrand = density_p * np.log(np.maximum(ratio, np.finfo(float).eps))

    grid_spacing = (grid_max - grid_min) / (n_grid - 1)
    kl_div = np.trapz(integrand, dx=grid_spacing)

    return max(0.0, kl_div)


def entropy_estimation(
    samples: np.ndarray,
    method: str = "histogram",
    bins: Optional[int] = None
) -> float:
    """Unified interface for entropy estimation (alias for differential_entropy)."""
    return differential_entropy(samples, method=method, bins=bins)


def copula_entropy(
    samples: np.ndarray,
    method: str = "histogram"
) -> float:
    """Calculate copula entropy (normalized entropy for dependence analysis).

    Copula entropy measures the dependence between variables while being
    invariant to monotonic transformations.

    Args:
        samples: 2D array (n_samples, n_variables)
        method: Estimation method

    Returns:
        Copula entropy value

    Raises:
        ValueError: If input is not 2D
    """
    samples = np.asarray(samples)
    if samples.ndim != 2:
        raise ValueError("Samples must be 2D array")

    n_vars = samples.shape[1]
    if n_vars < 2:
        raise ValueError("Need at least 2 variables for copula entropy")

    # Transform to copula space using empirical CDF
    copula_data = np.zeros_like(samples)

    for i in range(n_vars):
        # Rank transformation (empirical CDF)
        ranks = stats.rankdata(samples[:, i], method='average')
        copula_data[:, i] = ranks / (len(samples) + 1)

    # Calculate entropy of copula-transformed data
    entropy_val = 0.0
    for i in range(n_vars):
        h_i = differential_entropy(copula_data[:, i], method=method)
        entropy_val += h_i

    # Joint entropy
    h_joint = differential_entropy(copula_data, method=method)

    # Copula entropy: H_joint - sum(H_individual)
    copula_ent = h_joint - entropy_val

    return copula_ent


def transfer_entropy_continuous(
    x: np.ndarray,
    y: np.ndarray,
    lag: int = 1,
    method: str = "histogram"
) -> float:
    """Calculate transfer entropy for continuous time series.

    Args:
        x: Source time series
        y: Target time series
        lag: Time lag
        method: Entropy estimation method

    Returns:
        Transfer entropy T(X→Y)

    Raises:
        ValueError: If series have different lengths or lag is invalid
    """
    x = np.asarray(x).flatten()
    y = np.asarray(y).flatten()

    if len(x) != len(y):
        raise ValueError("Time series must have the same length")

    if lag < 1:
        raise ValueError("Lag must be >= 1")

    n = len(x)
    if n <= lag + 1:
        raise ValueError("Time series too short for given lag")

    # Transfer entropy: H(Y_{t+1} | Y_t) - H(Y_{t+1} | Y_t, X_t)
    # This is equivalent to: I(Y_{t+1}; X_t | Y_t)

    y_future = y[lag:]      # Y_{t+1}
    y_past = y[:-lag]       # Y_t
    x_past = x[:-lag]       # X_t

    # H(Y_{t+1} | Y_t)
    h_y_future_given_y_past = conditional_entropy_continuous(
        y_future, y_past, method=method
    )

    # H(Y_{t+1} | Y_t, X_t)
    h_y_future_given_y_past_x_past = conditional_entropy_continuous_3d(
        y_future, y_past, x_past, method=method
    )

    te = h_y_future_given_y_past - h_y_future_given_y_past_x_past
    return max(0.0, te)  # Ensure non-negative


def conditional_entropy_continuous(
    x: np.ndarray,
    y: np.ndarray,
    method: str = "histogram",
    bins: Optional[int] = None
) -> float:
    """Calculate conditional entropy H(X|Y) for continuous variables.

    Args:
        x: Samples from X
        y: Samples from Y
        method: Estimation method
        bins: Number of bins

    Returns:
        Conditional entropy estimate
    """
    # H(X|Y) = H(X,Y) - H(Y)
    xy = np.column_stack([x, y])
    h_xy = differential_entropy(xy, method=method, bins=bins)
    h_y = differential_entropy(y, method=method, bins=bins)

    return max(0.0, h_xy - h_y)


def conditional_entropy_continuous_3d(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    method: str = "histogram",
    bins: Optional[int] = None
) -> float:
    """Calculate conditional entropy H(X|Y,Z) for continuous variables.

    Args:
        x: Samples from X
        y: Samples from Y
        z: Samples from Z
        method: Estimation method
        bins: Number of bins

    Returns:
        Conditional entropy estimate
    """
    # H(X|Y,Z) = H(X,Y,Z) - H(Y,Z)
    xyz = np.column_stack([x, y, z])
    h_xyz = differential_entropy(xyz, method=method, bins=bins)

    yz = np.column_stack([y, z])
    h_yz = differential_entropy(yz, method=method, bins=bins)

    return max(0.0, h_xyz - h_yz)


def information_flow_network(
    time_series_data: np.ndarray,
    lag: int = 1,
    method: str = "histogram"
) -> np.ndarray:
    """Calculate information flow network from multivariate time series.

    Args:
        time_series_data: 2D array (n_variables, n_timepoints)
        lag: Time lag for transfer entropy
        method: Entropy estimation method

    Returns:
        Information flow matrix (n_variables x n_variables)

    Raises:
        ValueError: If input is not 2D
    """
    time_series_data = np.asarray(time_series_data)

    if time_series_data.ndim != 2:
        raise ValueError("Time series data must be 2D (n_variables x n_timepoints)")

    n_vars = time_series_data.shape[0]
    flow_matrix = np.zeros((n_vars, n_vars))

    for i in range(n_vars):
        for j in range(n_vars):
            if i != j:  # No self-flow
                te = transfer_entropy_continuous(
                    time_series_data[i], time_series_data[j],
                    lag=lag, method=method
                )
                flow_matrix[i, j] = te

    return flow_matrix
