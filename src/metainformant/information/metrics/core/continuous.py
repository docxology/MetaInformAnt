"""Continuous information theory measures for numerical data.

This module implements information-theoretic measures for continuous-valued
data including differential entropy, continuous mutual information, and
continuous divergence measures.
"""

from __future__ import annotations

import math
from typing import Optional

import numpy as np
from scipy import special, stats
from scipy.spatial import cKDTree

from metainformant.core.data import validation
from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def differential_entropy(
    samples: np.ndarray, method: str = "histogram", bins: Optional[int] = None
) -> float:
    """Calculate differential entropy of continuous data, in nats.

    A 1D input is treated as samples of a single continuous variable. A 2D
    input of shape ``(n_samples, n_features)`` is treated as joint samples of
    an ``n_features``-dimensional distribution and the JOINT differential
    entropy is estimated in the full ``n_features``-dimensional space; inputs
    are never flattened.

    Estimators (all return nats):
      - ``"histogram"``: equipartition histogram estimate
        ``H = -sum_i p_i * log(p_i / V)`` with bin probability masses ``p_i``
        and bin hypervolume ``V`` (Sturges' rule when ``bins`` is None).
      - ``"kde"``: leave-one-out Gaussian kernel density estimate (Scott
        bandwidth).
      - ``"knn"``: Kozachenko-Leonenko (KSG-style) k-nearest-neighbour
        estimator ``H = psi(n) - psi(k) + log(c_d) + (d/n) * sum_i log(eps_i)``
        with unit-ball volume ``c_d = pi^(d/2) / Gamma(d/2 + 1)`` and ``eps_i``
        the distance from sample ``i`` to its k-th nearest neighbour (k = 3).

    Args:
        samples: 1D array of samples, or 2D (n_samples, n_features) joint samples
        method: Estimation method ('histogram', 'kde', 'knn')
        bins: Number of bins per dimension for histogram method (auto if None)

    Returns:
        Differential entropy estimate in nats (can be negative for strongly
        concentrated distributions).

    Raises:
        ValueError: If samples is not 1D/2D or has insufficient data
    """
    validation.validate_type(samples, np.ndarray, "samples")

    arr = np.asarray(samples, dtype=float)
    if arr.ndim == 1:
        joint = arr.reshape(-1, 1)
    elif arr.ndim == 2:
        joint = arr
    else:
        raise ValueError(
            "samples must be a 1D array of samples or a 2D (n_samples, n_features) joint array"
        )

    n_samples = joint.shape[0]
    if n_samples < 10:
        raise ValueError("Need at least 10 samples for entropy estimation")

    if method == "histogram":
        return _histogram_entropy_nd(joint, bins)
    elif method == "kde":
        return _kde_entropy_nd(joint)
    elif method == "knn":
        return _knn_entropy_nd(joint)
    else:
        raise ValueError(f"Unknown method: {method}")


def _differential_entropy_histogram(
    samples: np.ndarray, bins: Optional[int] = None
) -> float:
    """Histogram differential entropy estimator for 1D data (nats).

    Equipartition histogram estimate ``H = -sum_i p_i * log(p_i / dx)`` where
    ``p_i`` are the bin probability masses and ``dx`` the bin width (Sturges'
    rule when ``bins`` is None). Expects a 1D sample array; joint (2D) data
    must go through :func:`differential_entropy` so the estimate is computed
    in the true joint space.
    """
    samples = np.asarray(samples, dtype=float)
    if samples.ndim != 1:
        raise ValueError(
            "_differential_entropy_histogram expects a 1D sample array; "
            "pass joint data as a 2D (n_samples, n_features) array to differential_entropy"
        )
    return _histogram_entropy_nd(samples.reshape(-1, 1), bins)


def _histogram_entropy_nd(joint: np.ndarray, bins: Optional[int] = None) -> float:
    """Equipartition histogram differential entropy in joint space (nats).

    ``joint`` has shape (n_samples, d); each row is one d-dimensional joint
    sample. With equal-width bins the estimator is
    ``H = -sum_i p_i * log(p_i / V)`` where ``p_i`` is the probability mass of
    bin i and ``V = prod_j dx_j`` is the bin hypervolume.
    """
    joint = np.asarray(joint, dtype=float)
    n_samples = joint.shape[0]

    if bins is None:
        # Sturges' rule for the number of bins per dimension
        bins = int(np.ceil(np.log2(n_samples) + 1))

    counts, edges = np.histogramdd(joint, bins=bins)
    masses = counts[counts > 0] / n_samples

    if masses.size == 0:
        return 0.0

    bin_volume = 1.0
    for dim_edges in edges:
        bin_volume *= dim_edges[1] - dim_edges[0]

    return float(-np.sum(masses * np.log(masses / bin_volume)))


def _differential_entropy_kde(samples: np.ndarray) -> float:
    """Leave-one-out Gaussian KDE differential entropy estimator for 1D data (nats).

    Expects a 1D sample array; joint (2D) data must go through
    :func:`differential_entropy` so the estimate is computed in the true
    joint space.
    """
    samples = np.asarray(samples, dtype=float)
    if samples.ndim != 1:
        raise ValueError(
            "_differential_entropy_kde expects a 1D sample array; "
            "pass joint data as a 2D (n_samples, n_features) array to differential_entropy"
        )
    return _kde_entropy_nd(samples.reshape(-1, 1))


def _kde_entropy_nd(joint: np.ndarray) -> float:
    """Leave-one-out Gaussian kernel-density differential entropy in joint space (nats).

    ``joint`` has shape (n_samples, d). The density at each sample is
    estimated with a product Gaussian kernel using Scott's per-dimension
    bandwidth ``h_j = std_j * n ** (-1 / (d + 4))``; the leave-one-out
    density (excluding the sample's own kernel) gives the cross-validation
    entropy estimate ``H = -(1/n) * sum_i log f_{-i}(x_i)``, which is well
    defined in any dimension.
    """
    joint = np.asarray(joint, dtype=float)
    n_samples, n_dims = joint.shape

    if n_samples < 2:
        raise ValueError("Need at least 2 samples for KDE entropy estimation")

    std = joint.std(axis=0, ddof=1)
    if np.any(std <= 0):
        raise ValueError("KDE entropy requires variation in every dimension")

    bandwidth = std * n_samples ** (-1.0 / (n_dims + 4))
    norm = 1.0 / ((2.0 * math.pi) ** (n_dims / 2.0) * float(np.prod(bandwidth)))

    diffs = (joint[:, None, :] - joint[None, :, :]) / bandwidth
    kernels = norm * np.exp(-0.5 * np.sum(diffs**2, axis=-1))
    density = kernels.mean(axis=1)

    # Leave-one-out density: f_{-i}(x_i) = (n * f(x_i) - phi_h(0)) / (n - 1)
    loo_density = (n_samples * density - norm) / (n_samples - 1)

    return float(-np.mean(np.log(loo_density)))


def _differential_entropy_knn(samples: np.ndarray, k: int = 3) -> float:
    """Kozachenko-Leonenko (KSG) k-NN differential entropy estimator for 1D data (nats).

    Expects a 1D sample array; joint (2D) data must go through
    :func:`differential_entropy` so the estimator runs in the true joint
    space.
    """
    samples = np.asarray(samples, dtype=float)
    if samples.ndim != 1:
        raise ValueError(
            "_differential_entropy_knn expects a 1D sample array; "
            "pass joint data as a 2D (n_samples, n_features) array to differential_entropy"
        )
    return _knn_entropy_nd(samples.reshape(-1, 1), k)


def _knn_entropy_nd(joint: np.ndarray, k: int = 3) -> float:
    """Kozachenko-Leonenko (KSG-style) k-nearest-neighbour entropy in joint space (nats).

    ``joint`` has shape (n_samples, d); each row is one d-dimensional joint
    sample. Implements the KL/KSG estimator::

        H = psi(n) - psi(k) + log(c_d) + (d/n) * sum_i log(eps_i)

    where ``psi`` is the digamma function, ``c_d = pi^(d/2) / Gamma(d/2 + 1)``
    is the volume of the d-dimensional unit ball, and ``eps_i`` is the
    Euclidean distance from sample i to its k-th nearest neighbour
    (k >= 1, excluding the sample itself). The data are assumed to have a
    non-degenerate d-dimensional distribution.

    Duplicate samples carry no distance information; when they occur a tiny
    deterministic jitter (uniform on ``[0, 1e-10] * data range``, fixed seed)
    is applied before the distance computation, following standard KSG
    practice.
    """
    joint = np.asarray(joint, dtype=float)
    n_samples, n_dims = joint.shape

    if k < 1:
        raise ValueError("k must be >= 1 for k-NN entropy estimation")
    if n_samples <= k:
        raise ValueError(f"Need more than {k} samples for k-NN entropy estimation")

    distances = cKDTree(joint).query(joint, k=k + 1)[0][:, k]  # k-th NN, self excluded

    if np.any(distances <= 0.0):
        # Duplicate samples: apply the standard tiny KSG jitter (deterministic).
        rng = np.random.default_rng(0)
        scale = 1e-10 * max(float(np.ptp(joint)), 1.0)
        jittered = joint + rng.uniform(0.0, scale, size=joint.shape)
        distances = cKDTree(jittered).query(jittered, k=k + 1)[0][:, k]

    log_unit_ball_volume = 0.5 * n_dims * math.log(math.pi) - math.lgamma(
        0.5 * n_dims + 1.0
    )
    entropy = (
        special.digamma(n_samples)
        - special.digamma(k)
        + log_unit_ball_volume
        + (n_dims / n_samples) * float(np.sum(np.log(distances)))
    )
    return float(entropy)


def mutual_information_continuous(
    x: np.ndarray, y: np.ndarray, method: str = "histogram", bins: Optional[int] = None
) -> float:
    """Calculate mutual information between two continuous variables, in nats.

    The joint entropy ``H(X,Y)`` is estimated in the true two-dimensional
    joint space (the input is never flattened to 1D).

    Args:
        x: Samples from first variable
        y: Samples from second variable
        method: Estimation method ('histogram', 'kde', 'knn')
        bins: Number of bins per dimension for histogram method (auto if None)

    Returns:
        Mutual information estimate in nats (clipped at 0; first-order
        estimation bias can push the raw estimate slightly negative)

    Raises:
        ValueError: If input arrays have different lengths
    """
    x_arr = np.asarray(x, dtype=float).ravel()
    y_arr = np.asarray(y, dtype=float).ravel()

    if len(x_arr) != len(y_arr):
        raise ValueError("Input arrays must have the same length")

    if len(x_arr) < 10:
        raise ValueError("Need at least 10 samples for MI estimation")

    # MI = H(X) + H(Y) - H(X,Y)
    h_x = differential_entropy(x_arr, method=method, bins=bins)
    h_y = differential_entropy(y_arr, method=method, bins=bins)

    # Joint entropy in the true 2D joint space
    xy = np.column_stack([x_arr, y_arr])
    h_xy = differential_entropy(xy, method=method, bins=bins)

    mi = h_x + h_y - h_xy
    return max(
        0.0, mi
    )  # MI >= 0 by definition; estimation bias can go slightly negative


def kl_divergence_continuous(
    p_samples: np.ndarray,
    q_samples: np.ndarray,
    method: str = "histogram",
    bins: Optional[int] = None,
) -> float:
    """Calculate KL divergence D_KL(P||Q) between two continuous distributions, in nats.

    Args:
        p_samples: Samples from distribution P
        q_samples: Samples from distribution Q
        method: Estimation method ('histogram', 'kde')
        bins: Number of bins for histogram method

    Returns:
        KL divergence estimate D_KL(P||Q) in nats

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
    p_samples: np.ndarray, q_samples: np.ndarray, bins: Optional[int] = None
) -> float:
    """Estimate KL divergence using histogram method."""
    if bins is None:
        # Use combined data range for bins
        all_samples = np.concatenate([p_samples, q_samples])
        bins = int(np.ceil(np.log2(len(all_samples)) + 1))

    # Create histograms
    combined_range = (
        min(p_samples.min(), q_samples.min()),
        max(p_samples.max(), q_samples.max()),
    )

    hist_p, bin_edges = np.histogram(
        p_samples, bins=bins, range=combined_range, density=True
    )
    hist_q, _ = np.histogram(q_samples, bins=bins, range=combined_range, density=True)

    # Avoid division by zero and log of zero
    hist_p = np.maximum(hist_p, np.finfo(float).eps)
    hist_q = np.maximum(hist_q, np.finfo(float).eps)

    # KL divergence: sum(p * log(p/q))
    kl_div = np.sum(hist_p * np.log(hist_p / hist_q))

    return max(0.0, float(kl_div))  # Ensure non-negative


def _kl_divergence_kde(p_samples: np.ndarray, q_samples: np.ndarray) -> float:
    """Estimate KL divergence using KDE method."""
    try:
        from sklearn.neighbors import KernelDensity
    except ImportError:
        raise ImportError("scikit-learn required for KDE KL divergence")

    # Fit KDEs
    kde_p = KernelDensity(bandwidth="scott", kernel="gaussian")
    kde_q = KernelDensity(bandwidth="scott", kernel="gaussian")

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
    kl_div = np.trapezoid(integrand, dx=grid_spacing)

    return max(0.0, float(kl_div))


def entropy_estimation(
    samples: np.ndarray, method: str = "histogram", bins: Optional[int] = None
) -> float:
    """Unified interface for entropy estimation in nats (alias for differential_entropy)."""
    return differential_entropy(samples, method=method, bins=bins)


def copula_entropy(samples: np.ndarray, method: str = "histogram") -> float:
    """Calculate copula entropy in nats (normalized entropy for dependence analysis).

    Copula entropy measures the dependence between variables while being
    invariant to monotonic transformations; it equals the negative multiinformation
    ``H_joint - sum(H_individual)`` of the copula-transformed data, estimated
    in the true joint space.

    Args:
        samples: 2D array (n_samples, n_variables)
        method: Estimation method

    Returns:
        Copula entropy value in nats (<= 0; 0 for independent variables)

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
        ranks = stats.rankdata(samples[:, i], method="average")
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
    x: np.ndarray, y: np.ndarray, lag: int = 1, method: str = "histogram"
) -> float:
    """Calculate transfer entropy for continuous time series, in nats.

    Args:
        x: Source time series
        y: Target time series
        lag: Time lag
        method: Entropy estimation method

    Returns:
        Transfer entropy T(X→Y) in nats (clipped at 0)

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

    y_future = y[lag:]  # Y_{t+1}
    y_past = y[:-lag]  # Y_t
    x_past = x[:-lag]  # X_t

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
    x: np.ndarray, y: np.ndarray, method: str = "histogram", bins: Optional[int] = None
) -> float:
    """Calculate conditional entropy H(X|Y) for continuous variables, in nats.

    ``H(X|Y) = H(X,Y) - H(Y)`` with the joint entropy estimated in the true
    two-dimensional joint space (the input is never flattened to 1D).

    Args:
        x: Samples from X
        y: Samples from Y
        method: Estimation method ('histogram', 'kde', 'knn')
        bins: Number of bins per dimension for histogram method (auto if None)

    Returns:
        Conditional entropy estimate in nats (clipped at 0)
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
    bins: Optional[int] = None,
) -> float:
    """Calculate conditional entropy H(X|Y,Z) for continuous variables, in nats.

    ``H(X|Y,Z) = H(X,Y,Z) - H(Y,Z)`` with the joint entropies estimated in
    the true three- and two-dimensional joint spaces (inputs are never
    flattened to 1D).

    Args:
        x: Samples from X
        y: Samples from Y
        z: Samples from Z
        method: Estimation method ('histogram', 'kde', 'knn')
        bins: Number of bins per dimension for histogram method (auto if None)

    Returns:
        Conditional entropy estimate in nats (clipped at 0)
    """
    # H(X|Y,Z) = H(X,Y,Z) - H(Y,Z)
    xyz = np.column_stack([x, y, z])
    h_xyz = differential_entropy(xyz, method=method, bins=bins)

    yz = np.column_stack([y, z])
    h_yz = differential_entropy(yz, method=method, bins=bins)

    return max(0.0, h_xyz - h_yz)


def information_flow_network(
    time_series_data: np.ndarray, lag: int = 1, method: str = "histogram"
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
                    time_series_data[i], time_series_data[j], lag=lag, method=method
                )
                flow_matrix[i, j] = te

    return flow_matrix
