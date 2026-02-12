"""Channel capacity and rate-distortion functions for information theory.

This module implements channel capacity computation via the Blahut-Arimoto algorithm,
rate-distortion function computation, the information bottleneck method, and analytical
capacity formulas for standard noisy channel models.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from metainformant.core.data import validation
from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

# Import shannon_entropy from syntactic module
try:
    from metainformant.information.metrics.core.syntactic import shannon_entropy
except ImportError:
    logger.warning("Could not import shannon_entropy from syntactic module")

    def shannon_entropy(probs: Sequence[float], base: float = 2.0) -> float:
        """Fallback Shannon entropy implementation."""
        p = np.array(probs, dtype=float)
        p = p[p > 0]
        return float(-np.sum(p * np.log(p) / np.log(base)))


def channel_capacity(
    transition_matrix: np.ndarray,
    max_iter: int = 100,
    tol: float = 1e-8,
) -> Dict[str, Any]:
    """Compute channel capacity using the Blahut-Arimoto algorithm.

    The Blahut-Arimoto algorithm iteratively computes the channel capacity
    C = max_{p(x)} I(X;Y) for a discrete memoryless channel defined by
    its transition probability matrix P(Y|X).

    Args:
        transition_matrix: Channel transition matrix P(Y|X) of shape (n_x, n_y).
            Each row must sum to 1.0 and all entries must be non-negative.
        max_iter: Maximum number of iterations for the algorithm.
        tol: Convergence tolerance for capacity estimate.

    Returns:
        Dictionary containing:
            - "capacity": Channel capacity in bits.
            - "optimal_input": Optimal input distribution as a list of floats.
            - "n_iterations": Number of iterations performed.

    Raises:
        ValueError: If transition_matrix has invalid shape, negative entries,
            or rows that do not sum to 1.
    """
    validation.validate_type(transition_matrix, np.ndarray, "transition_matrix")

    if transition_matrix.ndim != 2:
        raise ValueError("transition_matrix must be a 2D array")

    n_x, n_y = transition_matrix.shape

    if n_x == 0 or n_y == 0:
        raise ValueError("transition_matrix must have at least one row and one column")

    if np.any(transition_matrix < 0):
        raise ValueError("transition_matrix entries must be non-negative")

    row_sums = np.sum(transition_matrix, axis=1)
    if not np.allclose(row_sums, 1.0, atol=1e-6):
        raise ValueError(f"Each row of transition_matrix must sum to 1.0, " f"got row sums: {row_sums}")

    logger.debug(
        "Running Blahut-Arimoto on %dx%d channel matrix (max_iter=%d, tol=%g)",
        n_x,
        n_y,
        max_iter,
        tol,
    )

    # Initialize uniform input distribution
    p_x = np.ones(n_x) / n_x

    capacity = 0.0
    n_iterations = 0

    for iteration in range(max_iter):
        # Step 1: Compute output distribution p(y) = sum_x p(x) * P(y|x)
        p_y = p_x @ transition_matrix  # shape (n_y,)

        # Avoid division by zero in output distribution
        p_y = np.maximum(p_y, 1e-300)

        # Step 2: Compute the c(x) = prod_y [P(y|x) / p(y)]^{P(y|x)} for each x
        # Equivalently, log c(x) = sum_y P(y|x) * log(P(y|x) / p(y))
        # This is the KL divergence D(P(y|x) || p(y)) for each input x
        log_c = np.zeros(n_x)
        for i in range(n_x):
            for j in range(n_y):
                if transition_matrix[i, j] > 0:
                    log_c[i] += transition_matrix[i, j] * math.log2(transition_matrix[i, j] / p_y[j])

        # Step 3: Update input distribution
        # p_x(new) proportional to p_x(old) * 2^{log_c}
        c_values = np.power(2.0, log_c)
        p_x_new = p_x * c_values
        normalization = np.sum(p_x_new)

        if normalization > 0:
            p_x_new = p_x_new / normalization
        else:
            # Fallback to uniform if normalization fails
            p_x_new = np.ones(n_x) / n_x

        # Compute capacity estimate: C = log2(sum_x p_x * c_x)
        capacity_new = math.log2(normalization) if normalization > 0 else 0.0

        n_iterations = iteration + 1

        # Check convergence
        if abs(capacity_new - capacity) < tol:
            capacity = capacity_new
            p_x = p_x_new
            logger.debug(
                "Blahut-Arimoto converged after %d iterations (capacity=%.6f bits)",
                n_iterations,
                capacity,
            )
            break

        capacity = capacity_new
        p_x = p_x_new

    else:
        logger.debug(
            "Blahut-Arimoto reached max_iter=%d (capacity=%.6f bits)",
            max_iter,
            capacity,
        )

    return {
        "capacity": float(max(0.0, capacity)),
        "optimal_input": p_x.tolist(),
        "n_iterations": n_iterations,
    }


def rate_distortion(
    source_dist: Sequence[float],
    distortion_matrix: np.ndarray,
    max_rate: Optional[float] = None,
    n_points: int = 50,
) -> Dict[str, Any]:
    """Compute the rate-distortion function R(D) for a discrete source.

    Uses the Blahut-Arimoto algorithm for rate-distortion to compute the
    minimum rate required to represent a source within a given average
    distortion level.

    Args:
        source_dist: Source probability distribution P(X). Must sum to 1.0.
        distortion_matrix: Distortion measure d(x, x_hat) of shape
            (n_x, n_x_hat). Entry (i, j) is the distortion incurred when
            source symbol i is reproduced as symbol j.
        max_rate: Maximum rate to compute (in bits). If None, defaults to
            the source entropy H(X).
        n_points: Number of points to compute on the R(D) curve.

    Returns:
        Dictionary containing:
            - "rates": List of rate values (bits).
            - "distortions": List of corresponding distortion values.
            - "rd_curve": List of (rate, distortion) tuples.

    Raises:
        ValueError: If source_dist does not sum to 1 or dimensions mismatch.
    """
    validation.validate_type(source_dist, (list, tuple, np.ndarray), "source_dist")
    validation.validate_type(distortion_matrix, np.ndarray, "distortion_matrix")

    p_x = np.array(source_dist, dtype=float)

    if np.any(p_x < 0):
        raise ValueError("source_dist entries must be non-negative")
    if not np.isclose(np.sum(p_x), 1.0, atol=1e-6):
        raise ValueError("source_dist must sum to 1.0")

    if distortion_matrix.ndim != 2:
        raise ValueError("distortion_matrix must be a 2D array")

    n_x = len(p_x)
    n_x_hat = distortion_matrix.shape[1]

    if distortion_matrix.shape[0] != n_x:
        raise ValueError(
            f"distortion_matrix first dimension ({distortion_matrix.shape[0]}) "
            f"must match source_dist length ({n_x})"
        )

    if np.any(distortion_matrix < 0):
        raise ValueError("distortion_matrix entries must be non-negative")

    logger.debug(
        "Computing rate-distortion function for %d-symbol source with %d reproduction symbols",
        n_x,
        n_x_hat,
    )

    # Compute source entropy as upper bound on rate
    source_entropy = float(shannon_entropy(p_x.tolist(), base=2.0))
    if max_rate is None:
        max_rate = source_entropy

    # Maximum distortion (rate = 0): D_max = min_j sum_i p(x_i) * d(x_i, j)
    expected_distortions = np.array([np.sum(p_x * distortion_matrix[:, j]) for j in range(n_x_hat)])
    d_max = float(np.min(expected_distortions))

    # Minimum distortion (maximum rate): D_min = 0 if identity reproduction is possible
    d_min = 0.0

    # Sweep over beta (Lagrange multiplier) values to trace R(D) curve
    # Higher beta -> lower distortion, higher rate
    # Lower beta -> higher distortion, lower rate
    rates: List[float] = []
    distortions: List[float] = []

    # Use logarithmically spaced beta values for good coverage
    beta_values = np.concatenate(
        [
            np.array([0.01]),
            np.logspace(-1, 2, n_points - 2),
            np.array([1000.0]),
        ]
    )

    for beta in beta_values:
        # Blahut-Arimoto for rate-distortion
        # Initialize reproduction distribution q(x_hat|x) uniformly
        q_xhat_given_x = np.ones((n_x, n_x_hat)) / n_x_hat

        for _ in range(200):
            # Compute marginal reproduction distribution q(x_hat)
            q_xhat = np.zeros(n_x_hat)
            for j in range(n_x_hat):
                q_xhat[j] = np.sum(p_x * q_xhat_given_x[:, j])
            q_xhat = np.maximum(q_xhat, 1e-300)

            # Update conditional distribution
            q_new = np.zeros((n_x, n_x_hat))
            for i in range(n_x):
                for j in range(n_x_hat):
                    q_new[i, j] = q_xhat[j] * math.exp(-beta * distortion_matrix[i, j])
                row_sum = np.sum(q_new[i, :])
                if row_sum > 0:
                    q_new[i, :] /= row_sum
                else:
                    q_new[i, :] = 1.0 / n_x_hat

            # Check convergence
            if np.allclose(q_new, q_xhat_given_x, atol=1e-8):
                break
            q_xhat_given_x = q_new

        # Compute rate and distortion for this beta
        q_xhat = np.zeros(n_x_hat)
        for j in range(n_x_hat):
            q_xhat[j] = np.sum(p_x * q_xhat_given_x[:, j])
        q_xhat = np.maximum(q_xhat, 1e-300)

        rate = 0.0
        avg_distortion = 0.0
        for i in range(n_x):
            if p_x[i] <= 0:
                continue
            for j in range(n_x_hat):
                if q_xhat_given_x[i, j] > 0 and q_xhat[j] > 0:
                    rate += p_x[i] * q_xhat_given_x[i, j] * math.log2(q_xhat_given_x[i, j] / q_xhat[j])
                avg_distortion += p_x[i] * q_xhat_given_x[i, j] * distortion_matrix[i, j]

        rates.append(float(max(0.0, rate)))
        distortions.append(float(avg_distortion))

    # Sort by distortion for a proper R(D) curve
    sorted_pairs = sorted(zip(distortions, rates))
    distortions_sorted = [d for d, _ in sorted_pairs]
    rates_sorted = [r for _, r in sorted_pairs]

    # Remove duplicate distortion points (keep the one with lowest rate)
    final_distortions: List[float] = []
    final_rates: List[float] = []
    for d, r in zip(distortions_sorted, rates_sorted):
        if not final_distortions or abs(d - final_distortions[-1]) > 1e-10:
            final_distortions.append(d)
            final_rates.append(r)
        else:
            # Keep the minimum rate for this distortion level
            if r < final_rates[-1]:
                final_rates[-1] = r

    rd_curve = list(zip(final_rates, final_distortions))

    logger.debug(
        "Rate-distortion curve computed with %d points (D_min=%.4f, D_max=%.4f)",
        len(rd_curve),
        final_distortions[0] if final_distortions else 0.0,
        final_distortions[-1] if final_distortions else 0.0,
    )

    return {
        "rates": final_rates,
        "distortions": final_distortions,
        "rd_curve": rd_curve,
    }


def information_bottleneck(
    p_xy: np.ndarray,
    beta: float = 1.0,
    n_clusters: int = 2,
    max_iter: int = 100,
    tol: float = 1e-6,
    seed: int = 42,
) -> Dict[str, Any]:
    """Compute the information bottleneck compression.

    The information bottleneck method finds a compressed representation T of X
    that preserves as much information about Y as possible. It optimizes the
    Lagrangian: L = I(T;X) - beta * I(T;Y), finding the optimal tradeoff
    between compression and relevance.

    Args:
        p_xy: Joint distribution P(X,Y) as a 2D array of shape (n_x, n_y).
            Must be non-negative and sum to 1.0.
        beta: Tradeoff parameter controlling the balance between compression
            and information preservation. Higher beta preserves more
            information about Y at the cost of less compression.
        n_clusters: Number of clusters (compressed states) for T.
        max_iter: Maximum number of iterations.
        tol: Convergence tolerance.
        seed: Random seed for reproducible initialization.

    Returns:
        Dictionary containing:
            - "assignments": List of cluster assignments for each X value.
            - "I_T_X": Mutual information I(T;X) (compression cost).
            - "I_T_Y": Mutual information I(T;Y) (relevant information preserved).
            - "beta": The beta value used.

    Raises:
        ValueError: If p_xy is not a valid joint distribution or n_clusters
            exceeds the number of X values.
    """
    validation.validate_type(p_xy, np.ndarray, "p_xy")

    if p_xy.ndim != 2:
        raise ValueError("p_xy must be a 2D array")

    n_x, n_y = p_xy.shape

    if n_x == 0 or n_y == 0:
        raise ValueError("p_xy must have at least one row and one column")

    if np.any(p_xy < 0):
        raise ValueError("p_xy entries must be non-negative")

    total = np.sum(p_xy)
    if not np.isclose(total, 1.0, atol=1e-6):
        raise ValueError(f"p_xy must sum to 1.0, got {total}")

    if n_clusters > n_x:
        raise ValueError(f"n_clusters ({n_clusters}) cannot exceed number of X values ({n_x})")

    if n_clusters < 1:
        raise ValueError("n_clusters must be at least 1")

    logger.debug(
        "Running information bottleneck: n_x=%d, n_y=%d, n_clusters=%d, beta=%.3f",
        n_x,
        n_y,
        n_clusters,
        beta,
    )

    rng = np.random.RandomState(seed)

    # Compute marginals
    p_x = np.sum(p_xy, axis=1)  # P(X)
    p_y = np.sum(p_xy, axis=0)  # P(Y)

    # Compute P(Y|X) for each x
    p_y_given_x = np.zeros((n_x, n_y))
    for i in range(n_x):
        if p_x[i] > 0:
            p_y_given_x[i, :] = p_xy[i, :] / p_x[i]

    # Initialize P(T|X) randomly
    p_t_given_x = rng.dirichlet(np.ones(n_clusters), size=n_x)  # shape (n_x, n_t)

    for iteration in range(max_iter):
        # Compute P(T) = sum_x P(X=x) * P(T|X=x)
        p_t = np.zeros(n_clusters)
        for t in range(n_clusters):
            p_t[t] = np.sum(p_x * p_t_given_x[:, t])
        p_t = np.maximum(p_t, 1e-300)

        # Compute P(Y|T) = sum_x P(Y|X=x) * P(X=x) * P(T|X=x) / P(T)
        p_y_given_t = np.zeros((n_clusters, n_y))
        for t in range(n_clusters):
            if p_t[t] > 0:
                for i in range(n_x):
                    p_y_given_t[t, :] += p_x[i] * p_t_given_x[i, t] * p_y_given_x[i, :]
                p_y_given_t[t, :] /= p_t[t]
            else:
                p_y_given_t[t, :] = p_y  # fallback to marginal

        # Update P(T|X) using the IB update rule
        # P(T|X=x) proportional to P(T) * exp(-beta * D_KL(P(Y|X=x) || P(Y|T)))
        p_t_given_x_new = np.zeros((n_x, n_clusters))
        for i in range(n_x):
            for t in range(n_clusters):
                # Compute KL divergence D_KL(P(Y|X=x) || P(Y|T=t))
                kl = 0.0
                for j in range(n_y):
                    if p_y_given_x[i, j] > 0 and p_y_given_t[t, j] > 0:
                        kl += p_y_given_x[i, j] * math.log(p_y_given_x[i, j] / p_y_given_t[t, j])
                    elif p_y_given_x[i, j] > 0 and p_y_given_t[t, j] <= 0:
                        kl = float("inf")
                        break

                if math.isinf(kl):
                    p_t_given_x_new[i, t] = 0.0
                else:
                    p_t_given_x_new[i, t] = p_t[t] * math.exp(-beta * kl)

            # Normalize
            row_sum = np.sum(p_t_given_x_new[i, :])
            if row_sum > 0:
                p_t_given_x_new[i, :] /= row_sum
            else:
                # Uniform fallback
                p_t_given_x_new[i, :] = 1.0 / n_clusters

        # Check convergence
        diff = np.max(np.abs(p_t_given_x_new - p_t_given_x))
        p_t_given_x = p_t_given_x_new

        if diff < tol:
            logger.debug("Information bottleneck converged after %d iterations", iteration + 1)
            break

    # Compute hard assignments from soft P(T|X)
    assignments = np.argmax(p_t_given_x, axis=1).tolist()

    # Recompute final P(T) and P(Y|T)
    p_t = np.zeros(n_clusters)
    for t in range(n_clusters):
        p_t[t] = np.sum(p_x * p_t_given_x[:, t])
    p_t = np.maximum(p_t, 1e-300)

    p_y_given_t = np.zeros((n_clusters, n_y))
    for t in range(n_clusters):
        if p_t[t] > 1e-300:
            for i in range(n_x):
                p_y_given_t[t, :] += p_x[i] * p_t_given_x[i, t] * p_y_given_x[i, :]
            p_y_given_t[t, :] /= p_t[t]

    # Compute I(T;X) = sum_{x,t} P(X=x) * P(T=t|X=x) * log(P(T=t|X=x) / P(T=t))
    i_t_x = 0.0
    for i in range(n_x):
        if p_x[i] <= 0:
            continue
        for t in range(n_clusters):
            if p_t_given_x[i, t] > 0 and p_t[t] > 0:
                i_t_x += p_x[i] * p_t_given_x[i, t] * math.log2(p_t_given_x[i, t] / p_t[t])

    # Compute I(T;Y) = sum_{t,y} P(T=t) * P(Y=y|T=t) * log(P(Y=y|T=t) / P(Y=y))
    p_y = np.maximum(p_y, 1e-300)
    i_t_y = 0.0
    for t in range(n_clusters):
        if p_t[t] <= 1e-300:
            continue
        for j in range(n_y):
            if p_y_given_t[t, j] > 0 and p_y[j] > 0:
                i_t_y += p_t[t] * p_y_given_t[t, j] * math.log2(p_y_given_t[t, j] / p_y[j])

    return {
        "assignments": assignments,
        "I_T_X": float(max(0.0, i_t_x)),
        "I_T_Y": float(max(0.0, i_t_y)),
        "beta": float(beta),
    }


def channel_mutual_information(
    transition_matrix: np.ndarray,
    input_dist: Sequence[float],
) -> float:
    """Compute mutual information I(X;Y) for a channel with given input distribution.

    For a discrete memoryless channel with transition matrix P(Y|X) and input
    distribution P(X), computes I(X;Y) = H(Y) - H(Y|X).

    Args:
        transition_matrix: Channel transition matrix P(Y|X) of shape (n_x, n_y).
            Each row must sum to 1.0.
        input_dist: Input probability distribution P(X). Must sum to 1.0.

    Returns:
        Mutual information I(X;Y) in bits.

    Raises:
        ValueError: If dimensions mismatch, distributions are invalid,
            or matrix entries are negative.
    """
    validation.validate_type(transition_matrix, np.ndarray, "transition_matrix")
    validation.validate_type(input_dist, (list, tuple, np.ndarray), "input_dist")

    if transition_matrix.ndim != 2:
        raise ValueError("transition_matrix must be a 2D array")

    p_x = np.array(input_dist, dtype=float)
    n_x, n_y = transition_matrix.shape

    if len(p_x) != n_x:
        raise ValueError(f"input_dist length ({len(p_x)}) must match " f"transition_matrix rows ({n_x})")

    if np.any(p_x < 0):
        raise ValueError("input_dist entries must be non-negative")
    if not np.isclose(np.sum(p_x), 1.0, atol=1e-6):
        raise ValueError("input_dist must sum to 1.0")

    if np.any(transition_matrix < 0):
        raise ValueError("transition_matrix entries must be non-negative")

    row_sums = np.sum(transition_matrix, axis=1)
    if not np.allclose(row_sums, 1.0, atol=1e-6):
        raise ValueError("Each row of transition_matrix must sum to 1.0")

    # Compute output distribution P(Y) = sum_x P(X=x) * P(Y|X=x)
    p_y = p_x @ transition_matrix  # shape (n_y,)

    # H(Y) = -sum_y P(y) * log2(P(y))
    h_y = 0.0
    for j in range(n_y):
        if p_y[j] > 0:
            h_y -= p_y[j] * math.log2(p_y[j])

    # H(Y|X) = sum_x P(x) * H(Y|X=x) = -sum_x P(x) * sum_y P(y|x) * log2(P(y|x))
    h_y_given_x = 0.0
    for i in range(n_x):
        if p_x[i] <= 0:
            continue
        for j in range(n_y):
            if transition_matrix[i, j] > 0:
                h_y_given_x -= p_x[i] * transition_matrix[i, j] * math.log2(transition_matrix[i, j])

    # I(X;Y) = H(Y) - H(Y|X)
    mi = h_y - h_y_given_x
    return float(max(0.0, mi))


def noisy_channel_capacity(
    noise_level: float,
    channel_type: str = "binary_symmetric",
) -> float:
    """Compute analytical channel capacity for standard noisy channel models.

    Provides closed-form capacity expressions for well-known channel types.

    Args:
        noise_level: Channel noise parameter. Interpretation depends on
            channel_type:
            - "binary_symmetric": Crossover probability p (0 <= p <= 1).
            - "binary_erasure": Erasure probability epsilon (0 <= epsilon <= 1).
            - "awgn": Signal-to-noise ratio (SNR) in linear scale (>= 0).
        channel_type: Type of noisy channel. One of "binary_symmetric",
            "binary_erasure", or "awgn".

    Returns:
        Channel capacity in bits per channel use.

    Raises:
        ValueError: If noise_level is out of valid range or channel_type
            is unknown.
    """
    validation.validate_type(channel_type, str, "channel_type")

    if channel_type == "binary_symmetric":
        # Binary Symmetric Channel (BSC): C = 1 - H(p)
        p = noise_level
        if p < 0.0 or p > 1.0:
            raise ValueError(
                f"For binary_symmetric channel, noise_level (crossover probability) " f"must be in [0, 1], got {p}"
            )

        if p == 0.0 or p == 1.0:
            return 1.0

        # H(p) = -p*log2(p) - (1-p)*log2(1-p)
        h_p = -p * math.log2(p) - (1.0 - p) * math.log2(1.0 - p)
        return float(max(0.0, 1.0 - h_p))

    elif channel_type == "binary_erasure":
        # Binary Erasure Channel (BEC): C = 1 - epsilon
        epsilon = noise_level
        if epsilon < 0.0 or epsilon > 1.0:
            raise ValueError(
                f"For binary_erasure channel, noise_level (erasure probability) " f"must be in [0, 1], got {epsilon}"
            )

        return float(1.0 - epsilon)

    elif channel_type == "awgn":
        # Additive White Gaussian Noise (AWGN): C = 0.5 * log2(1 + SNR)
        snr = noise_level
        if snr < 0.0:
            raise ValueError(f"For awgn channel, noise_level (SNR) must be >= 0, got {snr}")

        return float(0.5 * math.log2(1.0 + snr))

    else:
        raise ValueError(
            f"Unknown channel_type: '{channel_type}'. " f"Must be one of 'binary_symmetric', 'binary_erasure', 'awgn'."
        )
