"""Information projection, divergences, channel capacity, and related measures.

This module implements information projection (m-projection via iterative
scaling), alpha-divergences, channel capacity (Blahut-Arimoto), rate-distortion
theory, the information bottleneck, entropy power inequality, and information
dimension (Grassberger-Procaccia).
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.spatial.distance import pdist

from metainformant.core.data import validation
from metainformant.core.utils import logging

from .fisher_rao import _validate_distribution

logger = logging.get_logger(__name__)

# ---------------------------------------------------------------------------
# Internal imports from sibling modules (deferred to avoid circular deps)
# ---------------------------------------------------------------------------

try:
    from metainformant.information.metrics.core.syntactic import kl_divergence, shannon_entropy
except ImportError:  # pragma: no cover
    kl_divergence = None  # type: ignore[assignment]
    shannon_entropy = None  # type: ignore[assignment]


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _safe_log(x: np.ndarray) -> np.ndarray:
    """Element-wise natural log that returns 0 for non-positive entries.

    Args:
        x: Input array.

    Returns:
        Array of same shape with log(x) where x > 0 and 0 elsewhere.
    """
    out = np.zeros_like(x, dtype=np.float64)
    mask = x > 0
    out[mask] = np.log(x[mask])
    return out


def _normalize_distribution(p: np.ndarray, name: str) -> np.ndarray:
    """Validate and normalise a probability vector, allowing small deviations.

    Unlike ``_validate_distribution`` this helper will silently renormalise
    a vector that is close to summing to 1 (within 1e-6) rather than raising
    an error, making it suitable for algorithm inputs that may be slightly
    off due to user rounding.

    Args:
        p: Array-like probability vector.
        name: Name used in error messages.

    Returns:
        Normalised float64 array.

    Raises:
        ValueError: If negative, fewer than 2 elements, or sums to zero.
    """
    p = np.asarray(p, dtype=np.float64).ravel()
    if p.size < 2:
        raise ValueError(f"{name} must have at least 2 elements, got {p.size}")
    if np.any(p < 0):
        raise ValueError(f"{name} contains negative values")
    total = p.sum()
    if total <= 0:
        raise ValueError(f"{name} sums to zero -- not a valid distribution")
    if not np.isclose(total, 1.0, atol=1e-6):
        logger.debug(
            "%s sums to %.8f -- renormalising to 1.0",
            name,
            total,
        )
    return p / total


def _validate_transition_matrix(matrix: np.ndarray) -> np.ndarray:
    """Return a validated copy of a channel transition matrix P(Y|X).

    Rows must sum to 1, entries must be non-negative, and the matrix must
    be at least 2x2.

    Args:
        matrix: Channel matrix P(Y|X) of shape ``(n_x, n_y)``.

    Returns:
        Row-normalised float64 copy of the matrix.

    Raises:
        ValueError: On any violation.
    """
    matrix = np.asarray(matrix, dtype=np.float64)
    if matrix.ndim != 2:
        raise ValueError(f"transition_matrix must be 2-dimensional, got {matrix.ndim}D")
    n_x, n_y = matrix.shape
    if n_x < 2 or n_y < 2:
        raise ValueError(f"transition_matrix must be at least 2x2, got {n_x}x{n_y}")
    if np.any(matrix < 0):
        raise ValueError("transition_matrix contains negative entries")
    row_sums = matrix.sum(axis=1)
    if not np.allclose(row_sums, 1.0, atol=1e-6):
        bad = int(np.argmax(np.abs(row_sums - 1.0)))
        raise ValueError(f"Row {bad} of transition_matrix sums to {row_sums[bad]:.8f}, " "expected 1.0")
    # Renormalise rows to machine precision
    return matrix / row_sums[:, np.newaxis]


# ---------------------------------------------------------------------------
# 3. Information Projection (m-projection via iterative scaling)
# ---------------------------------------------------------------------------


def information_projection(
    p: Sequence[float],
    constraint_set: List[Tuple[List[int], List[float]]],
    method: str = "iterative_scaling",
    max_iter: int = 100,
    tol: float = 1e-8,
) -> Dict[str, Any]:
    """Project a distribution onto the closest point in a constraint set (KL sense).

    Performs the *m-projection* (moment projection) of distribution *p* onto
    the set of distributions that satisfy a collection of marginal
    constraints.  This is equivalent to finding

        q* = argmin_{q in C} D_KL(q || p)

    where C is defined by the constraint set.  The iterative proportional
    fitting (Sinkhorn / iterative scaling) algorithm is used.

    Args:
        p: Reference distribution (must sum to 1).
        constraint_set: List of ``(index_set, target_marginal)`` tuples.
            Each *index_set* is a list of indices into *p* and
            *target_marginal* is the desired marginal probability mass for
            those indices (must sum to less than or equal to 1 and be
            non-negative).
        method: Projection algorithm.  Currently ``"iterative_scaling"`` is
            the only supported method.
        max_iter: Maximum number of iterations (must be >= 1).
        tol: Convergence tolerance on the maximum absolute change in the
            projected distribution between successive iterations.

    Returns:
        Dictionary with keys:

        * ``projected`` (*list[float]*) -- The projected distribution q*.
        * ``kl_divergence`` (*float*) -- D_KL(q* || p).
        * ``n_iterations`` (*int*) -- Number of iterations performed.

    Raises:
        ValueError: If the inputs are invalid or the method is unknown.
    """
    if method != "iterative_scaling":
        raise ValueError(f"Unknown method '{method}'; only 'iterative_scaling' is supported")
    if max_iter < 1:
        raise ValueError(f"max_iter must be >= 1, got {max_iter}")

    validation.validate_type(p, (list, tuple, np.ndarray), "p")
    validation.validate_type(constraint_set, list, "constraint_set")

    p_arr = _validate_distribution(np.asarray(p), "p")

    if not constraint_set:
        raise ValueError("constraint_set must not be empty")

    # Validate each constraint
    for idx, (index_set, target_marginal) in enumerate(constraint_set):
        if not index_set:
            raise ValueError(f"Constraint {idx}: index_set must not be empty")
        for i in index_set:
            if i < 0 or i >= p_arr.size:
                raise ValueError(f"Constraint {idx}: index {i} out of range [0, {p_arr.size})")
        target_arr = np.asarray(target_marginal, dtype=np.float64)
        if np.any(target_arr < 0):
            raise ValueError(f"Constraint {idx}: target_marginal contains negative values")

    # Iterative proportional fitting (iterative scaling)
    q = p_arr.copy()
    n_iterations = 0

    for iteration in range(1, max_iter + 1):
        q_prev = q.copy()

        for index_set, target_marginal in constraint_set:
            indices = np.array(index_set)
            target = np.asarray(target_marginal, dtype=np.float64)

            # Current marginal mass over the index set
            current_mass = q[indices].sum()

            if current_mass > 0 and target.sum() > 0:
                # Scale the entries in the index set to match the target mass
                # For simple marginal constraints: scale proportionally
                target_total = target.sum()
                scaling_factor = target_total / current_mass
                q[indices] *= scaling_factor

        # Renormalise to ensure valid distribution
        q_sum = q.sum()
        if q_sum > 0:
            q /= q_sum

        n_iterations = iteration

        # Check convergence
        max_change = float(np.max(np.abs(q - q_prev)))
        if max_change < tol:
            logger.debug(
                "Information projection converged after %d iterations " "(max_change=%.2e)",
                iteration,
                max_change,
            )
            break

    # Compute KL divergence D_KL(q || p)
    kl_div = 0.0
    for qi, pi in zip(q, p_arr):
        if qi > 0:
            if pi > 0:
                kl_div += qi * math.log(qi / pi)
            else:
                kl_div = float("inf")
                break

    result: Dict[str, Any] = {
        "projected": q.tolist(),
        "kl_divergence": kl_div,
        "n_iterations": n_iterations,
    }

    logger.debug(
        "Information projection: KL=%.6f, iterations=%d",
        kl_div,
        n_iterations,
    )

    return result


# ---------------------------------------------------------------------------
# 4. Statistical (Alpha) Divergence
# ---------------------------------------------------------------------------


def statistical_divergence(
    p: Sequence[float],
    q: Sequence[float],
    alpha: float = 0.5,
) -> float:
    """Compute the alpha-divergence between two discrete distributions.

    The alpha-divergence is a one-parameter family that interpolates between
    well-known divergences:

    * alpha -> 0: reverse KL divergence  D_KL(q || p)
    * alpha -> 1: forward KL divergence  D_KL(p || q)
    * alpha = 0.5: closely related to the squared Hellinger distance

    The general formula (for alpha not in {-1, +1}) is:

        D_alpha(p || q) = (4 / (1 - alpha^2))
                          * (1 - sum_i p_i^{(1+alpha)/2} * q_i^{(1-alpha)/2})

    Args:
        p: First probability distribution (must sum to 1).
        q: Second probability distribution (must sum to 1, same length as *p*).
        alpha: Divergence parameter.  Must satisfy ``|alpha| != 1`` for the
            general formula.  When ``alpha`` is close to 0 or 1 the function
            falls back to KL divergence computation.

    Returns:
        Alpha-divergence value (non-negative float).

    Raises:
        ValueError: If distributions have different lengths, contain negative
            values, do not sum to 1, or ``|alpha| == 1`` exactly.
    """
    validation.validate_type(p, (list, tuple, np.ndarray), "p")
    validation.validate_type(q, (list, tuple, np.ndarray), "q")

    p_arr = _validate_distribution(np.asarray(p), "p")
    q_arr = _validate_distribution(np.asarray(q), "q")

    if p_arr.size != q_arr.size:
        raise ValueError(f"Distributions must have the same length: " f"len(p)={p_arr.size}, len(q)={q_arr.size}")

    # Handle special cases where formula has 0/0 indeterminate form
    if abs(alpha - 1.0) < 1e-12:
        # Forward KL: D_KL(p || q)
        if kl_divergence is not None:
            return kl_divergence(p_arr, q_arr, base=math.e)
        raise ValueError("alpha=1 requires KL divergence but syntactic module is unavailable")

    if abs(alpha) < 1e-12:
        # Reverse KL: D_KL(q || p)
        if kl_divergence is not None:
            return kl_divergence(q_arr, p_arr, base=math.e)
        raise ValueError("alpha=0 requires KL divergence but syntactic module is unavailable")

    if abs(abs(alpha) - 1.0) < 1e-12:
        raise ValueError(f"|alpha| must not equal 1 for the general formula, got alpha={alpha}")

    # General alpha-divergence
    exp_p = (1.0 + alpha) / 2.0
    exp_q = (1.0 - alpha) / 2.0

    # Only consider indices where both p and q are positive
    # (0^positive_power = 0, which contributes 0 to the sum)
    integral = float(np.sum(np.power(p_arr, exp_p) * np.power(q_arr, exp_q)))

    divergence = (4.0 / (1.0 - alpha * alpha)) * (1.0 - integral)

    # Clamp to non-negative (numerical errors can produce tiny negatives)
    divergence = max(0.0, divergence)

    logger.debug(
        "Alpha-divergence (alpha=%.3f): %.6f",
        alpha,
        divergence,
    )

    return divergence


# ---------------------------------------------------------------------------
# 7. Channel Capacity (Blahut-Arimoto)
# ---------------------------------------------------------------------------


def channel_capacity(
    transition_matrix: np.ndarray,
    method: str = "blahut_arimoto",
    max_iter: int = 100,
    tol: float = 1e-8,
) -> Dict[str, Any]:
    """Compute channel capacity C = max_{p(x)} I(X;Y) for a discrete memoryless channel.

    The Blahut-Arimoto algorithm iteratively refines the input distribution
    p(x) to maximise the mutual information I(X;Y) for a fixed channel
    transition matrix P(Y|X).

    Args:
        transition_matrix: Channel matrix P(Y|X) of shape ``(n_x, n_y)``
            where ``transition_matrix[i, j] = P(Y=j | X=i)``.  Rows must
            sum to 1 and all entries must be non-negative.
        method: Optimisation algorithm.  Currently ``"blahut_arimoto"`` is
            the only supported method.
        max_iter: Maximum number of iterations (must be >= 1).
        tol: Convergence tolerance on the capacity estimate between
            successive iterations.

    Returns:
        Dictionary with keys:

        * ``capacity`` (*float*) -- Channel capacity in **nats**.
        * ``optimal_input`` (*np.ndarray*) -- Capacity-achieving input
          distribution p*(x), shape ``(n_x,)``.
        * ``n_iterations`` (*int*) -- Number of iterations performed.

    Raises:
        ValueError: If the transition matrix is invalid, the method is
            unknown, or *max_iter* < 1.
    """
    if method != "blahut_arimoto":
        raise ValueError(f"Unknown method '{method}'; only 'blahut_arimoto' is supported")
    if max_iter < 1:
        raise ValueError(f"max_iter must be >= 1, got {max_iter}")

    W = _validate_transition_matrix(transition_matrix)
    n_x, n_y = W.shape

    # Initialise with uniform input distribution
    p = np.full(n_x, 1.0 / n_x)

    capacity_prev = -np.inf

    for iteration in range(1, max_iter + 1):
        # E-step: compute q(x|y) proportional to p(x) W(y|x) for each y
        # q_xy[i, j] = p(x=i) * W(y=j | x=i)
        q_xy = p[:, np.newaxis] * W  # shape (n_x, n_y)

        # Marginal q(y) = sum_x p(x) W(y|x)
        q_y = q_xy.sum(axis=0)  # shape (n_y,)
        q_y = np.maximum(q_y, 1e-300)

        # M-step: update p(x)
        # For each x, compute  sum_y W(y|x) log[ W(y|x) / q(y) ]
        log_ratio = _safe_log(W) - _safe_log(q_y)[np.newaxis, :]  # (n_x, n_y)
        r = np.sum(W * log_ratio, axis=1)  # (n_x,)

        # Update: p_new(x) = p(x) exp(r(x)) / Z
        log_p_new = np.log(np.maximum(p, 1e-300)) + r
        log_p_new -= log_p_new.max()  # numerical stability
        p_new = np.exp(log_p_new)
        p_new /= p_new.sum()

        # Convergence check
        capacity_current = float(np.dot(p_new, r))

        if abs(capacity_current - capacity_prev) < tol:
            p = p_new
            logger.debug(
                "Blahut-Arimoto converged after %d iterations (C=%.8f nats)",
                iteration,
                capacity_current,
            )
            return {
                "capacity": float(capacity_current),
                "optimal_input": p,
                "n_iterations": iteration,
            }

        capacity_prev = capacity_current
        p = p_new

    logger.debug(
        "Blahut-Arimoto did not converge in %d iterations (C=%.8f nats)",
        max_iter,
        capacity_prev,
    )
    return {
        "capacity": float(capacity_prev),
        "optimal_input": p,
        "n_iterations": max_iter,
    }


# ---------------------------------------------------------------------------
# 8. Rate-Distortion Function
# ---------------------------------------------------------------------------


def rate_distortion_function(
    source_probs: Sequence[float],
    distortion_matrix: np.ndarray,
    max_rate: float = 5.0,
    n_points: int = 50,
) -> Dict[str, Any]:
    """Compute the rate-distortion curve R(D) for a discrete memoryless source.

    Uses the parametric (Lagrange-multiplier) representation: for each value
    of the slope parameter ``beta`` the rate R and expected distortion D are
    computed via the Blahut algorithm for rate-distortion.

    Args:
        source_probs: Source probability distribution p(x) over an alphabet
            of size *n*.  Must sum to 1 and be non-negative.
        distortion_matrix: Matrix ``d[i, j]`` giving the distortion of
            reproducing source symbol *i* as reproduction symbol *j*.
            Shape ``(n, m)`` where *n* is the source alphabet size.
        max_rate: Upper bound on the rate axis (nats).  The sweep of the
            Lagrange parameter is chosen so that R stays below this value.
        n_points: Number of (R, D) points to compute along the curve.

    Returns:
        Dictionary with keys:

        * ``rates`` (*list[float]*) -- Rate values (nats) in decreasing order.
        * ``distortions`` (*list[float]*) -- Corresponding distortion values
          in increasing order.
        * ``critical_distortion`` (*float*) -- D_max, the maximum expected
          distortion achievable at rate 0 (i.e. using the single best
          reproduction symbol).

    Raises:
        ValueError: If distributions or matrix dimensions are inconsistent.
    """
    p_x = _normalize_distribution(np.asarray(source_probs), "source_probs")
    D = np.asarray(distortion_matrix, dtype=np.float64)

    if D.ndim != 2:
        raise ValueError(f"distortion_matrix must be 2-dimensional, got {D.ndim}D")
    n_source, n_repr = D.shape
    if n_source != p_x.size:
        raise ValueError(f"distortion_matrix has {n_source} rows but source_probs has " f"{p_x.size} symbols")
    if np.any(D < 0):
        raise ValueError("distortion_matrix contains negative entries")
    if n_points < 2:
        raise ValueError(f"n_points must be >= 2, got {n_points}")

    # D_max: distortion when we ignore the source entirely and pick the
    # single reproduction symbol that minimises expected distortion.
    d_max = float(np.min(p_x @ D))

    # Sweep beta from small (high rate / low distortion) to large
    # (low rate / high distortion).
    betas = np.logspace(-2, 3, n_points)

    rates: List[float] = []
    distortions: List[float] = []

    for beta in betas:
        # Blahut iteration for rate-distortion at fixed beta
        # q(j|i) proportional to q(j) exp(-beta d(i,j))
        q_j = np.full(n_repr, 1.0 / n_repr)

        for _ in range(200):
            log_q_cond = np.log(np.maximum(q_j, 1e-300))[np.newaxis, :] - beta * D  # (n_source, n_repr)
            log_q_cond -= log_q_cond.max(axis=1, keepdims=True)
            q_cond = np.exp(log_q_cond)
            q_cond /= q_cond.sum(axis=1, keepdims=True)

            q_j_new = p_x @ q_cond  # (n_repr,)
            q_j_new = np.maximum(q_j_new, 1e-300)
            q_j_new /= q_j_new.sum()

            if np.allclose(q_j, q_j_new, atol=1e-10):
                q_j = q_j_new
                break
            q_j = q_j_new

        # Expected distortion at convergence
        expected_distortion = float(np.sum(p_x[:, np.newaxis] * q_cond * D))

        # Rate: R = sum_i p(i) sum_j q(j|i) log[ q(j|i) / q(j) ]
        log_ratio = _safe_log(q_cond) - _safe_log(q_j)[np.newaxis, :]
        rate = float(np.sum(p_x[:, np.newaxis] * q_cond * log_ratio))
        rate = max(0.0, rate)

        rates.append(rate)
        distortions.append(expected_distortion)

    # Sort by distortion ascending (rate descending)
    order = np.argsort(distortions)
    rates = [rates[i] for i in order]
    distortions = [distortions[i] for i in order]

    logger.debug(
        "Rate-distortion curve computed with %d points; D_max=%.6f",
        n_points,
        d_max,
    )

    return {
        "rates": rates,
        "distortions": distortions,
        "critical_distortion": d_max,
    }


# ---------------------------------------------------------------------------
# 9. Information Bottleneck
# ---------------------------------------------------------------------------


def information_bottleneck(
    joint_xy: np.ndarray,
    beta: float = 1.0,
    n_clusters: int = 2,
    max_iter: int = 100,
) -> Dict[str, Any]:
    """Compute Tishby's information bottleneck for a joint distribution P(X, Y).

    The information bottleneck seeks a compressed representation T of X
    (with at most *n_clusters* values) that maximises the functional

        L = I(T; Y) - (1/beta) I(T; X)

    The algorithm alternates between updating the encoder ``q(t|x)``, the
    marginal ``q(t)``, and the decoder ``q(y|t)`` until convergence.

    Args:
        joint_xy: Joint distribution P(X, Y) as a 2-D array of shape
            ``(|X|, |Y|)`` whose entries sum to 1.
        beta: Trade-off parameter.  Larger beta favours higher relevance
            ``I(T; Y)`` at the cost of less compression.
        n_clusters: Number of clusters (cardinality of T).  Must be >= 2
            and at most ``|X|``.
        max_iter: Maximum number of iterations.

    Returns:
        Dictionary with keys:

        * ``compression_rate`` (*float*) -- I(T; X) in nats.
        * ``relevance`` (*float*) -- I(T; Y) in nats.
        * ``cluster_assignments`` (*np.ndarray*) -- Hard assignment
          ``argmax_t q(t|x)`` for each x, shape ``(|X|,)``.
        * ``n_iterations`` (*int*) -- Number of iterations performed.

    Raises:
        ValueError: If inputs are invalid.
    """
    joint = np.asarray(joint_xy, dtype=np.float64)
    if joint.ndim != 2:
        raise ValueError(f"joint_xy must be 2-dimensional, got {joint.ndim}D")
    n_x, n_y = joint.shape
    if n_x < 2 or n_y < 2:
        raise ValueError(f"joint_xy must be at least 2x2, got {n_x}x{n_y}")
    if np.any(joint < 0):
        raise ValueError("joint_xy contains negative entries")
    total = joint.sum()
    if total <= 0:
        raise ValueError("joint_xy sums to zero")
    joint = joint / total  # normalise

    if beta <= 0:
        raise ValueError(f"beta must be positive, got {beta}")
    if n_clusters < 2:
        raise ValueError(f"n_clusters must be >= 2, got {n_clusters}")
    if n_clusters > n_x:
        raise ValueError(f"n_clusters ({n_clusters}) cannot exceed |X| ({n_x})")
    if max_iter < 1:
        raise ValueError(f"max_iter must be >= 1, got {max_iter}")

    # Marginals
    p_x = joint.sum(axis=1)  # (n_x,)
    p_y = joint.sum(axis=0)  # (n_y,)

    # P(Y|X)
    p_y_given_x = np.zeros_like(joint)
    nonzero_x = p_x > 0
    p_y_given_x[nonzero_x] = joint[nonzero_x] / p_x[nonzero_x, np.newaxis]

    n_t = n_clusters

    # Initialise q(t|x) with a deterministic-seeded random soft assignment
    rng = np.random.RandomState(42)
    q_t_given_x = rng.dirichlet(np.ones(n_t), size=n_x)  # (n_x, n_t)

    iteration = 0
    for iteration in range(1, max_iter + 1):
        # q(t) = sum_x p(x) q(t|x)
        q_t = p_x @ q_t_given_x  # (n_t,)
        q_t = np.maximum(q_t, 1e-300)

        # q(y|t) = sum_x p(x) q(t|x) p(y|x) / q(t)
        q_y_given_t = (q_t_given_x * p_x[:, np.newaxis]).T @ p_y_given_x
        q_y_given_t /= q_t[:, np.newaxis]
        q_y_given_t = np.maximum(q_y_given_t, 1e-300)
        q_y_given_t /= q_y_given_t.sum(axis=1, keepdims=True)

        # KL divergence D_KL( p(y|x) || q(y|t) ) for each (x, t)
        log_p = _safe_log(p_y_given_x)  # (n_x, n_y)
        log_q = _safe_log(q_y_given_t)  # (n_t, n_y)

        kl = np.zeros((n_x, n_t))
        for t in range(n_t):
            diff = log_p - log_q[t][np.newaxis, :]  # (n_x, n_y)
            kl[:, t] = np.sum(p_y_given_x * diff, axis=1)

        # Update q(t|x) proportional to q(t) exp(-beta * kl)
        log_q_t_given_x = np.log(np.maximum(q_t, 1e-300))[np.newaxis, :] - beta * kl  # (n_x, n_t)
        log_q_t_given_x -= log_q_t_given_x.max(axis=1, keepdims=True)
        q_t_given_x_new = np.exp(log_q_t_given_x)
        q_t_given_x_new /= q_t_given_x_new.sum(axis=1, keepdims=True)

        # Check convergence
        change = float(np.max(np.abs(q_t_given_x_new - q_t_given_x)))
        q_t_given_x = q_t_given_x_new

        if change < 1e-8:
            logger.debug(
                "Information bottleneck converged after %d iterations",
                iteration,
            )
            break

    # Compute I(T; X) and I(T; Y)
    q_t = p_x @ q_t_given_x  # (n_t,)
    q_t = np.maximum(q_t, 1e-300)

    # I(T; X) = sum_{x,t} p(x) q(t|x) log[ q(t|x) / q(t) ]
    log_ratio_tx = _safe_log(q_t_given_x) - _safe_log(q_t)[np.newaxis, :]
    i_tx = float(np.sum(p_x[:, np.newaxis] * q_t_given_x * log_ratio_tx))
    i_tx = max(0.0, i_tx)

    # q(y|t) at convergence
    q_y_given_t = (q_t_given_x * p_x[:, np.newaxis]).T @ p_y_given_x
    q_y_given_t /= q_t[:, np.newaxis]
    q_y_given_t = np.maximum(q_y_given_t, 1e-300)
    q_y_given_t /= q_y_given_t.sum(axis=1, keepdims=True)

    # I(T; Y) = sum_{t,y} q(t) q(y|t) log[ q(y|t) / p(y) ]
    log_ratio_ty = _safe_log(q_y_given_t) - _safe_log(p_y)[np.newaxis, :]
    i_ty = float(np.sum(q_t[:, np.newaxis] * q_y_given_t * log_ratio_ty))
    i_ty = max(0.0, i_ty)

    # Hard cluster assignments
    assignments = np.argmax(q_t_given_x, axis=1).astype(np.intp)

    return {
        "compression_rate": i_tx,
        "relevance": i_ty,
        "cluster_assignments": assignments,
        "n_iterations": iteration,
    }


# ---------------------------------------------------------------------------
# 10. Entropy Power Inequality
# ---------------------------------------------------------------------------


def entropy_power_inequality(
    variances: Sequence[float],
) -> Dict[str, Any]:
    """Compute entropy powers and verify the entropy power inequality (EPI).

    For independent continuous random variables X_1, ..., X_k the entropy
    power is defined as

        N(X) = (1 / (2 pi e)) exp(2 h(X))

    where h(X) is the differential entropy.  For a Gaussian with
    variance sigma^2, h(X) = 0.5 log(2 pi e sigma^2) and therefore
    N(X) = sigma^2.

    The EPI states that for independent X, Y:

        N(X + Y) >= N(X) + N(Y)

    with equality iff both X and Y are Gaussian.  This function assumes
    Gaussian inputs so that the inequality holds with equality.

    Args:
        variances: Sequence of variances (one per independent Gaussian
            random variable).  Each must be positive.

    Returns:
        Dictionary with keys:

        * ``entropy_powers`` (*list[float]*) -- Entropy power N(X_i) for
          each variable (equals the variance under the Gaussian
          assumption).
        * ``sum_entropy_power`` (*float*) -- N(X_1 + ... + X_k) computed
          from the variance of the sum (sigma_1^2 + ... + sigma_k^2) under
          independence.
        * ``epi_bound`` (*float*) -- Sum of individual entropy powers
          N(X_1) + ... + N(X_k), the lower bound from the EPI.
        * ``epi_satisfied`` (*bool*) -- Whether the inequality
          N(sum) >= sum(N(X_i)) holds (should always be True for Gaussians
          within numerical tolerance).

    Raises:
        ValueError: If any variance is non-positive or the sequence is empty.
    """
    variances_arr = np.asarray(variances, dtype=np.float64).ravel()
    if variances_arr.size == 0:
        raise ValueError("variances must contain at least one element")
    if np.any(variances_arr <= 0):
        raise ValueError(
            "All variances must be positive; got values: " + str(variances_arr[variances_arr <= 0].tolist())
        )

    # For Gaussian X_i with variance sigma_i^2:
    #   h(X_i) = 0.5 * ln(2 * pi * e * sigma_i^2)
    #   N(X_i) = (1/(2*pi*e)) * exp(2 * h(X_i)) = sigma_i^2
    entropy_powers = variances_arr.tolist()

    # Under independence the sum Z = X_1 + ... + X_k is Gaussian with
    # variance = sum(sigma_i^2), so N(Z) = sum(sigma_i^2).
    sum_variance = float(variances_arr.sum())
    sum_entropy_power = sum_variance  # N(Z)

    # EPI bound: N(Z) >= sum N(X_i)
    epi_bound = float(variances_arr.sum())

    # The inequality is tight for Gaussians, so check within tolerance
    epi_satisfied = bool(sum_entropy_power >= epi_bound - 1e-12)

    logger.debug(
        "Entropy power inequality: N(sum)=%.6f, bound=%.6f, satisfied=%s",
        sum_entropy_power,
        epi_bound,
        epi_satisfied,
    )

    return {
        "entropy_powers": entropy_powers,
        "sum_entropy_power": sum_entropy_power,
        "epi_bound": epi_bound,
        "epi_satisfied": epi_satisfied,
    }


# ---------------------------------------------------------------------------
# 11. Information Dimension (Grassberger-Procaccia)
# ---------------------------------------------------------------------------


def information_dimension(
    samples: np.ndarray,
    r_values: Optional[np.ndarray] = None,
    method: str = "correlation",
) -> Dict[str, Any]:
    """Estimate the information (Renyi) dimension from a point-cloud dataset.

    For the *correlation dimension* (d_2) the Grassberger-Procaccia algorithm
    computes the correlation integral

        C(r) = (2 / (N(N-1))) sum_{i<j} Theta(r - ||x_i - x_j||)

    where Theta is the Heaviside step function.  The correlation dimension is

        d_2 = lim_{r -> 0} log C(r) / log r

    which is estimated by a linear fit to log C(r) vs log r over a scaling
    region.

    Args:
        samples: Point-cloud data of shape ``(n_samples, n_dims)`` or 1-D
            array which is treated as ``(n_samples, 1)``.
        r_values: Radii at which to evaluate the correlation integral.  If
            ``None``, a logarithmically spaced grid is chosen automatically
            based on the pairwise distance distribution.
        method: Dimension estimation method.  Currently only
            ``"correlation"`` (Grassberger-Procaccia) is supported.

    Returns:
        Dictionary with keys:

        * ``dimension`` (*float*) -- Estimated correlation dimension d_2.
        * ``r_values`` (*np.ndarray*) -- Radii used.
        * ``correlation_integral`` (*np.ndarray*) -- C(r) at each radius.
        * ``r_squared`` (*float*) -- R^2 of the linear fit in log-log space
          (quality of the scaling region).

    Raises:
        ValueError: If samples are too few or the method is unknown.
    """
    if method != "correlation":
        raise ValueError(f"Unknown method '{method}'; only 'correlation' is supported")

    samples = np.asarray(samples, dtype=np.float64)
    if samples.ndim == 1:
        samples = samples.reshape(-1, 1)
    if samples.ndim != 2:
        raise ValueError(f"samples must be 1-D or 2-D, got {samples.ndim}-D")

    n_samples, _n_dims = samples.shape
    if n_samples < 10:
        raise ValueError(f"Need at least 10 samples for dimension estimation, got {n_samples}")

    # Pairwise distances (condensed form)
    dists = pdist(samples)
    if dists.size == 0:
        raise ValueError("Pairwise distance vector is empty")

    # Remove exact zeros (duplicate points)
    dists_nonzero = dists[dists > 0]
    if dists_nonzero.size == 0:
        raise ValueError("All pairwise distances are zero -- cannot estimate dimension")

    # Choose r_values if not provided
    if r_values is None:
        r_min = float(np.percentile(dists_nonzero, 1))
        r_max = float(np.percentile(dists_nonzero, 90))
        if r_min <= 0:
            r_min = float(dists_nonzero.min())
        if r_min >= r_max:
            r_max = float(dists_nonzero.max())
            r_min = float(dists_nonzero.min())
        if r_min <= 0 or r_max <= 0 or r_min >= r_max:
            raise ValueError("Cannot determine a valid range of radii from the data")
        r_values = np.logspace(np.log10(r_min), np.log10(r_max), 30)
    else:
        r_values = np.asarray(r_values, dtype=np.float64).ravel()
        if r_values.size < 2:
            raise ValueError("r_values must contain at least 2 elements")
        if np.any(r_values <= 0):
            raise ValueError("All r_values must be positive")

    # Correlation integral: C(r) = fraction of pairs with distance < r
    n_pairs = n_samples * (n_samples - 1) / 2.0
    corr_integral = np.array(
        [float((dists < r).sum()) / n_pairs for r in r_values],
        dtype=np.float64,
    )

    # Keep only points where C(r) > 0 for log-log fit
    valid = corr_integral > 0
    if valid.sum() < 2:
        logger.warning("Fewer than 2 valid C(r) points; returning dimension=0.0")
        return {
            "dimension": 0.0,
            "r_values": r_values,
            "correlation_integral": corr_integral,
            "r_squared": 0.0,
        }

    log_r = np.log(r_values[valid])
    log_c = np.log(corr_integral[valid])

    # OLS fit: log C = d * log r + intercept
    coeffs = np.polyfit(log_r, log_c, 1)
    dimension = float(coeffs[0])

    # R-squared
    predicted = np.polyval(coeffs, log_r)
    ss_res = float(np.sum((log_c - predicted) ** 2))
    ss_tot = float(np.sum((log_c - log_c.mean()) ** 2))
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

    logger.debug(
        "Correlation dimension estimate: d_2=%.4f (R^2=%.4f, %d points)",
        dimension,
        r_squared,
        int(valid.sum()),
    )

    return {
        "dimension": dimension,
        "r_values": r_values,
        "correlation_integral": corr_integral,
        "r_squared": r_squared,
    }
