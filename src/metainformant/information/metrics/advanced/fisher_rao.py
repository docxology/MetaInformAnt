"""Fisher-Rao distance, natural gradient, and related information geometry measures.

This module implements the Fisher-Rao geodesic distance on the statistical
manifold, the natural gradient via Fisher information, exponential family
entropy, Hellinger distance, and distribution validation utilities.
"""

from __future__ import annotations

import math
from typing import Any, Dict, Sequence

import numpy as np

from metainformant.core.data import validation
from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _validate_distribution(p: np.ndarray, name: str) -> np.ndarray:
    """Validate and return a probability distribution as a float64 array.

    Args:
        p: Array-like probability distribution.
        name: Name of the parameter (for error messages).

    Returns:
        Validated numpy array of dtype float64.

    Raises:
        ValueError: If the distribution contains negative values or does not
            sum to 1.
    """
    p = np.asarray(p, dtype=np.float64).ravel()
    if p.size == 0:
        raise ValueError(f"{name} must not be empty")
    if np.any(p < 0):
        raise ValueError(f"{name} contains negative values")
    if not np.isclose(np.sum(p), 1.0, atol=1e-6):
        raise ValueError(f"{name} must sum to 1.0, got {np.sum(p):.8f}")
    return p


# ---------------------------------------------------------------------------
# 1. Fisher-Rao Distance
# ---------------------------------------------------------------------------


def fisher_rao_distance(
    p: Sequence[float],
    q: Sequence[float],
) -> float:
    """Compute the Fisher-Rao geodesic distance between two discrete distributions.

    The Fisher-Rao distance is the geodesic distance on the statistical
    manifold equipped with the Fisher information metric.  For discrete
    distributions it equals

        d_FR(p, q) = 2 * arccos( sum_i sqrt(p_i * q_i) )

    which is twice the Bhattacharyya angle.  The result is non-negative and
    satisfies the triangle inequality, making it a true metric on the
    probability simplex.

    Args:
        p: First probability distribution.  Must sum to 1 and have the same
            length as *q*.
        q: Second probability distribution.  Must sum to 1 and have the same
            length as *p*.

    Returns:
        Fisher-Rao distance (non-negative float).  Returns 0.0 when the
        distributions are identical and pi when they have disjoint support.

    Raises:
        ValueError: If distributions have different lengths, contain negative
            values, or do not sum to 1.
    """
    validation.validate_type(p, (list, tuple, np.ndarray), "p")
    validation.validate_type(q, (list, tuple, np.ndarray), "q")

    p_arr = _validate_distribution(np.asarray(p), "p")
    q_arr = _validate_distribution(np.asarray(q), "q")

    if p_arr.size != q_arr.size:
        raise ValueError(f"Distributions must have the same length: " f"len(p)={p_arr.size}, len(q)={q_arr.size}")

    # Bhattacharyya coefficient BC(p, q) = sum sqrt(p_i * q_i)
    bc = float(np.sum(np.sqrt(p_arr * q_arr)))

    # Clamp to [0, 1] for numerical safety (sqrt products can slightly exceed 1)
    bc = max(0.0, min(1.0, bc))

    distance = 2.0 * math.acos(bc)

    logger.debug(
        "Fisher-Rao distance: %.6f (Bhattacharyya coefficient=%.6f)",
        distance,
        bc,
    )

    return distance


# ---------------------------------------------------------------------------
# 2. Natural Gradient
# ---------------------------------------------------------------------------


def natural_gradient(
    loss_gradient: np.ndarray,
    fisher_info_matrix: np.ndarray,
) -> np.ndarray:
    """Compute the natural gradient by pre-multiplying with the inverse Fisher information.

    The natural gradient rescales the ordinary (Euclidean) gradient by the
    inverse of the Fisher information matrix so that parameter updates are
    invariant to reparameterisation of the statistical model:

        natural_grad = F^{-1} @ grad

    For numerical stability the system ``F @ natural_grad = grad`` is solved
    directly via ``numpy.linalg.solve`` rather than explicitly inverting F.

    Args:
        loss_gradient: Ordinary gradient vector of shape ``(d,)``.
        fisher_info_matrix: Fisher information matrix of shape ``(d, d)``.
            Must be symmetric and positive semi-definite.

    Returns:
        Natural gradient vector of shape ``(d,)`` as a numpy array.

    Raises:
        ValueError: If dimensions are inconsistent or the Fisher information
            matrix is not square.
        numpy.linalg.LinAlgError: If the Fisher information matrix is
            singular and the system cannot be solved.
    """
    validation.validate_type(loss_gradient, (list, tuple, np.ndarray), "loss_gradient")
    validation.validate_type(fisher_info_matrix, (list, tuple, np.ndarray), "fisher_info_matrix")

    grad = np.asarray(loss_gradient, dtype=np.float64).ravel()
    fim = np.asarray(fisher_info_matrix, dtype=np.float64)

    if fim.ndim != 2:
        raise ValueError(f"fisher_info_matrix must be 2-dimensional, got {fim.ndim}D")
    if fim.shape[0] != fim.shape[1]:
        raise ValueError(f"fisher_info_matrix must be square, got shape {fim.shape}")
    if grad.size != fim.shape[0]:
        raise ValueError(
            f"Dimension mismatch: loss_gradient has {grad.size} elements but "
            f"fisher_info_matrix is {fim.shape[0]}x{fim.shape[1]}"
        )

    # Solve F @ nat_grad = grad  (more stable than inv(F) @ grad)
    try:
        nat_grad = np.linalg.solve(fim, grad)
    except np.linalg.LinAlgError:
        # Fall back to pseudo-inverse for singular / near-singular FIM
        logger.warning("Fisher information matrix is singular; using pseudo-inverse")
        nat_grad = np.linalg.lstsq(fim, grad, rcond=None)[0]

    logger.debug(
        "Natural gradient computed: ||grad||=%.6f, ||nat_grad||=%.6f",
        float(np.linalg.norm(grad)),
        float(np.linalg.norm(nat_grad)),
    )

    return nat_grad


# ---------------------------------------------------------------------------
# 5. Exponential Family Entropy
# ---------------------------------------------------------------------------


def exponential_family_entropy(
    natural_params: Sequence[float],
    sufficient_stats_expectation: Sequence[float],
    log_partition: float,
) -> float:
    """Compute the entropy of an exponential family distribution.

    For an exponential family distribution with density

        p(x; eta) = h(x) exp( eta^T T(x) - A(eta) )

    the entropy is given by the Legendre transform relation:

        H[p] = -eta^T E[T(x)] + A(eta)

    where eta is the natural parameter vector, E[T(x)] is the expectation of
    the sufficient statistics, and A(eta) is the log-partition function.

    Args:
        natural_params: Natural parameter vector eta of shape ``(d,)``.
        sufficient_stats_expectation: Expected sufficient statistics E[T(x)]
            of shape ``(d,)``.  Must have the same length as *natural_params*.
        log_partition: Value of the log-partition function A(eta) (scalar).

    Returns:
        Entropy value (in nats).

    Raises:
        ValueError: If the parameter vectors have different lengths.
    """
    validation.validate_type(natural_params, (list, tuple, np.ndarray), "natural_params")
    validation.validate_type(
        sufficient_stats_expectation,
        (list, tuple, np.ndarray),
        "sufficient_stats_expectation",
    )

    eta = np.asarray(natural_params, dtype=np.float64).ravel()
    e_t = np.asarray(sufficient_stats_expectation, dtype=np.float64).ravel()

    if eta.size != e_t.size:
        raise ValueError(
            f"natural_params and sufficient_stats_expectation must have the " f"same length: {eta.size} != {e_t.size}"
        )

    if eta.size == 0:
        raise ValueError("natural_params must not be empty")

    log_partition_val = float(log_partition)

    # H[p] = -eta^T E[T(x)] + A(eta)
    entropy = -float(np.dot(eta, e_t)) + log_partition_val

    logger.debug(
        "Exponential family entropy: %.6f nats " "(||eta||=%.4f, A(eta)=%.4f)",
        entropy,
        float(np.linalg.norm(eta)),
        log_partition_val,
    )

    return entropy


# ---------------------------------------------------------------------------
# 6. Hellinger Distance
# ---------------------------------------------------------------------------


def hellinger_distance(
    p: Sequence[float],
    q: Sequence[float],
) -> float:
    """Compute the Hellinger distance between two discrete distributions.

    The Hellinger distance is defined as

        H(p, q) = (1 / sqrt(2)) * sqrt( sum_i (sqrt(p_i) - sqrt(q_i))^2 )

    and takes values in [0, 1].  It equals 0 when the distributions are
    identical and 1 when they have disjoint support.  The Hellinger distance
    is related to the Bhattacharyya coefficient BC(p,q) via

        H(p, q) = sqrt(1 - BC(p, q))

    and to the Fisher-Rao distance (the geodesic distance on the probability
    simplex).

    Args:
        p: First probability distribution.  Must sum to 1 and have the same
            length as *q*.
        q: Second probability distribution.  Must sum to 1 and have the same
            length as *p*.

    Returns:
        Hellinger distance in the range [0, 1].

    Raises:
        ValueError: If distributions have different lengths, contain negative
            values, or do not sum to 1.
    """
    validation.validate_type(p, (list, tuple, np.ndarray), "p")
    validation.validate_type(q, (list, tuple, np.ndarray), "q")

    p_arr = _validate_distribution(np.asarray(p), "p")
    q_arr = _validate_distribution(np.asarray(q), "q")

    if p_arr.size != q_arr.size:
        raise ValueError(f"Distributions must have the same length: " f"len(p)={p_arr.size}, len(q)={q_arr.size}")

    # H(p, q) = (1/sqrt(2)) * || sqrt(p) - sqrt(q) ||_2
    diff = np.sqrt(p_arr) - np.sqrt(q_arr)
    distance = float(np.sqrt(np.sum(diff**2)) / math.sqrt(2.0))

    # Clamp to [0, 1] for numerical safety
    distance = max(0.0, min(1.0, distance))

    logger.debug("Hellinger distance: %.6f", distance)

    return distance
