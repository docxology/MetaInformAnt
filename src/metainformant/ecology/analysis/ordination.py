"""Ordination methods for ecological community analysis.

This module implements multivariate ordination techniques used to visualize and
analyze patterns in ecological community data. All algorithms are implemented
from first principles using only numpy for matrix operations -- no scipy or
scikit-bio dependencies required.

Supported methods:
    - PCoA (Principal Coordinates Analysis) via Gower's double-centering
    - NMDS (Non-metric Multidimensional Scaling) via Kruskal's stress minimization
    - CCA (Canonical Correspondence Analysis) via weighted averaging + eigenanalysis
    - Distance matrix computation (Bray-Curtis, Jaccard, Euclidean, Manhattan, Canberra)
    - Procrustes rotation for ordination comparison
"""

from __future__ import annotations

import math
import random
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from metainformant.core.data import validation
from metainformant.core.utils import errors
from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


# ---------------------------------------------------------------------------
# Internal helpers: linear algebra primitives (no scipy)
# ---------------------------------------------------------------------------


def _validate_square_symmetric(matrix: List[List[float]], name: str = "matrix") -> int:
    """Validate that a matrix is square and approximately symmetric.

    Args:
        matrix: The matrix to validate.
        name: Human-readable name for error messages.

    Returns:
        The dimension *n* of the n x n matrix.

    Raises:
        ValueError: If the matrix is empty, not square, or not symmetric.
    """
    validation.validate_not_empty(matrix, name)
    n = len(matrix)
    for i, row in enumerate(matrix):
        if len(row) != n:
            raise ValueError(f"{name} must be square: row {i} has length {len(row)}, expected {n}")
    # Symmetry check (tolerance for floating-point)
    for i in range(n):
        for j in range(i + 1, n):
            if abs(matrix[i][j] - matrix[j][i]) > 1e-8:
                raise ValueError(
                    f"{name} is not symmetric: matrix[{i}][{j}]={matrix[i][j]} != matrix[{j}][{i}]={matrix[j][i]}"
                )
    return n


def _eigen_symmetric(A: np.ndarray, n_components: int) -> Tuple[np.ndarray, np.ndarray]:
    """Compute the top *n_components* eigenvalues/vectors of a real symmetric matrix.

    Uses numpy's ``np.linalg.eigh`` which is guaranteed to handle symmetric
    matrices correctly and returns eigenvalues in ascending order.  We reverse
    to get the largest eigenvalues first.

    Args:
        A: Real symmetric matrix of shape (n, n).
        n_components: Number of leading eigen-pairs to return.

    Returns:
        Tuple of (eigenvalues, eigenvectors) where eigenvalues has shape
        (n_components,) and eigenvectors has shape (n, n_components).  The
        columns of eigenvectors correspond to the eigenvalues.
    """
    eigenvalues, eigenvectors = np.linalg.eigh(A)
    # eigh returns ascending order -- reverse for descending
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    return eigenvalues[:n_components], eigenvectors[:, :n_components]


# ---------------------------------------------------------------------------
# 1. PCoA -- Principal Coordinates Analysis
# ---------------------------------------------------------------------------


def pcoa(
    distance_matrix: List[List[float]],
    n_components: int = 2,
) -> Dict[str, Any]:
    """Principal Coordinates Analysis (classical multidimensional scaling).

    Embeds objects described by a distance matrix into a low-dimensional
    Euclidean space that best preserves pairwise distances.  Uses Gower's
    double-centering to convert the distance matrix into an inner-product
    matrix, then performs eigendecomposition.

    Args:
        distance_matrix: Square symmetric distance matrix of shape (n, n).
            Values on the diagonal should be zero.
        n_components: Number of dimensions in the embedding (default 2).

    Returns:
        Dictionary with keys:
            - ``coordinates``: List[List[float]] of shape (n, n_components).
            - ``eigenvalues``: List[float] of the top eigenvalues.
            - ``variance_explained``: List[float] fraction of variance per axis.

    Raises:
        ValueError: If the distance matrix is empty, not square, not symmetric,
            or n_components exceeds the matrix dimension.

    Example:
        >>> dm = [[0, 1, 2], [1, 0, 1.5], [2, 1.5, 0]]
        >>> result = pcoa(dm, n_components=2)
        >>> len(result["coordinates"])
        3
        >>> len(result["coordinates"][0])
        2
    """
    n = _validate_square_symmetric(distance_matrix, "distance_matrix")

    if n == 0:
        return {"coordinates": [], "eigenvalues": [], "variance_explained": []}

    if n == 1:
        return {
            "coordinates": [[0.0] * n_components],
            "eigenvalues": [0.0] * n_components,
            "variance_explained": [0.0] * n_components,
        }

    if n_components > n:
        raise ValueError(f"n_components ({n_components}) cannot exceed number of objects ({n})")

    D = np.array(distance_matrix, dtype=np.float64)

    # Gower's double-centering: B = -0.5 * J D^2 J
    # where J = I - (1/n) * 11^T is the centering matrix
    D_sq = D**2
    row_means = D_sq.mean(axis=1, keepdims=True)
    col_means = D_sq.mean(axis=0, keepdims=True)
    grand_mean = D_sq.mean()
    B = -0.5 * (D_sq - row_means - col_means + grand_mean)

    # Eigendecomposition of B
    eigenvalues, eigenvectors = _eigen_symmetric(B, n_components)

    # Coordinates = eigenvectors * sqrt(max(eigenvalue, 0))
    # Negative eigenvalues are clamped to zero (may occur with non-Euclidean distances)
    coords = np.zeros((n, n_components), dtype=np.float64)
    for k in range(n_components):
        lam = max(eigenvalues[k], 0.0)
        coords[:, k] = eigenvectors[:, k] * math.sqrt(lam)

    # Variance explained
    positive_eigenvalues = eigenvalues[eigenvalues > 0]
    total_variance = float(positive_eigenvalues.sum()) if len(positive_eigenvalues) > 0 else 1.0
    variance_explained = [max(float(ev), 0.0) / total_variance for ev in eigenvalues]

    logger.info(
        "PCoA: %d objects -> %d components, variance explained: %s",
        n,
        n_components,
        [f"{v:.3f}" for v in variance_explained],
    )

    return {
        "coordinates": coords.tolist(),
        "eigenvalues": [float(ev) for ev in eigenvalues],
        "variance_explained": variance_explained,
    }


# ---------------------------------------------------------------------------
# 2. NMDS -- Non-metric Multidimensional Scaling
# ---------------------------------------------------------------------------


def _isotonic_regression(values: np.ndarray, weights: Optional[np.ndarray] = None) -> np.ndarray:
    """Pool adjacent violators algorithm for monotone (isotonic) regression.

    Fits a non-decreasing step function to *values* that minimizes the
    weighted sum of squared residuals.

    Args:
        values: 1-D array of target values.
        weights: Optional 1-D array of weights (default: uniform).

    Returns:
        1-D array of fitted monotone non-decreasing values.
    """
    n = len(values)
    if n == 0:
        return values.copy()
    if weights is None:
        weights = np.ones(n, dtype=np.float64)

    result = values.astype(np.float64).copy()
    w = weights.astype(np.float64).copy()

    # Pool adjacent violators
    i = 0
    while i < n - 1:
        if result[i] > result[i + 1]:
            # Merge blocks
            pooled = (result[i] * w[i] + result[i + 1] * w[i + 1]) / (w[i] + w[i + 1])
            result[i] = pooled
            result[i + 1] = pooled
            w[i] = w[i] + w[i + 1]
            w[i + 1] = w[i]
            # Walk back
            j = i
            while j > 0 and result[j - 1] > result[j]:
                pooled = (result[j - 1] * w[j - 1] + result[j] * w[j]) / (w[j - 1] + w[j])
                result[j - 1] = pooled
                result[j] = pooled
                w[j - 1] = w[j - 1] + w[j]
                w[j] = w[j - 1]
                j -= 1
        i += 1

    # Forward pass to propagate merged values
    for i in range(1, n):
        if result[i] < result[i - 1]:
            result[i] = result[i - 1]

    return result


def _kruskal_stress(
    d_orig: np.ndarray,
    d_config: np.ndarray,
    d_hat: np.ndarray,
) -> float:
    """Compute Kruskal's stress-1 statistic.

    stress = sqrt( sum((d_hat - d_config)^2) / sum(d_config^2) )

    Args:
        d_orig: Original distances (upper triangle, flattened).
        d_config: Distances in the current configuration (same shape).
        d_hat: Disparities from isotonic regression (same shape).

    Returns:
        Stress-1 value (0 = perfect, higher = worse).
    """
    numerator = float(np.sum((d_hat - d_config) ** 2))
    denominator = float(np.sum(d_config**2))
    if denominator < 1e-15:
        return 0.0
    return math.sqrt(numerator / denominator)


def _pairwise_distances_config(X: np.ndarray) -> np.ndarray:
    """Compute pairwise Euclidean distances for an (n, p) configuration.

    Returns a 1-D array corresponding to the upper triangle (i < j),
    in row-major order.

    Args:
        X: Configuration matrix of shape (n, p).

    Returns:
        1-D array of pairwise distances.
    """
    n = X.shape[0]
    dists = []
    for i in range(n):
        for j in range(i + 1, n):
            dists.append(float(np.sqrt(np.sum((X[i] - X[j]) ** 2))))
    return np.array(dists, dtype=np.float64)


def nmds(
    distance_matrix: List[List[float]],
    n_components: int = 2,
    max_iter: int = 300,
    n_init: int = 4,
    tol: float = 1e-7,
    seed: Optional[int] = None,
) -> Dict[str, Any]:
    """Non-metric Multidimensional Scaling using Kruskal's algorithm.

    Finds a low-dimensional configuration whose rank-order of pairwise
    distances best matches the rank-order of the original distance matrix.
    Multiple random initialisations are tried and the result with the lowest
    stress is returned.

    Args:
        distance_matrix: Square symmetric distance matrix of shape (n, n).
        n_components: Dimensionality of the embedding (default 2).
        max_iter: Maximum gradient descent iterations per initialisation.
        n_init: Number of random initialisations to try.
        tol: Convergence tolerance on stress improvement.
        seed: Optional random seed for reproducibility.

    Returns:
        Dictionary with keys:
            - ``coordinates``: List[List[float]] of shape (n, n_components).
            - ``stress``: Final Kruskal stress-1 value.
            - ``n_iter``: Number of iterations used in the best run.

    Raises:
        ValueError: If the distance matrix is invalid or n < 2.

    Example:
        >>> dm = [[0, 1, 2], [1, 0, 1.5], [2, 1.5, 0]]
        >>> result = nmds(dm, n_components=2, n_init=2, seed=42)
        >>> result["stress"] < 0.5
        True
    """
    n = _validate_square_symmetric(distance_matrix, "distance_matrix")

    if n < 2:
        coords = [[0.0] * n_components] * n
        return {"coordinates": coords, "stress": 0.0, "n_iter": 0}

    if n_components >= n:
        raise ValueError(f"n_components ({n_components}) must be less than n ({n})")

    rng = random.Random(seed)
    D = np.array(distance_matrix, dtype=np.float64)

    # Extract upper triangle distances
    d_orig_list: List[float] = []
    for i in range(n):
        for j in range(i + 1, n):
            d_orig_list.append(D[i, j])
    d_orig = np.array(d_orig_list, dtype=np.float64)

    # Sort indices for isotonic regression (by original distance rank)
    rank_order = np.argsort(d_orig)

    best_X: Optional[np.ndarray] = None
    best_stress = float("inf")
    best_n_iter = 0

    for init_idx in range(n_init):
        # Random initialisation (centered)
        X = np.array(
            [[rng.gauss(0, 1) for _ in range(n_components)] for _ in range(n)],
            dtype=np.float64,
        )
        X -= X.mean(axis=0)

        stress_prev = float("inf")
        step_size = 0.05

        for iteration in range(1, max_iter + 1):
            # Current configuration distances
            d_config = _pairwise_distances_config(X)

            # Isotonic regression on configuration distances in rank order of original distances
            d_sorted = d_config[rank_order]
            d_hat_sorted = _isotonic_regression(d_sorted)
            d_hat = np.empty_like(d_config)
            d_hat[rank_order] = d_hat_sorted

            # Compute stress
            stress = _kruskal_stress(d_orig, d_config, d_hat)

            # Check convergence
            if abs(stress_prev - stress) < tol:
                break
            stress_prev = stress

            # Gradient descent: move points to reduce (d_config - d_hat)
            # Derivative of stress w.r.t. X
            grad = np.zeros_like(X)
            idx = 0
            for i in range(n):
                for j in range(i + 1, n):
                    dc = d_config[idx]
                    if dc < 1e-12:
                        idx += 1
                        continue
                    ratio = (dc - d_hat[idx]) / dc
                    diff = X[i] - X[j]
                    grad[i] += ratio * diff
                    grad[j] -= ratio * diff
                    idx += 1

            # Normalise gradient and step
            grad_norm = float(np.sqrt(np.sum(grad**2)))
            if grad_norm > 1e-12:
                X -= step_size * (grad / grad_norm) * np.std(d_config)

            # Re-center
            X -= X.mean(axis=0)

        # Recompute final stress
        d_config = _pairwise_distances_config(X)
        d_sorted = d_config[rank_order]
        d_hat_sorted = _isotonic_regression(d_sorted)
        d_hat = np.empty_like(d_config)
        d_hat[rank_order] = d_hat_sorted
        stress = _kruskal_stress(d_orig, d_config, d_hat)

        if stress < best_stress:
            best_stress = stress
            best_X = X.copy()
            best_n_iter = iteration

    if best_X is None:
        best_X = np.zeros((n, n_components), dtype=np.float64)

    logger.info(
        "NMDS: %d objects -> %d components, stress=%.4f, iterations=%d",
        n,
        n_components,
        best_stress,
        best_n_iter,
    )

    return {
        "coordinates": best_X.tolist(),
        "stress": float(best_stress),
        "n_iter": best_n_iter,
    }


# ---------------------------------------------------------------------------
# 3. CCA -- Canonical Correspondence Analysis
# ---------------------------------------------------------------------------


def cca(
    species_matrix: List[List[float]],
    environmental_matrix: List[List[float]],
    n_components: Optional[int] = None,
) -> Dict[str, Any]:
    """Canonical Correspondence Analysis.

    CCA is a constrained ordination method that relates community composition
    to measured environmental variables.  It performs weighted averaging of the
    species data followed by eigenanalysis of the environmental-constrained
    chi-square distances.

    Args:
        species_matrix: Sites x species abundance matrix of shape (n, p).
            All values must be non-negative.
        environmental_matrix: Sites x environmental-variable matrix of shape
            (n, q).  Should be quantitative (centred internally).
        n_components: Number of ordination axes to return.  Defaults to
            min(q, p-1, n-1).

    Returns:
        Dictionary with keys:
            - ``site_scores``: List[List[float]] shape (n, n_components).
            - ``species_scores``: List[List[float]] shape (p, n_components).
            - ``eigenvalues``: List[float] of constrained eigenvalues.
            - ``variance_explained``: List[float] fraction of inertia per axis.

    Raises:
        ValueError: If dimensions are inconsistent, matrices are empty, or
            species abundances contain negative values.

    Example:
        >>> sp = [[10, 0, 5], [0, 8, 2], [3, 3, 3]]
        >>> env = [[1.0, 2.0], [3.0, 1.0], [2.0, 3.0]]
        >>> result = cca(sp, env)
        >>> len(result["site_scores"])
        3
    """
    validation.validate_not_empty(species_matrix, "species_matrix")
    validation.validate_not_empty(environmental_matrix, "environmental_matrix")

    n_sites = len(species_matrix)
    n_species = len(species_matrix[0])
    n_env = len(environmental_matrix[0])

    if len(environmental_matrix) != n_sites:
        raise ValueError(
            f"Row count mismatch: species_matrix has {n_sites} rows, "
            f"environmental_matrix has {len(environmental_matrix)} rows"
        )

    Y = np.array(species_matrix, dtype=np.float64)
    X = np.array(environmental_matrix, dtype=np.float64)

    if np.any(Y < 0):
        raise ValueError("species_matrix must contain only non-negative values")

    # Total abundance
    grand_total = Y.sum()
    if grand_total < 1e-15:
        raise ValueError("species_matrix total abundance is zero")

    # Relative frequencies
    P = Y / grand_total  # (n, p) proportional abundance matrix

    # Row and column weights
    r = P.sum(axis=1)  # site weights (n,)
    c = P.sum(axis=0)  # species weights (p,)

    # Avoid division by zero for sites/species with zero total
    r_safe = np.where(r > 0, r, 1.0)
    c_safe = np.where(c > 0, c, 1.0)

    # Chi-square standardised residuals: Q_{ij} = (P_{ij} - r_i*c_j) / sqrt(r_i * c_j)
    expected = np.outer(r, c)
    Q = (P - expected) / np.sqrt(np.where(expected > 0, expected, 1.0))

    # Weighted environmental matrix: centre and weight by sqrt(r)
    # Centre X with row weights
    X_mean = (r[:, np.newaxis] * X).sum(axis=0) / r.sum()
    X_centered = X - X_mean

    # Weight rows by sqrt(site weights)
    sqrt_r = np.sqrt(r_safe)
    X_w = X_centered * sqrt_r[:, np.newaxis]  # (n, q) weighted-centred env

    # Projection: constrained Q via environmental regression
    # Z = X_w (X_w^T X_w)^{-1} X_w^T Q  (projection of Q onto column space of X_w)
    XtX = X_w.T @ X_w
    try:
        XtX_inv = np.linalg.inv(XtX)
    except np.linalg.LinAlgError:
        # Singular -- use pseudoinverse
        XtX_inv = np.linalg.pinv(XtX)

    hat_matrix = X_w @ XtX_inv @ X_w.T  # (n, n)
    Q_hat = hat_matrix @ Q  # (n, p) constrained residuals

    # Eigenanalysis of Q_hat^T Q_hat (p x p) for species scores
    # or equivalently SVD of Q_hat
    M = Q_hat.T @ Q_hat  # (p, p)

    max_components = min(n_env, n_species - 1, n_sites - 1)
    max_components = max(max_components, 1)
    if n_components is None:
        n_components = max_components
    else:
        n_components = min(n_components, max_components)

    eigenvalues_all, V = np.linalg.eigh(M)
    idx = np.argsort(eigenvalues_all)[::-1]
    eigenvalues_all = eigenvalues_all[idx]
    V = V[:, idx]

    eigenvalues = eigenvalues_all[:n_components]
    V_k = V[:, :n_components]

    # Species scores: scale by 1/sqrt(c) and eigenvectors
    species_scores = np.zeros((n_species, n_components), dtype=np.float64)
    for k in range(n_components):
        lam = max(float(eigenvalues[k]), 0.0)
        if lam > 1e-15:
            species_scores[:, k] = V_k[:, k] / np.sqrt(c_safe)

    # Site scores: weighted average of species scores
    # u_i = (1/r_i) * sum_j P_{ij} * v_j  (weighted average)
    site_scores = np.zeros((n_sites, n_components), dtype=np.float64)
    for k in range(n_components):
        lam = max(float(eigenvalues[k]), 0.0)
        if lam > 1e-15:
            for i in range(n_sites):
                if r[i] > 0:
                    site_scores[i, k] = np.sum(P[i, :] * species_scores[:, k]) / r[i]

    # Variance explained (as fraction of total inertia)
    total_inertia = float(np.sum(eigenvalues_all[eigenvalues_all > 0]))
    if total_inertia < 1e-15:
        total_inertia = 1.0
    variance_explained = [max(float(ev), 0.0) / total_inertia for ev in eigenvalues]

    logger.info(
        "CCA: %d sites x %d species, %d env vars -> %d axes, variance explained: %s",
        n_sites,
        n_species,
        n_env,
        n_components,
        [f"{v:.3f}" for v in variance_explained],
    )

    return {
        "site_scores": site_scores.tolist(),
        "species_scores": species_scores.tolist(),
        "eigenvalues": [max(float(ev), 0.0) for ev in eigenvalues],
        "variance_explained": variance_explained,
    }


# ---------------------------------------------------------------------------
# 4. Distance matrix computation
# ---------------------------------------------------------------------------

_DISTANCE_METHODS = ("bray_curtis", "jaccard", "euclidean", "manhattan", "canberra")


def _bray_curtis(a: np.ndarray, b: np.ndarray) -> float:
    """Bray-Curtis dissimilarity between two abundance vectors."""
    numerator = float(np.sum(np.abs(a - b)))
    denominator = float(np.sum(a + b))
    if denominator < 1e-15:
        return 0.0
    return numerator / denominator


def _jaccard(a: np.ndarray, b: np.ndarray) -> float:
    """Jaccard distance (1 - Jaccard index) based on presence/absence."""
    present_a = a > 0
    present_b = b > 0
    intersection = int(np.sum(present_a & present_b))
    union = int(np.sum(present_a | present_b))
    if union == 0:
        return 0.0
    return 1.0 - intersection / union


def _euclidean(a: np.ndarray, b: np.ndarray) -> float:
    """Euclidean distance."""
    return float(np.sqrt(np.sum((a - b) ** 2)))


def _manhattan(a: np.ndarray, b: np.ndarray) -> float:
    """Manhattan (city-block) distance."""
    return float(np.sum(np.abs(a - b)))


def _canberra(a: np.ndarray, b: np.ndarray) -> float:
    """Canberra distance."""
    total = 0.0
    for ai, bi in zip(a, b):
        denom = abs(float(ai)) + abs(float(bi))
        if denom > 1e-15:
            total += abs(float(ai) - float(bi)) / denom
    return total


_DISTANCE_FUNCS = {
    "bray_curtis": _bray_curtis,
    "jaccard": _jaccard,
    "euclidean": _euclidean,
    "manhattan": _manhattan,
    "canberra": _canberra,
}


def distance_matrix(
    communities: List[List[float]],
    method: str = "bray_curtis",
) -> List[List[float]]:
    """Compute a pairwise distance matrix from community abundance vectors.

    Args:
        communities: List of abundance vectors, each of the same length.
            Shape is (n_sites, n_species).
        method: Distance metric.  One of ``"bray_curtis"``, ``"jaccard"``,
            ``"euclidean"``, ``"manhattan"``, ``"canberra"``.

    Returns:
        Symmetric distance matrix of shape (n, n) as nested lists.

    Raises:
        ValueError: If communities is empty, vectors have unequal length,
            or an unsupported method is given.

    Example:
        >>> comms = [[10, 0, 5], [0, 8, 2], [3, 3, 3]]
        >>> dm = distance_matrix(comms, method="bray_curtis")
        >>> dm[0][0]
        0.0
        >>> dm[0][1] == dm[1][0]
        True
    """
    validation.validate_not_empty(communities, "communities")

    if method not in _DISTANCE_FUNCS:
        raise ValueError(f"Unsupported distance method: {method!r}. Choose from {list(_DISTANCE_FUNCS.keys())}")

    n = len(communities)

    if n == 0:
        return []
    if n == 1:
        return [[0.0]]

    # Validate equal lengths
    expected_len = len(communities[0])
    for i, comm in enumerate(communities):
        if len(comm) != expected_len:
            raise ValueError(
                f"All community vectors must have the same length. "
                f"Community 0 has length {expected_len}, community {i} has length {len(comm)}"
            )

    func = _DISTANCE_FUNCS[method]
    arr = [np.array(c, dtype=np.float64) for c in communities]
    dm: List[List[float]] = [[0.0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i + 1, n):
            d = func(arr[i], arr[j])
            dm[i][j] = d
            dm[j][i] = d

    logger.info("Distance matrix: %d communities, method=%s", n, method)
    return dm


# ---------------------------------------------------------------------------
# 5. Procrustes rotation
# ---------------------------------------------------------------------------


def procrustes(
    coords1: List[List[float]],
    coords2: List[List[float]],
) -> Dict[str, Any]:
    """Procrustes rotation to optimally superimpose two ordination configurations.

    Translates, uniformly scales, and rotates *coords2* to best match
    *coords1* in the least-squares sense.  This is the orthogonal Procrustes
    analysis (no reflection allowed by default, but reflection is included
    when it reduces the residual).

    Args:
        coords1: Reference configuration of shape (n, p).
        coords2: Target configuration of shape (n, p) to be transformed.

    Returns:
        Dictionary with keys:
            - ``transformed_coords``: List[List[float]] -- coords2 after
              translation, scaling, and rotation to match coords1.
            - ``m2``: float -- sum of squared differences (Procrustes
              statistic) after transformation.
            - ``correlation``: float -- Procrustes correlation (1 - m2)
              normalised to [0, 1].

    Raises:
        ValueError: If configurations have different shapes, are empty, or
            have fewer than 2 points.

    Example:
        >>> c1 = [[0, 0], [1, 0], [0, 1]]
        >>> c2 = [[0, 0], [0, 1], [-1, 0]]
        >>> result = procrustes(c1, c2)
        >>> result["correlation"] > 0.9
        True
    """
    validation.validate_not_empty(coords1, "coords1")
    validation.validate_not_empty(coords2, "coords2")

    n1 = len(coords1)
    n2 = len(coords2)
    if n1 != n2:
        raise ValueError(f"Configurations must have the same number of points: {n1} != {n2}")

    p1 = len(coords1[0])
    p2 = len(coords2[0])
    if p1 != p2:
        raise ValueError(f"Configurations must have the same dimensionality: {p1} != {p2}")

    n = n1
    p = p1

    if n < 2:
        raise ValueError(f"Need at least 2 points for Procrustes analysis, got {n}")

    X = np.array(coords1, dtype=np.float64)
    Y = np.array(coords2, dtype=np.float64)

    # Step 1: Translate both to centroid at origin
    X_centroid = X.mean(axis=0)
    Y_centroid = Y.mean(axis=0)
    X_c = X - X_centroid
    Y_c = Y - Y_centroid

    # Step 2: Scale X to unit sum-of-squares (Frobenius norm)
    ss_X = float(np.sum(X_c**2))
    ss_Y = float(np.sum(Y_c**2))

    if ss_X < 1e-15 or ss_Y < 1e-15:
        # Degenerate: all points coincide
        return {
            "transformed_coords": X_c.tolist(),
            "m2": 0.0,
            "correlation": 1.0,
        }

    norm_X = math.sqrt(ss_X)
    norm_Y = math.sqrt(ss_Y)
    X_n = X_c / norm_X
    Y_n = Y_c / norm_Y

    # Step 3: Optimal rotation via SVD
    # Find R that minimises || X_n - Y_n @ R ||^2
    # Solution: R = V @ U^T where U S V^T = SVD(Y_n^T @ X_n)
    M = Y_n.T @ X_n  # (p, p)
    U, S, Vt = np.linalg.svd(M)
    V = Vt.T

    # Check for reflection: det(V @ U^T) should be +1 for proper rotation
    # If det < 0, flip sign of last column of V
    d = np.linalg.det(V @ U.T)
    if d < 0:
        V[:, -1] *= -1

    R = V @ U.T  # (p, p) rotation matrix

    # Step 4: Apply rotation and scale Y to match X
    Y_rotated = Y_n @ R  # (n, p)

    # Scale back to X's scale
    Y_transformed = Y_rotated * norm_X + X_centroid

    # Procrustes statistic m2 = sum of squared differences after optimal superimposition
    # Computed on the normalised coordinates
    m2 = float(np.sum((X_n - Y_rotated) ** 2))

    # Correlation: related to the trace of S
    # Procrustes correlation = sqrt(1 - m2) when properly normalised
    trace_S = float(np.sum(S))
    # m2 = 2 * (1 - trace_S) when both are normalised to unit SS
    correlation = max(0.0, min(1.0, trace_S))

    logger.info("Procrustes: %d points, %d dims, m2=%.6f, correlation=%.4f", n, p, m2, correlation)

    return {
        "transformed_coords": Y_transformed.tolist(),
        "m2": m2,
        "correlation": correlation,
    }
