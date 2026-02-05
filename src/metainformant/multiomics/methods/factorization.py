"""Matrix factorization and network fusion methods for multi-omic integration.

Implements core integration algorithms including joint non-negative matrix
factorization, simplified MOFA (Multi-Omics Factor Analysis), CP tensor
decomposition, Similarity Network Fusion, and regularized Canonical
Correlation Analysis. All methods operate on plain Python structures with
optional numpy acceleration.

Algorithms:
    - Joint NMF: multiplicative update rules across concatenated omic layers.
    - MOFA simple: EM-based Bayesian factor model with ARD sparsity priors.
    - Tensor decomposition: alternating least squares CP decomposition.
    - SNF: iterative diffusion to fuse patient similarity networks.
    - CCA: SVD-based regularized canonical correlation for two views.
"""

from __future__ import annotations

import math
import random
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Optional scientific dependencies
try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None  # type: ignore[assignment]

try:
    from scipy import linalg as sp_linalg

    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    sp_linalg = None  # type: ignore[assignment]


# ---------------------------------------------------------------------------
# Pure-Python helpers (used when numpy is unavailable)
# ---------------------------------------------------------------------------


def _py_zeros(rows: int, cols: int) -> list[list[float]]:
    """Create a rows x cols matrix of zeros in pure Python."""
    return [[0.0] * cols for _ in range(rows)]


def _py_random_matrix(rows: int, cols: int) -> list[list[float]]:
    """Create a rows x cols matrix of random non-negative values."""
    return [[random.random() * 0.5 + 0.01 for _ in range(cols)] for _ in range(rows)]


def _py_matmul(A: list[list[float]], B: list[list[float]]) -> list[list[float]]:
    """Pure-Python matrix multiplication."""
    rows_a = len(A)
    cols_a = len(A[0])
    cols_b = len(B[0])
    result = _py_zeros(rows_a, cols_b)
    for i in range(rows_a):
        for k in range(cols_a):
            a_ik = A[i][k]
            for j in range(cols_b):
                result[i][j] += a_ik * B[k][j]
    return result


def _py_transpose(M: list[list[float]]) -> list[list[float]]:
    """Pure-Python matrix transpose."""
    rows = len(M)
    cols = len(M[0])
    return [[M[i][j] for i in range(rows)] for j in range(cols)]


def _py_frobenius_norm(M: list[list[float]]) -> float:
    """Frobenius norm of a matrix."""
    total = 0.0
    for row in M:
        for val in row:
            total += val * val
    return math.sqrt(total)


def _py_matrix_subtract(A: list[list[float]], B: list[list[float]]) -> list[list[float]]:
    """Element-wise A - B."""
    rows = len(A)
    cols = len(A[0])
    return [[A[i][j] - B[i][j] for j in range(cols)] for i in range(rows)]


# ---------------------------------------------------------------------------
# Joint NMF
# ---------------------------------------------------------------------------


def joint_nmf(
    data_matrices: dict[str, Any],
    k: int = 10,
    max_iter: int = 200,
    tol: float = 1e-4,
) -> dict[str, Any]:
    """Joint Non-negative Matrix Factorization across multiple omic layers.

    Factorises a set of non-negative data matrices that share the same
    sample dimension into a shared basis matrix W and per-omic coefficient
    matrices H_i.  Uses multiplicative update rules that guarantee
    non-negativity and monotonic decrease of the reconstruction error.

    Each omic matrix V_i (samples x features_i) is approximated as W @ H_i
    where W (samples x k) is shared.

    Args:
        data_matrices: Mapping from omic name to data matrix.  Each value is
            either a numpy ndarray or a list-of-lists with shape
            (n_samples, n_features_i).  All matrices must share the same
            number of rows (samples).
        k: Number of latent factors.
        max_iter: Maximum number of multiplicative-update iterations.
        tol: Convergence tolerance on the relative change of
            reconstruction error between successive iterations.

    Returns:
        Dictionary with keys:
            W: Shared factor matrix (n_samples x k).
            H_dict: Per-omic coefficient matrices {name: H_i}.
            reconstruction_error: Final Frobenius reconstruction error.
            n_iter: Number of iterations performed.
            converged: Whether the algorithm converged within *tol*.

    Raises:
        ValueError: If matrices have incompatible sample dimensions, if k < 1,
            or if any matrix contains negative values.
    """
    if k < 1:
        raise ValueError(f"k must be >= 1, got {k}")

    if not data_matrices:
        raise ValueError("data_matrices must contain at least one entry")

    logger.info("Running joint NMF with k=%d, max_iter=%d", k, max_iter)

    # --- Convert to numpy or keep as list-of-lists -----------------------
    use_np = HAS_NUMPY
    views: dict[str, Any] = {}
    n_samples: int | None = None

    for name, mat in data_matrices.items():
        if use_np:
            arr = np.asarray(mat, dtype=np.float64)
            if arr.ndim != 2:
                raise ValueError(f"Matrix '{name}' must be 2-D, got {arr.ndim}-D")
            if np.any(arr < 0):
                raise ValueError(f"Matrix '{name}' contains negative values; NMF requires non-negative data")
            if n_samples is None:
                n_samples = arr.shape[0]
            elif arr.shape[0] != n_samples:
                raise ValueError(
                    f"Sample-dimension mismatch: expected {n_samples}, " f"got {arr.shape[0]} for '{name}'"
                )
            views[name] = arr
        else:
            # list-of-lists path
            if not isinstance(mat, list) or not mat:
                raise ValueError(f"Matrix '{name}' must be a non-empty list-of-lists")
            rows = len(mat)
            if n_samples is None:
                n_samples = rows
            elif rows != n_samples:
                raise ValueError(f"Sample-dimension mismatch: expected {n_samples}, " f"got {rows} for '{name}'")
            for row in mat:
                for val in row:
                    if val < 0:
                        raise ValueError(f"Matrix '{name}' contains negative values; " "NMF requires non-negative data")
            views[name] = mat

    assert n_samples is not None  # at least one matrix

    # --- Initialise W and H_i -------------------------------------------
    eps = 1e-10  # small constant to avoid division by zero

    if use_np:
        rng = np.random.RandomState(42)
        W = rng.rand(n_samples, k).astype(np.float64) * 0.5 + 0.01
        H_dict: dict[str, Any] = {}
        for name, V in views.items():
            n_feat = V.shape[1]
            H_dict[name] = rng.rand(k, n_feat).astype(np.float64) * 0.5 + 0.01
    else:
        random.seed(42)
        W = _py_random_matrix(n_samples, k)
        H_dict = {}
        for name, V in views.items():
            n_feat = len(V[0])
            H_dict[name] = _py_random_matrix(k, n_feat)

    # --- Multiplicative updates ------------------------------------------
    prev_error = float("inf")
    converged = False
    n_iter = 0

    for iteration in range(1, max_iter + 1):
        n_iter = iteration

        if use_np:
            # --- Update each H_i ---
            for name in views:
                V = views[name]
                H = H_dict[name]
                numerator = W.T @ V  # (k, f_i)
                denominator = (W.T @ W) @ H + eps  # (k, f_i)
                H_dict[name] = H * (numerator / denominator)

            # --- Update W (aggregate gradients across all views) ---
            W_num = np.zeros_like(W)
            W_den = np.zeros_like(W)
            for name in views:
                V = views[name]
                H = H_dict[name]
                W_num += V @ H.T  # (n, k)
                W_den += (W @ H) @ H.T + eps  # (n, k)
            W = W * (W_num / W_den)

            # --- Compute reconstruction error ---
            total_error = 0.0
            for name in views:
                diff = views[name] - W @ H_dict[name]
                total_error += float(np.sum(diff**2))
            total_error = math.sqrt(total_error)
        else:
            # Pure Python path
            for name in views:
                V = views[name]
                H = H_dict[name]
                Wt = _py_transpose(W)
                numerator = _py_matmul(Wt, V)
                WtW = _py_matmul(Wt, W)
                denominator_mat = _py_matmul(WtW, H)
                k_dim = len(H)
                f_dim = len(H[0])
                new_H = [
                    [H[r][c] * (numerator[r][c] / (denominator_mat[r][c] + eps)) for c in range(f_dim)]
                    for r in range(k_dim)
                ]
                H_dict[name] = new_H

            W_num = _py_zeros(n_samples, k)
            W_den = _py_zeros(n_samples, k)
            for name in views:
                V = views[name]
                H = H_dict[name]
                Ht = _py_transpose(H)
                vht = _py_matmul(V, Ht)
                wh = _py_matmul(W, H)
                whht = _py_matmul(wh, Ht)
                for i in range(n_samples):
                    for j in range(k):
                        W_num[i][j] += vht[i][j]
                        W_den[i][j] += whht[i][j] + eps
            W = [[W[i][j] * (W_num[i][j] / W_den[i][j]) for j in range(k)] for i in range(n_samples)]

            total_error = 0.0
            for name in views:
                V = views[name]
                WH = _py_matmul(W, H_dict[name])
                diff = _py_matrix_subtract(V, WH)
                total_error += sum(v * v for row in diff for v in row)
            total_error = math.sqrt(total_error)

        # Convergence check
        rel_change = abs(prev_error - total_error) / (prev_error + eps)
        if iteration > 1 and rel_change < tol:
            converged = True
            logger.info("Joint NMF converged at iteration %d (rel_change=%.2e)", iteration, rel_change)
            break
        prev_error = total_error

        if iteration % 50 == 0:
            logger.debug("Joint NMF iteration %d, error=%.6f", iteration, total_error)

    if not converged:
        logger.warning("Joint NMF did not converge within %d iterations", max_iter)

    # Convert results to lists if numpy was used (for JSON serialisability)
    result_W = W.tolist() if use_np else W
    result_H: dict[str, Any] = {}
    for name, H in H_dict.items():
        result_H[name] = H.tolist() if use_np else H

    return {
        "W": result_W,
        "H_dict": result_H,
        "reconstruction_error": total_error,
        "n_iter": n_iter,
        "converged": converged,
    }


# ---------------------------------------------------------------------------
# Simplified MOFA (Multi-Omics Factor Analysis)
# ---------------------------------------------------------------------------


def mofa_simple(
    data_matrices: dict[str, Any],
    k: int = 10,
    max_iter: int = 100,
    tol: float = 1e-3,
) -> dict[str, Any]:
    """Simplified MOFA+ (Multi-Omics Factor Analysis).

    Implements a Bayesian factor model with Automatic Relevance Determination
    (ARD) priors that encourage sparsity across views.  Uses an EM algorithm
    where the E-step computes posterior factor expectations and the M-step
    updates weights and precision parameters.

    Each omic matrix Y_m (samples x features_m) is modelled as
        Y_m ~ Z @ W_m^T + noise
    where Z (samples x k) are shared latent factors and W_m (features_m x k)
    are per-view weights.  ARD priors on columns of W_m allow entire factors
    to be switched off for irrelevant views.

    Args:
        data_matrices: Mapping from omic name to data matrix (n_samples x
            n_features_m).  numpy arrays or list-of-lists.
        k: Number of latent factors.
        max_iter: Maximum EM iterations.
        tol: Convergence tolerance on the change in ELBO (evidence lower
            bound) between iterations.

    Returns:
        Dictionary with keys:
            factors: Latent factor matrix Z (n_samples x k).
            weights_per_view: Per-omic weight matrices {name: W_m}.
            variance_explained: Fraction of variance explained per factor
                per view {name: list[float]}.
            active_factors: List of factor indices with non-negligible
                contribution (ARD precision below threshold).

    Raises:
        ValueError: If matrices have incompatible sample dimensions or k < 1.
    """
    if k < 1:
        raise ValueError(f"k must be >= 1, got {k}")
    if not data_matrices:
        raise ValueError("data_matrices must contain at least one entry")

    logger.info("Running simplified MOFA with k=%d, max_iter=%d", k, max_iter)

    if not HAS_NUMPY:
        raise ImportError("numpy is required for MOFA. Install with: uv pip install numpy")

    # --- Prepare data -----------------------------------------------------
    views: dict[str, Any] = {}
    n_samples: int | None = None

    for name, mat in data_matrices.items():
        arr = np.asarray(mat, dtype=np.float64)
        if arr.ndim != 2:
            raise ValueError(f"Matrix '{name}' must be 2-D")
        if n_samples is None:
            n_samples = arr.shape[0]
        elif arr.shape[0] != n_samples:
            raise ValueError(f"Sample-dimension mismatch: expected {n_samples}, " f"got {arr.shape[0]} for '{name}'")
        # Center each feature
        arr = arr - arr.mean(axis=0, keepdims=True)
        views[name] = arr

    assert n_samples is not None

    rng = np.random.RandomState(42)

    # --- Initialise parameters -------------------------------------------
    Z = rng.randn(n_samples, k).astype(np.float64)  # Factors
    weights: dict[str, Any] = {}
    tau: dict[str, Any] = {}  # Noise precision per view
    alpha: dict[str, Any] = {}  # ARD precision per factor per view

    for name, Y in views.items():
        d_m = Y.shape[1]
        weights[name] = rng.randn(d_m, k).astype(np.float64) * 0.1
        tau[name] = np.ones(d_m, dtype=np.float64)
        alpha[name] = np.ones(k, dtype=np.float64)

    # --- EM iterations ----------------------------------------------------
    prev_elbo = -float("inf")
    eps = 1e-10

    for iteration in range(1, max_iter + 1):
        # ---- E-step: update posterior of Z ----
        # Posterior precision: Lambda = I + sum_m W_m^T diag(tau_m) W_m
        Lambda = np.eye(k, dtype=np.float64)
        rhs = np.zeros((n_samples, k), dtype=np.float64)

        for name, Y in views.items():
            W_m = weights[name]  # (d_m, k)
            tau_m = tau[name]  # (d_m,)
            # W_m^T diag(tau_m) W_m
            WtTW = (W_m * tau_m[:, None]).T @ W_m  # (k, k)
            Lambda += WtTW
            # Y_m diag(tau_m) W_m
            rhs += (Y * tau_m[None, :]) @ W_m  # (n, k)

        # Posterior mean of Z
        Lambda_inv = np.linalg.inv(Lambda)
        Z = rhs @ Lambda_inv.T  # (n, k)

        # E[Z^T Z] = n * Lambda_inv + Z^T Z
        ZtZ = Z.T @ Z + n_samples * Lambda_inv  # (k, k)

        # ---- M-step: update W_m, tau_m, alpha_m ----
        for name, Y in views.items():
            W_m = weights[name]
            d_m = Y.shape[1]
            alpha_m = alpha[name]  # (k,)

            # Update W_m: each row w_j ~ N(mu_j, Sigma_j)
            # Sigma_j = (tau_j * ZtZ + diag(alpha_m))^{-1}
            # mu_j = tau_j * Sigma_j * Z^T y_j
            for j in range(d_m):
                Sigma_j_inv = tau[name][j] * ZtZ + np.diag(alpha_m)
                Sigma_j = np.linalg.inv(Sigma_j_inv)
                mu_j = tau[name][j] * Sigma_j @ (Z.T @ Y[:, j])
                W_m[j, :] = mu_j

            weights[name] = W_m

            # Update tau_m (noise precision per feature)
            residual = Y - Z @ W_m.T  # (n, d_m)
            for j in range(d_m):
                ss = float(np.sum(residual[:, j] ** 2))
                # Add trace correction
                w_j = W_m[j, :]
                trace_corr = float(w_j @ (n_samples * Lambda_inv) @ w_j)
                tau[name][j] = n_samples / (ss + trace_corr + eps)

            # Update alpha_m (ARD precision per factor)
            for f in range(k):
                col_sq = float(np.sum(W_m[:, f] ** 2))
                alpha[name][f] = d_m / (col_sq + eps)

        # ---- Compute approximate ELBO for convergence --------------------
        elbo = 0.0
        for name, Y in views.items():
            residual = Y - Z @ weights[name].T
            d_m = Y.shape[1]
            for j in range(d_m):
                elbo += 0.5 * n_samples * math.log(tau[name][j] + eps)
                elbo -= 0.5 * tau[name][j] * float(np.sum(residual[:, j] ** 2))
        # Regularisation terms (simplified)
        elbo -= 0.5 * float(np.sum(Z**2))

        if iteration > 1 and abs(elbo - prev_elbo) / (abs(prev_elbo) + eps) < tol:
            logger.info("MOFA converged at iteration %d", iteration)
            break
        prev_elbo = elbo

        if iteration % 20 == 0:
            logger.debug("MOFA iteration %d, ELBO=%.4f", iteration, elbo)

    # --- Compute variance explained per factor per view -------------------
    variance_explained: dict[str, list[float]] = {}
    for name, Y in views.items():
        total_var = float(np.sum(Y**2))
        ve_per_factor: list[float] = []
        for f in range(k):
            reconstruction_f = np.outer(Z[:, f], weights[name][:, f])
            ve = float(np.sum(reconstruction_f**2)) / (total_var + eps)
            ve_per_factor.append(ve)
        variance_explained[name] = ve_per_factor

    # --- Determine active factors (ARD threshold) -------------------------
    # A factor is active if at least one view has low ARD precision
    ard_threshold = 1e4
    active_factors: list[int] = []
    for f in range(k):
        for name in views:
            if alpha[name][f] < ard_threshold:
                active_factors.append(f)
                break

    logger.info(
        "MOFA complete: %d/%d active factors across %d views",
        len(active_factors),
        k,
        len(views),
    )

    return {
        "factors": Z.tolist(),
        "weights_per_view": {name: w.tolist() for name, w in weights.items()},
        "variance_explained": variance_explained,
        "active_factors": active_factors,
    }


# ---------------------------------------------------------------------------
# Tensor Decomposition (CP / CANDECOMP-PARAFAC)
# ---------------------------------------------------------------------------


def tensor_decomposition(
    tensor: list[list[list[float]]] | Any,
    rank: int = 5,
    method: str = "cp",
    max_iter: int = 100,
) -> dict[str, Any]:
    """CP (CANDECOMP/PARAFAC) tensor decomposition via alternating least squares.

    Decomposes a three-way tensor X (I x J x K) into a sum of *rank*
    rank-one components:  X ~ sum_r a_r outer b_r outer c_r.

    The ALS algorithm alternately fixes two factor matrices and solves for
    the third in a least-squares sense.

    Args:
        tensor: A 3-D tensor as a nested list or numpy array with shape
            (I, J, K).
        rank: Number of components (rank of the decomposition).
        method: Decomposition method.  Currently only ``"cp"`` is supported.
        max_iter: Maximum ALS iterations.

    Returns:
        Dictionary with keys:
            factors: List of three factor matrices [A, B, C] where A is
                (I x rank), B is (J x rank), C is (K x rank).
            fit: Proportion of variance explained (1 - ||X - X_hat||^2 / ||X||^2).
            core_consistency: Core consistency diagnostic (percentage, 0-100).

    Raises:
        ValueError: If tensor is not 3-D, rank < 1, or method is unsupported.
    """
    if method != "cp":
        raise ValueError(f"Unsupported decomposition method: '{method}'. Use 'cp'.")
    if rank < 1:
        raise ValueError(f"rank must be >= 1, got {rank}")

    logger.info("Running CP tensor decomposition with rank=%d, max_iter=%d", rank, max_iter)

    if not HAS_NUMPY:
        raise ImportError("numpy is required for tensor decomposition. " "Install with: uv pip install numpy")

    T = np.asarray(tensor, dtype=np.float64)
    if T.ndim != 3:
        raise ValueError(f"Tensor must be 3-D, got {T.ndim}-D")

    I, J, K = T.shape
    rng = np.random.RandomState(42)

    # Initialise factor matrices
    A = rng.randn(I, rank).astype(np.float64)
    B = rng.randn(J, rank).astype(np.float64)
    C = rng.randn(K, rank).astype(np.float64)

    # Unfold the tensor along each mode
    # Mode-0 unfolding: (I, J*K)
    T0 = T.reshape(I, J * K)
    # Mode-1 unfolding: (J, I*K)
    T1 = T.transpose(1, 0, 2).reshape(J, I * K)
    # Mode-2 unfolding: (K, I*J)
    T2 = T.transpose(2, 0, 1).reshape(K, I * J)

    T_norm_sq = float(np.sum(T**2))
    eps = 1e-10

    for iteration in range(1, max_iter + 1):
        # Update A: T_(0) ~ A (C khatri-rao B)^T
        CB = np.zeros((J * K, rank), dtype=np.float64)
        for r in range(rank):
            CB[:, r] = np.kron(C[:, r], B[:, r])
        A = np.linalg.lstsq(CB.T @ CB + eps * np.eye(rank), CB.T @ T0.T, rcond=None)[0].T

        # Update B: T_(1) ~ B (C khatri-rao A)^T
        CA = np.zeros((I * K, rank), dtype=np.float64)
        for r in range(rank):
            CA[:, r] = np.kron(C[:, r], A[:, r])
        B = np.linalg.lstsq(CA.T @ CA + eps * np.eye(rank), CA.T @ T1.T, rcond=None)[0].T

        # Update C: T_(2) ~ C (B khatri-rao A)^T
        BA = np.zeros((I * J, rank), dtype=np.float64)
        for r in range(rank):
            BA[:, r] = np.kron(B[:, r], A[:, r])
        C = np.linalg.lstsq(BA.T @ BA + eps * np.eye(rank), BA.T @ T2.T, rcond=None)[0].T

        if iteration % 20 == 0:
            logger.debug("CP-ALS iteration %d", iteration)

    # Reconstruct tensor to measure fit
    T_hat = np.zeros_like(T)
    for r in range(rank):
        T_hat += np.einsum("i,j,k->ijk", A[:, r], B[:, r], C[:, r])

    residual_sq = float(np.sum((T - T_hat) ** 2))
    fit = 1.0 - residual_sq / (T_norm_sq + eps)

    # Core consistency diagnostic (Tucker core should be superdiagonal)
    # Simplified: compare reconstructed core to ideal superdiagonal
    try:
        Q_A, _ = np.linalg.qr(A)
        Q_B, _ = np.linalg.qr(B)
        Q_C, _ = np.linalg.qr(C)
        # Project tensor into factor space
        core = np.einsum("ijk,ia,jb,kc->abc", T, Q_A[:, :rank], Q_B[:, :rank], Q_C[:, :rank])
        ideal_core = np.zeros_like(core)
        for r in range(min(rank, core.shape[0], core.shape[1], core.shape[2])):
            ideal_core[r, r, r] = core[r, r, r]
        core_norm = float(np.sum(core**2))
        core_consistency = 100.0 * (1.0 - float(np.sum((core - ideal_core) ** 2)) / (core_norm + eps))
        core_consistency = max(0.0, min(100.0, core_consistency))
    except Exception:
        core_consistency = float("nan")
        logger.warning("Could not compute core consistency diagnostic")

    logger.info("CP decomposition complete: fit=%.4f, core_consistency=%.1f%%", fit, core_consistency)

    return {
        "factors": [A.tolist(), B.tolist(), C.tolist()],
        "fit": fit,
        "core_consistency": core_consistency,
    }


# ---------------------------------------------------------------------------
# Similarity Network Fusion (SNF)
# ---------------------------------------------------------------------------


def similarity_network_fusion(
    networks: list[list[list[float]]] | list[Any],
    k_neighbors: int = 20,
    n_iter: int = 20,
    alpha: float = 0.5,
) -> dict[str, Any]:
    """Similarity Network Fusion (SNF) for multi-omic patient networks.

    Fuses multiple patient similarity networks from different omic layers
    into a single integrated network.  The algorithm iteratively updates
    each network by diffusing information through a shared neighbourhood
    structure.

    Steps:
        1. Compute normalised weight matrices (kernel) from each similarity.
        2. Compute local neighbourhood kernel (top-k neighbours).
        3. Iteratively update: P_m = S_m @ (avg of other P_j) @ S_m^T.
        4. Final fused network = average of all converged P_m.

    Args:
        networks: List of symmetric similarity matrices (one per omic layer).
            Each matrix is n_patients x n_patients as a numpy array or
            list-of-lists.
        k_neighbors: Number of neighbours for local kernel construction.
        n_iter: Number of diffusion iterations.
        alpha: Scaling parameter for the kernel normalisation (controls
            kernel bandwidth).  Higher values produce sparser kernels.

    Returns:
        Dictionary with keys:
            fused_network: Fused similarity matrix (n x n).
            cluster_labels: Cluster labels from spectral clustering of the
                fused network (uses 2 clusters by default).
            silhouette_score: Silhouette score of the clustering.

    Raises:
        ValueError: If fewer than 2 networks are provided or matrices have
            different shapes.
    """
    if len(networks) < 2:
        raise ValueError("SNF requires at least 2 networks to fuse")

    logger.info("Running SNF with %d networks, k=%d, n_iter=%d", len(networks), k_neighbors, n_iter)

    if not HAS_NUMPY:
        raise ImportError("numpy is required for SNF. Install with: uv pip install numpy")

    # Convert to numpy
    net_arrays: list[Any] = []
    n: int | None = None
    for i, net in enumerate(networks):
        arr = np.asarray(net, dtype=np.float64)
        if arr.ndim != 2 or arr.shape[0] != arr.shape[1]:
            raise ValueError(f"Network {i} must be a square matrix")
        if n is None:
            n = arr.shape[0]
        elif arr.shape[0] != n:
            raise ValueError(f"Network {i} has {arr.shape[0]} nodes, expected {n}")
        net_arrays.append(arr)

    assert n is not None
    eps = 1e-10
    k_neighbors = min(k_neighbors, n - 1)

    # --- Step 1: normalised full kernel ---
    def _normalise_kernel(W: Any) -> Any:
        """Row-normalise a similarity matrix to create a transition matrix."""
        row_sums = W.sum(axis=1, keepdims=True)
        return W / (row_sums + eps)

    # --- Step 2: local (sparse) kernel based on k nearest neighbours ---
    def _local_kernel(W: Any, k: int) -> Any:
        """Construct sparse local affinity using k nearest neighbours."""
        S = np.zeros_like(W)
        for i_node in range(W.shape[0]):
            row = W[i_node, :].copy()
            row[i_node] = 0  # exclude self
            if k >= len(row):
                neighbours = np.arange(len(row))
            else:
                neighbours = np.argpartition(row, -k)[-k:]
            S[i_node, neighbours] = row[neighbours]
            S[i_node, i_node] = 0  # no self-loops in local kernel
        # Symmetrise
        S = (S + S.T) / 2.0
        return _normalise_kernel(S)

    # Build normalised kernels
    P_list = [_normalise_kernel(W) for W in net_arrays]
    S_list = [_local_kernel(W, k_neighbors) for W in net_arrays]

    n_views = len(P_list)

    # --- Step 3: iterative diffusion ---
    for iteration in range(n_iter):
        P_new: list[Any] = []
        for m in range(n_views):
            # Average of all other P matrices
            other_avg = np.zeros((n, n), dtype=np.float64)
            for j in range(n_views):
                if j != m:
                    other_avg += P_list[j]
            other_avg /= n_views - 1
            # Diffusion: P_m = S_m @ other_avg @ S_m^T
            updated = S_list[m] @ other_avg @ S_list[m].T
            # Re-normalise
            updated = _normalise_kernel(updated)
            P_new.append(updated)
        P_list = P_new

    # --- Step 4: fuse by averaging ---
    fused = np.zeros((n, n), dtype=np.float64)
    for P in P_list:
        fused += P
    fused /= n_views
    # Symmetrise
    fused = (fused + fused.T) / 2.0

    # --- Spectral clustering on fused network ---
    cluster_labels, silhouette = _spectral_cluster_from_similarity(fused, n_clusters=2)

    logger.info("SNF complete: silhouette=%.4f", silhouette)

    return {
        "fused_network": fused.tolist(),
        "cluster_labels": cluster_labels,
        "silhouette_score": silhouette,
    }


def _spectral_cluster_from_similarity(similarity: Any, n_clusters: int = 2) -> tuple[list[int], float]:
    """Spectral clustering from a similarity matrix.

    Args:
        similarity: Symmetric similarity matrix (n x n).
        n_clusters: Number of clusters.

    Returns:
        Tuple of (cluster labels list, silhouette score).
    """
    n = similarity.shape[0]
    if n < n_clusters:
        return list(range(n)), 0.0

    # Laplacian
    D = np.diag(similarity.sum(axis=1))
    L = D - similarity
    D_inv_sqrt = np.diag(1.0 / (np.sqrt(np.diag(D)) + 1e-10))
    L_norm = D_inv_sqrt @ L @ D_inv_sqrt

    # Eigendecomposition
    eigenvalues, eigenvectors = np.linalg.eigh(L_norm)
    # Take eigenvectors corresponding to smallest eigenvalues
    idx = np.argsort(eigenvalues)[:n_clusters]
    embedding = eigenvectors[:, idx]

    # Normalise rows
    row_norms = np.linalg.norm(embedding, axis=1, keepdims=True)
    embedding = embedding / (row_norms + 1e-10)

    # K-means on the embedding
    labels = _simple_kmeans(embedding, n_clusters)

    # Silhouette score
    silhouette = _compute_silhouette(similarity, labels)

    return labels, silhouette


def _simple_kmeans(X: Any, k: int, max_iter: int = 50) -> list[int]:
    """Simple k-means clustering.

    Args:
        X: Data matrix (n x d).
        k: Number of clusters.
        max_iter: Maximum iterations.

    Returns:
        List of cluster labels.
    """
    n = X.shape[0]
    rng = np.random.RandomState(42)
    indices = rng.choice(n, size=k, replace=False)
    centroids = X[indices].copy()

    labels = np.zeros(n, dtype=int)
    for _ in range(max_iter):
        # Assign
        for i in range(n):
            dists = [float(np.sum((X[i] - centroids[c]) ** 2)) for c in range(k)]
            labels[i] = int(np.argmin(dists))
        # Update centroids
        new_centroids = np.zeros_like(centroids)
        for c in range(k):
            members = X[labels == c]
            if len(members) > 0:
                new_centroids[c] = members.mean(axis=0)
            else:
                new_centroids[c] = centroids[c]
        if np.allclose(centroids, new_centroids):
            break
        centroids = new_centroids

    return labels.tolist()


def _compute_silhouette(similarity: Any, labels: list[int]) -> float:
    """Compute silhouette score from similarity matrix and labels.

    Args:
        similarity: Symmetric similarity matrix.
        labels: Cluster labels.

    Returns:
        Mean silhouette score.
    """
    n = len(labels)
    if n < 2:
        return 0.0

    unique_labels = list(set(labels))
    if len(unique_labels) < 2:
        return 0.0

    # Convert similarity to distance
    distance = 1.0 - similarity
    np.fill_diagonal(distance, 0.0)

    silhouettes: list[float] = []
    for i in range(n):
        label_i = labels[i]
        # a(i): mean distance to same-cluster points
        same = [j for j in range(n) if labels[j] == label_i and j != i]
        if not same:
            silhouettes.append(0.0)
            continue
        a_i = float(np.mean(distance[i, same]))

        # b(i): min mean distance to other cluster
        b_i = float("inf")
        for other_label in unique_labels:
            if other_label == label_i:
                continue
            other = [j for j in range(n) if labels[j] == other_label]
            if other:
                b_i = min(b_i, float(np.mean(distance[i, other])))

        denom = max(a_i, b_i)
        if denom > 0:
            silhouettes.append((b_i - a_i) / denom)
        else:
            silhouettes.append(0.0)

    return float(np.mean(silhouettes))


# ---------------------------------------------------------------------------
# Regularised Canonical Correlation Analysis
# ---------------------------------------------------------------------------


def canonical_correlation(
    X: Any,
    Y: Any,
    n_components: int = 2,
    regularization: float = 0.1,
) -> dict[str, Any]:
    """Regularised Canonical Correlation Analysis for two omic matrices.

    Finds linear combinations of variables in X and Y that are maximally
    correlated.  Regularisation (ridge penalty) is applied to the
    covariance matrices to handle high-dimensional / collinear data.

    SVD-based implementation:
        1. Compute within- and between-view covariance matrices.
        2. Add regularisation to within-view covariances.
        3. Solve the generalised eigenvalue problem via SVD of the
           whitened cross-covariance matrix.

    Args:
        X: First omic data matrix (n_samples x p) as numpy array or
            list-of-lists.
        Y: Second omic data matrix (n_samples x q) as numpy array or
            list-of-lists.
        n_components: Number of canonical components to extract.
        regularization: Ridge regularisation parameter (lambda >= 0).

    Returns:
        Dictionary with keys:
            x_scores: Canonical scores for X (n_samples x n_components).
            y_scores: Canonical scores for Y (n_samples x n_components).
            correlations: Canonical correlations per component.
            x_loadings: Loading matrix for X (p x n_components).
            y_loadings: Loading matrix for Y (q x n_components).

    Raises:
        ValueError: If X and Y have different numbers of samples, or
            n_components < 1.
        ImportError: If numpy is not available.
    """
    if n_components < 1:
        raise ValueError(f"n_components must be >= 1, got {n_components}")

    logger.info(
        "Running regularised CCA with n_components=%d, reg=%.4f",
        n_components,
        regularization,
    )

    if not HAS_NUMPY:
        raise ImportError("numpy is required for CCA. Install with: uv pip install numpy")

    X_arr = np.asarray(X, dtype=np.float64)
    Y_arr = np.asarray(Y, dtype=np.float64)

    if X_arr.ndim != 2 or Y_arr.ndim != 2:
        raise ValueError("X and Y must be 2-D matrices")
    if X_arr.shape[0] != Y_arr.shape[0]:
        raise ValueError(
            f"Sample-dimension mismatch: X has {X_arr.shape[0]} samples, " f"Y has {Y_arr.shape[0]} samples"
        )

    n = X_arr.shape[0]
    p = X_arr.shape[1]
    q = Y_arr.shape[1]
    n_components = min(n_components, p, q, n)

    # Center
    X_c = X_arr - X_arr.mean(axis=0, keepdims=True)
    Y_c = Y_arr - Y_arr.mean(axis=0, keepdims=True)

    # Covariance matrices
    Cxx = (X_c.T @ X_c) / (n - 1) + regularization * np.eye(p)
    Cyy = (Y_c.T @ Y_c) / (n - 1) + regularization * np.eye(q)
    Cxy = (X_c.T @ Y_c) / (n - 1)

    # Whiten X and Y
    Ux, Sx, _ = np.linalg.svd(Cxx)
    Uy, Sy, _ = np.linalg.svd(Cyy)

    # Regularised inverse square root
    Sx_inv_sqrt = np.diag(1.0 / np.sqrt(Sx + 1e-10))
    Sy_inv_sqrt = np.diag(1.0 / np.sqrt(Sy + 1e-10))

    Cxx_inv_sqrt = Ux @ Sx_inv_sqrt @ Ux.T
    Cyy_inv_sqrt = Uy @ Sy_inv_sqrt @ Uy.T

    # SVD of whitened cross-covariance
    M = Cxx_inv_sqrt @ Cxy @ Cyy_inv_sqrt
    U, S, Vt = np.linalg.svd(M, full_matrices=False)

    # Canonical directions
    A = Cxx_inv_sqrt @ U[:, :n_components]  # (p, n_components)
    B = Cyy_inv_sqrt @ Vt[:n_components, :].T  # (q, n_components)

    # Canonical correlations
    correlations = S[:n_components].tolist()

    # Scores
    x_scores = (X_c @ A).tolist()
    y_scores = (Y_c @ B).tolist()

    # Loadings: correlation of original variables with canonical scores
    x_scores_arr = X_c @ A
    y_scores_arr = Y_c @ B

    x_loadings = np.zeros((p, n_components))
    for j in range(p):
        for c in range(n_components):
            x_loadings[j, c] = float(np.corrcoef(X_c[:, j], x_scores_arr[:, c])[0, 1])

    y_loadings = np.zeros((q, n_components))
    for j in range(q):
        for c in range(n_components):
            y_loadings[j, c] = float(np.corrcoef(Y_c[:, j], y_scores_arr[:, c])[0, 1])

    logger.info("CCA complete: correlations=%s", correlations)

    return {
        "x_scores": x_scores,
        "y_scores": y_scores,
        "correlations": correlations,
        "x_loadings": x_loadings.tolist(),
        "y_loadings": y_loadings.tolist(),
    }
