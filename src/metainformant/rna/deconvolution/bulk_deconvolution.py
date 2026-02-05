"""Cell type deconvolution from bulk RNA-seq expression data.

This module provides algorithms for estimating the proportions of different
cell types in bulk RNA-seq samples using reference expression signatures.
Includes NNLS (non-negative least squares) deconvolution with a pure Python
projected gradient descent implementation, SVR-based deconvolution following
the CIBERSORT methodology, signature matrix construction from reference
profiles, marker gene selection, result validation, and batch processing.

All core algorithms have pure Python fallbacks; numpy/scipy/sklearn are
used when available for performance.

Main Functions:
    Deconvolution:
        - deconvolve_nnls: Non-negative least squares deconvolution
        - deconvolve_svr: SVR-based deconvolution (CIBERSORT-style)
        - batch_deconvolve: Deconvolve multiple samples

    Reference Construction:
        - build_signature_matrix: Build cell-type signature matrix
        - select_marker_genes: Select informative marker genes

    Validation:
        - validate_deconvolution: Assess deconvolution accuracy

Example:
    >>> from metainformant.rna.deconvolution import bulk_deconvolution
    >>> mixture = [10.5, 3.2, 8.1, 1.5, 6.7]
    >>> signature = [[8.0, 2.0], [2.5, 1.0], [7.0, 3.0], [1.0, 0.5], [5.0, 4.0]]
    >>> result = bulk_deconvolution.deconvolve_nnls(mixture, signature)
"""

from __future__ import annotations

import math
from collections import defaultdict
from typing import Any, Dict, List, Literal, Optional, Sequence, Tuple, Union

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Optional dependency handling
try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None  # type: ignore[assignment]

try:
    from scipy import stats as scipy_stats
    from scipy.optimize import nnls as scipy_nnls

    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    scipy_stats = None  # type: ignore[assignment]
    scipy_nnls = None  # type: ignore[assignment]

try:
    from sklearn.svm import NuSVR

    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False
    NuSVR = None  # type: ignore[assignment]


# =============================================================================
# NNLS Deconvolution
# =============================================================================


def deconvolve_nnls(
    mixture: list[float] | Any,
    signature_matrix: list[list[float]] | Any,
    normalize: bool = True,
) -> dict:
    """Estimate cell type proportions using non-negative least squares.

    Solves the deconvolution problem: mixture ~ signature_matrix @ proportions,
    subject to proportions >= 0. Uses scipy.optimize.nnls when available,
    with a pure Python projected gradient descent fallback.

    The proportions are optionally normalized to sum to 1, representing
    fractional composition of the bulk sample.

    Args:
        mixture: Expression profile of the bulk sample. Either a list of
            floats or a numpy array of shape (n_genes,).
        signature_matrix: Reference expression profiles for each cell type.
            Either a list of lists (n_genes x n_cell_types) or a numpy
            array of shape (n_genes, n_cell_types). Each column represents
            a cell type's expression signature.
        normalize: Whether to normalize estimated proportions to sum to 1.0.

    Returns:
        Dictionary with keys:
            - proportions (dict): Mapping of cell type index (int) to
              estimated proportion (float). If normalize=True, sums to 1.0.
            - residual (float): Residual norm of the NNLS fit.
            - r_squared (float): Coefficient of determination (R-squared).
            - rmse (float): Root mean squared error of the fit.
            - method (str): Algorithm used ("scipy_nnls" or "projected_gradient").
            - n_genes (int): Number of genes used.
            - n_cell_types (int): Number of cell types.

    Raises:
        ValueError: If dimensions are incompatible or inputs are empty.

    Example:
        >>> mixture = [10.0, 3.0, 8.0, 1.5, 6.0]
        >>> signature = [[8.0, 2.0], [2.5, 1.0], [7.0, 3.0], [1.0, 0.5], [5.0, 4.0]]
        >>> result = deconvolve_nnls(mixture, signature)
        >>> len(result["proportions"])
        2
    """
    # Convert to numpy arrays if available
    mix_arr, sig_arr, n_genes, n_types = _validate_and_convert(mixture, signature_matrix)

    # Solve NNLS
    if HAS_NUMPY and HAS_SCIPY:
        method = "scipy_nnls"
        proportions_raw, residual = scipy_nnls(sig_arr, mix_arr)
    elif HAS_NUMPY:
        method = "projected_gradient"
        proportions_raw, residual = _nnls_projected_gradient_numpy(sig_arr, mix_arr)
    else:
        method = "projected_gradient"
        proportions_raw, residual = _nnls_projected_gradient_pure(
            _to_nested_list(signature_matrix, n_genes, n_types),
            _to_flat_list(mixture, n_genes),
        )

    # Convert raw solution to list
    if HAS_NUMPY:
        props_list = proportions_raw.tolist()
    else:
        props_list = list(proportions_raw)

    # Normalize proportions
    if normalize:
        total = sum(props_list)
        if total > 0:
            props_list = [p / total for p in props_list]

    # Compute fit quality metrics
    predicted = _matrix_vector_multiply(
        _to_nested_list(signature_matrix, n_genes, n_types) if not HAS_NUMPY else sig_arr,
        props_list if not HAS_NUMPY else proportions_raw,
    )

    mix_list = _to_flat_list(mixture, n_genes)
    r_squared = _compute_r_squared(mix_list, predicted)
    rmse = _compute_rmse(mix_list, predicted)

    proportions_dict = {i: round(p, 6) for i, p in enumerate(props_list)}

    return {
        "proportions": proportions_dict,
        "residual": round(float(residual), 6),
        "r_squared": round(r_squared, 6),
        "rmse": round(rmse, 6),
        "method": method,
        "n_genes": n_genes,
        "n_cell_types": n_types,
    }


def _validate_and_convert(
    mixture: Any,
    signature_matrix: Any,
) -> tuple[Any, Any, int, int]:
    """Validate and convert inputs to working format.

    Args:
        mixture: Mixture expression vector.
        signature_matrix: Signature matrix.

    Returns:
        Tuple of (mixture_array, signature_array, n_genes, n_cell_types).

    Raises:
        ValueError: If dimensions are incompatible.
    """
    if HAS_NUMPY:
        mix_arr = np.asarray(mixture, dtype=float).ravel()
        sig_arr = np.asarray(signature_matrix, dtype=float)

        if sig_arr.ndim != 2:
            raise ValueError(f"Signature matrix must be 2D, got {sig_arr.ndim}D")

        n_genes = sig_arr.shape[0]
        n_types = sig_arr.shape[1]

        if len(mix_arr) != n_genes:
            raise ValueError(f"Mixture length ({len(mix_arr)}) doesn't match signature " f"rows ({n_genes})")

        return mix_arr, sig_arr, n_genes, n_types

    # Pure Python path
    if isinstance(mixture, (list, tuple)):
        mix_list = [float(x) for x in mixture]
    else:
        mix_list = [float(x) for x in list(mixture)]

    if isinstance(signature_matrix, (list, tuple)):
        sig_list = [[float(x) for x in row] for row in signature_matrix]
    else:
        sig_list = [[float(x) for x in row] for row in list(signature_matrix)]

    n_genes = len(sig_list)
    if n_genes == 0:
        raise ValueError("Signature matrix is empty")

    n_types = len(sig_list[0])
    if n_types == 0:
        raise ValueError("Signature matrix has no cell types (columns)")

    if len(mix_list) != n_genes:
        raise ValueError(f"Mixture length ({len(mix_list)}) doesn't match signature " f"rows ({n_genes})")

    return mix_list, sig_list, n_genes, n_types


def _nnls_projected_gradient_numpy(
    A: Any,
    b: Any,
    max_iterations: int = 5000,
    tolerance: float = 1e-8,
) -> tuple[Any, float]:
    """NNLS solver using projected gradient descent with numpy.

    Args:
        A: Signature matrix (n_genes x n_types).
        b: Mixture vector (n_genes,).
        max_iterations: Maximum iterations.
        tolerance: Convergence tolerance.

    Returns:
        Tuple of (solution vector, residual).
    """
    n_types = A.shape[1]

    # Initialize with least squares solution, then project
    x = np.zeros(n_types)

    # Precompute A^T A and A^T b
    AtA = A.T @ A
    Atb = A.T @ b

    # Step size from Lipschitz constant
    L = np.linalg.norm(AtA, ord=2)
    step_size = 1.0 / L if L > 0 else 0.01

    for _ in range(max_iterations):
        gradient = AtA @ x - Atb
        x_new = np.maximum(0, x - step_size * gradient)

        if np.max(np.abs(x_new - x)) < tolerance:
            x = x_new
            break

        x = x_new

    residual = float(np.linalg.norm(b - A @ x))
    return x, residual


def _nnls_projected_gradient_pure(
    A: list[list[float]],
    b: list[float],
    max_iterations: int = 5000,
    tolerance: float = 1e-8,
) -> tuple[list[float], float]:
    """NNLS solver using projected gradient descent, pure Python.

    Args:
        A: Signature matrix (n_genes x n_types).
        b: Mixture vector (n_genes,).
        max_iterations: Maximum iterations.
        tolerance: Convergence tolerance.

    Returns:
        Tuple of (solution list, residual).
    """
    n_genes = len(A)
    n_types = len(A[0])

    # Compute A^T A
    AtA = [[0.0] * n_types for _ in range(n_types)]
    for i in range(n_types):
        for j in range(n_types):
            for k in range(n_genes):
                AtA[i][j] += A[k][i] * A[k][j]

    # Compute A^T b
    Atb = [0.0] * n_types
    for i in range(n_types):
        for k in range(n_genes):
            Atb[i] += A[k][i] * b[k]

    # Estimate step size (1/spectral_norm approximation)
    max_diag = max(AtA[i][i] for i in range(n_types))
    step_size = 1.0 / (max_diag * n_types) if max_diag > 0 else 0.01

    x = [0.0] * n_types

    for _ in range(max_iterations):
        # Gradient: AtA @ x - Atb
        gradient = [0.0] * n_types
        for i in range(n_types):
            for j in range(n_types):
                gradient[i] += AtA[i][j] * x[j]
            gradient[i] -= Atb[i]

        # Projected gradient step
        x_new = [max(0.0, x[i] - step_size * gradient[i]) for i in range(n_types)]

        # Check convergence
        max_change = max(abs(x_new[i] - x[i]) for i in range(n_types))
        x = x_new

        if max_change < tolerance:
            break

    # Compute residual
    residual = 0.0
    for k in range(n_genes):
        pred = sum(A[k][j] * x[j] for j in range(n_types))
        residual += (b[k] - pred) ** 2
    residual = math.sqrt(residual)

    return x, residual


def _to_flat_list(data: Any, expected_len: int) -> list[float]:
    """Convert data to a flat list of floats."""
    if HAS_NUMPY and isinstance(data, np.ndarray):
        return data.ravel().tolist()
    return [float(x) for x in data]


def _to_nested_list(data: Any, n_rows: int, n_cols: int) -> list[list[float]]:
    """Convert data to nested list of floats."""
    if HAS_NUMPY and isinstance(data, np.ndarray):
        return data.tolist()
    return [[float(x) for x in row] for row in data]


def _matrix_vector_multiply(matrix: Any, vector: Any) -> list[float]:
    """Multiply matrix by vector, return list."""
    if HAS_NUMPY:
        m = np.asarray(matrix, dtype=float)
        v = np.asarray(vector, dtype=float)
        return (m @ v).tolist()

    n_rows = len(matrix)
    n_cols = len(matrix[0]) if n_rows > 0 else 0
    result = []
    for i in range(n_rows):
        val = sum(matrix[i][j] * vector[j] for j in range(n_cols))
        result.append(val)
    return result


def _compute_r_squared(observed: list[float], predicted: list[float]) -> float:
    """Compute R-squared (coefficient of determination).

    Args:
        observed: Observed values.
        predicted: Predicted values.

    Returns:
        R-squared value. Can be negative for very poor fits.
    """
    n = len(observed)
    if n == 0:
        return 0.0

    mean_obs = sum(observed) / n
    ss_res = sum((o - p) ** 2 for o, p in zip(observed, predicted))
    ss_tot = sum((o - mean_obs) ** 2 for o in observed)

    if ss_tot < 1e-15:
        return 1.0 if ss_res < 1e-15 else 0.0

    return 1.0 - ss_res / ss_tot


def _compute_rmse(observed: list[float], predicted: list[float]) -> float:
    """Compute root mean squared error.

    Args:
        observed: Observed values.
        predicted: Predicted values.

    Returns:
        RMSE value.
    """
    n = len(observed)
    if n == 0:
        return 0.0

    mse = sum((o - p) ** 2 for o, p in zip(observed, predicted)) / n
    return math.sqrt(mse)


# =============================================================================
# SVR-based Deconvolution (CIBERSORT-style)
# =============================================================================


def deconvolve_svr(
    mixture: list[float] | Any,
    signature_matrix: list[list[float]] | Any,
    nu: float = 0.5,
    normalize: bool = True,
) -> dict:
    """Estimate cell type proportions using support vector regression.

    Implements a CIBERSORT-style deconvolution approach using nu-SVR.
    Each cell type's proportion is estimated by training a separate SVR
    model. When sklearn is not available, falls back to NNLS deconvolution.

    Args:
        mixture: Expression profile of the bulk sample. Either a list of
            floats or a numpy array of shape (n_genes,).
        signature_matrix: Reference expression profiles for each cell type.
            Either a list of lists (n_genes x n_cell_types) or a numpy
            array of shape (n_genes, n_cell_types).
        nu: Nu parameter for NuSVR, controlling the fraction of support
            vectors. Must be in (0, 1]. Lower values = fewer support vectors.
        normalize: Whether to normalize proportions to sum to 1.0.

    Returns:
        Dictionary with keys:
            - proportions (dict): Mapping of cell type index to proportion.
            - residual (float): Residual norm of the fit.
            - r_squared (float): Coefficient of determination.
            - rmse (float): Root mean squared error.
            - method (str): "svr_cibersort" or "nnls_fallback".
            - correlation (float): Pearson correlation between mixture and
              reconstructed signal.
            - p_value (float): Empirical p-value from permutation test
              (only if n_genes > 20, otherwise 0.0).
            - n_genes (int): Number of genes used.
            - n_cell_types (int): Number of cell types.

    Raises:
        ValueError: If dimensions are incompatible or nu is out of range.

    Example:
        >>> mixture = [10.0, 3.0, 8.0, 1.5, 6.0]
        >>> signature = [[8.0, 2.0], [2.5, 1.0], [7.0, 3.0], [1.0, 0.5], [5.0, 4.0]]
        >>> result = deconvolve_svr(mixture, signature, nu=0.5)
    """
    if not (0 < nu <= 1):
        raise ValueError(f"nu must be in (0, 1], got {nu}")

    mix_arr, sig_arr, n_genes, n_types = _validate_and_convert(mixture, signature_matrix)

    if not HAS_SKLEARN or not HAS_NUMPY:
        logger.warning("sklearn not available, falling back to NNLS deconvolution")
        result = deconvolve_nnls(mixture, signature_matrix, normalize=normalize)
        result["method"] = "nnls_fallback"
        result["correlation"] = 0.0
        result["p_value"] = 1.0
        return result

    # Ensure numpy arrays
    mix_np = np.asarray(mix_arr, dtype=float).ravel()
    sig_np = np.asarray(sig_arr, dtype=float)

    # CIBERSORT approach: fit SVR with signature as features, mixture as target
    # Train: each gene is a sample, cell types are features
    svr = NuSVR(nu=nu, kernel="linear", C=1.0)

    try:
        svr.fit(sig_np, mix_np)
    except Exception as e:
        logger.warning(f"SVR fitting failed: {e}, falling back to NNLS")
        result = deconvolve_nnls(mixture, signature_matrix, normalize=normalize)
        result["method"] = "nnls_fallback"
        result["correlation"] = 0.0
        result["p_value"] = 1.0
        return result

    # Extract coefficients (cell type weights)
    coefs = svr.coef_.ravel()

    # Enforce non-negativity
    coefs = np.maximum(0, coefs)

    # Normalize if requested
    coefs_list = coefs.tolist()
    if normalize:
        total = sum(coefs_list)
        if total > 0:
            coefs_list = [c / total for c in coefs_list]

    # Reconstruct mixture from estimated proportions
    reconstructed = (sig_np @ np.array(coefs_list)).tolist()

    # Fit quality
    mix_list = mix_np.tolist()
    r_squared = _compute_r_squared(mix_list, reconstructed)
    rmse = _compute_rmse(mix_list, reconstructed)
    residual = float(np.linalg.norm(mix_np - sig_np @ np.array(coefs_list)))

    # Pearson correlation between mixture and reconstructed
    correlation = _pearson_corr_numpy(mix_np, np.array(reconstructed))

    # Empirical p-value via permutation (CIBERSORT-style)
    p_value = 0.0
    if n_genes > 20:
        p_value = _permutation_pvalue(mix_np, sig_np, correlation, n_permutations=100)

    proportions_dict = {i: round(p, 6) for i, p in enumerate(coefs_list)}

    return {
        "proportions": proportions_dict,
        "residual": round(residual, 6),
        "r_squared": round(r_squared, 6),
        "rmse": round(rmse, 6),
        "method": "svr_cibersort",
        "correlation": round(correlation, 6),
        "p_value": round(p_value, 6),
        "n_genes": n_genes,
        "n_cell_types": n_types,
    }


def _pearson_corr_numpy(x: Any, y: Any) -> float:
    """Compute Pearson correlation using numpy.

    Args:
        x: First array.
        y: Second array.

    Returns:
        Correlation coefficient.
    """
    if len(x) < 2:
        return 0.0

    corr_matrix = np.corrcoef(x, y)
    r = float(corr_matrix[0, 1])
    return r if not np.isnan(r) else 0.0


def _permutation_pvalue(
    mixture: Any,
    signature: Any,
    observed_corr: float,
    n_permutations: int = 100,
) -> float:
    """Compute empirical p-value by permuting the mixture vector.

    Args:
        mixture: Original mixture vector.
        signature: Signature matrix.
        observed_corr: Observed correlation.
        n_permutations: Number of permutations.

    Returns:
        Empirical p-value.
    """
    rng = np.random.default_rng(42)
    count_better = 0

    for _ in range(n_permutations):
        perm_mix = rng.permutation(mixture)
        try:
            svr = NuSVR(nu=0.5, kernel="linear", C=1.0)
            svr.fit(signature, perm_mix)
            coefs = np.maximum(0, svr.coef_.ravel())
            total = coefs.sum()
            if total > 0:
                coefs = coefs / total
            reconstructed = signature @ coefs
            perm_corr = _pearson_corr_numpy(perm_mix, reconstructed)
            if perm_corr >= observed_corr:
                count_better += 1
        except Exception:
            pass

    return (count_better + 1) / (n_permutations + 1)


# =============================================================================
# Signature Matrix Construction
# =============================================================================


def build_signature_matrix(
    expression_profiles: dict,
    cell_types: list[str],
    n_markers: int = 50,
    method: str = "fold_change",
) -> dict:
    """Build a cell-type signature matrix from reference expression profiles.

    Selects informative marker genes for each cell type and constructs
    a signature matrix suitable for deconvolution algorithms.

    Args:
        expression_profiles: Dict mapping cell type labels to expression dicts.
            Structure: {cell_type: {gene_id: expression_value, ...}, ...}
            Each cell type can have multiple samples stored as lists:
            {cell_type: {gene_id: [sample1_expr, sample2_expr, ...], ...}, ...}
        cell_types: List of cell type labels to include in the signature.
            Must be keys in expression_profiles.
        n_markers: Number of marker genes to select per cell type.
        method: Marker selection method:
            - "fold_change": Select genes with highest fold change (default)
            - "specificity": Select genes with highest specificity index
            - "t_test": Select by t-test (requires list-valued expression)

    Returns:
        Dictionary with keys:
            - matrix (list[list[float]]): Signature matrix (n_genes x n_cell_types).
              Rows are selected marker genes, columns are cell types.
            - gene_ids (list[str]): Gene identifiers for each row.
            - cell_types (list[str]): Cell type labels for each column.
            - markers_per_type (dict): Mapping of cell type to its marker genes.
            - n_genes (int): Total number of genes in the matrix.
            - n_cell_types (int): Number of cell types.

    Raises:
        ValueError: If expression_profiles is empty, cell_types not found,
            or no marker genes can be selected.

    Example:
        >>> profiles = {
        ...     "T_cell": {"CD3": 100, "CD4": 80, "CD19": 2, "HBA1": 1},
        ...     "B_cell": {"CD3": 5, "CD4": 2, "CD19": 95, "HBA1": 1},
        ... }
        >>> result = build_signature_matrix(profiles, ["T_cell", "B_cell"], n_markers=2)
        >>> result["n_cell_types"]
        2
    """
    if not expression_profiles:
        raise ValueError("Expression profiles dictionary is empty")

    for ct in cell_types:
        if ct not in expression_profiles:
            raise ValueError(f"Cell type '{ct}' not found in expression_profiles")

    # Select marker genes
    all_markers = select_marker_genes(expression_profiles, cell_types, n_genes=n_markers, method=method)

    if not all_markers:
        raise ValueError("No marker genes could be selected")

    # Build the signature matrix using mean expression per cell type
    # for the selected marker genes
    unique_markers = list(dict.fromkeys(all_markers))  # Preserve order, remove dupes

    matrix: list[list[float]] = []
    markers_per_type: dict[str, list[str]] = defaultdict(list)

    for gene_id in unique_markers:
        row: list[float] = []
        for ct in cell_types:
            expr = expression_profiles[ct].get(gene_id, 0)
            if isinstance(expr, (list, tuple)):
                mean_expr = sum(expr) / len(expr) if expr else 0.0
            else:
                mean_expr = float(expr)
            row.append(mean_expr)
        matrix.append(row)

    # Track which markers belong to which cell type
    # (determined by which type they have highest expression in)
    for gene_id in unique_markers:
        best_type = ""
        best_expr = -1.0
        for ct in cell_types:
            expr = expression_profiles[ct].get(gene_id, 0)
            if isinstance(expr, (list, tuple)):
                mean_expr = sum(expr) / len(expr) if expr else 0.0
            else:
                mean_expr = float(expr)
            if mean_expr > best_expr:
                best_expr = mean_expr
                best_type = ct
        markers_per_type[best_type].append(gene_id)

    logger.info(f"Built signature matrix: {len(unique_markers)} genes x " f"{len(cell_types)} cell types")

    return {
        "matrix": matrix,
        "gene_ids": unique_markers,
        "cell_types": cell_types,
        "markers_per_type": dict(markers_per_type),
        "n_genes": len(unique_markers),
        "n_cell_types": len(cell_types),
    }


# =============================================================================
# Marker Gene Selection
# =============================================================================


def select_marker_genes(
    expression: dict,
    labels: list[str],
    n_genes: int = 50,
    method: str = "fold_change",
) -> list[str]:
    """Select marker genes per cell type from reference expression profiles.

    Identifies genes that are specifically upregulated in one cell type
    compared to all others, making them informative markers for deconvolution.

    Args:
        expression: Dict mapping cell type labels to expression dicts.
            Structure: {cell_type: {gene_id: expression_value, ...}, ...}
        labels: List of cell type labels to consider.
        n_genes: Number of top marker genes to select per cell type.
        method: Selection method:
            - "fold_change": Select by maximum fold change vs other types.
              Genes with highest ratio of target type expression to mean
              of other types are selected.
            - "specificity": Select by specificity index (proportion of
              total expression attributed to the target type).
            - "t_test": Select by t-test p-value between target type and
              others (requires multiple samples per type as lists).

    Returns:
        Flat list of selected marker gene IDs. May contain up to
        n_genes * len(labels) genes (fewer if there are overlaps or
        insufficient genes).

    Raises:
        ValueError: If expression is empty or labels not found.

    Example:
        >>> expr = {
        ...     "type_a": {"g1": 100, "g2": 5, "g3": 50},
        ...     "type_b": {"g1": 3, "g2": 90, "g3": 40},
        ... }
        >>> markers = select_marker_genes(expr, ["type_a", "type_b"], n_genes=1)
        >>> "g1" in markers  # Top marker for type_a
        True
    """
    if not expression:
        raise ValueError("Expression dictionary is empty")

    for label in labels:
        if label not in expression:
            raise ValueError(f"Label '{label}' not found in expression")

    # Get all gene IDs across all cell types
    all_genes: set[str] = set()
    for label in labels:
        all_genes.update(expression[label].keys())

    if not all_genes:
        return []

    selected_markers: list[str] = []

    for target_type in labels:
        # Compute scores for each gene
        gene_scores: list[tuple[str, float]] = []

        for gene_id in all_genes:
            target_expr = _get_mean_expression(expression[target_type], gene_id)

            # Mean expression in other cell types
            other_exprs = [_get_mean_expression(expression[ct], gene_id) for ct in labels if ct != target_type]

            if method == "fold_change":
                other_mean = sum(other_exprs) / len(other_exprs) if other_exprs else 0
                # Fold change with pseudocount
                fc = (target_expr + 1.0) / (other_mean + 1.0)
                gene_scores.append((gene_id, fc))

            elif method == "specificity":
                total = target_expr + sum(other_exprs)
                specificity = target_expr / total if total > 0 else 0.0
                gene_scores.append((gene_id, specificity))

            elif method == "t_test":
                # Need list values for t-test
                target_vals = _get_expression_list(expression[target_type], gene_id)
                other_vals = []
                for ct in labels:
                    if ct != target_type:
                        other_vals.extend(_get_expression_list(expression[ct], gene_id))

                if len(target_vals) >= 2 and len(other_vals) >= 2:
                    t_stat = _simple_t_test(target_vals, other_vals)
                    gene_scores.append((gene_id, t_stat))
                else:
                    # Fall back to fold change
                    other_mean = sum(other_exprs) / len(other_exprs) if other_exprs else 0
                    fc = (target_expr + 1.0) / (other_mean + 1.0)
                    gene_scores.append((gene_id, fc))

            else:
                raise ValueError(f"Unknown method: {method}. Valid: fold_change, specificity, t_test")

        # Sort by score descending and select top n
        gene_scores.sort(key=lambda x: x[1], reverse=True)
        top_genes = [g for g, _ in gene_scores[:n_genes]]
        selected_markers.extend(top_genes)

    # Remove duplicates while preserving order
    seen: set[str] = set()
    unique_markers: list[str] = []
    for g in selected_markers:
        if g not in seen:
            seen.add(g)
            unique_markers.append(g)

    logger.info(
        f"Selected {len(unique_markers)} unique marker genes "
        f"({n_genes} per type, {len(labels)} types, method={method})"
    )

    return unique_markers


def _get_mean_expression(type_expr: dict, gene_id: str) -> float:
    """Get mean expression value for a gene in a cell type.

    Args:
        type_expr: Expression dict for one cell type.
        gene_id: Gene identifier.

    Returns:
        Mean expression value (float).
    """
    val = type_expr.get(gene_id, 0)
    if isinstance(val, (list, tuple)):
        return sum(val) / len(val) if val else 0.0
    return float(val)


def _get_expression_list(type_expr: dict, gene_id: str) -> list[float]:
    """Get expression values as a list.

    Args:
        type_expr: Expression dict for one cell type.
        gene_id: Gene identifier.

    Returns:
        List of expression values.
    """
    val = type_expr.get(gene_id, 0)
    if isinstance(val, (list, tuple)):
        return [float(v) for v in val]
    return [float(val)]


def _simple_t_test(group1: list[float], group2: list[float]) -> float:
    """Simple two-sample t-test statistic (Welch's).

    Args:
        group1: Values from group 1.
        group2: Values from group 2.

    Returns:
        t-statistic (positive means group1 > group2).
    """
    n1 = len(group1)
    n2 = len(group2)

    if n1 < 2 or n2 < 2:
        return 0.0

    mean1 = sum(group1) / n1
    mean2 = sum(group2) / n2

    var1 = sum((x - mean1) ** 2 for x in group1) / (n1 - 1)
    var2 = sum((x - mean2) ** 2 for x in group2) / (n2 - 1)

    se = math.sqrt(var1 / n1 + var2 / n2)

    if se < 1e-15:
        return 0.0

    return (mean1 - mean2) / se


# =============================================================================
# Deconvolution Validation
# =============================================================================


def validate_deconvolution(
    estimated: dict,
    ground_truth: dict | None = None,
) -> dict:
    """Validate deconvolution results with optional ground truth comparison.

    Assesses the quality of deconvolution estimates through internal
    consistency checks and, when ground truth is available, accuracy metrics.

    Args:
        estimated: Estimated proportions dict, mapping cell type identifiers
            (str or int) to proportion values (float).
        ground_truth: Optional ground truth proportions dict with same keys
            as estimated. If provided, accuracy metrics are computed.

    Returns:
        Dictionary with keys:
            - sum_proportions (float): Sum of all estimated proportions
              (should be close to 1.0 if normalized)
            - n_cell_types (int): Number of cell types estimated
            - n_nonzero (int): Number of cell types with non-zero proportion
            - max_proportion (float): Largest estimated proportion
            - max_proportion_type (str): Cell type with largest proportion
            - min_proportion (float): Smallest non-zero proportion
            - entropy (float): Shannon entropy of the proportion distribution
            - valid (bool): Whether proportions pass basic sanity checks
              (all non-negative, sum within [0.5, 1.5])
            If ground_truth is provided, also includes:
            - correlation (float): Pearson correlation with ground truth
            - rmse (float): Root mean squared error vs ground truth
            - mae (float): Mean absolute error vs ground truth
            - max_error (float): Maximum absolute error for any cell type
            - max_error_type (str): Cell type with maximum error

    Example:
        >>> estimated = {0: 0.6, 1: 0.3, 2: 0.1}
        >>> result = validate_deconvolution(estimated)
        >>> result["valid"]
        True
        >>> result["sum_proportions"]
        1.0
    """
    if not estimated:
        return {
            "sum_proportions": 0.0,
            "n_cell_types": 0,
            "n_nonzero": 0,
            "max_proportion": 0.0,
            "max_proportion_type": "",
            "min_proportion": 0.0,
            "entropy": 0.0,
            "valid": False,
        }

    values = list(estimated.values())
    keys = list(estimated.keys())

    sum_props = sum(values)
    n_types = len(values)
    n_nonzero = sum(1 for v in values if v > 0)

    # Find max/min
    max_val = max(values)
    max_type = str(keys[values.index(max_val)])
    nonzero_vals = [v for v in values if v > 0]
    min_val = min(nonzero_vals) if nonzero_vals else 0.0

    # Shannon entropy
    if sum_props > 0:
        proportions = [v / sum_props for v in values]
        entropy = -sum(p * math.log2(p) for p in proportions if p > 0)
    else:
        entropy = 0.0

    # Validity checks
    all_nonneg = all(v >= 0 for v in values)
    sum_reasonable = 0.5 <= sum_props <= 1.5  # Allow some tolerance
    valid = all_nonneg and sum_reasonable

    result: dict[str, Any] = {
        "sum_proportions": round(sum_props, 6),
        "n_cell_types": n_types,
        "n_nonzero": n_nonzero,
        "max_proportion": round(max_val, 6),
        "max_proportion_type": max_type,
        "min_proportion": round(min_val, 6),
        "entropy": round(entropy, 6),
        "valid": valid,
    }

    # Ground truth comparison
    if ground_truth is not None:
        all_keys = sorted(set(list(estimated.keys()) + list(ground_truth.keys())))
        est_vals = [estimated.get(k, 0.0) for k in all_keys]
        gt_vals = [ground_truth.get(k, 0.0) for k in all_keys]

        # Correlation
        n = len(all_keys)
        if n >= 2:
            mean_est = sum(est_vals) / n
            mean_gt = sum(gt_vals) / n
            cov = sum((e - mean_est) * (g - mean_gt) for e, g in zip(est_vals, gt_vals))
            var_est = sum((e - mean_est) ** 2 for e in est_vals)
            var_gt = sum((g - mean_gt) ** 2 for g in gt_vals)
            denom = math.sqrt(var_est * var_gt)
            correlation = cov / denom if denom > 1e-15 else 0.0
        else:
            correlation = 0.0

        # Error metrics
        errors = [abs(e - g) for e, g in zip(est_vals, gt_vals)]
        squared_errors = [(e - g) ** 2 for e, g in zip(est_vals, gt_vals)]

        rmse = math.sqrt(sum(squared_errors) / n) if n > 0 else 0.0
        mae = sum(errors) / n if n > 0 else 0.0
        max_error = max(errors) if errors else 0.0
        max_error_idx = errors.index(max_error) if errors else 0
        max_error_type = str(all_keys[max_error_idx]) if all_keys else ""

        result.update(
            {
                "correlation": round(correlation, 6),
                "rmse": round(rmse, 6),
                "mae": round(mae, 6),
                "max_error": round(max_error, 6),
                "max_error_type": max_error_type,
            }
        )

    return result


# =============================================================================
# Batch Deconvolution
# =============================================================================


def batch_deconvolve(
    samples: dict,
    signature_matrix: dict,
    method: str = "nnls",
) -> dict:
    """Deconvolve multiple bulk RNA-seq samples in batch.

    Applies deconvolution to each sample independently using the same
    signature matrix, and returns a matrix of cell type proportions.

    Args:
        samples: Dictionary mapping sample IDs to expression vectors.
            Structure: {sample_id: [gene1_expr, gene2_expr, ...], ...}
            or {sample_id: {gene_id: expr_value, ...}, ...}
            Gene order must match the signature matrix rows.
        signature_matrix: Pre-built signature matrix dict from
            build_signature_matrix() with keys "matrix", "gene_ids",
            "cell_types". Alternatively, a raw matrix (list of lists)
            can be passed.
        method: Deconvolution method:
            - "nnls": Non-negative least squares (default)
            - "svr": SVR-based (CIBERSORT-style)

    Returns:
        Dictionary with keys:
            - proportions (dict): Nested dict mapping sample_id to cell
              type proportions {sample_id: {cell_type: proportion, ...}}.
            - summary (dict): Summary statistics per cell type across all
              samples (mean, std, min, max).
            - n_samples (int): Number of samples processed.
            - n_cell_types (int): Number of cell types.
            - method (str): Method used.
            - failed_samples (list[str]): Sample IDs that failed deconvolution.

    Raises:
        ValueError: If samples or signature_matrix is empty.

    Example:
        >>> samples = {
        ...     "sample1": [10.0, 3.0, 8.0],
        ...     "sample2": [5.0, 7.0, 4.0],
        ... }
        >>> sig = {"matrix": [[8, 2], [2, 6], [7, 3]], "gene_ids": ["g1", "g2", "g3"],
        ...        "cell_types": ["type_a", "type_b"]}
        >>> result = batch_deconvolve(samples, sig, method="nnls")
        >>> result["n_samples"]
        2
    """
    if not samples:
        raise ValueError("No samples provided")

    # Extract signature matrix
    if isinstance(signature_matrix, dict) and "matrix" in signature_matrix:
        sig_mat = signature_matrix["matrix"]
        cell_type_labels = signature_matrix.get("cell_types", None)
        gene_ids = signature_matrix.get("gene_ids", None)
    else:
        sig_mat = signature_matrix
        cell_type_labels = None
        gene_ids = None

    # Determine deconvolution function
    if method == "nnls":
        deconvolve_fn = deconvolve_nnls
    elif method == "svr":
        deconvolve_fn = deconvolve_svr
    else:
        raise ValueError(f"Unknown method: {method}. Valid: nnls, svr")

    proportions: dict[str, dict] = {}
    failed_samples: list[str] = []

    for sample_id, sample_expr in samples.items():
        # Convert dict-based expression to ordered list if needed
        if isinstance(sample_expr, dict) and gene_ids is not None:
            mixture = [sample_expr.get(gid, 0.0) for gid in gene_ids]
        else:
            mixture = list(sample_expr)

        try:
            result = deconvolve_fn(mixture, sig_mat)

            # Map integer indices to cell type labels if available
            sample_props: dict[str, float] = {}
            for idx, prop in result["proportions"].items():
                if cell_type_labels and int(idx) < len(cell_type_labels):
                    label = cell_type_labels[int(idx)]
                else:
                    label = str(idx)
                sample_props[label] = prop

            proportions[sample_id] = sample_props

        except Exception as e:
            logger.warning(f"Deconvolution failed for sample {sample_id}: {e}")
            failed_samples.append(sample_id)

    # Compute summary statistics per cell type
    if proportions:
        all_ct_labels = set()
        for sp in proportions.values():
            all_ct_labels.update(sp.keys())

        summary: dict[str, dict[str, float]] = {}
        for ct in sorted(all_ct_labels):
            ct_values = [proportions[sid].get(ct, 0.0) for sid in proportions]
            summary[ct] = {
                "mean": round(sum(ct_values) / len(ct_values), 6),
                "std": round(
                    math.sqrt(
                        sum((v - sum(ct_values) / len(ct_values)) ** 2 for v in ct_values) / max(len(ct_values) - 1, 1)
                    ),
                    6,
                ),
                "min": round(min(ct_values), 6),
                "max": round(max(ct_values), 6),
            }
        n_types = len(all_ct_labels)
    else:
        summary = {}
        n_types = 0

    logger.info(
        f"Batch deconvolution complete: {len(proportions)}/{len(samples)} samples "
        f"processed ({len(failed_samples)} failed)"
    )

    return {
        "proportions": proportions,
        "summary": summary,
        "n_samples": len(proportions),
        "n_cell_types": n_types,
        "method": method,
        "failed_samples": failed_samples,
    }
