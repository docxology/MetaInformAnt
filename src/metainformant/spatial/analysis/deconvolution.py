"""Cell type deconvolution for spatial transcriptomics.

Deconvolution estimates the proportion of different cell types within each
spatial spot (especially relevant for Visium where each spot captures multiple
cells). Implements Non-Negative Least Squares (NNLS) and NMF-based approaches.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Optional dependencies
try:
    import numpy as np
    from numpy.typing import NDArray
except ImportError:
    np = None  # type: ignore[assignment]
    NDArray = None  # type: ignore[assignment,misc]

try:
    from scipy import sparse as sp_sparse
    from scipy.optimize import nnls
except ImportError:
    sp_sparse = None  # type: ignore[assignment]
    nnls = None  # type: ignore[assignment]

try:
    from sklearn.decomposition import NMF
    from sklearn.preprocessing import normalize
except ImportError:
    NMF = None  # type: ignore[assignment,misc]
    normalize = None  # type: ignore[assignment]


@dataclass
class DeconvolutionResult:
    """Result of cell type deconvolution.

    Attributes:
        weights: Raw deconvolution weights matrix (n_spots x n_types).
        fractions: Normalized cell type fractions (n_spots x n_types), rows sum to 1.
        cell_type_names: List of cell type names.
        residuals: Per-spot fitting residuals.
        method: Deconvolution method used.
        metadata: Additional result metadata.
    """

    weights: Any  # np.ndarray (n_spots, n_types)
    fractions: Any  # np.ndarray (n_spots, n_types)
    cell_type_names: list[str]
    residuals: Any  # np.ndarray (n_spots,)
    method: str
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def n_spots(self) -> int:
        """Number of spatial spots."""
        if hasattr(self.weights, "shape"):
            return int(self.weights.shape[0])
        return 0

    @property
    def n_types(self) -> int:
        """Number of cell types."""
        return len(self.cell_type_names)


def create_reference_profiles(
    scrna_data: Any,
    cell_type_labels: Any,
    *,
    gene_names: list[str] | None = None,
    method: str = "mean",
) -> tuple[Any, list[str], list[str]]:
    """Build cell type reference expression profiles from scRNA-seq data.

    For each cell type, computes the average (or median) expression profile
    across all cells of that type.

    Args:
        scrna_data: scRNA-seq expression matrix (n_cells x n_genes), dense or sparse.
        cell_type_labels: Array or list of cell type labels (length n_cells).
        gene_names: List of gene names (length n_genes). If None, uses indices.
        method: Aggregation method ("mean" or "median").

    Returns:
        Tuple of (reference_profiles, cell_type_names, gene_names) where
        reference_profiles is (n_types x n_genes).
    """
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")

    # Convert to dense
    if sp_sparse is not None and sp_sparse.issparse(scrna_data):
        data = scrna_data.toarray()
    else:
        data = np.asarray(scrna_data, dtype=np.float64)

    labels = np.asarray(cell_type_labels)
    unique_types = sorted(set(labels.tolist()))

    n_genes = data.shape[1]
    if gene_names is None:
        gene_names = [f"gene_{i}" for i in range(n_genes)]

    profiles = np.zeros((len(unique_types), n_genes), dtype=np.float64)

    for idx, ct in enumerate(unique_types):
        mask = labels == ct
        subset = data[mask, :]

        if method == "median":
            profiles[idx, :] = np.median(subset, axis=0)
        else:
            profiles[idx, :] = np.mean(subset, axis=0)

    type_names = [str(t) for t in unique_types]
    logger.info(
        f"Created reference profiles: {len(type_names)} cell types, "
        f"{n_genes} genes, method={method}"
    )
    return profiles, type_names, gene_names


def nnls_deconvolution(
    bulk_expression: Any,
    reference_signatures: Any,
) -> tuple[Any, Any]:
    """Non-negative least squares deconvolution for a single spot or batch.

    Solves: min ||bulk - reference^T * x||^2  subject to x >= 0

    For each spot, finds the non-negative weights that best reconstruct the
    observed expression from the reference cell type signatures.

    Args:
        bulk_expression: Expression vector(s). If 1D (n_genes,), single spot.
            If 2D (n_spots x n_genes), batch of spots.
        reference_signatures: Reference profiles matrix (n_types x n_genes).

    Returns:
        Tuple of (weights, residuals). weights has shape matching input:
        (n_types,) for single spot or (n_spots, n_types) for batch.
        residuals is the per-spot/vector fitting residual norm.

    Raises:
        ImportError: If scipy is not installed.
    """
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")
    if nnls is None:
        raise ImportError("scipy.optimize.nnls is required: uv pip install scipy")

    ref = np.asarray(reference_signatures, dtype=np.float64)  # (n_types, n_genes)
    bulk = np.asarray(bulk_expression, dtype=np.float64)

    if bulk.ndim == 1:
        # Single spot: solve ref^T @ x = bulk
        weights, residual = nnls(ref.T, bulk)
        return weights, residual

    # Batch of spots
    n_spots = bulk.shape[0]
    n_types = ref.shape[0]
    all_weights = np.zeros((n_spots, n_types), dtype=np.float64)
    all_residuals = np.zeros(n_spots, dtype=np.float64)

    ref_t = ref.T  # (n_genes, n_types)

    for i in range(n_spots):
        w, r = nnls(ref_t, bulk[i, :])
        all_weights[i, :] = w
        all_residuals[i] = r

    return all_weights, all_residuals


def estimate_cell_fractions(deconvolution_result: DeconvolutionResult | Any) -> Any:
    """Normalize deconvolution weights to cell type fractions (proportions summing to 1).

    Args:
        deconvolution_result: Either a DeconvolutionResult object or a raw weights
            matrix (n_spots x n_types).

    Returns:
        Normalized fractions array (n_spots x n_types) where each row sums to 1.
    """
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")

    if isinstance(deconvolution_result, DeconvolutionResult):
        weights = np.asarray(deconvolution_result.weights, dtype=np.float64)
    else:
        weights = np.asarray(deconvolution_result, dtype=np.float64)

    if weights.ndim == 1:
        total = weights.sum()
        if total > 0:
            return weights / total
        return weights

    # Per-row normalization
    row_sums = weights.sum(axis=1, keepdims=True)
    # Avoid division by zero
    row_sums[row_sums == 0] = 1.0
    fractions = weights / row_sums

    return fractions


def enrichment_score(
    observed: Any,
    expected: Any,
) -> Any:
    """Compute cell type enrichment score per spot.

    Enrichment is the log2 fold change of observed fractions vs expected
    (background) fractions, with a pseudocount for numerical stability.

    enrichment = log2((observed + epsilon) / (expected + epsilon))

    Args:
        observed: Observed cell type fractions (n_spots x n_types) or (n_types,).
        expected: Expected (background) cell type fractions (n_types,).

    Returns:
        Enrichment scores with same shape as observed.
    """
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")

    obs = np.asarray(observed, dtype=np.float64)
    exp = np.asarray(expected, dtype=np.float64)

    epsilon = 1e-10
    scores = np.log2((obs + epsilon) / (exp + epsilon))

    return scores


def deconvolve_spots(
    spatial_expression: Any,
    reference_profiles: Any,
    method: str = "nnls",
    *,
    cell_type_names: list[str] | None = None,
    gene_names: list[str] | None = None,
    spatial_gene_names: list[str] | None = None,
    reference_gene_names: list[str] | None = None,
    alpha: float = 0.0,
) -> DeconvolutionResult:
    """Deconvolve spatial spots to estimate cell type composition.

    Supports multiple deconvolution methods:
    - "nnls": Non-negative least squares (fast, robust).
    - "nmf": Non-negative matrix factorization based approach.

    When gene_names are provided for both spatial and reference data, the function
    automatically intersects to shared genes.

    Args:
        spatial_expression: Spatial expression matrix (n_spots x n_genes), dense or sparse.
        reference_profiles: Reference cell type profiles (n_types x n_genes).
        method: Deconvolution method ("nnls" or "nmf").
        cell_type_names: Names of cell types (length n_types).
        gene_names: Deprecated, use spatial_gene_names and reference_gene_names.
        spatial_gene_names: Gene names for spatial data.
        reference_gene_names: Gene names for reference profiles.
        alpha: Regularization parameter (for NMF).

    Returns:
        DeconvolutionResult with weights, fractions, and residuals.
    """
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")

    # Convert to dense
    if sp_sparse is not None and sp_sparse.issparse(spatial_expression):
        spatial = spatial_expression.toarray()
    else:
        spatial = np.asarray(spatial_expression, dtype=np.float64)

    ref = np.asarray(reference_profiles, dtype=np.float64)
    n_types = ref.shape[0]

    # Handle gene name intersection
    if spatial_gene_names is not None and reference_gene_names is not None:
        spatial_set = set(spatial_gene_names)
        ref_set = set(reference_gene_names)
        shared_genes = sorted(spatial_set & ref_set)

        if len(shared_genes) == 0:
            raise ValueError("No shared genes between spatial and reference data")

        spatial_indices = [spatial_gene_names.index(g) for g in shared_genes]
        ref_indices = [reference_gene_names.index(g) for g in shared_genes]

        spatial = spatial[:, spatial_indices]
        ref = ref[:, ref_indices]
        logger.info(f"Using {len(shared_genes)} shared genes for deconvolution")

    if cell_type_names is None:
        cell_type_names = [f"type_{i}" for i in range(n_types)]

    if method == "nnls":
        weights, residuals = nnls_deconvolution(spatial, ref)
    elif method == "nmf":
        weights, residuals = _nmf_deconvolution(spatial, ref, alpha=alpha)
    else:
        raise ValueError(f"Unknown deconvolution method: {method}. Use 'nnls' or 'nmf'.")

    fractions = estimate_cell_fractions(weights)

    result = DeconvolutionResult(
        weights=weights,
        fractions=fractions,
        cell_type_names=cell_type_names,
        residuals=residuals,
        method=method,
        metadata={
            "n_spots": spatial.shape[0],
            "n_genes_used": spatial.shape[1],
            "n_types": n_types,
            "alpha": alpha,
        },
    )

    logger.info(
        f"Deconvolution ({method}): {result.n_spots} spots, "
        f"{result.n_types} cell types, mean residual={np.mean(residuals):.4f}"
    )
    return result


def _nmf_deconvolution(
    spatial: Any,
    reference: Any,
    alpha: float = 0.0,
) -> tuple[Any, Any]:
    """NMF-based deconvolution.

    Uses the reference profiles as a fixed basis (W) and solves for the
    coefficient matrix (H) that reconstructs the spatial expression.

    spatial ~= H @ reference (where H is n_spots x n_types)

    Args:
        spatial: Expression matrix (n_spots x n_genes).
        reference: Reference profiles (n_types x n_genes).
        alpha: Regularization strength.

    Returns:
        Tuple of (weights, residuals).
    """
    if NMF is None:
        raise ImportError("scikit-learn is required for NMF deconvolution: uv pip install scikit-learn")

    n_spots = spatial.shape[0]
    n_types = reference.shape[0]

    # Ensure non-negative
    spatial_nn = np.maximum(spatial, 0)
    ref_nn = np.maximum(reference, 0)

    # For each spot, solve via NNLS against reference (more stable than NMF projection)
    # This is equivalent to NMF with fixed W
    weights = np.zeros((n_spots, n_types), dtype=np.float64)
    residuals = np.zeros(n_spots, dtype=np.float64)

    ref_t = ref_nn.T  # (n_genes, n_types)

    if alpha > 0:
        # Add L2 regularization: augment system with sqrt(alpha) * I
        n_genes = ref_t.shape[0]
        reg_block = np.sqrt(alpha) * np.eye(n_types)
        ref_aug = np.vstack([ref_t, reg_block])
    else:
        ref_aug = ref_t

    for i in range(n_spots):
        if alpha > 0:
            target = np.concatenate([spatial_nn[i, :], np.zeros(n_types)])
        else:
            target = spatial_nn[i, :]

        w, r = nnls(ref_aug, target)
        weights[i, :] = w
        # Compute actual residual (without regularization term)
        actual_resid = np.linalg.norm(spatial_nn[i, :] - ref_nn.T @ w)
        residuals[i] = actual_resid

    return weights, residuals
