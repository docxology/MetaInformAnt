"""Map scRNA-seq data to spatial transcriptomics (label transfer and imputation).

Provides methods to transfer cell type annotations from scRNA-seq reference
data to spatial spots, and to impute spatially-unmeasured genes using scRNA-seq
expression profiles.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal

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
    from scipy.stats import pearsonr, spearmanr
except ImportError:
    sp_sparse = None  # type: ignore[assignment]
    pearsonr = None  # type: ignore[assignment]
    spearmanr = None  # type: ignore[assignment]

try:
    from sklearn.neighbors import NearestNeighbors
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
except ImportError:
    NearestNeighbors = None  # type: ignore[assignment,misc]
    StandardScaler = None  # type: ignore[assignment,misc]
    PCA = None  # type: ignore[assignment,misc]


@dataclass
class MappingResult:
    """Result of scRNA-seq to spatial mapping.

    Attributes:
        predicted_labels: Predicted cell type label per spatial spot (length n_spots).
        prediction_scores: Confidence score per spot (length n_spots).
        label_probabilities: Probability matrix (n_spots x n_types) for each cell type.
        cell_type_names: List of cell type names.
        method: Mapping method used.
        metadata: Additional result metadata.
    """

    predicted_labels: list[str]
    prediction_scores: Any  # np.ndarray (n_spots,)
    label_probabilities: Any  # np.ndarray (n_spots, n_types)
    cell_type_names: list[str]
    method: str
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass
class ImputationResult:
    """Result of spatial gene imputation.

    Attributes:
        imputed_expression: Imputed expression matrix (n_spots x n_imputed_genes).
        gene_names: Names of imputed genes.
        confidence: Per-gene confidence scores.
        method: Imputation method used.
    """

    imputed_expression: Any  # np.ndarray (n_spots, n_genes)
    gene_names: list[str]
    confidence: Any  # np.ndarray (n_genes,)
    method: str


def _to_dense(matrix: Any) -> Any:
    """Convert sparse matrix to dense numpy array."""
    if sp_sparse is not None and sp_sparse.issparse(matrix):
        return matrix.toarray()
    return np.asarray(matrix, dtype=np.float64)


def _intersect_genes(
    spatial_genes: list[str],
    ref_genes: list[str],
) -> tuple[list[int], list[int], list[str]]:
    """Find shared genes between spatial and reference datasets.

    Returns:
        Tuple of (spatial_indices, ref_indices, shared_gene_names).
    """
    spatial_set = {g: i for i, g in enumerate(spatial_genes)}
    ref_set = {g: i for i, g in enumerate(ref_genes)}

    shared = sorted(set(spatial_genes) & set(ref_genes))
    spatial_idx = [spatial_set[g] for g in shared]
    ref_idx = [ref_set[g] for g in shared]

    return spatial_idx, ref_idx, shared


def correlation_mapping(
    spatial_expression: Any,
    scrna_profiles: Any,
    *,
    cell_type_names: list[str] | None = None,
    correlation_method: Literal["pearson", "spearman"] = "pearson",
) -> MappingResult:
    """Correlation-based mapping of spatial spots to cell type profiles.

    For each spatial spot, computes the correlation between its expression
    vector and each reference cell type profile, then assigns the most
    correlated type.

    Args:
        spatial_expression: Spatial expression matrix (n_spots x n_genes).
        scrna_profiles: Reference cell type profiles (n_types x n_genes).
            Each row is the mean expression for one cell type.
        cell_type_names: Names of cell types. If None, uses indices.
        correlation_method: "pearson" or "spearman".

    Returns:
        MappingResult with predicted labels and probabilities.
    """
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")
    if pearsonr is None:
        raise ImportError("scipy.stats is required: uv pip install scipy")

    spatial = _to_dense(spatial_expression)
    profiles = np.asarray(scrna_profiles, dtype=np.float64)

    n_spots = spatial.shape[0]
    n_types = profiles.shape[0]

    if cell_type_names is None:
        cell_type_names = [f"type_{i}" for i in range(n_types)]

    # Compute correlation between each spot and each profile
    corr_matrix = np.zeros((n_spots, n_types), dtype=np.float64)

    corr_fn = pearsonr if correlation_method == "pearson" else spearmanr

    for i in range(n_spots):
        spot_expr = spatial[i, :]
        for j in range(n_types):
            ref_expr = profiles[j, :]
            # Handle constant vectors
            if np.std(spot_expr) == 0 or np.std(ref_expr) == 0:
                corr_matrix[i, j] = 0.0
            else:
                r, _ = corr_fn(spot_expr, ref_expr)
                corr_matrix[i, j] = r

    # Convert correlations to probabilities via softmax
    # Shift for numerical stability
    corr_shifted = corr_matrix - corr_matrix.max(axis=1, keepdims=True)
    exp_corr = np.exp(corr_shifted * 5.0)  # temperature scaling
    probabilities = exp_corr / exp_corr.sum(axis=1, keepdims=True)

    # Assign labels
    best_idx = np.argmax(corr_matrix, axis=1)
    predicted_labels = [cell_type_names[int(idx)] for idx in best_idx]
    prediction_scores = np.array([corr_matrix[i, best_idx[i]] for i in range(n_spots)])

    logger.info(
        f"Correlation mapping ({correlation_method}): {n_spots} spots, "
        f"{n_types} types, mean score={prediction_scores.mean():.3f}"
    )

    return MappingResult(
        predicted_labels=predicted_labels,
        prediction_scores=prediction_scores,
        label_probabilities=probabilities,
        cell_type_names=cell_type_names,
        method=f"correlation_{correlation_method}",
    )


def anchor_based_transfer(
    spatial_data: Any,
    scrna_data: Any,
    anchors: Any,
    *,
    scrna_labels: Any | None = None,
    cell_type_names: list[str] | None = None,
    spatial_genes: list[str] | None = None,
    scrna_genes: list[str] | None = None,
    n_pcs: int = 30,
    k_anchor: int = 5,
    seed: int = 42,
) -> MappingResult:
    """Anchor-based label transfer from scRNA-seq to spatial data.

    Inspired by Seurat's anchor-based integration approach:
    1. Project both datasets into shared PCA space using shared genes.
    2. Find mutual nearest neighbors (MNNs) as anchors between datasets.
    3. Transfer labels through anchor-weighted voting.

    If pre-computed anchors are provided (as index pairs), uses those directly.
    Otherwise, computes MNN anchors.

    Args:
        spatial_data: Spatial expression matrix (n_spots x n_spatial_genes).
        scrna_data: scRNA-seq expression matrix (n_cells x n_scrna_genes).
        anchors: Either:
            - Array of (spatial_idx, scrna_idx) pairs, shape (n_anchors, 2).
            - None to auto-compute anchors via MNN.
        scrna_labels: Cell type labels for scRNA-seq cells (length n_cells).
            Required for label transfer.
        cell_type_names: List of possible cell types.
        spatial_genes: Gene names for spatial data.
        scrna_genes: Gene names for scRNA-seq data.
        n_pcs: Number of PCA components for embedding.
        k_anchor: Number of nearest neighbors for MNN anchor finding.
        seed: Random seed.

    Returns:
        MappingResult with transferred labels.
    """
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")
    if NearestNeighbors is None or PCA is None or StandardScaler is None:
        raise ImportError("scikit-learn is required: uv pip install scikit-learn")

    spatial = _to_dense(spatial_data)
    scrna = _to_dense(scrna_data)

    n_spots = spatial.shape[0]
    n_cells = scrna.shape[0]

    # Intersect genes if names provided
    if spatial_genes is not None and scrna_genes is not None:
        s_idx, r_idx, shared = _intersect_genes(spatial_genes, scrna_genes)
        if len(shared) == 0:
            raise ValueError("No shared genes between spatial and scRNA-seq data")
        spatial = spatial[:, s_idx]
        scrna = scrna[:, r_idx]
        logger.info(f"Using {len(shared)} shared genes for anchor-based transfer")

    # Standardize
    scaler = StandardScaler()
    combined = np.vstack([spatial, scrna])
    combined_scaled = scaler.fit_transform(combined)

    # PCA
    n_components = min(n_pcs, combined_scaled.shape[1], combined_scaled.shape[0] - 1)
    if n_components < 1:
        n_components = 1
    pca = PCA(n_components=n_components, random_state=seed)
    combined_pca = pca.fit_transform(combined_scaled)

    spatial_pca = combined_pca[:n_spots]
    scrna_pca = combined_pca[n_spots:]

    # Find or use anchors
    if anchors is None:
        # Compute MNN anchors
        k = min(k_anchor, n_cells, n_spots)

        # Spatial -> scRNA neighbors
        nn_scrna = NearestNeighbors(n_neighbors=k, metric="euclidean")
        nn_scrna.fit(scrna_pca)
        _, s_to_r = nn_scrna.kneighbors(spatial_pca)

        # scRNA -> spatial neighbors
        nn_spatial = NearestNeighbors(n_neighbors=k, metric="euclidean")
        nn_spatial.fit(spatial_pca)
        _, r_to_s = nn_spatial.kneighbors(scrna_pca)

        # Find mutual nearest neighbors
        anchor_pairs: list[tuple[int, int]] = []
        for si in range(n_spots):
            for ri_idx in range(k):
                ri = int(s_to_r[si, ri_idx])
                # Check if spatial spot is also a neighbor of this scRNA cell
                if si in r_to_s[ri, :]:
                    anchor_pairs.append((si, ri))

        anchors = np.array(anchor_pairs, dtype=np.int32) if anchor_pairs else np.empty((0, 2), dtype=np.int32)
        logger.info(f"Found {len(anchor_pairs)} MNN anchors")
    else:
        anchors = np.asarray(anchors, dtype=np.int32)

    if scrna_labels is None:
        raise ValueError("scrna_labels must be provided for label transfer")

    labels = np.asarray(scrna_labels)
    if cell_type_names is None:
        cell_type_names = sorted(set(labels.tolist()))

    type_to_idx = {t: i for i, t in enumerate(cell_type_names)}
    n_types = len(cell_type_names)

    # Transfer labels through anchors
    # For each spatial spot, find its anchors and vote-transfer labels
    # Weight anchors by distance in PCA space
    label_probs = np.zeros((n_spots, n_types), dtype=np.float64)

    if len(anchors) > 0:
        # Build spatial -> anchor mapping
        spot_anchors: dict[int, list[tuple[int, float]]] = {}
        for ai in range(anchors.shape[0]):
            si, ri = int(anchors[ai, 0]), int(anchors[ai, 1])
            d = float(np.linalg.norm(spatial_pca[si] - scrna_pca[ri]))
            weight = 1.0 / (1.0 + d)
            if si not in spot_anchors:
                spot_anchors[si] = []
            spot_anchors[si].append((ri, weight))

        for si in range(n_spots):
            if si in spot_anchors:
                for ri, weight in spot_anchors[si]:
                    ct = labels[ri]
                    if ct in type_to_idx:
                        label_probs[si, type_to_idx[ct]] += weight
            else:
                # No direct anchors: use KNN in PCA space to scRNA-seq
                nn = NearestNeighbors(n_neighbors=min(5, n_cells), metric="euclidean")
                nn.fit(scrna_pca)
                dists, indices = nn.kneighbors(spatial_pca[si : si + 1])
                for j_idx in range(indices.shape[1]):
                    ri = int(indices[0, j_idx])
                    d = float(dists[0, j_idx])
                    weight = 1.0 / (1.0 + d)
                    ct = labels[ri]
                    if ct in type_to_idx:
                        label_probs[si, type_to_idx[ct]] += weight
    else:
        # Fallback: direct KNN transfer
        nn = NearestNeighbors(n_neighbors=min(k_anchor, n_cells), metric="euclidean")
        nn.fit(scrna_pca)
        dists, indices = nn.kneighbors(spatial_pca)
        for si in range(n_spots):
            for j_idx in range(indices.shape[1]):
                ri = int(indices[si, j_idx])
                d = float(dists[si, j_idx])
                weight = 1.0 / (1.0 + d)
                ct = labels[ri]
                if ct in type_to_idx:
                    label_probs[si, type_to_idx[ct]] += weight

    # Normalize probabilities
    row_sums = label_probs.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    label_probs = label_probs / row_sums

    # Assign labels
    best_idx = np.argmax(label_probs, axis=1)
    predicted_labels = [cell_type_names[int(idx)] for idx in best_idx]
    prediction_scores = np.array([label_probs[i, best_idx[i]] for i in range(n_spots)])

    logger.info(
        f"Anchor-based transfer: {n_spots} spots, {len(anchors)} anchors, "
        f"mean confidence={prediction_scores.mean():.3f}"
    )

    return MappingResult(
        predicted_labels=predicted_labels,
        prediction_scores=prediction_scores,
        label_probabilities=label_probs,
        cell_type_names=cell_type_names,
        method="anchor_based",
        metadata={
            "n_anchors": len(anchors),
            "n_pcs": n_components,
            "k_anchor": k_anchor,
        },
    )


def map_scrna_to_spatial(
    scrna_data: Any,
    spatial_data: Any,
    method: Literal["correlation", "anchor"] = "correlation",
    *,
    scrna_labels: Any | None = None,
    scrna_genes: list[str] | None = None,
    spatial_genes: list[str] | None = None,
    cell_type_names: list[str] | None = None,
    **kwargs: Any,
) -> MappingResult:
    """Map scRNA-seq cell type annotations to spatial spots.

    Unified interface for label transfer from scRNA-seq reference to
    spatial transcriptomics data.

    Args:
        scrna_data: scRNA-seq expression matrix (n_cells x n_genes).
        spatial_data: Spatial expression matrix (n_spots x n_genes).
        method: Mapping method ("correlation" or "anchor").
        scrna_labels: Cell type labels for scRNA-seq cells.
        scrna_genes: Gene names for scRNA-seq data.
        spatial_genes: Gene names for spatial data.
        cell_type_names: List of cell type names.
        **kwargs: Additional arguments passed to the specific method.

    Returns:
        MappingResult with predicted spatial labels.
    """
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")

    scrna = _to_dense(scrna_data)
    spatial = _to_dense(spatial_data)

    # Intersect genes
    if scrna_genes is not None and spatial_genes is not None:
        s_idx, r_idx, shared = _intersect_genes(spatial_genes, scrna_genes)
        spatial = spatial[:, s_idx]
        scrna = scrna[:, r_idx]
        logger.info(f"Using {len(shared)} shared genes")

    if method == "correlation":
        if scrna_labels is None:
            raise ValueError("scrna_labels required for correlation mapping")

        # Build profiles from scRNA-seq data
        from metainformant.spatial.analysis.deconvolution import create_reference_profiles

        labels = np.asarray(scrna_labels)
        profiles, type_names, _ = create_reference_profiles(scrna, labels)

        if cell_type_names is None:
            cell_type_names = type_names

        return correlation_mapping(
            spatial,
            profiles,
            cell_type_names=cell_type_names,
            correlation_method=kwargs.get("correlation_method", "pearson"),
        )

    elif method == "anchor":
        return anchor_based_transfer(
            spatial,
            scrna,
            anchors=kwargs.get("anchors", None),
            scrna_labels=scrna_labels,
            cell_type_names=cell_type_names,
            n_pcs=kwargs.get("n_pcs", 30),
            k_anchor=kwargs.get("k_anchor", 5),
            seed=kwargs.get("seed", 42),
        )

    else:
        raise ValueError(f"Unknown mapping method: {method}. Use 'correlation' or 'anchor'.")


def impute_spatial_genes(
    spatial_data: Any,
    scrna_data: Any,
    genes: list[str],
    *,
    spatial_genes: list[str] | None = None,
    scrna_genes: list[str] | None = None,
    n_neighbors: int = 10,
    n_pcs: int = 30,
    seed: int = 42,
) -> ImputationResult:
    """Impute unmeasured genes in spatial data using scRNA-seq reference.

    For genes not present in the spatial panel but measured in scRNA-seq:
    1. Project both datasets into shared PCA space (using overlapping genes).
    2. For each spatial spot, find K nearest scRNA-seq cells in PCA space.
    3. Impute target gene expression as weighted average of neighbor values.

    Args:
        spatial_data: Spatial expression matrix (n_spots x n_spatial_genes).
        scrna_data: scRNA-seq expression matrix (n_cells x n_scrna_genes).
        genes: List of gene names to impute (must be present in scrna_genes).
        spatial_genes: Gene names for spatial data.
        scrna_genes: Gene names for scRNA-seq data.
        n_neighbors: Number of scRNA-seq neighbors for imputation.
        n_pcs: Number of PCA components.
        seed: Random seed.

    Returns:
        ImputationResult with imputed expression for requested genes.
    """
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")
    if NearestNeighbors is None or PCA is None or StandardScaler is None:
        raise ImportError("scikit-learn is required: uv pip install scikit-learn")

    spatial = _to_dense(spatial_data)
    scrna = _to_dense(scrna_data)

    n_spots = spatial.shape[0]
    n_cells = scrna.shape[0]

    if spatial_genes is None or scrna_genes is None:
        raise ValueError("Both spatial_genes and scrna_genes must be provided for imputation")

    # Find shared genes for PCA projection
    s_idx, r_idx, shared = _intersect_genes(spatial_genes, scrna_genes)
    if len(shared) == 0:
        raise ValueError("No shared genes between spatial and scRNA-seq data")

    spatial_shared = spatial[:, s_idx]
    scrna_shared = scrna[:, r_idx]

    # Find target gene indices in scRNA-seq data
    scrna_gene_map = {g: i for i, g in enumerate(scrna_genes)}
    target_indices: list[int] = []
    target_names: list[str] = []
    for g in genes:
        if g in scrna_gene_map:
            target_indices.append(scrna_gene_map[g])
            target_names.append(g)
        else:
            logger.warning(f"Gene '{g}' not found in scRNA-seq data; skipping")

    if len(target_indices) == 0:
        raise ValueError("None of the requested genes found in scRNA-seq data")

    # Joint PCA on shared genes
    combined = np.vstack([spatial_shared, scrna_shared])
    scaler = StandardScaler()
    combined_scaled = scaler.fit_transform(combined)

    n_components = min(n_pcs, combined_scaled.shape[1], combined_scaled.shape[0] - 1)
    if n_components < 1:
        n_components = 1
    pca = PCA(n_components=n_components, random_state=seed)
    combined_pca = pca.fit_transform(combined_scaled)

    spatial_pca = combined_pca[:n_spots]
    scrna_pca = combined_pca[n_spots:]

    # Find KNN in PCA space
    k = min(n_neighbors, n_cells)
    nn = NearestNeighbors(n_neighbors=k, metric="euclidean")
    nn.fit(scrna_pca)
    distances, indices = nn.kneighbors(spatial_pca)

    # Impute target genes via weighted averaging
    n_targets = len(target_indices)
    imputed = np.zeros((n_spots, n_targets), dtype=np.float64)
    confidence = np.zeros(n_targets, dtype=np.float64)

    # Inverse distance weights
    weights = 1.0 / (distances + 1e-10)
    weights = weights / weights.sum(axis=1, keepdims=True)

    scrna_targets = scrna[:, target_indices]  # (n_cells, n_targets)

    for si in range(n_spots):
        nbr_idx = indices[si, :]
        nbr_weights = weights[si, :]
        imputed[si, :] = nbr_weights @ scrna_targets[nbr_idx, :]

    # Confidence: average correlation of imputed values across neighbors
    # Higher agreement among neighbors = higher confidence
    for gi in range(n_targets):
        neighbor_values = scrna_targets[indices, gi]  # (n_spots, k)
        # Coefficient of variation across neighbors
        neighbor_std = neighbor_values.std(axis=1)
        neighbor_mean = np.abs(neighbor_values.mean(axis=1)) + 1e-10
        cv = neighbor_std / neighbor_mean
        # Lower CV = higher confidence; map to [0, 1]
        confidence[gi] = float(1.0 / (1.0 + np.mean(cv)))

    logger.info(f"Imputed {n_targets} genes for {n_spots} spots, " f"mean confidence={confidence.mean():.3f}")

    return ImputationResult(
        imputed_expression=imputed,
        gene_names=target_names,
        confidence=confidence,
        method="knn_weighted",
    )
