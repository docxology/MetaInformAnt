"""Single-cell data integration and batch correction methods.

This module provides functions for integrating multiple single-cell datasets
by correcting for batch effects and technical variation. Methods include
mutual nearest neighbors (MNN), BBKNN, Harmony, and Scanorama-like approaches.

Only :func:`bbknn_integration` follows its published algorithm (a batch-balanced k-nearest-neighbor graph).
:func:`harmony_integration`, :func:`scanorama_integration`, :func:`mnn_integration` and
:func:`combat_integration` are simplified location-scale or nearest-neighbor corrections named after the
method they approximate; their docstrings say what they actually compute. None of them are ports of the
reference implementations, so their output is not comparable to those packages'.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from metainformant.core.data import validation
from metainformant.core.utils import errors, logging
from metainformant.singlecell.data.preprocessing import SingleCellData

# Optional scientific dependencies
try:
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False
    StandardScaler = None
    PCA = None

logger = logging.get_logger(__name__)


def _to_dense(x: Any) -> Any:
    """Return a dense view/copy of ``x`` (sparse matrices are converted)."""
    return x.toarray() if hasattr(x, "toarray") else x


def _require_same_genes(data_list: List[SingleCellData]) -> None:
    """Refuse datasets whose genes differ: the integrations stack the matrices column by column.

    Raises:
        ValidationError: If a dataset has another number of genes than the first, or other gene names
    """
    first = data_list[0]
    for i, data in enumerate(data_list[1:], 1):
        if data.n_vars != first.n_vars:
            raise errors.ValidationError(
                f"data_list[{i}] has {data.n_vars} genes but data_list[0] has {first.n_vars}: "
                "the datasets must have the same genes in the same order"
            )
        if first.var is not None and data.var is not None and not data.var.index.equals(first.var.index):
            raise errors.ValidationError(
                f"data_list[{i}] has other gene names than data_list[0]: "
                "the datasets must have the same genes in the same order"
            )


def bbknn_integration(data: SingleCellData, batch_key: str, n_neighbors: int = 15) -> SingleCellData:
    """Build a batch-balanced k-nearest-neighbor graph (BBKNN, Polanski et al. 2020).

    The cells are embedded in a PCA space (of the standardized expression matrix, at most 50 components) and every
    cell is connected to its ``n_neighbors`` nearest cells *in each batch* separately, its own included in its own
    batch, so that a cell's neighborhood draws equally on every batch instead of on the batch it came from. The
    expression matrix is not changed: BBKNN corrects the neighbor graph, not the data.

    Args:
        data: SingleCellData object with batch labels
        batch_key: Column name containing batch labels
        n_neighbors: Number of neighbors per cell in each batch (fewer for a batch with fewer cells)

    Returns:
        SingleCellData with the graph in ``uns["bbknn"]`` (``indices`` and ``distances``, one row per cell and the
        neighbors of every batch side by side, ``connectivities`` as a symmetric sparse matrix, ``batch_key``,
        ``n_batches``, ``n_neighbors`` and ``batch_adjacency_shape``) and in ``uns["neighbors"]`` (the same
        arrays with ``params``), where the diffusion map, DPT and PAGA read it. ``obs["bbknn_1"]`` and
        ``obs["bbknn_2"]`` are a t-SNE of the PCA space, for plotting only; that embedding is not batch corrected.

    Raises:
        ValidationError: If data is not SingleCellData or batch_key is not in data.obs
        ImportError: If scikit-learn not available
    """
    if not HAS_SKLEARN:
        raise ImportError(
            "scikit-learn is required for BBKNN integration. " "Install with: uv pip install scikit-learn"
        )

    validation.validate_type(data, SingleCellData, "data")

    if data.obs is None or batch_key not in data.obs.columns:
        raise errors.ValidationError(f"Batch key '{batch_key}' not found in data.obs")

    logger.info(f"Performing BBKNN integration using {batch_key} for batches")

    # Create copy
    result = data.copy()

    # Get expression matrix and batch labels
    X = _to_dense(data.X)
    batch_labels = data.obs[batch_key].values
    unique_batches = np.unique(batch_labels)

    logger.info(f"Found {len(unique_batches)} batches: {unique_batches}")

    # PCA embedding in which the neighbors are searched
    pca = PCA(n_components=min(50, X.shape[1]), random_state=42)
    X_pca = pca.fit_transform(StandardScaler().fit_transform(X))

    n_cells = X.shape[0]

    # Batch-balanced graph: the nearest cells of every batch, found in that batch alone
    from scipy import sparse
    from sklearn.neighbors import NearestNeighbors

    index_blocks: List[np.ndarray] = []
    distance_blocks: List[np.ndarray] = []
    for batch in unique_batches:
        members = np.flatnonzero(batch_labels == batch)
        k_batch = min(n_neighbors, len(members))
        batch_distances, batch_indices = NearestNeighbors(n_neighbors=k_batch).fit(X_pca[members]).kneighbors(X_pca)
        index_blocks.append(members[batch_indices])
        distance_blocks.append(batch_distances)
    indices = np.hstack(index_blocks)
    distances = np.hstack(distance_blocks)

    # Gaussian-kernel connectivities, symmetrized by averaging, as compute_neighbors does
    positive = distances[distances > 0]
    sigma = float(np.median(positive)) if positive.size else 1.0
    weights = np.exp(-(distances.reshape(-1) ** 2) / (2 * sigma**2))
    rows = np.repeat(np.arange(n_cells), indices.shape[1])
    directed = sparse.csr_matrix((weights, (rows, indices.reshape(-1))), shape=(n_cells, n_cells))
    connectivities = sparse.csr_matrix((directed + directed.T) / 2)

    result.uns["bbknn"] = {
        "batch_key": batch_key,
        "n_neighbors": n_neighbors,
        "n_batches": len(unique_batches),
        "indices": indices,
        "distances": distances,
        "connectivities": connectivities,
        "batch_adjacency_shape": connectivities.shape,
    }
    result.uns["neighbors"] = {
        "indices": indices,
        "distances": distances,
        "connectivities": connectivities,
        "params": {
            "n_neighbors": n_neighbors,
            "metric": "euclidean",
            "use_rep": "X_pca",
            "method": "bbknn",
            "batch_key": batch_key,
        },
    }

    # For plotting, a 2D embedding of the (uncorrected) PCA space; it is not part of the graph
    from sklearn.manifold import TSNE

    result.obs = result.obs.copy() if result.obs is not None else pd.DataFrame()
    if n_cells >= 4:
        perplexity = min(30.0, (n_cells - 1) / 3.0)
        X_embed = TSNE(n_components=2, random_state=42, perplexity=perplexity).fit_transform(X_pca)
        result.obs["bbknn_1"] = X_embed[:, 0]
        result.obs["bbknn_2"] = X_embed[:, 1]
    else:
        logger.warning("Fewer than 4 cells: no t-SNE embedding for plotting is added")

    logger.info("BBKNN integration completed")
    return result


def harmony_integration(data: SingleCellData, batch_key: str, n_components: int = 50) -> SingleCellData:
    """Correct batch effects in a PCA space by centering each batch (a simplified stand-in for Harmony).

    The standardized expression matrix is embedded with PCA (at most 50 components), the mean of each batch is
    subtracted from its cells, and the centered coordinates are re-embedded with PCA. Real Harmony (Korsunsky et al.
    2019) iterates soft clustering and per-cluster linear corrections and can mix batches that differ in composition;
    this one only aligns batch centroids, so it is not Harmony's output.

    Args:
        data: SingleCellData object with batch labels
        batch_key: Column name containing batch labels
        n_components: Number of components of the re-embedding

    Returns:
        SingleCellData with batch-centered PCA coordinates

    Raises:
        ValidationError: If data is not SingleCellData or batch_key is not in data.obs
        ImportError: If scikit-learn not available
    """
    if not HAS_SKLEARN:
        raise ImportError(
            "scikit-learn is required for Harmony integration. " "Install with: uv pip install scikit-learn"
        )

    validation.validate_type(data, SingleCellData, "data")

    if data.obs is None or batch_key not in data.obs.columns:
        raise errors.ValidationError(f"Batch key '{batch_key}' not found in data.obs")

    logger.info(f"Performing Harmony integration using {batch_key} for batches")

    # Create copy
    result = data.copy()

    # Get expression matrix and batch labels
    X = _to_dense(data.X)
    batch_labels = data.obs[batch_key].values

    # Harmony-like batch correction (simplified implementation)
    # Real Harmony uses a more sophisticated approach

    # Step 1: PCA on original data
    pca = PCA(n_components=min(50, X.shape[1]), random_state=42)
    X_pca = pca.fit_transform(StandardScaler().fit_transform(X))

    # Step 2: Regress out batch effects
    unique_batches = np.unique(batch_labels)
    X_corrected = X_pca.copy()

    # Simple batch effect correction: center each batch
    for batch in unique_batches:
        batch_mask = batch_labels == batch
        if np.sum(batch_mask) > 0:
            batch_mean = np.mean(X_pca[batch_mask], axis=0)
            X_corrected[batch_mask] -= batch_mean

    # Step 3: Re-embed the corrected data
    pca_corrected = PCA(n_components=n_components, random_state=42)
    X_harmony = pca_corrected.fit_transform(X_corrected)

    # Store Harmony results
    result.obs = result.obs.copy() if result.obs is not None else pd.DataFrame()
    for i in range(min(n_components, 10)):  # Store first 10 components
        result.obs[f"harmony_{i+1}"] = X_harmony[:, i]

    result.uns["harmony"] = {
        "batch_key": batch_key,
        "n_components": n_components,
        "n_batches": len(unique_batches),
        "explained_variance_ratio": pca.explained_variance_ratio_[:n_components].tolist(),
    }

    logger.info(f"Harmony integration completed: {n_components} components")
    return result


def scanorama_integration(data_list: List[SingleCellData], batch_key: str) -> SingleCellData:
    """Concatenate datasets, standardize each one per gene and embed them with PCA (Scanorama-like).

    This is a simplified integration, not Scanorama (Hie et al. 2019): there is no panorama stitching by mutual
    nearest neighbors of the datasets' embeddings. The datasets must have the same genes in the same order: nothing
    is aligned by gene name, so datasets that disagree are refused.

    Args:
        data_list: List of SingleCellData objects to integrate
        batch_key: Column name that will contain batch labels

    Returns:
        Integrated SingleCellData object (its columns are principal components)

    Raises:
        ValidationError: If data_list is empty, contains an object that is not SingleCellData, or the datasets do not
            have the same genes in the same order
        ImportError: If scikit-learn not available
    """
    if not HAS_SKLEARN:
        raise ImportError(
            "scikit-learn is required for Scanorama integration. " "Install with: uv pip install scikit-learn"
        )

    if not data_list:
        raise errors.ValidationError("data_list cannot be empty")

    for i, data in enumerate(data_list):
        if not isinstance(data, SingleCellData):
            raise errors.ValidationError(f"data_list[{i}] is not a SingleCellData object")
    _require_same_genes(data_list)

    logger.info(f"Performing Scanorama integration of {len(data_list)} datasets")

    # Concatenate all datasets
    all_X = []
    all_obs = []
    batch_labels = []

    for i, data in enumerate(data_list):
        X = _to_dense(data.X)
        all_X.append(X)

        obs = data.obs.copy() if data.obs is not None else pd.DataFrame(index=range(data.n_obs))
        all_obs.append(obs)

        batch_labels.extend([f"batch_{i}"] * data.n_obs)

    # Combine matrices
    X_combined = np.vstack(all_X)
    obs_combined = pd.concat(all_obs, ignore_index=True)
    obs_combined[batch_key] = batch_labels

    # Scanorama-like correction (simplified)
    # Real Scanorama uses linear corrections in a shared low-dimensional space

    # Step 1: Find shared highly variable genes (simplified)
    # For this implementation, assume all genes are shared

    # Step 2: Compute dataset-specific corrections
    unique_batches = np.unique(batch_labels)
    X_corrected = X_combined.copy()

    for batch in unique_batches:
        batch_mask = np.array(batch_labels) == batch
        if np.sum(batch_mask) > 1:  # Need at least 2 cells
            batch_data = X_combined[batch_mask]

            # Simple correction: center and scale within batch
            batch_mean = np.mean(batch_data, axis=0)
            batch_std = np.std(batch_data, axis=0)
            batch_std = np.where(batch_std == 0, 1, batch_std)  # Avoid division by zero

            X_corrected[batch_mask] = (batch_data - batch_mean) / batch_std

    # Step 3: Joint embedding
    pca = PCA(n_components=min(50, X_corrected.shape[1]), random_state=42)
    X_integrated = pca.fit_transform(X_corrected)

    # Create integrated SingleCellData
    integrated_data = SingleCellData(
        X=X_integrated, obs=obs_combined, var=pd.DataFrame(index=[f"PC{i+1}" for i in range(X_integrated.shape[1])])
    )

    # Add integration metadata
    integrated_data.uns["scanorama"] = {
        "n_datasets": len(data_list),
        "batch_key": batch_key,
        "total_cells": X_combined.shape[0],
        "total_genes": X_combined.shape[1],
        "explained_variance_ratio": pca.explained_variance_ratio_.tolist(),
    }

    logger.info("Scanorama integration completed")
    return integrated_data


def mnn_integration(data_list: List[SingleCellData], batch_key: str) -> SingleCellData:
    """Integrate datasets pairwise with a simplified mutual-nearest-neighbor correction.

    Mutual nearest neighbors between each new dataset and the integrated data define a correction vector that is
    applied to the new dataset; it follows Haghverdi et al. (2018) in outline, not in detail (no cosine
    normalization, no smoothing of the correction vectors), so it is not their output.

    Args:
        data_list: List of SingleCellData objects to integrate
        batch_key: Column name that will contain batch labels

    Returns:
        Integrated SingleCellData object

    Raises:
        ValidationError: If data_list has fewer than 2 datasets, contains an object that is not SingleCellData,
            the datasets do not have the same genes in the same order, or a dataset has no mutual nearest neighbor
            among the datasets before it (it would otherwise be left out of the result)
        ImportError: If scikit-learn not available
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn is required for MNN integration. " "Install with: uv pip install scikit-learn")

    if len(data_list) < 2:
        raise errors.ValidationError("MNN integration requires at least 2 datasets")

    for i, data in enumerate(data_list):
        if not isinstance(data, SingleCellData):
            raise errors.ValidationError(f"data_list[{i}] is not a SingleCellData object")
    _require_same_genes(data_list)

    logger.info(f"Performing MNN integration of {len(data_list)} datasets")

    # Start with the first dataset
    integrated_data = data_list[0].copy()
    integrated_data.obs = integrated_data.obs.copy() if integrated_data.obs is not None else pd.DataFrame()
    integrated_data.obs[batch_key] = "batch_0"

    # Iteratively integrate additional datasets
    for i, data in enumerate(data_list[1:], 1):
        logger.info(f"Integrating dataset {i}")

        # Get expression matrices
        X1 = _to_dense(integrated_data.X)
        X2 = _to_dense(data.X)

        # Find mutual nearest neighbors
        mnn_pairs = _find_mutual_nearest_neighbors(X1, X2, k=20)

        if len(mnn_pairs) == 0:
            # Skipping the dataset would drop its cells from the result without a trace
            raise errors.ValidationError(
                f"No mutual nearest neighbors between dataset {i} and the datasets before it: the batches cannot be "
                "aligned, and leaving the dataset out would lose its cells"
            )

        # Correct batch effects using MNN pairs
        X1_corrected, X2_corrected = _correct_batch_effects_mnn(X1, X2, mnn_pairs)

        # Update integrated data
        integrated_data.X = X1_corrected

        # Add the new dataset
        obs2 = data.obs.copy() if data.obs is not None else pd.DataFrame(index=range(data.n_obs))
        obs2[batch_key] = f"batch_{i}"

        # Concatenate observations
        integrated_data.obs = pd.concat([integrated_data.obs, obs2], ignore_index=True)

        # Concatenate expression matrices
        integrated_data.X = np.vstack([integrated_data.X, X2_corrected])

    # Final dimensionality reduction for integrated data
    X_integrated = integrated_data.X
    pca = PCA(n_components=min(50, X_integrated.shape[1]), random_state=42)
    X_pca = pca.fit_transform(X_integrated)

    # Update the integrated data with PCA coordinates
    integrated_data.X = X_pca
    integrated_data.var = pd.DataFrame(index=[f"PC{i+1}" for i in range(X_pca.shape[1])])

    # Add integration metadata
    integrated_data.uns["mnn"] = {
        "n_datasets": len(data_list),
        "batch_key": batch_key,
        "total_cells": X_integrated.shape[0],
        "explained_variance_ratio": pca.explained_variance_ratio_.tolist(),
    }

    logger.info("MNN integration completed")
    return integrated_data


def combat_integration(data: SingleCellData, batch_key: str, covariates: Optional[List[str]] = None) -> SingleCellData:
    """Correct batch effects gene by gene with a location-scale adjustment in the spirit of ComBat.

    Every gene of each non-reference batch is standardized with its own mean and variance and rescaled to the mean
    and variance of the first batch. Real ComBat (Johnson et al. 2007) shrinks the batch estimates toward a common
    prior by empirical Bayes and can protect covariates; this function does neither, so it is not ComBat's output.
    A batch with a single cell cannot be estimated and is left uncorrected.

    Args:
        data: SingleCellData object with batch labels
        batch_key: Column name containing batch labels
        covariates: Not supported. Covariates cannot be preserved by this adjustment, so a non-empty list is
            refused rather than accepted and ignored.

    Returns:
        Batch-corrected SingleCellData

    Raises:
        ValidationError: If data is not SingleCellData, batch_key is not in data.obs, or the first batch has fewer
            than 2 cells
        NotImplementedError: If ``covariates`` is not empty
    """
    validation.validate_type(data, SingleCellData, "data")

    if covariates:
        raise NotImplementedError(
            "combat_integration() cannot preserve covariates: the location-scale adjustment does not model them, "
            f"so {list(covariates)} would be silently ignored. Run ComBat itself (for example scanpy.pp.combat with "
            "covariates) when covariates matter."
        )

    if data.obs is None or batch_key not in data.obs.columns:
        raise errors.ValidationError(f"Batch key '{batch_key}' not found in data.obs")

    logger.info(f"Performing ComBat integration using {batch_key} for batches")

    # Create copy
    result = data.copy()

    # Get expression matrix and batch labels
    X = _to_dense(data.X)
    batch_labels = data.obs[batch_key].values

    # Simplified ComBat-like correction
    # Real ComBat uses empirical Bayes for parameter estimation

    unique_batches = np.unique(batch_labels)
    X_corrected = X.copy()

    # For each gene, estimate batch effects and correct
    for gene_idx in range(X.shape[1]):
        gene_expr = X[:, gene_idx]

        # Estimate batch-specific parameters
        batch_means = {}
        batch_vars = {}

        for batch in unique_batches:
            batch_mask = batch_labels == batch
            if np.sum(batch_mask) > 1:
                batch_expr = gene_expr[batch_mask]
                batch_means[batch] = np.mean(batch_expr)
                batch_vars[batch] = np.var(batch_expr)

        # Use reference batch (first batch) as baseline
        ref_batch = unique_batches[0]
        if ref_batch not in batch_means:
            raise errors.ValidationError(f"The reference batch {ref_batch!r} has fewer than 2 cells")
        ref_mean = batch_means[ref_batch]
        ref_var = batch_vars[ref_batch]

        # Correct each batch to reference (a batch with one cell has no estimate and is left as it is)
        for batch in unique_batches:
            if batch != ref_batch and batch in batch_vars:
                batch_mask = batch_labels == batch
                batch_expr = gene_expr[batch_mask]

                # Simple standardization and recentering
                if batch_vars[batch] > 0:
                    standardized = (batch_expr - batch_means[batch]) / np.sqrt(batch_vars[batch])
                    corrected = standardized * np.sqrt(ref_var) + ref_mean
                    X_corrected[batch_mask, gene_idx] = corrected

    result.X = X_corrected

    # Store ComBat metadata
    result.uns["combat"] = {
        "batch_key": batch_key,
        "n_batches": len(unique_batches),
        "covariates": [],
        "method": "location_scale",
    }

    logger.info("ComBat integration completed")
    return result


def _find_mutual_nearest_neighbors(X1: np.ndarray, X2: np.ndarray, k: int = 20) -> List[Tuple[int, int]]:
    """Find mutual nearest neighbors between two datasets."""
    from sklearn.neighbors import NearestNeighbors

    # Find kNN from X1 to X2
    nbrs1to2 = NearestNeighbors(n_neighbors=k, metric="euclidean")
    nbrs1to2.fit(X2)
    distances1to2, indices1to2 = nbrs1to2.kneighbors(X1)

    # Find kNN from X2 to X1
    nbrs2to1 = NearestNeighbors(n_neighbors=k, metric="euclidean")
    nbrs2to1.fit(X1)
    distances2to1, indices2to1 = nbrs2to1.kneighbors(X2)

    # Find mutual pairs
    mnn_pairs = []

    for i in range(X1.shape[0]):
        for j in indices1to2[i]:
            # Check if i is in j's nearest neighbors
            if i in indices2to1[j]:
                mnn_pairs.append((i, j))

    return mnn_pairs


def _correct_batch_effects_mnn(
    X1: np.ndarray, X2: np.ndarray, mnn_pairs: List[Tuple[int, int]]
) -> Tuple[np.ndarray, np.ndarray]:
    """Correct batch effects using MNN pairs."""
    if not mnn_pairs:
        return X1, X2

    X1_corrected = X1.copy()
    X2_corrected = X2.copy()

    # For each MNN pair, compute correction vector
    for i, j in mnn_pairs:
        correction_vector = X1[i] - X2[j]

        # Apply correction to nearby cells in X2
        # Simplified: just correct this pair
        X2_corrected[j] = X2[j] + correction_vector * 0.5

    return X1_corrected, X2_corrected


def integrate_multiple_batches(
    data_list: List[SingleCellData], integration_method: str = "scanorama", **kwargs: Any
) -> SingleCellData:
    """Integrate multiple single-cell datasets using specified method.

    Args:
        data_list: List of SingleCellData objects to integrate
        integration_method: Integration method ("scanorama", "mnn", "harmony", "bbknn")
        **kwargs: Additional arguments for the integration method

    Returns:
        Integrated SingleCellData object

    Raises:
        ValidationError: If integration_method is unsupported, or the method rejects data_list (for example when it is
            empty)
    """
    valid_methods = ["scanorama", "mnn", "harmony", "bbknn"]
    if integration_method not in valid_methods:
        raise errors.ValidationError(f"Unsupported integration method: {integration_method}")

    logger.info(f"Integrating {len(data_list)} datasets using {integration_method}")

    if integration_method == "scanorama":
        return scanorama_integration(data_list, kwargs.pop("batch_key", "batch"), **kwargs)
    elif integration_method == "mnn":
        return mnn_integration(data_list, kwargs.pop("batch_key", "batch"), **kwargs)
    elif integration_method == "harmony":
        # For harmony, we need to concatenate first
        if not data_list:
            raise errors.ValidationError("No datasets provided")

        # Concatenate all datasets
        combined_data = scanorama_integration(data_list, batch_key="batch")
        return harmony_integration(combined_data, batch_key="batch", **kwargs)
    elif integration_method == "bbknn":
        # For BBKNN, we need to concatenate first
        if not data_list:
            raise errors.ValidationError("No datasets provided")

        combined_data = scanorama_integration(data_list, batch_key="batch")
        return bbknn_integration(combined_data, batch_key="batch", **kwargs)
    else:  # pragma: no cover - guarded by the valid_methods check above
        raise errors.ValidationError(f"Unsupported integration method: {integration_method}")


def evaluate_integration_quality(
    integrated_data: SingleCellData, batch_key: str, cluster_key: Optional[str] = None
) -> Dict[str, Any]:
    """Evaluate the quality of data integration.

    Args:
        integrated_data: Integrated SingleCellData object
        batch_key: Column name containing batch labels
        cluster_key: Optional column name containing cluster labels

    Returns:
        Dictionary with integration quality metrics

    Raises:
        ValidationError: If integrated_data is not SingleCellData or batch_key is not in integrated_data.obs
    """
    validation.validate_type(integrated_data, SingleCellData, "integrated_data")

    if integrated_data.obs is None or batch_key not in integrated_data.obs.columns:
        raise errors.ValidationError(f"Batch key '{batch_key}' not found in integrated_data.obs")

    logger.info("Evaluating integration quality")

    batch_labels = integrated_data.obs[batch_key].values
    X = _to_dense(integrated_data.X)

    unique_batches = np.unique(batch_labels)
    n_batches = len(unique_batches)

    metrics = {
        "n_batches": n_batches,
        "batch_sizes": [np.sum(batch_labels == batch) for batch in unique_batches],
    }

    # Batch mixing score (simplified)
    # Compute kNN and see batch composition
    from sklearn.neighbors import NearestNeighbors

    nbrs = NearestNeighbors(n_neighbors=min(50, X.shape[0] - 1))
    nbrs.fit(X)
    distances, indices = nbrs.kneighbors(X)

    batch_mixing_scores = []

    for i in range(X.shape[0]):
        neighbors = indices[i]
        neighbor_batches = batch_labels[neighbors]
        own_batch = batch_labels[i]

        # Fraction of neighbors from different batches
        different_batch_fraction = np.mean(neighbor_batches != own_batch)
        batch_mixing_scores.append(different_batch_fraction)

    metrics["mean_batch_mixing_score"] = np.mean(batch_mixing_scores)

    # Within-batch vs between-batch distances
    within_batch_distances = []
    between_batch_distances = []

    for i in range(min(1000, X.shape[0])):  # Sample for performance
        for j in range(i + 1, min(i + 11, X.shape[0])):  # Sample pairs
            dist = distances[i, j - i - 1] if j - i - 1 < distances.shape[1] else 0
            if batch_labels[i] == batch_labels[j]:
                within_batch_distances.append(dist)
            else:
                between_batch_distances.append(dist)

    if within_batch_distances and between_batch_distances:
        metrics["mean_within_batch_distance"] = np.mean(within_batch_distances)
        metrics["mean_between_batch_distance"] = np.mean(between_batch_distances)
        metrics["batch_distance_ratio"] = np.mean(between_batch_distances) / np.mean(within_batch_distances)

    # Cluster-batch association if clusters available
    if cluster_key and cluster_key in integrated_data.obs.columns:
        cluster_labels = integrated_data.obs[cluster_key].values

        # Compute contingency table
        from scipy.stats import chi2_contingency

        contingency = pd.crosstab(cluster_labels, batch_labels)
        chi2, p_value, dof, expected = chi2_contingency(contingency)

        metrics["cluster_batch_chi2"] = chi2
        metrics["cluster_batch_p_value"] = p_value
        metrics["cluster_batch_independence"] = p_value > 0.05  # True if independent

    logger.info("Integration quality evaluation completed")
    return metrics
