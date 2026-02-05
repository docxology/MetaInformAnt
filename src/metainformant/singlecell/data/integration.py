"""Single-cell data integration and batch correction methods.

This module provides functions for integrating multiple single-cell datasets
by correcting for batch effects and technical variation. Methods include
mutual nearest neighbors (MNN), BBKNN, Harmony, and Scanorama-like approaches.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from metainformant.core import errors, logging, validation

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

# Import our SingleCellData
from .preprocessing import SingleCellData


def bbknn_integration(data: SingleCellData, batch_key: str, n_neighbors: int = 15) -> SingleCellData:
    """Perform BBKNN (Batch Balanced K-Nearest Neighbors) integration.

    Args:
        data: SingleCellData object with batch labels
        batch_key: Column name containing batch labels
        n_neighbors: Number of neighbors for kNN graph

    Returns:
        SingleCellData with integrated coordinates

    Raises:
        TypeError: If data is not SingleCellData
        ValueError: If batch_key not found
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
    X = data.X.toarray() if hasattr(data.X, "toarray") else data.X
    batch_labels = data.obs[batch_key].values
    unique_batches = np.unique(batch_labels)

    logger.info(f"Found {len(unique_batches)} batches: {unique_batches}")

    # For BBKNN, we need to construct a batch-aware kNN graph
    # This is a simplified implementation

    # First, compute PCA for dimensionality reduction
    pca = PCA(n_components=min(50, X.shape[1]), random_state=42)
    X_pca = pca.fit_transform(StandardScaler().fit_transform(X))

    # Compute batch-balanced kNN graph
    n_cells = X.shape[0]
    batch_adjacency = np.zeros((n_cells, n_cells))

    # For each batch pair, compute cross-batch neighbors
    for i, batch1 in enumerate(unique_batches):
        for j, batch2 in enumerate(unique_batches):
            if i <= j:  # Include self-batch and upper triangle
                batch1_mask = batch_labels == batch1
                batch2_mask = batch_labels == batch2

                batch1_indices = np.where(batch1_mask)[0]
                batch2_indices = np.where(batch2_mask)[0]

                if len(batch1_indices) == 0 or len(batch2_indices) == 0:
                    continue

                # Compute distances between batches
                X_batch1 = X_pca[batch1_indices]
                X_batch2 = X_pca[batch2_indices]

                # Use kNN to find mutual nearest neighbors
                from sklearn.neighbors import NearestNeighbors

                # Find neighbors from batch1 to batch2
                if len(batch2_indices) >= n_neighbors:
                    nbrs = NearestNeighbors(n_neighbors=n_neighbors, metric="euclidean")
                    nbrs.fit(X_batch2)
                    distances, indices = nbrs.kneighbors(X_batch1)

                    # Convert local indices to global indices
                    for local_i, global_i in enumerate(batch1_indices):
                        for k in range(n_neighbors):
                            global_j = batch2_indices[indices[local_i, k]]
                            batch_adjacency[global_i, global_j] = 1.0 / (distances[local_i, k] + 1e-6)

                # Find neighbors from batch2 to batch1 (for undirected graph)
                if i != j and len(batch1_indices) >= n_neighbors:  # Don't duplicate self-batch
                    nbrs = NearestNeighbors(n_neighbors=n_neighbors, metric="euclidean")
                    nbrs.fit(X_batch1)
                    distances, indices = nbrs.kneighbors(X_batch2)

                    for local_i, global_i in enumerate(batch2_indices):
                        for k in range(n_neighbors):
                            global_j = batch1_indices[indices[local_i, k]]
                            batch_adjacency[global_i, global_j] = 1.0 / (distances[local_i, k] + 1e-6)

    # Store the integrated adjacency matrix
    result.uns["bbknn"] = {
        "batch_key": batch_key,
        "n_neighbors": n_neighbors,
        "n_batches": len(unique_batches),
        "batch_adjacency_shape": batch_adjacency.shape,
    }

    # For visualization, we can compute a 2D embedding from the graph
    # This is a simplified approach
    try:
        from sklearn.manifold import TSNE

        X_embed = TSNE(n_components=2, random_state=42, perplexity=min(30, n_cells - 1)).fit_transform(X_pca)
        result.obs["bbknn_1"] = X_embed[:, 0]
        result.obs["bbknn_2"] = X_embed[:, 1]
    except Exception as e:
        logger.warning(f"Could not compute BBKNN embedding: {e}")

    logger.info("BBKNN integration completed")
    return result


def harmony_integration(data: SingleCellData, batch_key: str, n_components: int = 50) -> SingleCellData:
    """Perform Harmony integration for batch correction.

    Args:
        data: SingleCellData object with batch labels
        batch_key: Column name containing batch labels
        n_components: Number of Harmony components

    Returns:
        SingleCellData with Harmony-corrected coordinates

    Raises:
        TypeError: If data is not SingleCellData
        ValueError: If batch_key not found
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
    X = data.X.toarray() if hasattr(data.X, "toarray") else data.X
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
    """Perform Scanorama-like integration of multiple datasets.

    Args:
        data_list: List of SingleCellData objects to integrate
        batch_key: Column name that will contain batch labels

    Returns:
        Integrated SingleCellData object

    Raises:
        TypeError: If data_list contains non-SingleCellData objects
        ValueError: If data_list is empty
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

    logger.info(f"Performing Scanorama integration of {len(data_list)} datasets")

    # Concatenate all datasets
    all_X = []
    all_obs = []
    batch_labels = []

    for i, data in enumerate(data_list):
        X = data.X.toarray() if hasattr(data.X, "toarray") else data.X
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
    """Perform Mutual Nearest Neighbors (MNN) integration.

    Args:
        data_list: List of SingleCellData objects to integrate
        batch_key: Column name that will contain batch labels

    Returns:
        Integrated SingleCellData object

    Raises:
        TypeError: If data_list contains non-SingleCellData objects
        ValueError: If data_list has fewer than 2 datasets
        ImportError: If scikit-learn not available
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn is required for MNN integration. " "Install with: uv pip install scikit-learn")

    if len(data_list) < 2:
        raise errors.ValidationError("MNN integration requires at least 2 datasets")

    for i, data in enumerate(data_list):
        if not isinstance(data, SingleCellData):
            raise errors.ValidationError(f"data_list[{i}] is not a SingleCellData object")

    logger.info(f"Performing MNN integration of {len(data_list)} datasets")

    # Start with the first dataset
    integrated_data = data_list[0].copy()
    integrated_data.obs = integrated_data.obs.copy() if integrated_data.obs is not None else pd.DataFrame()
    integrated_data.obs[batch_key] = "batch_0"

    # Iteratively integrate additional datasets
    for i, data in enumerate(data_list[1:], 1):
        logger.info(f"Integrating dataset {i}")

        # Get expression matrices
        X1 = integrated_data.X.toarray() if hasattr(integrated_data.X, "toarray") else integrated_data.X
        X2 = data.X.toarray() if hasattr(data.X, "toarray") else data.X

        # Find mutual nearest neighbors
        mnn_pairs = _find_mutual_nearest_neighbors(X1, X2, k=20)

        if len(mnn_pairs) == 0:
            logger.warning(f"No MNN pairs found between integrated data and dataset {i}")
            continue

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
    """Perform ComBat batch effect correction.

    Args:
        data: SingleCellData object with batch labels
        batch_key: Column name containing batch labels
        covariates: Optional list of covariate columns to preserve

    Returns:
        Batch-corrected SingleCellData

    Raises:
        TypeError: If data is not SingleCellData
        ValueError: If batch_key not found
    """
    validation.validate_type(data, SingleCellData, "data")

    if data.obs is None or batch_key not in data.obs.columns:
        raise errors.ValidationError(f"Batch key '{batch_key}' not found in data.obs")

    logger.info(f"Performing ComBat integration using {batch_key} for batches")

    # Create copy
    result = data.copy()

    # Get expression matrix and batch labels
    X = data.X.toarray() if hasattr(data.X, "toarray") else data.X
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
        ref_mean = batch_means[ref_batch]
        ref_var = batch_vars[ref_batch]

        # Correct each batch to reference
        for batch in unique_batches:
            if batch != ref_batch:
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
        "covariates": covariates or [],
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
    data_list: List[SingleCellData], integration_method: str = "scanorama", **kwargs
) -> SingleCellData:
    """Integrate multiple single-cell datasets using specified method.

    Args:
        data_list: List of SingleCellData objects to integrate
        integration_method: Integration method ("scanorama", "mnn", "harmony")
        **kwargs: Additional arguments for the integration method

    Returns:
        Integrated SingleCellData object

    Raises:
        ValueError: If integration_method is unsupported
    """
    valid_methods = ["scanorama", "mnn", "harmony", "bbknn"]
    if integration_method not in valid_methods:
        raise errors.ValidationError(f"Unsupported integration method: {integration_method}")

    logger.info(f"Integrating {len(data_list)} datasets using {integration_method}")

    if integration_method == "scanorama":
        return scanorama_integration(data_list, **kwargs)
    elif integration_method == "mnn":
        return mnn_integration(data_list, **kwargs)
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
    else:
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
        TypeError: If integrated_data is not SingleCellData
        ValueError: If required columns not found
    """
    validation.validate_type(integrated_data, SingleCellData, "integrated_data")

    if integrated_data.obs is None or batch_key not in integrated_data.obs.columns:
        raise errors.ValidationError(f"Batch key '{batch_key}' not found in integrated_data.obs")

    logger.info("Evaluating integration quality")

    batch_labels = integrated_data.obs[batch_key].values
    X = integrated_data.X.toarray() if hasattr(integrated_data.X, "toarray") else integrated_data.X

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
