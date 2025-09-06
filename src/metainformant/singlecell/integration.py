"""Integration methods for single-cell data.

This module provides methods for integrating multiple single-cell datasets,
including batch correction and harmonization techniques.
All implementations use real computational methods without mocking.
"""

from __future__ import annotations

import warnings
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

try:
    from scipy import sparse

    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    sparse = None

try:
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False
    PCA = None
    StandardScaler = None

from .preprocessing import SingleCellData


def batch_correction(
    data: SingleCellData,
    batch_key: str,
    method: str = "combat",
    n_components: int = 50,
) -> SingleCellData:
    """Perform batch correction on single-cell data.

    Args:
        data: SingleCellData object
        batch_key: Column name in obs indicating batch
        method: Batch correction method ('combat', 'mnn', 'scaling')
        n_components: Number of components for dimensionality-based methods

    Returns:
        SingleCellData with batch-corrected expression
    """
    if batch_key not in data.obs.columns:
        raise ValueError(f"Batch key '{batch_key}' not found in obs")

    data = data.copy()

    if method == "scaling":
        # Simple scaling-based batch correction
        X_corrected = _scaling_batch_correction(data.X, data.obs[batch_key])
        data.X = X_corrected

    elif method == "combat":
        # ComBat-style batch correction (simplified implementation)
        X_corrected = _combat_batch_correction(data.X, data.obs[batch_key])
        data.X = X_corrected

    else:
        raise ValueError(f"Unknown batch correction method: {method}")

    data.uns["batch_correction"] = {
        "method": method,
        "batch_key": batch_key,
    }

    return data


def _scaling_batch_correction(X: np.ndarray, batch_labels: pd.Series) -> np.ndarray:
    """Simple scaling-based batch correction."""
    if sparse.issparse(X):
        X = X.toarray()

    X_corrected = X.copy()
    batches = batch_labels.unique()

    for batch in batches:
        batch_mask = batch_labels == batch
        batch_data = X[batch_mask, :]

        # Center and scale each batch
        batch_mean = np.mean(batch_data, axis=0)
        batch_std = np.std(batch_data, axis=0)
        batch_std[batch_std == 0] = 1  # Avoid division by zero

        # Global statistics
        global_mean = np.mean(X, axis=0)
        global_std = np.std(X, axis=0)
        global_std[global_std == 0] = 1

        # Apply correction
        X_corrected[batch_mask, :] = (batch_data - batch_mean) / batch_std * global_std + global_mean

    return X_corrected


def _combat_batch_correction(X: np.ndarray, batch_labels: pd.Series) -> np.ndarray:
    """Simplified ComBat-style batch correction."""
    if sparse.issparse(X):
        X = X.toarray()

    X_corrected = X.copy()
    batches = batch_labels.unique()

    if len(batches) <= 1:
        return X_corrected

    # Estimate batch effects
    for gene_idx in range(X.shape[1]):
        gene_expr = X[:, gene_idx]

        # Fit model: expression = overall_mean + batch_effect + error
        overall_mean = np.mean(gene_expr)

        batch_effects = {}
        for batch in batches:
            batch_mask = batch_labels == batch
            batch_mean = np.mean(gene_expr[batch_mask])
            batch_effects[batch] = batch_mean - overall_mean

        # Remove batch effects
        for batch in batches:
            batch_mask = batch_labels == batch
            X_corrected[batch_mask, gene_idx] -= batch_effects[batch]

    return X_corrected


def integrate_datasets(
    datasets: List[SingleCellData],
    method: str = "concat",
    batch_keys: Optional[List[str]] = None,
) -> SingleCellData:
    """Integrate multiple single-cell datasets.

    Args:
        datasets: List of SingleCellData objects to integrate
        method: Integration method ('concat', 'intersect', 'union')
        batch_keys: Batch identifiers for each dataset

    Returns:
        Integrated SingleCellData object
    """
    if len(datasets) == 0:
        raise ValueError("No datasets provided")

    if len(datasets) == 1:
        return datasets[0].copy()

    if batch_keys is None:
        batch_keys = [f"batch_{i}" for i in range(len(datasets))]

    if method == "concat":
        return _concatenate_datasets(datasets, batch_keys)
    elif method == "intersect":
        return _intersect_datasets(datasets, batch_keys)
    elif method == "union":
        return _union_datasets(datasets, batch_keys)
    else:
        raise ValueError(f"Unknown integration method: {method}")


def _concatenate_datasets(datasets: List[SingleCellData], batch_keys: List[str]) -> SingleCellData:
    """Concatenate datasets using common genes."""
    # Find common genes
    common_genes = None
    for data in datasets:
        if common_genes is None:
            common_genes = set(data.var.index)
        else:
            common_genes = common_genes.intersection(set(data.var.index))

    common_genes = sorted(list(common_genes))

    if len(common_genes) == 0:
        raise ValueError("No common genes found between datasets")

    # Extract common genes from each dataset
    X_list = []
    obs_list = []

    for i, data in enumerate(datasets):
        # Get indices of common genes
        gene_indices = [data.var.index.get_loc(gene) for gene in common_genes]

        if sparse.issparse(data.X):
            X_subset = data.X[:, gene_indices]
        else:
            X_subset = data.X[:, gene_indices]

        X_list.append(X_subset)

        # Add batch information
        obs_subset = data.obs.copy()
        obs_subset["batch"] = batch_keys[i]
        obs_list.append(obs_subset)

    # Concatenate
    if sparse.issparse(X_list[0]):
        X_integrated = sparse.vstack(X_list)
    else:
        X_integrated = np.vstack(X_list)

    obs_integrated = pd.concat(obs_list, ignore_index=True)

    # Create var DataFrame for common genes
    var_integrated = datasets[0].var.loc[common_genes].copy()

    return SingleCellData(
        X=X_integrated,
        obs=obs_integrated,
        var=var_integrated,
    )


def _intersect_datasets(datasets: List[SingleCellData], batch_keys: List[str]) -> SingleCellData:
    """Same as concatenate for now - could be enhanced with more sophisticated intersection."""
    return _concatenate_datasets(datasets, batch_keys)


def _union_datasets(datasets: List[SingleCellData], batch_keys: List[str]) -> SingleCellData:
    """Union of all genes across datasets (fill missing with zeros)."""
    # Find all genes
    all_genes = set()
    for data in datasets:
        all_genes.update(data.var.index)

    all_genes = sorted(list(all_genes))

    # Create expression matrices with all genes
    X_list = []
    obs_list = []

    for i, data in enumerate(datasets):
        # Create matrix for all genes
        if sparse.issparse(data.X):
            X_full = sparse.lil_matrix((data.n_obs, len(all_genes)))
        else:
            X_full = np.zeros((data.n_obs, len(all_genes)))

        # Fill in existing genes
        for j, gene in enumerate(all_genes):
            if gene in data.var.index:
                gene_idx = data.var.index.get_loc(gene)
                if sparse.issparse(data.X):
                    X_full[:, j] = data.X[:, gene_idx]
                else:
                    X_full[:, j] = data.X[:, gene_idx]

        if sparse.issparse(X_full):
            X_full = X_full.tocsr()

        X_list.append(X_full)

        # Add batch information
        obs_subset = data.obs.copy()
        obs_subset["batch"] = batch_keys[i]
        obs_list.append(obs_subset)

    # Concatenate
    if sparse.issparse(X_list[0]):
        X_integrated = sparse.vstack(X_list)
    else:
        X_integrated = np.vstack(X_list)

    obs_integrated = pd.concat(obs_list, ignore_index=True)

    # Create var DataFrame for all genes
    var_integrated = pd.DataFrame(index=all_genes)
    var_integrated["gene_name"] = all_genes

    return SingleCellData(
        X=X_integrated,
        obs=obs_integrated,
        var=var_integrated,
    )


def harmony_integration(
    data: SingleCellData,
    batch_key: str,
    n_components: int = 50,
    theta: float = 2.0,
    max_iter: int = 20,
) -> SingleCellData:
    """Harmony-inspired batch integration (simplified implementation).

    Args:
        data: SingleCellData with batch information
        batch_key: Column name indicating batch
        n_components: Number of harmony components
        theta: Diversity clustering penalty
        max_iter: Maximum iterations

    Returns:
        SingleCellData with harmony-corrected coordinates
    """
    if batch_key not in data.obs.columns:
        raise ValueError(f"Batch key '{batch_key}' not found in obs")

    data = data.copy()

    # Use PCA coordinates if available, otherwise compute them
    if "X_pca" in data.obsm:
        X_embed = data.obsm["X_pca"][:, :n_components]
    else:
        # Compute PCA
        if sparse.issparse(data.X):
            X_pca_input = data.X.toarray()
        else:
            X_pca_input = data.X

        pca = PCA(n_components=n_components)
        X_embed = pca.fit_transform(X_pca_input)

    # Get batch information
    batch_labels = data.obs[batch_key]
    unique_batches = batch_labels.unique()
    n_batches = len(unique_batches)

    if n_batches <= 1:
        warnings.warn("Only one batch found, no correction needed")
        data.obsm["X_harmony"] = X_embed
        return data

    # Simplified harmony correction
    X_corrected = X_embed.copy()

    for iteration in range(max_iter):
        # For each batch, adjust embeddings to reduce batch effects
        for batch in unique_batches:
            batch_mask = batch_labels == batch
            if batch_mask.sum() == 0:
                continue

            batch_cells = X_embed[batch_mask, :]
            other_cells = X_embed[~batch_mask, :]

            if len(other_cells) == 0:
                continue

            # Compute batch center and global center
            batch_center = np.mean(batch_cells, axis=0)
            global_center = np.mean(X_embed, axis=0)

            # Correction factor (simplified)
            correction = (global_center - batch_center) * 0.1  # Small correction

            # Apply correction
            X_corrected[batch_mask, :] += correction

        X_embed = X_corrected.copy()

    # Store corrected embeddings
    data.obsm["X_harmony"] = X_corrected
    data.uns["harmony"] = {
        "batch_key": batch_key,
        "n_components": n_components,
        "theta": theta,
        "max_iter": max_iter,
    }

    return data
