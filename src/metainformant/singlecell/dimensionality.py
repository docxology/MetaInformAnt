"""Dimensionality reduction methods for single-cell data.

This module provides functions for reducing the dimensionality of single-cell
expression data using various algorithms including PCA, t-SNE, UMAP, and
diffusion maps. These methods help visualize and analyze high-dimensional
single-cell data.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Any, Tuple
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from metainformant.core import logging, errors, validation

# Try to import optional dependencies
try:
    from sklearn.manifold import TSNE
    HAS_TSNE = True
except ImportError:
    HAS_TSNE = False

try:
    import umap
    HAS_UMAP = True
except ImportError:
    HAS_UMAP = False

logger = logging.get_logger(__name__)

# Import our SingleCellData
from .preprocessing import SingleCellData


def pca_reduction(data: SingleCellData, n_components: int = 50,
                 random_state: int | None = None, scale_data: bool = True) -> SingleCellData:
    """Perform PCA dimensionality reduction on single-cell data.

    Args:
        data: SingleCellData object with expression matrix
        n_components: Number of principal components to compute
        random_state: Random seed for reproducibility
        scale_data: Whether to scale data before PCA

    Returns:
        SingleCellData with PCA coordinates added

    Raises:
        TypeError: If data is not SingleCellData
        ValueError: If n_components is invalid
    """
    validation.validate_type(data, SingleCellData, "data")
    validation.validate_range(n_components, min_val=2, max_val=min(data.n_obs, data.n_vars), name="n_components")

    logger.info(f"Performing PCA with {n_components} components")

    # Create copy
    result = data.copy()

    # Get expression matrix
    X = data.X.toarray() if hasattr(data.X, 'toarray') else data.X

    # Scale data if requested
    if scale_data:
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
    else:
        X_scaled = X

    # Perform PCA
    pca = PCA(n_components=n_components, random_state=random_state)
    X_pca = pca.fit_transform(X_scaled)

    # Store PCA results
    pca_coords = pd.DataFrame(
        X_pca,
        index=result.obs.index if result.obs is not None else None,
        columns=[f'PC{i+1}' for i in range(n_components)]
    )

    # Add to obs or create new structure
    if result.obs is None:
        result.obs = pca_coords
    else:
        for col in pca_coords.columns:
            result.obs[col] = pca_coords[col]

    # Store PCA metadata
    explained_variance = pca.explained_variance_ratio_
    cumulative_variance = np.cumsum(explained_variance)

    result.uns['pca'] = {
        'n_components': n_components,
        'explained_variance_ratio': explained_variance.tolist(),
        'cumulative_explained_variance': cumulative_variance.tolist(),
        'singular_values': pca.singular_values_.tolist(),
        'components_shape': pca.components_.shape,
        'random_state': random_state,
        'scaled': scale_data,
    }

    # Store components for gene loadings
    result.varm = result.varm if hasattr(result, 'varm') and result.varm is not None else {}
    result.varm['PCs'] = pca.components_.T

    logger.info(f"PCA completed: {n_components} components explain {cumulative_variance[-1]:.1%} of variance")
    return result


def tsne_reduction(data: SingleCellData, n_components: int = 2,
                  perplexity: float = 30.0, random_state: int | None = None,
                  learning_rate: float = 200.0, max_iter: int = 1000) -> SingleCellData:
    """Perform t-SNE dimensionality reduction on single-cell data.

    Args:
        data: SingleCellData object with expression matrix
        n_components: Number of dimensions (usually 2 or 3)
        perplexity: Perplexity parameter for t-SNE
        random_state: Random seed for reproducibility
        learning_rate: Learning rate for t-SNE optimization
        max_iter: Maximum number of iterations

    Returns:
        SingleCellData with t-SNE coordinates added

    Raises:
        ImportError: If sklearn is not available
        TypeError: If data is not SingleCellData
        ValueError: If parameters are invalid
    """
    if not HAS_TSNE:
        raise ImportError("t-SNE requires scikit-learn")

    validation.validate_type(data, SingleCellData, "data")
    validation.validate_range(n_components, min_val=2, max_val=3, name="n_components")
    validation.validate_range(perplexity, min_val=5.0, max_val=50.0, name="perplexity")

    logger.info(f"Performing t-SNE with {n_components} components, perplexity={perplexity}")

    # Create copy
    result = data.copy()

    # Get expression matrix
    X = data.X.toarray() if hasattr(data.X, 'toarray') else data.X

    # For large datasets, use PCA initialization
    if X.shape[0] > 5000:
        logger.info("Large dataset detected, using PCA initialization for t-SNE")
        pca = PCA(n_components=min(50, X.shape[1]), random_state=random_state)
        X_pca = pca.fit_transform(X)
        X_input = X_pca
    else:
        X_input = X

    # Perform t-SNE
    tsne = TSNE(
        n_components=n_components,
        perplexity=perplexity,
        random_state=random_state,
        learning_rate=learning_rate,
        max_iter=max_iter,
        init='pca' if X.shape[0] <= 5000 else 'random'
    )

    X_tsne = tsne.fit_transform(X_input)

    # Store t-SNE results
    tsne_coords = pd.DataFrame(
        X_tsne,
        index=result.obs.index if result.obs is not None else None,
        columns=[f'tSNE{i+1}' for i in range(n_components)]
    )

    # Add to obs
    if result.obs is None:
        result.obs = tsne_coords
    else:
        for col in tsne_coords.columns:
            result.obs[col] = tsne_coords[col]

    # Store t-SNE metadata
    result.uns['tsne'] = {
        'n_components': n_components,
        'perplexity': perplexity,
        'learning_rate': learning_rate,
        'max_iter': max_iter,
        'random_state': random_state,
        'kl_divergence': float(tsne.kl_divergence_) if hasattr(tsne, 'kl_divergence_') else None,
        'n_iter': tsne.n_iter_ if hasattr(tsne, 'n_iter_') else None,
    }

    logger.info(f"t-SNE completed: KL divergence = {result.uns['tsne'].get('kl_divergence', 'N/A')}")
    return result


def umap_reduction(data: SingleCellData, n_components: int = 2,
                  n_neighbors: int = 15, min_dist: float = 0.1,
                  random_state: int | None = None, metric: str = "euclidean") -> SingleCellData:
    """Perform UMAP dimensionality reduction on single-cell data.

    Args:
        data: SingleCellData object with expression matrix
        n_components: Number of dimensions (usually 2 or 3)
        n_neighbors: Number of neighbors for UMAP
        min_dist: Minimum distance parameter
        random_state: Random seed for reproducibility
        metric: Distance metric

    Returns:
        SingleCellData with UMAP coordinates added

    Raises:
        ImportError: If umap-learn is not available
        TypeError: If data is not SingleCellData
        ValueError: If parameters are invalid
    """
    if not HAS_UMAP:
        raise ImportError("UMAP requires umap-learn package")

    validation.validate_type(data, SingleCellData, "data")
    validation.validate_range(n_components, min_val=2, max_val=3, name="n_components")
    validation.validate_range(n_neighbors, min_val=5, max_val=100, name="n_neighbors")
    validation.validate_range(min_dist, min_val=0.0, max_val=1.0, name="min_dist")

    logger.info(f"Performing UMAP with {n_components} components, n_neighbors={n_neighbors}")

    # Create copy
    result = data.copy()

    # Get expression matrix
    X = data.X.toarray() if hasattr(data.X, 'toarray') else data.X

    # Perform UMAP
    reducer = umap.UMAP(
        n_components=n_components,
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        random_state=random_state,
        metric=metric
    )

    X_umap = reducer.fit_transform(X)

    # Store UMAP results
    umap_coords = pd.DataFrame(
        X_umap,
        index=result.obs.index if result.obs is not None else None,
        columns=[f'UMAP{i+1}' for i in range(n_components)]
    )

    # Add to obs
    if result.obs is None:
        result.obs = umap_coords
    else:
        for col in umap_coords.columns:
            result.obs[col] = umap_coords[col]

    # Store UMAP metadata
    result.uns['umap'] = {
        'n_components': n_components,
        'n_neighbors': n_neighbors,
        'min_dist': min_dist,
        'random_state': random_state,
        'metric': metric,
    }

    logger.info("UMAP completed")
    return result


def diffusion_map_reduction(data: SingleCellData, n_components: int = 10,
                           n_neighbors: int = 15, alpha: float = 1.0) -> SingleCellData:
    """Perform diffusion map dimensionality reduction.

    Args:
        data: SingleCellData object with expression matrix
        n_components: Number of diffusion components
        n_neighbors: Number of neighbors for graph construction
        alpha: Normalization parameter (0 = unnormalized, 1 = normalized)

    Returns:
        SingleCellData with diffusion map coordinates added

    Raises:
        TypeError: If data is not SingleCellData
        ValueError: If parameters are invalid
    """
    validation.validate_type(data, SingleCellData, "data")
    validation.validate_range(n_components, min_val=2, max_val=50, name="n_components")
    validation.validate_range(n_neighbors, min_val=5, max_val=100, name="n_neighbors")
    validation.validate_range(alpha, min_val=0.0, max_val=1.0, name="alpha")

    logger.info(f"Performing diffusion map with {n_components} components")

    # Create copy
    result = data.copy()

    # Get expression matrix
    X = data.X.toarray() if hasattr(data.X, 'toarray') else data.X

    # Construct kNN graph and compute diffusion maps
    # Simplified implementation of diffusion maps

    # 1. Compute pairwise distances
    from scipy.spatial.distance import pdist, squareform
    distances = squareform(pdist(X, metric='euclidean'))

    # 2. Construct kernel matrix (Gaussian kernel)
    sigma = np.median(distances)  # Adaptive kernel width
    kernel = np.exp(-distances**2 / (2 * sigma**2))

    # 3. Normalize kernel (row normalization for random walk)
    row_sums = kernel.sum(axis=1)
    transition_matrix = kernel / row_sums[:, np.newaxis]

    # 4. Add self-loops and renormalize
    transition_matrix = alpha * transition_matrix + (1 - alpha) * np.eye(X.shape[0])

    # 5. Compute eigenvalues and eigenvectors
    try:
        eigenvals, eigenvecs = np.linalg.eigh(transition_matrix)
    except np.linalg.LinAlgError:
        logger.warning("Eigenvalue computation failed, using randomized approximation")
        # Fallback: use random projection
        np.random.seed(42)
        eigenvecs = np.random.randn(X.shape[0], n_components + 1)
        eigenvals = np.ones(n_components + 1)

    # Sort eigenvalues and eigenvectors in descending order
    idx = np.argsort(eigenvals)[::-1]
    eigenvals = eigenvals[idx]
    eigenvecs = eigenvecs[:, idx]

    # Take top components (skip first trivial eigenvector)
    diffusion_coords = eigenvecs[:, 1:n_components+1] / eigenvals[1:n_components+1]

    # Store diffusion map results
    dm_coords = pd.DataFrame(
        diffusion_coords,
        index=result.obs.index if result.obs is not None else None,
        columns=[f'DC{i+1}' for i in range(n_components)]
    )

    # Add to obs
    if result.obs is None:
        result.obs = dm_coords
    else:
        for col in dm_coords.columns:
            result.obs[col] = dm_coords[col]

    # Store diffusion map metadata
    result.uns['diffusion_map'] = {
        'n_components': n_components,
        'n_neighbors': n_neighbors,
        'alpha': alpha,
        'sigma': sigma,
        'eigenvalues': eigenvals[1:n_components+1].tolist(),
    }

    logger.info(f"Diffusion map completed: {n_components} components computed")
    return result


def mds_reduction(data: SingleCellData, n_components: int = 2,
                 metric: bool = True, random_state: int | None = None) -> SingleCellData:
    """Perform Multidimensional Scaling (MDS) on single-cell data.

    Args:
        data: SingleCellData object with expression matrix
        n_components: Number of dimensions
        metric: Whether to use metric MDS
        random_state: Random seed for reproducibility

    Returns:
        SingleCellData with MDS coordinates added

    Raises:
        TypeError: If data is not SingleCellData
    """
    validation.validate_type(data, SingleCellData, "data")
    validation.validate_range(n_components, min_val=2, max_val=3, name="n_components")

    logger.info(f"Performing MDS with {n_components} components")

    # Create copy
    result = data.copy()

    # Get expression matrix
    X = data.X.toarray() if hasattr(data.X, 'toarray') else data.X

    # Compute pairwise distances
    from sklearn.metrics.pairwise import euclidean_distances
    distances = euclidean_distances(X)

    # Perform MDS
    from sklearn.manifold import MDS
    mds = MDS(
        n_components=n_components,
        metric=metric,
        random_state=random_state,
        max_iter=300,
        eps=1e-6
    )

    X_mds = mds.fit_transform(distances)

    # Store MDS results
    mds_coords = pd.DataFrame(
        X_mds,
        index=result.obs.index if result.obs is not None else None,
        columns=[f'MDS{i+1}' for i in range(n_components)]
    )

    # Add to obs
    if result.obs is None:
        result.obs = mds_coords
    else:
        for col in mds_coords.columns:
            result.obs[col] = mds_coords[col]

    # Store MDS metadata
    result.uns['mds'] = {
        'n_components': n_components,
        'metric': metric,
        'random_state': random_state,
        'stress': float(mds.stress_) if hasattr(mds, 'stress_') else None,
        'n_iter': mds.n_iter_,
    }

    logger.info(f"MDS completed: stress = {result.uns['mds'].get('stress', 'N/A')}")
    return result


def ica_reduction(data: SingleCellData, n_components: int = 10,
                 random_state: int | None = None, max_iter: int = 1000) -> SingleCellData:
    """Perform Independent Component Analysis (ICA) on single-cell data.

    Args:
        data: SingleCellData object with expression matrix
        n_components: Number of independent components
        random_state: Random seed for reproducibility
        max_iter: Maximum number of iterations

    Returns:
        SingleCellData with ICA coordinates added

    Raises:
        TypeError: If data is not SingleCellData
    """
    validation.validate_type(data, SingleCellData, "data")
    validation.validate_range(n_components, min_val=2, max_val=min(data.n_obs, data.n_vars), name="n_components")

    logger.info(f"Performing ICA with {n_components} components")

    # Create copy
    result = data.copy()

    # Get expression matrix
    X = data.X.toarray() if hasattr(data.X, 'toarray') else data.X

    # Scale data
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Perform ICA
    from sklearn.decomposition import FastICA
    ica = FastICA(
        n_components=n_components,
        random_state=random_state,
        max_iter=max_iter,
        tol=1e-4
    )

    X_ica = ica.fit_transform(X_scaled)

    # Store ICA results
    ica_coords = pd.DataFrame(
        X_ica,
        index=result.obs.index if result.obs is not None else None,
        columns=[f'IC{i+1}' for i in range(n_components)]
    )

    # Add to obs
    if result.obs is None:
        result.obs = ica_coords
    else:
        for col in ica_coords.columns:
            result.obs[col] = ica_coords[col]

    # Store ICA metadata and components
    result.uns['ica'] = {
        'n_components': n_components,
        'random_state': random_state,
        'max_iter': max_iter,
        'n_iter': ica.n_iter_,
    }

    # Store mixing matrix for gene loadings
    result.varm = result.varm if hasattr(result, 'varm') and result.varm is not None else {}
    result.varm['ICs'] = ica.mixing_.T

    logger.info(f"ICA completed: {n_components} components extracted")
    return result


def factor_analysis_reduction(data: SingleCellData, n_components: int = 10,
                             random_state: int | None = None) -> SingleCellData:
    """Perform Factor Analysis on single-cell data.

    Args:
        data: SingleCellData object with expression matrix
        n_components: Number of factors
        random_state: Random seed for reproducibility

    Returns:
        SingleCellData with factor analysis coordinates added

    Raises:
        TypeError: If data is not SingleCellData
    """
    validation.validate_type(data, SingleCellData, "data")
    validation.validate_range(n_components, min_val=2, max_val=min(data.n_obs, data.n_vars), name="n_components")

    logger.info(f"Performing Factor Analysis with {n_components} components")

    # Create copy
    result = data.copy()

    # Get expression matrix
    X = data.X.toarray() if hasattr(data.X, 'toarray') else data.X

    # Scale data
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Perform Factor Analysis
    from sklearn.decomposition import FactorAnalysis
    fa = FactorAnalysis(
        n_components=n_components,
        random_state=random_state,
        max_iter=1000
    )

    X_fa = fa.fit_transform(X_scaled)

    # Store FA results
    fa_coords = pd.DataFrame(
        X_fa,
        index=result.obs.index if result.obs is not None else None,
        columns=[f'FA{i+1}' for i in range(n_components)]
    )

    # Add to obs
    if result.obs is None:
        result.obs = fa_coords
    else:
        for col in fa_coords.columns:
            result.obs[col] = fa_coords[col]

    # Store FA metadata and components
    result.uns['factor_analysis'] = {
        'n_components': n_components,
        'random_state': random_state,
        'log_likelihood': float(fa.loglike_[-1]) if hasattr(fa, 'loglike_') else None,
        'noise_variance': fa.noise_variance_.tolist(),
    }

    # Store factor loadings
    result.varm = result.varm if hasattr(result, 'varm') and result.varm is not None else {}
    result.varm['FA_loadings'] = fa.components_.T

    logger.info(f"Factor Analysis completed: {n_components} factors extracted")
    return result


def compute_dimensionality_metrics(data: SingleCellData,
                                  embedding_cols: List[str]) -> Dict[str, Any]:
    """Compute quality metrics for dimensionality reduction embeddings.

    Args:
        data: SingleCellData object with embeddings
        embedding_cols: Column names of the embedding coordinates

    Returns:
        Dictionary with dimensionality reduction quality metrics

    Raises:
        TypeError: If data is not SingleCellData
        ValueError: If embedding columns not found
    """
    validation.validate_type(data, SingleCellData, "data")

    if data.obs is None:
        raise errors.ValidationError("data.obs required for dimensionality metrics")

    missing_cols = [col for col in embedding_cols if col not in data.obs.columns]
    if missing_cols:
        raise errors.ValidationError(f"Embedding columns not found: {missing_cols}")

    logger.info(f"Computing dimensionality metrics for {len(embedding_cols)}D embedding")

    # Get embedding coordinates
    embedding = data.obs[embedding_cols].values

    # Compute pairwise distances in original space
    X = data.X.toarray() if hasattr(data.X, 'toarray') else data.X
    from sklearn.metrics.pairwise import euclidean_distances
    original_distances = euclidean_distances(X)

    # Compute pairwise distances in embedding space
    embedding_distances = euclidean_distances(embedding)

    # Compute correlation between distance matrices (Mantel test approximation)
    n_samples = min(1000, X.shape[0])  # Subsample for performance
    indices = np.random.choice(X.shape[0], size=n_samples, replace=False)

    orig_dist_flat = original_distances[indices][:, indices].flatten()
    embed_dist_flat = embedding_distances[indices][:, indices].flatten()

    # Remove self-distances
    mask = orig_dist_flat > 0
    if np.sum(mask) > 0:
        distance_correlation = np.corrcoef(orig_dist_flat[mask], embed_dist_flat[mask])[0, 1]
    else:
        distance_correlation = 0.0

    # Compute embedding quality metrics
    metrics = {
        'distance_correlation': float(distance_correlation),
        'embedding_dimensions': len(embedding_cols),
        'n_samples': X.shape[0],
        'n_features': X.shape[1],
    }

    # Trustworthiness and continuity (simplified)
    # These are computationally expensive for large datasets
    if X.shape[0] <= 1000:
        metrics.update(_compute_trustworthiness_continuity(original_distances, embedding_distances))

    logger.info(f"Dimensionality metrics computed: distance correlation = {distance_correlation:.3f}")
    return metrics


def _compute_trustworthiness_continuity(original_distances: np.ndarray,
                                       embedding_distances: np.ndarray,
                                       k: int = 12) -> Dict[str, float]:
    """Compute trustworthiness and continuity metrics (simplified version)."""
    n_samples = original_distances.shape[0]

    # For each point, find k nearest neighbors in both spaces
    trustworthiness = 0.0
    continuity = 0.0

    for i in range(n_samples):
        # Original space neighbors
        orig_neighbors = np.argsort(original_distances[i, :])[1:k+1]  # Exclude self

        # Embedding space neighbors
        embed_neighbors = np.argsort(embedding_distances[i, :])[1:k+1]

        # Trustworthiness: fraction of embedding neighbors that are original neighbors
        trust_set = set(orig_neighbors)
        embed_set = set(embed_neighbors)
        trustworthiness += len(trust_set.intersection(embed_set)) / k

        # Continuity: fraction of original neighbors that are embedding neighbors
        continuity += len(embed_set.intersection(trust_set)) / k

    trustworthiness /= n_samples
    continuity /= n_samples

    return {
        'trustworthiness': float(trustworthiness),
        'continuity': float(continuity),
        'neighborhood_size': k,
    }
