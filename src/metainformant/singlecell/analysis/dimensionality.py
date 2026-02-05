"""Dimensionality reduction methods for single-cell data.

This module provides functions for reducing the dimensionality of single-cell
expression data using various algorithms including PCA, t-SNE, UMAP, and
diffusion maps. These methods help visualize and analyze high-dimensional
single-cell data.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from metainformant.core import errors, logging, validation

# Try to import optional dependencies
try:
    from sklearn.decomposition import PCA
    from sklearn.manifold import TSNE
    from sklearn.preprocessing import StandardScaler

    HAS_SKLEARN = True
    HAS_TSNE = True
except ImportError:
    HAS_SKLEARN = False
    HAS_TSNE = False
    PCA = None
    StandardScaler = None
    TSNE = None

try:
    import umap

    HAS_UMAP = True
except ImportError:
    HAS_UMAP = False

logger = logging.get_logger(__name__)

# Import our SingleCellData
from metainformant.singlecell.data.preprocessing import SingleCellData


def pca_reduction(
    data: SingleCellData,
    n_components: int = 50,
    random_state: int | None = None,
    scale_data: bool = True,
    *,
    use_hvgs: bool = False,
) -> SingleCellData:
    """Perform PCA dimensionality reduction on single-cell data.

    Args:
        data: SingleCellData object with expression matrix
        n_components: Number of principal components to compute
        random_state: Random seed for reproducibility
        scale_data: Whether to scale data before PCA
        use_hvgs: Whether using only highly variable genes

    Returns:
        SingleCellData with PCA coordinates added

    Raises:
        TypeError: If data is not SingleCellData
        ValueError: If n_components is invalid
        ImportError: If scikit-learn not available
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn is required for PCA. " "Install with: uv pip install scikit-learn")

    validation.validate_type(data, SingleCellData, "data")
    validation.validate_range(n_components, min_val=2, max_val=min(data.n_obs, data.n_vars), name="n_components")

    logger.info(f"Performing PCA with {n_components} components")

    # Create copy
    result = data.copy()

    # Get expression matrix
    X = data.X.toarray() if hasattr(data.X, "toarray") else data.X

    # Scale data if requested
    if scale_data:
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
    else:
        X_scaled = X

    # Perform PCA
    pca = PCA(n_components=n_components, random_state=random_state)
    X_pca = pca.fit_transform(X_scaled)

    # Store PCA results in obsm (standard location)
    result.obsm = result.obsm if hasattr(result, "obsm") and result.obsm is not None else {}
    result.obsm["X_pca"] = X_pca

    # Also add PCA coordinates to obs for convenience
    pca_coords = pd.DataFrame(
        X_pca,
        index=result.obs.index if result.obs is not None else None,
        columns=[f"PC{i+1}" for i in range(n_components)],
    )
    if result.obs is None:
        result.obs = pca_coords
    else:
        for col in pca_coords.columns:
            result.obs[col] = pca_coords[col]

    # Store PCA metadata
    explained_variance = pca.explained_variance_ratio_
    cumulative_variance = np.cumsum(explained_variance)

    result.uns["pca"] = {
        "n_components": n_components,
        "explained_variance_ratio": explained_variance.tolist(),
        "variance_ratio": explained_variance.tolist(),  # Alias
        "variance": explained_variance.tolist(),  # Alias
        "cumulative_explained_variance": cumulative_variance.tolist(),
        "cumulative_variance": cumulative_variance.tolist(),  # Alias
        "singular_values": pca.singular_values_.tolist(),
        "components_shape": pca.components_.shape,
        "random_state": random_state,
        "scaled": scale_data,
        "use_hvgs": use_hvgs,
    }

    # Store components for gene loadings
    result.varm = result.varm if hasattr(result, "varm") and result.varm is not None else {}
    result.varm["PCs"] = pca.components_.T

    logger.info(f"PCA completed: {n_components} components explain {cumulative_variance[-1]:.1%} of variance")
    return result


def tsne_reduction(
    data: SingleCellData,
    n_components: int = 2,
    perplexity: float = 30.0,
    random_state: int | None = None,
    learning_rate: float = 200.0,
    max_iter: int = 1000,
) -> SingleCellData:
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

    # Create copy
    result = data.copy()

    # Get expression matrix
    X = data.X.toarray() if hasattr(data.X, "toarray") else data.X

    # Automatically adjust perplexity if too high for dataset size
    n_samples = X.shape[0]
    actual_perplexity = perplexity
    if perplexity >= n_samples:
        actual_perplexity = max(1.0, n_samples - 1)
        logger.warning(f"perplexity reduced from {perplexity} to {actual_perplexity} (must be < n_samples={n_samples})")

    logger.info(f"Performing t-SNE with {n_components} components, perplexity={actual_perplexity}")

    # Ensure max_iter is at least 250 (sklearn's minimum)
    actual_max_iter = max(250, max_iter)
    if actual_max_iter != max_iter:
        logger.warning(f"max_iter increased from {max_iter} to {actual_max_iter} (sklearn minimum)")

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
        perplexity=actual_perplexity,
        random_state=random_state,
        learning_rate=learning_rate,
        max_iter=actual_max_iter,
        init="pca" if X.shape[0] <= 5000 else "random",
    )

    X_tsne = tsne.fit_transform(X_input)

    # Store t-SNE results
    tsne_coords = pd.DataFrame(
        X_tsne,
        index=result.obs.index if result.obs is not None else None,
        columns=[f"tSNE{i+1}" for i in range(n_components)],
    )

    # Add to obs
    if result.obs is None:
        result.obs = tsne_coords
    else:
        for col in tsne_coords.columns:
            result.obs[col] = tsne_coords[col]

    # Store in obsm as well (standard location)
    result.obsm = result.obsm if hasattr(result, "obsm") and result.obsm is not None else {}
    result.obsm["X_tsne"] = X_tsne

    # Store t-SNE metadata
    result.uns["tsne"] = {
        "n_components": n_components,
        "params": {
            "n_components": n_components,
            "perplexity": actual_perplexity,
            "learning_rate": learning_rate,
            "max_iter": actual_max_iter,
            "n_iter": actual_max_iter,  # Alias for max_iter (consistent with param name)
        },
        "perplexity": actual_perplexity,
        "learning_rate": learning_rate,
        "max_iter": actual_max_iter,
        "random_state": random_state,
        "kl_divergence": float(tsne.kl_divergence_) if hasattr(tsne, "kl_divergence_") else None,
        "n_iter": tsne.n_iter_ if hasattr(tsne, "n_iter_") else None,
    }

    logger.info(f"t-SNE completed: KL divergence = {result.uns['tsne'].get('kl_divergence', 'N/A')}")
    return result


def umap_reduction(
    data: SingleCellData,
    n_components: int = 2,
    n_neighbors: int = 15,
    min_dist: float = 0.1,
    random_state: int | None = None,
    metric: str = "euclidean",
    *,
    n_epochs: int | None = None,
) -> SingleCellData:
    """Perform UMAP dimensionality reduction on single-cell data.

    Args:
        data: SingleCellData object with expression matrix
        n_components: Number of dimensions (usually 2 or 3)
        n_neighbors: Number of neighbors for UMAP
        min_dist: Minimum distance parameter
        random_state: Random seed for reproducibility
        metric: Distance metric
        n_epochs: Number of training epochs (optional)

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
    X = data.X.toarray() if hasattr(data.X, "toarray") else data.X

    # Perform UMAP
    umap_kwargs = {
        "n_components": n_components,
        "n_neighbors": n_neighbors,
        "min_dist": min_dist,
        "random_state": random_state,
        "metric": metric,
    }
    if n_epochs is not None:
        umap_kwargs["n_epochs"] = n_epochs

    reducer = umap.UMAP(**umap_kwargs)

    X_umap = reducer.fit_transform(X)

    # Store UMAP results
    umap_coords = pd.DataFrame(
        X_umap,
        index=result.obs.index if result.obs is not None else None,
        columns=[f"UMAP{i+1}" for i in range(n_components)],
    )

    # Add to obs
    if result.obs is None:
        result.obs = umap_coords
    else:
        for col in umap_coords.columns:
            result.obs[col] = umap_coords[col]

    # Store UMAP metadata
    result.uns["umap"] = {
        "n_components": n_components,
        "n_neighbors": n_neighbors,
        "min_dist": min_dist,
        "random_state": random_state,
        "metric": metric,
    }

    logger.info("UMAP completed")
    return result


def diffusion_map_reduction(
    data: SingleCellData, n_components: int = 10, n_neighbors: int = 15, alpha: float = 1.0
) -> SingleCellData:
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
    X = data.X.toarray() if hasattr(data.X, "toarray") else data.X

    # Construct kNN graph and compute diffusion maps
    # Simplified implementation of diffusion maps

    # 1. Compute pairwise distances
    from scipy.spatial.distance import pdist, squareform

    distances = squareform(pdist(X, metric="euclidean"))

    # 2. Construct kernel matrix (Gaussian kernel)
    sigma = np.median(distances)  # Adaptive kernel width
    kernel = np.exp(-(distances**2) / (2 * sigma**2))

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
    diffusion_coords = eigenvecs[:, 1 : n_components + 1] / eigenvals[1 : n_components + 1]

    # Store diffusion map results
    dm_coords = pd.DataFrame(
        diffusion_coords,
        index=result.obs.index if result.obs is not None else None,
        columns=[f"DC{i+1}" for i in range(n_components)],
    )

    # Add to obs
    if result.obs is None:
        result.obs = dm_coords
    else:
        for col in dm_coords.columns:
            result.obs[col] = dm_coords[col]

    # Store diffusion map metadata
    result.uns["diffusion_map"] = {
        "n_components": n_components,
        "n_neighbors": n_neighbors,
        "alpha": alpha,
        "sigma": sigma,
        "eigenvalues": eigenvals[1 : n_components + 1].tolist(),
    }

    logger.info(f"Diffusion map completed: {n_components} components computed")
    return result


def mds_reduction(
    data: SingleCellData, n_components: int = 2, metric: bool = True, random_state: int | None = None
) -> SingleCellData:
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
    X = data.X.toarray() if hasattr(data.X, "toarray") else data.X

    # Compute pairwise distances
    from sklearn.metrics.pairwise import euclidean_distances

    distances = euclidean_distances(X)

    # Perform MDS
    from sklearn.manifold import MDS

    mds = MDS(n_components=n_components, metric=metric, random_state=random_state, max_iter=300, eps=1e-6)

    X_mds = mds.fit_transform(distances)

    # Store MDS results
    mds_coords = pd.DataFrame(
        X_mds,
        index=result.obs.index if result.obs is not None else None,
        columns=[f"MDS{i+1}" for i in range(n_components)],
    )

    # Add to obs
    if result.obs is None:
        result.obs = mds_coords
    else:
        for col in mds_coords.columns:
            result.obs[col] = mds_coords[col]

    # Store MDS metadata
    result.uns["mds"] = {
        "n_components": n_components,
        "metric": metric,
        "random_state": random_state,
        "stress": float(mds.stress_) if hasattr(mds, "stress_") else None,
        "n_iter": mds.n_iter_,
    }

    logger.info(f"MDS completed: stress = {result.uns['mds'].get('stress', 'N/A')}")
    return result


def ica_reduction(
    data: SingleCellData, n_components: int = 10, random_state: int | None = None, max_iter: int = 1000
) -> SingleCellData:
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
    X = data.X.toarray() if hasattr(data.X, "toarray") else data.X

    # Scale data
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Perform ICA
    from sklearn.decomposition import FastICA

    ica = FastICA(n_components=n_components, random_state=random_state, max_iter=max_iter, tol=1e-4)

    X_ica = ica.fit_transform(X_scaled)

    # Store ICA results
    ica_coords = pd.DataFrame(
        X_ica,
        index=result.obs.index if result.obs is not None else None,
        columns=[f"IC{i+1}" for i in range(n_components)],
    )

    # Add to obs
    if result.obs is None:
        result.obs = ica_coords
    else:
        for col in ica_coords.columns:
            result.obs[col] = ica_coords[col]

    # Store ICA metadata and components
    result.uns["ica"] = {
        "n_components": n_components,
        "random_state": random_state,
        "max_iter": max_iter,
        "n_iter": ica.n_iter_,
    }

    # Store mixing matrix for gene loadings
    result.varm = result.varm if hasattr(result, "varm") and result.varm is not None else {}
    result.varm["ICs"] = ica.mixing_.T

    logger.info(f"ICA completed: {n_components} components extracted")
    return result


def factor_analysis_reduction(
    data: SingleCellData, n_components: int = 10, random_state: int | None = None
) -> SingleCellData:
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
    X = data.X.toarray() if hasattr(data.X, "toarray") else data.X

    # Scale data
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Perform Factor Analysis
    from sklearn.decomposition import FactorAnalysis

    fa = FactorAnalysis(n_components=n_components, random_state=random_state, max_iter=1000)

    X_fa = fa.fit_transform(X_scaled)

    # Store FA results
    fa_coords = pd.DataFrame(
        X_fa,
        index=result.obs.index if result.obs is not None else None,
        columns=[f"FA{i+1}" for i in range(n_components)],
    )

    # Add to obs
    if result.obs is None:
        result.obs = fa_coords
    else:
        for col in fa_coords.columns:
            result.obs[col] = fa_coords[col]

    # Store FA metadata and components
    result.uns["factor_analysis"] = {
        "n_components": n_components,
        "random_state": random_state,
        "log_likelihood": float(fa.loglike_[-1]) if hasattr(fa, "loglike_") else None,
        "noise_variance": fa.noise_variance_.tolist(),
    }

    # Store factor loadings
    result.varm = result.varm if hasattr(result, "varm") and result.varm is not None else {}
    result.varm["FA_loadings"] = fa.components_.T

    logger.info(f"Factor Analysis completed: {n_components} factors extracted")
    return result


def compute_dimensionality_metrics(data: SingleCellData, embedding_cols: List[str]) -> Dict[str, Any]:
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
    X = data.X.toarray() if hasattr(data.X, "toarray") else data.X
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
        "distance_correlation": float(distance_correlation),
        "embedding_dimensions": len(embedding_cols),
        "n_samples": X.shape[0],
        "n_features": X.shape[1],
    }

    # Trustworthiness and continuity (simplified)
    # These are computationally expensive for large datasets
    if X.shape[0] <= 1000:
        metrics.update(_compute_trustworthiness_continuity(original_distances, embedding_distances))

    logger.info(f"Dimensionality metrics computed: distance correlation = {distance_correlation:.3f}")
    return metrics


def _compute_trustworthiness_continuity(
    original_distances: np.ndarray, embedding_distances: np.ndarray, k: int = 12
) -> Dict[str, float]:
    """Compute trustworthiness and continuity metrics (simplified version)."""
    n_samples = original_distances.shape[0]

    # For each point, find k nearest neighbors in both spaces
    trustworthiness = 0.0
    continuity = 0.0

    for i in range(n_samples):
        # Original space neighbors
        orig_neighbors = np.argsort(original_distances[i, :])[1 : k + 1]  # Exclude self

        # Embedding space neighbors
        embed_neighbors = np.argsort(embedding_distances[i, :])[1 : k + 1]

        # Trustworthiness: fraction of embedding neighbors that are original neighbors
        trust_set = set(orig_neighbors)
        embed_set = set(embed_neighbors)
        trustworthiness += len(trust_set.intersection(embed_set)) / k

        # Continuity: fraction of original neighbors that are embedding neighbors
        continuity += len(embed_set.intersection(trust_set)) / k

    trustworthiness /= n_samples
    continuity /= n_samples

    return {
        "trustworthiness": float(trustworthiness),
        "continuity": float(continuity),
        "neighborhood_size": k,
    }


# Check for sklearn availability
try:
    import sklearn

    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False


def run_pca(
    data: SingleCellData, n_components: int = 50, random_state: int | None = None, scale_data: bool = True
) -> SingleCellData:
    """Alias for pca_reduction - run PCA dimensionality reduction.

    Args:
        data: SingleCellData object
        n_components: Number of components
        random_state: Random seed
        scale_data: Whether to scale data

    Returns:
        SingleCellData with PCA results
    """
    return pca_reduction(data, n_components, random_state, scale_data)


def run_tsne(
    data: SingleCellData,
    n_components: int = 2,
    perplexity: float = 30.0,
    random_state: int | None = None,
    learning_rate: float = 200.0,
    max_iter: int = 1000,
) -> SingleCellData:
    """Alias for tsne_reduction - run t-SNE dimensionality reduction.

    Args:
        data: SingleCellData object
        n_components: Number of components
        perplexity: t-SNE perplexity
        random_state: Random seed
        learning_rate: Learning rate
        max_iter: Maximum iterations

    Returns:
        SingleCellData with t-SNE results
    """
    return tsne_reduction(data, n_components, perplexity, random_state, learning_rate, max_iter)


def run_umap(
    data: SingleCellData,
    n_components: int = 2,
    n_neighbors: int = 15,
    min_dist: float = 0.1,
    random_state: int | None = None,
    metric: str = "euclidean",
) -> SingleCellData:
    """Alias for umap_reduction - run UMAP dimensionality reduction.

    Args:
        data: SingleCellData object
        n_components: Number of components
        n_neighbors: Number of neighbors
        min_dist: Minimum distance
        random_state: Random seed
        metric: Distance metric

    Returns:
        SingleCellData with UMAP results
    """
    return umap_reduction(data, n_components, n_neighbors, min_dist, random_state, metric)


def compute_pca(
    data: SingleCellData,
    n_components: int = 50,
    random_state: int | None = None,
    scale_data: bool = True,
    *,
    use_hvgs: bool = False,
) -> SingleCellData:
    """Compute PCA dimensionality reduction (alias for pca_reduction).

    Args:
        data: SingleCellData object
        n_components: Number of components
        random_state: Random seed
        scale_data: Whether to scale data
        use_hvgs: Whether to use only highly variable genes

    Returns:
        SingleCellData with PCA results
    """
    # Note: use_hvgs is just a flag - the data should already have HVGs marked
    # This parameter exists for API compatibility
    return pca_reduction(data, n_components, random_state, scale_data, use_hvgs=use_hvgs)


def compute_tsne(
    data: SingleCellData,
    n_components: int = 2,
    perplexity: float = 30.0,
    random_state: int | None = None,
    learning_rate: float = 200.0,
    max_iter: int = 1000,
    *,
    n_iter: int | None = None,
) -> SingleCellData:
    """Compute t-SNE dimensionality reduction (alias for tsne_reduction).

    Args:
        data: SingleCellData object
        n_components: Number of components
        perplexity: t-SNE perplexity
        random_state: Random seed
        learning_rate: Learning rate
        max_iter: Maximum iterations
        n_iter: Alias for max_iter

    Returns:
        SingleCellData with t-SNE results
    """
    # Handle alias
    if n_iter is not None:
        max_iter = n_iter
    return tsne_reduction(data, n_components, perplexity, random_state, learning_rate, max_iter)


def compute_umap(
    data: SingleCellData,
    n_components: int = 2,
    n_neighbors: int = 15,
    min_dist: float = 0.1,
    random_state: int | None = None,
    metric: str = "euclidean",
    *,
    n_epochs: int | None = None,
) -> SingleCellData:
    """Compute UMAP dimensionality reduction (alias for umap_reduction).

    Args:
        data: SingleCellData object
        n_components: Number of components
        n_neighbors: Number of neighbors
        min_dist: Minimum distance
        random_state: Random seed
        metric: Distance metric
        n_epochs: Number of training epochs (optional)

    Returns:
        SingleCellData with UMAP results
    """
    return umap_reduction(data, n_components, n_neighbors, min_dist, random_state, metric, n_epochs=n_epochs)


def compute_neighbors(
    data: SingleCellData,
    n_neighbors: int = 15,
    metric: str = "euclidean",
    use_rep: str = "X_pca",
    *,
    n_pcs: int | None = None,
) -> SingleCellData:
    """Compute neighbor graph for single-cell data.

    Args:
        data: SingleCellData object
        n_neighbors: Number of neighbors
        metric: Distance metric
        use_rep: Representation to use for neighbor calculation
        n_pcs: Number of principal components to use (optional)

    Returns:
        SingleCellData with neighbor graph
    """
    if not SKLEARN_AVAILABLE:
        raise ImportError("scikit-learn required for neighbor computation")

    from sklearn.neighbors import NearestNeighbors

    # Get the representation to use
    if use_rep == "X_pca" and hasattr(data, "obsm") and "X_pca" in data.obsm:
        X = data.obsm["X_pca"]
        # Limit to n_pcs if specified
        if n_pcs is not None and X.shape[1] > n_pcs:
            X = X[:, :n_pcs]
    elif hasattr(data, "X"):
        X = data.X
        if hasattr(X, "toarray"):
            X = X.toarray()
    else:
        raise ValueError(f"Representation {use_rep} not found in data")

    # Compute neighbors
    nbrs = NearestNeighbors(n_neighbors=n_neighbors, metric=metric)
    nbrs.fit(X)

    distances, indices = nbrs.kneighbors(X)

    # Build connectivity matrix (sparse)
    from scipy import sparse

    n_cells = X.shape[0]
    rows = np.repeat(np.arange(n_cells), n_neighbors)
    cols = indices.flatten()
    # Connectivity: convert distances to similarities using Gaussian kernel
    # sigma = median distance
    median_dist = np.median(distances[distances > 0]) if np.any(distances > 0) else 1.0
    data_flat = np.exp(-(distances.flatten() ** 2) / (2 * median_dist**2))
    connectivities = sparse.csr_matrix((data_flat, (rows, cols)), shape=(n_cells, n_cells))

    # Symmetrize the connectivity matrix (average of A and A.T)
    connectivities = (connectivities + connectivities.T) / 2

    # Store results
    data.uns["neighbors"] = {
        "indices": indices,
        "distances": distances,
        "connectivities": connectivities,
        "params": {"n_neighbors": n_neighbors, "metric": metric, "use_rep": use_rep},
    }

    logger.info(f"Computed neighbor graph with {n_neighbors} neighbors for {len(X)} cells")
    return data


def compute_diffusion_map(
    data: SingleCellData, n_components: int = 10, n_neighbors: int = 15, alpha: float = 1.0
) -> SingleCellData:
    """Compute diffusion map for single-cell data.

    Args:
        data: SingleCellData object
        n_components: Number of diffusion components
        n_neighbors: Number of neighbors for graph construction
        alpha: Normalization parameter

    Returns:
        SingleCellData with diffusion map coordinates
    """
    if not SKLEARN_AVAILABLE:
        raise ImportError("scikit-learn required for diffusion maps")

    # First compute neighbors if not already done
    if "neighbors" not in data.uns:
        data = compute_neighbors(data, n_neighbors)

    # Get neighbor information
    indices = data.uns["neighbors"]["indices"]
    distances = data.uns["neighbors"]["distances"]

    n_cells = len(data)
    # Create adjacency matrix from neighbors
    from scipy import sparse

    rows, cols = [], []
    data_weights = []

    for i in range(n_cells):
        for j_idx, j in enumerate(indices[i]):
            if i != j:  # Exclude self
                rows.append(i)
                cols.append(j)
                # Use distance-based weights
                weight = 1.0 / (1.0 + distances[i, j_idx])
                data_weights.append(weight)

    # Create symmetric adjacency matrix
    adjacency = sparse.csr_matrix((data_weights + data_weights, (rows + cols, cols + rows)), shape=(n_cells, n_cells))

    # Normalize adjacency matrix (create diffusion operator)
    degrees = np.array(adjacency.sum(axis=1)).flatten()
    degrees[degrees == 0] = 1  # Avoid division by zero

    if alpha == 0:
        # Unnormalized Laplacian
        diffusion_operator = adjacency
    elif alpha == 1:
        # Normalized Laplacian
        D_inv_sqrt = sparse.diags(1.0 / np.sqrt(degrees))
        diffusion_operator = D_inv_sqrt @ adjacency @ D_inv_sqrt
    else:
        # Generalized normalization
        D_alpha = sparse.diags(degrees ** (-alpha))
        diffusion_operator = D_alpha @ adjacency @ D_alpha

    # Compute eigenvectors (diffusion components)
    from scipy.sparse.linalg import eigsh

    eigenvalues, eigenvectors = eigsh(diffusion_operator, k=n_components + 1, which="LM", sigma=1.0)

    # Sort by eigenvalue magnitude (largest first)
    idx = np.argsort(np.abs(eigenvalues))[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    # Store results (skip first eigenvector which is constant)
    data.obsm["X_diffmap"] = eigenvectors[:, 1 : n_components + 1]
    data.uns["diffmap"] = {
        "eigenvalues": eigenvalues[1 : n_components + 1],
        "params": {"n_components": n_components, "n_neighbors": n_neighbors, "alpha": alpha},
    }

    logger.info(f"Computed diffusion map with {n_components} components")
    return data


def select_hvgs(
    data: SingleCellData,
    n_top_genes: int = 2000,
    flavor: str | None = None,
    *,
    method: str | None = None,
) -> SingleCellData:
    """Select highly variable genes from single-cell data.

    Args:
        data: SingleCellData object
        n_top_genes: Number of top variable genes to select
        flavor: Method for HVG selection ('seurat', 'cell_ranger', 'seurat_v3', 'variance')
        method: Alias for flavor

    Returns:
        SingleCellData with highly variable genes marked
    """
    # Handle parameter alias
    if method is not None and flavor is None:
        flavor = method
    if flavor is None:
        flavor = "seurat"

    if not hasattr(data, "X") or data.X is None:
        raise ValueError("Data must contain expression matrix X")

    X = data.X
    if hasattr(X, "toarray"):  # sparse matrix
        X = X.toarray()

    n_cells, n_genes = X.shape

    if flavor == "seurat":
        # Seurat-style HVG selection
        # Calculate mean and variance for each gene
        gene_means = np.mean(X, axis=0)
        gene_vars = np.var(X, axis=0)

        # Fit curve: variance = a * mean^b
        # Use only genes with mean > 0 for fitting
        nonzero_mask = gene_means > 0
        if np.sum(nonzero_mask) < 10:
            # Fallback: just select by variance
            top_indices = np.argsort(gene_vars)[::-1][:n_top_genes]
        else:
            log_means = np.log10(gene_means[nonzero_mask])
            log_vars = np.log10(gene_vars[nonzero_mask])

            # Simple linear fit
            coeffs = np.polyfit(log_means, log_vars, 1)

            # Only apply to genes with non-zero mean
            predicted_vars = np.zeros(n_genes)
            predicted_vars[nonzero_mask] = 10 ** (coeffs[0] * np.log10(gene_means[nonzero_mask]) + coeffs[1])

            # For zero-mean genes, set predicted variance to small value
            predicted_vars[~nonzero_mask] = 1e-10

            # Calculate standardized variance
            standardized_vars = gene_vars / np.maximum(predicted_vars, 1e-10)

            # Select top genes by standardized variance (only consider valid genes)
            valid_mask = gene_means > 0
            valid_std_vars = standardized_vars[valid_mask]
            sorted_indices = np.argsort(valid_std_vars)[::-1]

            # Map back to original gene indices
            valid_gene_indices = np.where(valid_mask)[0]
            n_valid = min(n_top_genes, len(valid_gene_indices))
            top_indices = valid_gene_indices[sorted_indices[:n_valid]]

            # If we need more, add zero-mean genes by variance
            if n_valid < n_top_genes:
                zero_mean_indices = np.where(~valid_mask)[0]
                zero_mean_vars = gene_vars[~valid_mask]
                remaining = n_top_genes - n_valid
                additional = zero_mean_indices[np.argsort(zero_mean_vars)[::-1][:remaining]]
                top_indices = np.concatenate([top_indices, additional])

    elif flavor == "cell_ranger":
        # Cell Ranger-style HVG selection
        # Use coefficient of variation squared
        gene_means = np.mean(X, axis=0)
        gene_vars = np.var(X, axis=0)

        # Avoid division by zero
        gene_means = np.maximum(gene_means, 1e-10)

        cv_squared = gene_vars / (gene_means**2)

        # Select top genes by CV²
        top_indices = np.argsort(cv_squared)[::-1][:n_top_genes]

    elif flavor == "variance":
        # Simple variance-based selection
        gene_vars = np.var(X, axis=0)
        top_indices = np.argsort(gene_vars)[::-1][:n_top_genes]

    else:
        raise ValueError(f"Unknown HVG selection method: {flavor}. Use 'seurat', 'cell_ranger', or 'variance'")

    # Mark highly variable genes
    highly_variable = np.zeros(n_genes, dtype=bool)
    highly_variable[top_indices] = True

    # Calculate gene statistics for all methods
    gene_means_all = np.mean(X, axis=0)
    gene_vars_all = np.var(X, axis=0)
    # Dispersions as coefficient of variation
    gene_means_safe = np.maximum(gene_means_all, 1e-10)
    dispersions = gene_vars_all / (gene_means_safe**2)
    # Normalized variance (z-score of variance)
    var_mean = np.mean(gene_vars_all)
    var_std = np.std(gene_vars_all)
    if var_std > 0:
        variances_norm = (gene_vars_all - var_mean) / var_std
    else:
        variances_norm = np.zeros_like(gene_vars_all)

    if hasattr(data, "var") and data.var is not None:
        data.var["highly_variable"] = highly_variable
        data.var["means"] = gene_means_all
        data.var["variances"] = gene_vars_all
        data.var["variances_norm"] = variances_norm
        data.var["dispersions"] = dispersions
    else:
        # Create var dataframe if it doesn't exist
        import pandas as pd

        var_df = pd.DataFrame(index=range(n_genes))
        var_df["highly_variable"] = highly_variable
        var_df["means"] = gene_means_all
        var_df["variances"] = gene_vars_all
        var_df["variances_norm"] = variances_norm
        var_df["dispersions"] = dispersions
        data.var = var_df

    data.uns["hvg"] = {"flavor": flavor, "n_top_genes": n_top_genes, "n_selected": len(top_indices)}

    logger.info(f"Selected {len(top_indices)} highly variable genes using {flavor} method")
    return data
