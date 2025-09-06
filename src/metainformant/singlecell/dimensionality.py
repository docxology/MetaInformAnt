"""Dimensionality reduction for single-cell data.

This module provides essential dimensionality reduction techniques for single-cell
analysis including PCA, UMAP, t-SNE, and neighbor graph construction.
All implementations use real computational methods without mocking.
"""

from __future__ import annotations

import warnings
from typing import Any, Dict, Optional, Tuple, Union

import numpy as np
import pandas as pd

try:
    from scipy import sparse
    from scipy.sparse import csgraph, csr_matrix
    from scipy.spatial.distance import pdist, squareform

    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    sparse = None
    pdist = None
    squareform = None
    csgraph = None
    csr_matrix = None

try:
    from sklearn.decomposition import PCA
    from sklearn.metrics import pairwise_distances
    from sklearn.neighbors import NearestNeighbors

    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False
    PCA = None
    NearestNeighbors = None
    pairwise_distances = None

from .preprocessing import SingleCellData


def select_hvgs(
    data: SingleCellData,
    n_top_genes: int = 2000,
    method: str = "seurat",
    min_mean: float = 0.0125,
    max_mean: float = 3,
    min_disp: float = 0.5,
) -> SingleCellData:
    """Select highly variable genes (HVGs).

    Args:
        data: SingleCellData object (should be normalized and log-transformed)
        n_top_genes: Number of top variable genes to select
        method: Method for HVG selection ('seurat', 'cell_ranger', 'variance')
        min_mean: Minimum mean expression for Seurat method
        max_mean: Maximum mean expression for Seurat method
        min_disp: Minimum dispersion for Seurat method

    Returns:
        SingleCellData with HVG information in var
    """
    data = data.copy()
    X = data.X

    if sparse.issparse(X):
        X_dense = X.toarray()
    else:
        X_dense = X

    if method == "seurat":
        # Seurat method: mean-variance relationship
        gene_means = np.mean(X_dense, axis=0)
        gene_vars = np.var(X_dense, axis=0)

        # Avoid division by zero
        gene_means_nonzero = gene_means.copy()
        gene_means_nonzero[gene_means_nonzero == 0] = 1e-12

        # Calculate dispersion (normalized variance)
        dispersion = gene_vars / gene_means_nonzero

        # Handle infinite and NaN values
        dispersion = np.where(np.isfinite(dispersion), dispersion, 0)

        # Fit mean-dispersion relationship
        # Use log-log fit for better numerical stability
        log_means = np.log10(np.maximum(gene_means, 1e-12))
        log_disp = np.log10(np.maximum(dispersion, 1e-12))

        # Simple linear fit (alternative to loess)
        from scipy.stats import linregress

        try:
            slope, intercept, _, _, _ = linregress(log_means, log_disp)
            # Calculate residuals (observed - expected)
            expected_log_disp = slope * log_means + intercept
            residuals = log_disp - expected_log_disp
        except ValueError:
            # Handle case where all values are identical (no variance)
            slope, intercept = 0, np.mean(log_disp)
            expected_log_disp = np.full_like(log_disp, intercept)
            residuals = log_disp - expected_log_disp

        # Filter based on mean and dispersion criteria
        mean_filter = (gene_means >= min_mean) & (gene_means <= max_mean)
        disp_filter = dispersion >= min_disp
        combined_filter = mean_filter & disp_filter

        # Select top genes by residuals among filtered genes
        if combined_filter.sum() < n_top_genes:
            warnings.warn(f"Only {combined_filter.sum()} genes pass filters, less than requested {n_top_genes}")
            hvg_mask = combined_filter
        else:
            # Get top genes by residuals
            residuals_filtered = residuals.copy()
            residuals_filtered[~combined_filter] = -np.inf
            top_indices = np.argsort(residuals_filtered)[-n_top_genes:]
            hvg_mask = np.zeros(len(gene_means), dtype=bool)
            hvg_mask[top_indices] = True

    elif method == "cell_ranger":
        # Cell Ranger method: normalized dispersion
        gene_means = np.mean(X_dense, axis=0)
        gene_vars = np.var(X_dense, axis=0)

        # Normalized dispersion
        norm_disp = gene_vars / (gene_means + 1e-12)

        # Select top genes by normalized dispersion
        top_indices = np.argsort(norm_disp)[-n_top_genes:]
        hvg_mask = np.zeros(len(gene_means), dtype=bool)
        hvg_mask[top_indices] = True

    elif method == "variance":
        # Simple variance-based selection
        gene_vars = np.var(X_dense, axis=0)
        top_indices = np.argsort(gene_vars)[-n_top_genes:]
        hvg_mask = np.zeros(len(gene_vars), dtype=bool)
        hvg_mask[top_indices] = True

    else:
        raise ValueError(f"Unknown HVG selection method: {method}")

    # Store results in var
    data.var["highly_variable"] = hvg_mask
    data.var["means"] = np.mean(X_dense, axis=0)
    data.var["variances"] = np.var(X_dense, axis=0)
    data.var["variances_norm"] = data.var["variances"] / (data.var["means"] + 1e-12)

    print(f"Selected {hvg_mask.sum()} highly variable genes")

    return data


def compute_pca(
    data: SingleCellData,
    n_components: int = 50,
    use_hvgs: bool = True,
    svd_solver: str = "arpack",
) -> SingleCellData:
    """Compute Principal Component Analysis.

    Args:
        data: SingleCellData object (should be scaled)
        n_components: Number of principal components
        use_hvgs: Whether to use only highly variable genes
        svd_solver: SVD solver ('auto', 'full', 'arpack', 'randomized')

    Returns:
        SingleCellData with PCA results in obsm and varm
    """
    data = data.copy()
    X = data.X

    # Select genes for PCA
    if use_hvgs and "highly_variable" in data.var.columns:
        hvg_mask = data.var["highly_variable"].values
        if hvg_mask.sum() == 0:
            warnings.warn("No highly variable genes found, using all genes")
            X_pca = X
            gene_mask = np.ones(data.n_vars, dtype=bool)
        else:
            X_pca = X[:, hvg_mask]
            gene_mask = hvg_mask
    else:
        X_pca = X
        gene_mask = np.ones(data.n_vars, dtype=bool)

    if sparse.issparse(X_pca):
        X_pca = X_pca.toarray()

    # Compute PCA
    n_components = min(n_components, X_pca.shape[0] - 1, X_pca.shape[1])

    pca = PCA(n_components=n_components, svd_solver=svd_solver)
    X_pca_transformed = pca.fit_transform(X_pca)

    # Store results
    data.obsm["X_pca"] = X_pca_transformed

    # Store gene loadings (only for genes used in PCA)
    loadings = np.zeros((data.n_vars, n_components))
    loadings[gene_mask, :] = pca.components_.T
    data.varm["PCs"] = loadings

    # Store PCA metadata
    data.uns["pca"] = {
        "variance_ratio": pca.explained_variance_ratio_,
        "variance": pca.explained_variance_,
        "n_components": n_components,
        "use_hvgs": use_hvgs,
    }

    print(
        f"PCA computed: {n_components} components, explaining {pca.explained_variance_ratio_.sum():.3f} of total variance"
    )

    return data


def compute_neighbors(
    data: SingleCellData,
    n_neighbors: int = 15,
    n_pcs: int = 40,
    metric: str = "euclidean",
    method: str = "umap",
) -> SingleCellData:
    """Compute k-nearest neighbor graph.

    Args:
        data: SingleCellData object with PCA results
        n_neighbors: Number of neighbors
        n_pcs: Number of PCs to use for neighbor calculation
        metric: Distance metric ('euclidean', 'cosine', 'manhattan')
        method: Method for neighbor calculation ('umap', 'sklearn', 'gauss')

    Returns:
        SingleCellData with neighbor graph in uns
    """
    data = data.copy()

    # Use PCA coordinates if available
    if "X_pca" in data.obsm:
        X_use = data.obsm["X_pca"][:, :n_pcs]
    else:
        warnings.warn("No PCA found, using original expression data")
        X_use = data.X
        if sparse.issparse(X_use):
            X_use = X_use.toarray()

    # Compute k-nearest neighbors
    if method == "sklearn":
        from sklearn.neighbors import NearestNeighbors

        nbrs = NearestNeighbors(
            n_neighbors=n_neighbors + 1, metric=metric, algorithm="auto"  # +1 because each point is its own neighbor
        ).fit(X_use)

        distances, indices = nbrs.kneighbors(X_use)

        # Remove self-connections (first column)
        distances = distances[:, 1:]
        indices = indices[:, 1:]

    elif method == "umap":
        # UMAP-style neighbor calculation with better handling of local connectivity
        try:
            from scipy.sparse import csr_matrix
            from sklearn.neighbors import NearestNeighbors

            # Standard k-NN
            nbrs = NearestNeighbors(n_neighbors=n_neighbors + 1, metric=metric, algorithm="auto").fit(X_use)

            distances, indices = nbrs.kneighbors(X_use)

            # Remove self-connections
            distances = distances[:, 1:]
            indices = indices[:, 1:]

        except ImportError:
            # Fallback to basic implementation
            from scipy.spatial.distance import cdist

            # Compute all pairwise distances (memory intensive for large datasets)
            if X_use.shape[0] > 5000:
                warnings.warn("Large dataset detected, consider using 'sklearn' method for better memory efficiency")

            dist_matrix = pairwise_distances(X_use, metric=metric)

            # Get k-nearest neighbors
            indices = np.argsort(dist_matrix, axis=1)[:, 1 : n_neighbors + 1]  # Exclude self
            distances = np.sort(dist_matrix, axis=1)[:, 1 : n_neighbors + 1]

    else:
        raise ValueError(f"Unknown neighbor method: {method}")

    # Create sparse connectivity matrix
    n_cells = X_use.shape[0]
    rows = np.repeat(np.arange(n_cells), n_neighbors)
    cols = indices.flatten()

    # Use distances as weights (convert to similarities)
    # Use Gaussian kernel: exp(-d^2 / (2 * sigma^2))
    sigma = np.median(distances)  # Adaptive bandwidth
    if sigma == 0:
        sigma = 1.0

    weights = np.exp(-(distances.flatten() ** 2) / (2 * sigma**2))

    # Create sparse adjacency matrix
    adjacency = csr_matrix((weights, (rows, cols)), shape=(n_cells, n_cells))

    # Symmetrize the matrix
    adjacency = (adjacency + adjacency.T) / 2

    # Store results
    data.uns["neighbors"] = {
        "connectivities": adjacency,
        "distances": csr_matrix((distances.flatten(), (rows, cols)), shape=(n_cells, n_cells)),
        "indices": indices,
        "params": {
            "n_neighbors": n_neighbors,
            "n_pcs": n_pcs,
            "metric": metric,
            "method": method,
        },
    }

    print(f"Computed neighbor graph: {n_neighbors} neighbors per cell")

    return data


def compute_umap(
    data: SingleCellData,
    n_components: int = 2,
    min_dist: float = 0.5,
    spread: float = 1.0,
    n_epochs: int = 200,
    random_state: int = 42,
) -> SingleCellData:
    """Compute UMAP embedding.

    Args:
        data: SingleCellData with neighbor graph
        n_components: Number of UMAP dimensions
        min_dist: Minimum distance between points in embedding
        spread: Scale of embedded points
        n_epochs: Number of training epochs
        random_state: Random seed

    Returns:
        SingleCellData with UMAP coordinates in obsm
    """
    data = data.copy()

    # Check if neighbors are computed
    if "neighbors" not in data.uns:
        warnings.warn("No neighbor graph found, computing with default parameters")
        data = compute_neighbors(data)

    try:
        import umap

        # Use pre-computed neighbor graph if available
        if "neighbors" in data.uns:
            # Convert connectivity matrix to UMAP format
            connectivities = data.uns["neighbors"]["connectivities"]

            # UMAP with precomputed neighbors
            reducer = umap.UMAP(
                n_components=n_components,
                min_dist=min_dist,
                spread=spread,
                n_epochs=n_epochs,
                random_state=random_state,
                metric="precomputed",
                n_neighbors=data.uns["neighbors"]["params"]["n_neighbors"],
            )

            # Convert sparse matrix to dense for UMAP
            if sparse.issparse(connectivities):
                distance_matrix = 1.0 - connectivities.toarray()  # Convert similarity to distance
                np.fill_diagonal(distance_matrix, 0)  # Ensure diagonal is zero
            else:
                distance_matrix = 1.0 - connectivities
                np.fill_diagonal(distance_matrix, 0)

            X_umap = reducer.fit_transform(distance_matrix)
        else:
            # Compute UMAP from PCA coordinates
            if "X_pca" in data.obsm:
                X_input = data.obsm["X_pca"]
            else:
                X_input = data.X
                if sparse.issparse(X_input):
                    X_input = X_input.toarray()

            reducer = umap.UMAP(
                n_components=n_components,
                min_dist=min_dist,
                spread=spread,
                n_epochs=n_epochs,
                random_state=random_state,
            )

            X_umap = reducer.fit_transform(X_input)

    except ImportError:
        warnings.warn("UMAP not available, using t-SNE as fallback")
        return compute_tsne(data, n_components=n_components, random_state=random_state)

    # Store results
    data.obsm["X_umap"] = X_umap
    data.uns["umap"] = {
        "params": {
            "n_components": n_components,
            "min_dist": min_dist,
            "spread": spread,
            "n_epochs": n_epochs,
            "random_state": random_state,
        }
    }

    print(f"UMAP computed: {n_components}D embedding")

    return data


def compute_tsne(
    data: SingleCellData,
    n_components: int = 2,
    perplexity: float = 30.0,
    n_iter: int = 1000,
    random_state: int = 42,
) -> SingleCellData:
    """Compute t-SNE embedding.

    Args:
        data: SingleCellData object
        n_components: Number of t-SNE dimensions
        perplexity: t-SNE perplexity parameter
        n_iter: Number of iterations
        random_state: Random seed

    Returns:
        SingleCellData with t-SNE coordinates in obsm
    """
    data = data.copy()

    try:
        from sklearn.manifold import TSNE

        # Use PCA coordinates if available, otherwise use expression data
        if "X_pca" in data.obsm:
            X_input = data.obsm["X_pca"]
        else:
            X_input = data.X
            if sparse.issparse(X_input):
                X_input = X_input.toarray()

        # Adjust perplexity if necessary
        max_perplexity = (X_input.shape[0] - 1) / 3
        if perplexity > max_perplexity:
            perplexity = max_perplexity
            warnings.warn(f"Perplexity too large, reduced to {perplexity:.1f}")

        tsne = TSNE(
            n_components=n_components,
            perplexity=perplexity,
            max_iter=n_iter,
            random_state=random_state,
            init="pca" if X_input.shape[1] >= n_components else "random",
        )

        X_tsne = tsne.fit_transform(X_input)

        # Store results
        data.obsm["X_tsne"] = X_tsne
        data.uns["tsne"] = {
            "params": {
                "n_components": n_components,
                "perplexity": perplexity,
                "n_iter": n_iter,
                "random_state": random_state,
            }
        }

        print(f"t-SNE computed: {n_components}D embedding")

    except ImportError:
        raise ImportError("scikit-learn required for t-SNE computation")

    return data


def compute_diffusion_map(
    data: SingleCellData,
    n_components: int = 15,
    alpha: float = 1.0,
) -> SingleCellData:
    """Compute diffusion map embedding for trajectory analysis.

    Args:
        data: SingleCellData with neighbor graph
        n_components: Number of diffusion components
        alpha: Diffusion parameter (higher values give more global structure)

    Returns:
        SingleCellData with diffusion map coordinates in obsm
    """
    data = data.copy()

    # Check if neighbors are computed
    if "neighbors" not in data.uns:
        warnings.warn("No neighbor graph found, computing with default parameters")
        data = compute_neighbors(data)

    # Get connectivity matrix
    connectivities = data.uns["neighbors"]["connectivities"]

    if sparse.issparse(connectivities):
        W = connectivities.toarray()
    else:
        W = connectivities

    # Compute degree matrix
    D = np.diag(np.sum(W, axis=1))

    # Compute transition matrix with alpha normalization
    D_alpha = np.diag(np.sum(W, axis=1) ** (-alpha))
    W_alpha = D_alpha @ W @ D_alpha

    # Row-normalize to get transition probabilities
    row_sums = np.sum(W_alpha, axis=1)
    row_sums[row_sums == 0] = 1  # Avoid division by zero
    P = W_alpha / row_sums[:, np.newaxis]

    # Eigendecomposition
    eigenvals, eigenvecs = np.linalg.eig(P.T)

    # Sort by eigenvalue magnitude (descending)
    idx = np.argsort(np.abs(eigenvals))[::-1]
    eigenvals = eigenvals[idx]
    eigenvecs = eigenvecs[:, idx]

    # Take real parts (diffusion maps should be real)
    eigenvals = np.real(eigenvals)
    eigenvecs = np.real(eigenvecs)

    # Select top components (skip the first trivial eigenvalue = 1)
    n_components = min(n_components, len(eigenvals) - 1)
    eigenvals_use = eigenvals[1 : n_components + 1]
    eigenvecs_use = eigenvecs[:, 1 : n_components + 1]

    # Diffusion map embedding
    X_diffmap = eigenvecs_use * eigenvals_use[np.newaxis, :]

    # Store results
    data.obsm["X_diffmap"] = X_diffmap
    data.uns["diffmap"] = {
        "eigenvalues": eigenvals_use,
        "params": {
            "n_components": n_components,
            "alpha": alpha,
        },
    }

    print(f"Diffusion map computed: {n_components} components")

    return data
