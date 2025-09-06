"""Dimensionality reduction methods for biological data."""

from __future__ import annotations

import warnings
from typing import Any, Dict, List, Optional, Tuple

import numpy as np


def reduce_dimensions_pca(
    X: np.ndarray, n_components: int = 50, standardize: bool = True
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Principal Component Analysis for dimensionality reduction.

    Args:
        X: Data matrix (samples x features)
        n_components: Number of components to keep
        standardize: Whether to standardize features

    Returns:
        Tuple of (transformed_data, components, explained_variance)
    """
    if standardize:
        X_centered = X - np.mean(X, axis=0, keepdims=True)
        X_std = np.std(X, axis=0, keepdims=True) + 1e-8
        X_scaled = X_centered / X_std
    else:
        X_scaled = X - np.mean(X, axis=0, keepdims=True)

    # Compute covariance matrix
    cov_matrix = np.cov(X_scaled.T)

    # Eigendecomposition
    eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)

    # Sort by eigenvalues (descending)
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    # Select top components
    n_components = min(n_components, len(eigenvalues))
    selected_components = eigenvectors[:, :n_components]
    selected_eigenvalues = eigenvalues[:n_components]

    # Transform data
    X_transformed = X_scaled @ selected_components

    # Explained variance ratio
    explained_variance = selected_eigenvalues / np.sum(eigenvalues)

    return X_transformed, selected_components, explained_variance


def reduce_dimensions_umap(
    X: np.ndarray,
    n_components: int = 2,
    n_neighbors: int = 15,
    min_dist: float = 0.1,
    random_state: Optional[int] = None,
) -> np.ndarray:
    """Simplified UMAP-like dimensionality reduction.

    Args:
        X: Data matrix
        n_components: Number of output dimensions
        n_neighbors: Number of neighbors for local structure
        min_dist: Minimum distance between points in low-dim space
        random_state: Random seed

    Returns:
        Low-dimensional embedding
    """
    if random_state is not None:
        np.random.seed(random_state)

    n_samples, n_features = X.shape

    # Start with PCA initialization
    X_pca, _, _ = reduce_dimensions_pca(X, n_components=n_components)

    # Simple optimization (very simplified UMAP)
    embedding = X_pca.copy()

    # Build k-nearest neighbor graph
    distances = np.zeros((n_samples, n_samples))
    for i in range(n_samples):
        for j in range(i + 1, n_samples):
            dist = np.linalg.norm(X[i] - X[j])
            distances[i, j] = dist
            distances[j, i] = dist

    # Simple force-based layout adjustment
    learning_rate = 1.0
    for iteration in range(100):
        forces = np.zeros_like(embedding)

        for i in range(n_samples):
            # Find k nearest neighbors in original space
            neighbor_indices = np.argsort(distances[i])[: n_neighbors + 1][1:]  # Exclude self

            for j in neighbor_indices:
                # Attractive force for neighbors
                diff = embedding[j] - embedding[i]
                dist = np.linalg.norm(diff) + 1e-8

                # Target distance based on original distance
                target_dist = max(min_dist, distances[i, j] / np.max(distances) * 5.0)

                force = (target_dist - dist) / dist * diff * 0.1
                forces[i] += force

        # Update positions
        embedding += learning_rate * forces
        learning_rate *= 0.99  # Decay learning rate

    return embedding


def reduce_dimensions_tsne(
    X: np.ndarray, n_components: int = 2, perplexity: float = 30.0, random_state: Optional[int] = None
) -> np.ndarray:
    """Simplified t-SNE dimensionality reduction.

    Args:
        X: Data matrix
        n_components: Number of output dimensions
        perplexity: t-SNE perplexity parameter
        random_state: Random seed

    Returns:
        Low-dimensional embedding
    """
    if random_state is not None:
        np.random.seed(random_state)

    n_samples = X.shape[0]

    # Adjust perplexity if too large for dataset
    perplexity = min(perplexity, (n_samples - 1) / 3.0)

    # Initialize embedding randomly
    embedding = np.random.randn(n_samples, n_components) * 0.01

    # Compute pairwise distances in high-dimensional space
    distances_hd = np.zeros((n_samples, n_samples))
    for i in range(n_samples):
        for j in range(i + 1, n_samples):
            dist = np.linalg.norm(X[i] - X[j])
            distances_hd[i, j] = dist
            distances_hd[j, i] = dist

    # Convert distances to probabilities (simplified)
    # In real t-SNE, this involves binary search for optimal sigma
    sigma = np.std(distances_hd) / 2.0
    P = np.exp(-(distances_hd**2) / (2 * sigma**2))
    np.fill_diagonal(P, 0)
    P = P / np.sum(P)  # Normalize

    # Optimization loop (simplified)
    learning_rate = 200.0
    momentum = 0.8
    velocity = np.zeros_like(embedding)

    for iteration in range(300):
        # Compute pairwise distances in low-dimensional space
        distances_ld = np.zeros((n_samples, n_samples))
        for i in range(n_samples):
            for j in range(i + 1, n_samples):
                dist = np.linalg.norm(embedding[i] - embedding[j])
                distances_ld[i, j] = dist
                distances_ld[j, i] = dist

        # Convert to probabilities using Student-t distribution
        Q = 1.0 / (1.0 + distances_ld**2)
        np.fill_diagonal(Q, 0)
        Q = Q / (np.sum(Q) + 1e-12)

        # Compute gradients
        gradients = np.zeros_like(embedding)
        for i in range(n_samples):
            for j in range(n_samples):
                if i != j:
                    factor = (P[i, j] - Q[i, j]) * Q[i, j]
                    gradients[i] += factor * (embedding[i] - embedding[j])

        # Update embedding with momentum
        velocity = momentum * velocity - learning_rate * gradients
        embedding += velocity

        # Decay learning rate
        if iteration == 100:
            learning_rate = 50.0

    return embedding


def biological_embedding(
    X: np.ndarray, method: str = "pca", n_components: int = 2, feature_weights: Optional[np.ndarray] = None, **kwargs
) -> Dict[str, Any]:
    """Create biological-aware dimensionality reduction.

    Args:
        X: Data matrix
        method: Reduction method ("pca", "umap", "tsne")
        n_components: Number of output dimensions
        feature_weights: Weights for features (e.g., biological importance)
        **kwargs: Method-specific parameters

    Returns:
        Dictionary with embedding and metadata
    """
    # Apply feature weights if provided
    if feature_weights is not None:
        if len(feature_weights) != X.shape[1]:
            raise ValueError("Feature weights length must match number of features")
        X_weighted = X * feature_weights
    else:
        X_weighted = X

    # Apply dimensionality reduction
    if method == "pca":
        embedding, components, explained_var = reduce_dimensions_pca(X_weighted, n_components=n_components, **kwargs)
        metadata = {
            "components": components,
            "explained_variance": explained_var,
            "cumulative_variance": np.cumsum(explained_var),
        }

    elif method == "umap":
        embedding = reduce_dimensions_umap(X_weighted, n_components=n_components, **kwargs)
        metadata = {"method": "umap"}

    elif method == "tsne":
        embedding = reduce_dimensions_tsne(X_weighted, n_components=n_components, **kwargs)
        metadata = {"method": "tsne"}

    else:
        raise ValueError(f"Unknown method: {method}")

    return {
        "embedding": embedding,
        "method": method,
        "n_components": n_components,
        "feature_weights_used": feature_weights is not None,
        **metadata,
    }
