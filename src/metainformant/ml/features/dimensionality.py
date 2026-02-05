"""Dimensionality reduction utilities for METAINFORMANT.

This module provides dimensionality reduction methods specifically designed
for biological data analysis, including PCA, ICA, UMAP, and t-SNE.
"""

from __future__ import annotations

import numpy as np
from typing import Any, Dict, List, Optional, Tuple, Union

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Optional imports for dimensionality reduction
try:
    from sklearn.decomposition import PCA, FastICA
    from sklearn.preprocessing import StandardScaler

    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False
    logger.warning("scikit-learn not available, dimensionality reduction disabled")

# Optional imports for advanced methods
try:
    import umap

    HAS_UMAP = True
except ImportError:
    HAS_UMAP = False
    logger.warning("umap-learn not available, UMAP disabled")

try:
    from sklearn.manifold import TSNE

    HAS_TSNE = True
except ImportError:
    HAS_TSNE = False
    logger.warning("sklearn.manifold not available, t-SNE disabled")


def pca_reduction(
    X: np.ndarray,
    n_components: int | None = None,
    scale_data: bool = True,
    random_state: int | None = None,
    **kwargs: Any,
) -> Tuple[np.ndarray, PCA]:
    """Perform PCA dimensionality reduction.

    Args:
        X: Input data matrix (n_samples, n_features)
        n_components: Number of components to keep
        scale_data: Whether to standardize data before PCA
        random_state: Random state for reproducibility
        **kwargs: Additional parameters for PCA

    Returns:
        Tuple of (transformed_data, pca_model)

    Raises:
        ImportError: If scikit-learn not available
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn required for PCA")

    # Scale data if requested
    if scale_data:
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
    else:
        X_scaled = X

    # Perform PCA
    pca = PCA(n_components=n_components, random_state=random_state, **kwargs)
    X_pca = pca.fit_transform(X_scaled)

    # Log explained variance
    explained_var = pca.explained_variance_ratio_
    cumulative_var = np.cumsum(explained_var)

    logger.info(
        f"PCA reduction: {X.shape[1]} → {X_pca.shape[1]} dimensions, " f"explained variance: {cumulative_var[-1]:.3f}"
    )

    return X_pca, pca


def ica_reduction(
    X: np.ndarray,
    n_components: int | None = None,
    scale_data: bool = True,
    random_state: int | None = None,
    **kwargs: Any,
) -> Tuple[np.ndarray, FastICA]:
    """Perform ICA dimensionality reduction.

    Args:
        X: Input data matrix (n_samples, n_features)
        n_components: Number of components to extract
        scale_data: Whether to standardize data before ICA
        random_state: Random state for reproducibility
        **kwargs: Additional parameters for FastICA

    Returns:
        Tuple of (transformed_data, ica_model)

    Raises:
        ImportError: If scikit-learn not available
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn required for ICA")

    # Scale data if requested
    if scale_data:
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
    else:
        X_scaled = X

    # Perform ICA
    ica = FastICA(n_components=n_components, random_state=random_state, **kwargs)
    X_ica = ica.fit_transform(X_scaled)

    logger.info(f"ICA reduction: {X.shape[1]} → {X_ica.shape[1]} dimensions")

    return X_ica, ica


def umap_reduction(
    X: np.ndarray,
    n_components: int = 2,
    n_neighbors: int = 15,
    min_dist: float = 0.1,
    metric: str = "euclidean",
    random_state: int | None = None,
    **kwargs: Any,
) -> np.ndarray:
    """Perform UMAP dimensionality reduction.

    Args:
        X: Input data matrix (n_samples, n_features)
        n_components: Number of dimensions to reduce to
        n_neighbors: Number of neighbors for UMAP
        min_dist: Minimum distance between points in embedding
        metric: Distance metric to use
        random_state: Random state for reproducibility
        **kwargs: Additional parameters for UMAP

    Returns:
        Reduced dimensionality data

    Raises:
        ImportError: If umap-learn not available
    """
    if not HAS_UMAP:
        raise ImportError("umap-learn required for UMAP")

    # Perform UMAP
    reducer = umap.UMAP(
        n_components=n_components,
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        metric=metric,
        random_state=random_state,
        **kwargs,
    )

    X_umap = reducer.fit_transform(X)

    logger.info(
        f"UMAP reduction: {X.shape[1]} → {n_components} dimensions " f"(n_neighbors={n_neighbors}, min_dist={min_dist})"
    )

    return X_umap


def tsne_reduction(
    X: np.ndarray,
    n_components: int = 2,
    perplexity: float = 30.0,
    learning_rate: float = 200.0,
    random_state: int | None = None,
    **kwargs: Any,
) -> np.ndarray:
    """Perform t-SNE dimensionality reduction.

    Args:
        X: Input data matrix (n_samples, n_features)
        n_components: Number of dimensions to reduce to
        perplexity: Perplexity parameter for t-SNE
        learning_rate: Learning rate for optimization
        random_state: Random state for reproducibility
        **kwargs: Additional parameters for t-SNE

    Returns:
        Reduced dimensionality data

    Raises:
        ImportError: If sklearn.manifold not available
    """
    if not HAS_TSNE:
        raise ImportError("sklearn.manifold required for t-SNE")

    # Perform t-SNE
    tsne = TSNE(
        n_components=n_components,
        perplexity=perplexity,
        learning_rate=learning_rate,
        random_state=random_state,
        **kwargs,
    )

    X_tsne = tsne.fit_transform(X)

    logger.info(
        f"t-SNE reduction: {X.shape[1]} → {n_components} dimensions "
        f"(perplexity={perplexity}, learning_rate={learning_rate})"
    )

    return X_tsne


def compare_dimensionality_methods(
    X: np.ndarray,
    methods: List[str] = None,
    n_components: int = 2,
    random_state: int | None = None,
) -> Dict[str, Any]:
    """Compare different dimensionality reduction methods.

    Args:
        X: Input data matrix
        methods: List of methods to compare
        n_components: Number of components for reduction
        random_state: Random state for reproducibility

    Returns:
        Dictionary with comparison results
    """
    if methods is None:
        methods = ["pca", "ica"]
        if HAS_UMAP:
            methods.append("umap")
        if HAS_TSNE:
            methods.append("tsne")

    results = {}
    embeddings = {}

    for method in methods:
        try:
            if method == "pca":
                X_reduced, model = pca_reduction(X, n_components=n_components, random_state=random_state)
                explained_var = model.explained_variance_ratio_.sum()

            elif method == "ica":
                X_reduced, model = ica_reduction(X, n_components=n_components, random_state=random_state)
                explained_var = None

            elif method == "umap":
                if not HAS_UMAP:
                    raise ImportError("umap-learn required")
                X_reduced = umap_reduction(X, n_components=n_components, random_state=random_state)
                explained_var = None

            elif method == "tsne":
                if not HAS_TSNE:
                    raise ImportError("sklearn.manifold required")
                X_reduced = tsne_reduction(X, n_components=n_components, random_state=random_state)
                explained_var = None

            else:
                raise ValueError(f"Unknown method: {method}")

            embeddings[method] = X_reduced
            results[method] = {
                "success": True,
                "shape": X_reduced.shape,
                "explained_variance": explained_var,
            }

        except Exception as e:
            logger.error(f"Method {method} failed: {e}")
            results[method] = {
                "success": False,
                "error": str(e),
            }

    return {
        "comparison": results,
        "embeddings": embeddings,
        "input_shape": X.shape,
        "n_components": n_components,
    }


def optimize_dimensionality_parameters(
    X: np.ndarray,
    method: str = "pca",
    param_grid: Dict[str, List] = None,
    cv_metric: str = "reconstruction_error",
    random_state: int | None = None,
) -> Dict[str, Any]:
    """Optimize parameters for dimensionality reduction.

    Args:
        X: Input data matrix
        method: Dimensionality reduction method
        param_grid: Parameter grid to search
        cv_metric: Metric for cross-validation
        random_state: Random state for reproducibility

    Returns:
        Dictionary with optimization results
    """
    if param_grid is None:
        if method == "pca":
            param_grid = {"n_components": [2, 5, 10, 20, 50]}
        elif method == "umap":
            param_grid = {
                "n_neighbors": [5, 15, 30, 50],
                "min_dist": [0.1, 0.25, 0.5],
            }
        elif method == "tsne":
            param_grid = {
                "perplexity": [5, 15, 30, 50],
                "learning_rate": [10, 100, 1000],
            }

    best_params = None
    best_score = float("inf") if cv_metric == "reconstruction_error" else float("-inf")
    results = []

    # Simple grid search (could be improved with proper CV)
    from itertools import product

    param_names = list(param_grid.keys())
    param_values = list(param_grid.values())

    for param_combination in product(*param_values):
        params = dict(zip(param_names, param_combination))

        try:
            if method == "pca":
                X_reduced, model = pca_reduction(X, random_state=random_state, **params)
                if cv_metric == "reconstruction_error":
                    # Calculate reconstruction error
                    X_reconstructed = model.inverse_transform(X_reduced)
                    score = np.mean((X - X_reconstructed) ** 2)
                else:
                    score = model.explained_variance_ratio_.sum()

            elif method == "umap" and HAS_UMAP:
                X_reduced = umap_reduction(X, random_state=random_state, **params)
                # UMAP doesn't have a direct reconstruction, use embedding quality proxy
                score = -np.mean([np.std(X_reduced[:, i]) for i in range(X_reduced.shape[1])])

            elif method == "tsne" and HAS_TSNE:
                X_reduced = tsne_reduction(X, random_state=random_state, **params)
                # t-SNE doesn't have reconstruction, use embedding spread
                score = np.mean([np.std(X_reduced[:, i]) for i in range(X_reduced.shape[1])])

            else:
                continue

            results.append(
                {
                    "params": params,
                    "score": score,
                }
            )

            # Update best
            if cv_metric == "reconstruction_error":
                if score < best_score:
                    best_score = score
                    best_params = params
            else:
                if score > best_score:
                    best_score = score
                    best_params = params

        except Exception as e:
            logger.warning(f"Parameter combination {params} failed: {e}")

    return {
        "method": method,
        "best_params": best_params,
        "best_score": best_score,
        "cv_metric": cv_metric,
        "all_results": results,
        "n_combinations_tested": len(results),
    }


def biological_dimensionality_analysis(
    X: np.ndarray,
    y: Optional[np.ndarray] = None,
    methods: List[str] = None,
    n_components: int = 2,
    random_state: int | None = None,
) -> Dict[str, Any]:
    """Comprehensive dimensionality reduction analysis for biological data.

    Args:
        X: Input data matrix
        y: Optional labels for supervised analysis
        methods: List of methods to apply
        n_components: Number of components to reduce to
        random_state: Random state for reproducibility

    Returns:
        Dictionary with comprehensive analysis results
    """
    if methods is None:
        methods = ["pca"]
        if HAS_UMAP:
            methods.append("umap")

    results = {
        "input_analysis": {
            "shape": X.shape,
            "data_type": str(X.dtype),
            "has_labels": y is not None,
        },
        "methods": {},
        "comparison": {},
    }

    # Analyze each method
    for method in methods:
        try:
            if method == "pca":
                X_reduced, model = pca_reduction(X, n_components=n_components, random_state=random_state)
                results["methods"][method] = {
                    "embedding": X_reduced,
                    "explained_variance": model.explained_variance_ratio_.tolist(),
                    "cumulative_variance": np.cumsum(model.explained_variance_ratio_).tolist(),
                    "components": model.components_,
                }

            elif method == "ica":
                X_reduced, model = ica_reduction(X, n_components=n_components, random_state=random_state)
                results["methods"][method] = {
                    "embedding": X_reduced,
                    "mixing_matrix": model.mixing_,
                    "unmixing_matrix": model.unmixing_,
                }

            elif method == "umap" and HAS_UMAP:
                X_reduced = umap_reduction(X, n_components=n_components, random_state=random_state)
                results["methods"][method] = {
                    "embedding": X_reduced,
                }

            elif method == "tsne" and HAS_TSNE:
                X_reduced = tsne_reduction(X, n_components=n_components, random_state=random_state)
                results["methods"][method] = {
                    "embedding": X_reduced,
                }

            # Calculate embedding statistics
            embedding = results["methods"][method]["embedding"]
            results["methods"][method]["statistics"] = {
                "mean": embedding.mean(axis=0).tolist(),
                "std": embedding.std(axis=0).tolist(),
                "range": (embedding.min(axis=0), embedding.max(axis=0)),
            }

        except Exception as e:
            logger.error(f"Method {method} failed: {e}")
            results["methods"][method] = {"error": str(e)}

    # Compare methods if multiple available
    if len([m for m in results["methods"] if "embedding" in results["methods"][m]]) > 1:
        comparison = compare_dimensionality_methods(
            X, methods=methods, n_components=n_components, random_state=random_state
        )
        results["comparison"] = comparison

    logger.info(
        f"Biological dimensionality analysis completed: " f"{len(results['methods'])} methods tested on {X.shape} data"
    )

    return results


def biological_embedding(
    sequences: List[Any] | None = None,
    embedding_dim: int = 100,
    method: str = "pca",
    X: np.ndarray | None = None,
    n_components: int | None = None,
    labels: np.ndarray | None = None,
    **kwargs: Any,
) -> Dict[str, Any]:
    """Create biological embeddings from sequences or features.

    Args:
        sequences: List of sequences or feature vectors
        embedding_dim: Dimension of the embedding
        method: Embedding method ('pca', 'umap', 'tsne')
        X: Feature matrix (alternative to sequences)
        n_components: Number of components (alternative to embedding_dim)
        labels: Optional sample labels
        **kwargs: Additional arguments for the embedding method

    Returns:
        Dictionary with embedding results
    """
    # Handle alternative parameter names
    if X is not None:
        data = X
    elif sequences is not None:
        if not sequences:
            raise ValueError("No sequences provided")
        if isinstance(sequences[0], str):
            features = []
            for seq in sequences:
                features.append(_sequence_to_features(seq))
            data = np.array(features)
        else:
            data = np.array(sequences)
    else:
        raise ValueError("Must provide either sequences or X")

    if n_components is not None:
        embedding_dim = n_components

    if data.shape[0] == 0:
        raise ValueError("No valid features extracted")

    # Apply dimensionality reduction
    result: Dict[str, Any] = {"method": method}

    if method == "pca":
        X_reduced, components, explained_var = reduce_dimensions_pca(data, n_components=embedding_dim, **kwargs)
        result["embedding"] = X_reduced
        result["explained_variance"] = explained_var
        result["components"] = components
    elif method == "umap" and HAS_UMAP:
        X_reduced = reduce_dimensions_umap(data, n_components=embedding_dim, **kwargs)
        result["embedding"] = X_reduced
        result["explained_variance"] = None
    elif method == "tsne" and HAS_TSNE:
        X_reduced = reduce_dimensions_tsne(data, n_components=embedding_dim, **kwargs)
        result["embedding"] = X_reduced
        result["explained_variance"] = None
    else:
        raise ValueError(f"Unsupported embedding method: {method}")

    if labels is not None:
        result["labels"] = labels

    return result


def reduce_dimensions_pca(
    X: np.ndarray,
    n_components: int | None = None,
    scale_data: bool = True,
    standardize: bool | None = None,
    **kwargs: Any,
) -> tuple:
    """Reduce dimensions using PCA.

    Args:
        X: Input data matrix
        n_components: Number of components to keep
        scale_data: Whether to scale data before PCA
        standardize: Alias for scale_data
        **kwargs: Additional arguments for PCA

    Returns:
        Tuple of (X_reduced, components, explained_variance_ratio)
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn required for PCA")

    if standardize is not None:
        scale_data = standardize

    if n_components is None:
        n_components = min(X.shape[0], X.shape[1], 50)

    X_reduced, pca_obj = pca_reduction(X, n_components=n_components, scale_data=scale_data, **kwargs)

    # Return components and explained variance alongside reduced data
    components = pca_obj.components_.T if hasattr(pca_obj, "components_") else np.eye(X.shape[1], n_components)
    explained_var = (
        pca_obj.explained_variance_ratio_ if hasattr(pca_obj, "explained_variance_ratio_") else np.zeros(n_components)
    )

    return X_reduced, components, explained_var


def reduce_dimensions_tsne(X: np.ndarray, n_components: int = 2, perplexity: float = 30.0, **kwargs: Any) -> np.ndarray:
    """Reduce dimensions using t-SNE.

    Args:
        X: Input data matrix
        n_components: Number of components to keep
        perplexity: t-SNE perplexity parameter
        **kwargs: Additional arguments for t-SNE

    Returns:
        Reduced dimensionality data
    """
    if not HAS_TSNE:
        raise ImportError("scikit-learn required for t-SNE")

    return tsne_reduction(X, n_components=n_components, perplexity=perplexity, **kwargs)


def reduce_dimensions_umap(X: np.ndarray, n_components: int = 2, n_neighbors: int = 15, **kwargs: Any) -> np.ndarray:
    """Reduce dimensions using UMAP.

    Args:
        X: Input data matrix
        n_components: Number of components to keep
        n_neighbors: UMAP n_neighbors parameter
        **kwargs: Additional arguments for UMAP

    Returns:
        Reduced dimensionality data
    """
    if not HAS_UMAP:
        raise ImportError("umap-learn required for UMAP")

    X_reduced = umap_reduction(X, n_components=n_components, n_neighbors=n_neighbors, **kwargs)
    return X_reduced


def _sequence_to_features(sequence: str) -> List[float]:
    """Convert a sequence to basic feature vector.

    This is a simplified implementation for demonstration.
    In practice, you'd use more sophisticated sequence embeddings.
    """
    if not sequence:
        return [0.0] * 20  # Basic amino acid features

    # Simple k-mer frequencies (k=1)
    features = []
    sequence = sequence.upper()

    # Nucleotide/amino acid counts
    bases = "ATCG" if all(c in "ATCGN" for c in sequence) else "ACDEFGHIKLMNPQRSTVWY"
    for base in bases:
        features.append(sequence.count(base) / len(sequence))

    # Sequence length
    features.append(len(sequence))

    # GC content (if DNA)
    if all(c in "ATCGN" for c in sequence):
        gc_count = sequence.count("G") + sequence.count("C")
        features.append(gc_count / len(sequence) if sequence else 0)

    return features
