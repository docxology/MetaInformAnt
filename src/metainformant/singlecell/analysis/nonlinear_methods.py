"""Nonlinear dimensionality reduction methods for single-cell data.

This module provides t-SNE, UMAP, diffusion maps, MDS, and neighbor
graph computation for single-cell expression data, along with embedding
quality metrics.
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.linalg import eigsh

from metainformant.core.data import validation
from metainformant.core.utils import errors, logging
from metainformant.singlecell.data.preprocessing import SingleCellData

from .pca_methods import _add_embedding_to_obs

# Try to import optional dependencies
try:
    from sklearn.decomposition import PCA

    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False
    PCA = None

try:
    from sklearn.manifold import TSNE

    HAS_TSNE = True
except ImportError:
    HAS_TSNE = False
    TSNE = None

try:
    import umap

    HAS_UMAP = True
except ImportError:
    HAS_UMAP = False

SKLEARN_AVAILABLE = HAS_SKLEARN

logger = logging.get_logger(__name__)

# Eigenproblems up to this many cells use the dense symmetric solver (exact, and fast at this size);
# larger ones ask ARPACK for just the leading eigenpairs.
_DENSE_EIGEN_LIMIT = 1000


def _leading_eigenpairs(operator: Any, k: int) -> tuple[np.ndarray, np.ndarray]:
    """The ``k`` algebraically largest eigenpairs of a symmetric matrix (dense or sparse), largest eigenvalue first.

    Which solver runs depends on the size of the problem alone: the dense symmetric solver for up to
    ``_DENSE_EIGEN_LIMIT`` rows (and whenever ``k`` is close to the number of rows, where ARPACK is not
    applicable), ARPACK's ``eigsh`` for the leading pairs of larger matrices. Nothing relies on a solver
    raising to select the other one.
    """
    n = operator.shape[0]
    if n <= _DENSE_EIGEN_LIMIT or k >= n - 1:
        # (scipy-stubs narrow issparse() to a method-less SparseABC; every sparse format provides toarray())
        dense = operator.toarray() if hasattr(operator, "toarray") else np.asarray(operator)
        values, vectors = np.linalg.eigh(dense)
    else:
        # Fixed start vector: ARPACK's own default is drawn from internal state that changes between calls
        values, vectors = eigsh(operator, k=k, which="LA", v0=np.random.default_rng(0).standard_normal(n))
    order = np.argsort(values)[::-1][:k]
    return values[order], vectors[:, order]


def _canonical_signs(vectors: np.ndarray) -> np.ndarray:
    """Flip each column so that its entry of largest magnitude is positive.

    An eigenvector is only defined up to its sign and LAPACK/ARPACK builds differ in the one they return, so
    the sign is fixed by this convention (the first such entry if several tie) instead of left to the solver.
    """
    if vectors.shape[1] == 0:
        return vectors
    pivot = np.argmax(np.abs(vectors), axis=0)
    signs = np.sign(vectors[pivot, np.arange(vectors.shape[1])])
    signs[signs == 0] = 1.0
    return vectors * signs


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
        ValidationError: If data is not SingleCellData
        ValueError: If a parameter is invalid (raised by scikit-learn as ``InvalidParameterError``)
    """
    if not HAS_TSNE:
        raise ImportError("t-SNE requires scikit-learn")

    validation.validate_type(data, SingleCellData, "data")

    # Create copy
    result = data.copy()

    # Get expression matrix
    X = data.X.toarray() if hasattr(data.X, "toarray") else data.X

    # Automatically adjust perplexity for dataset size
    # sklearn requires: perplexity < (n_samples - 1) / 3
    n_samples = X.shape[0]
    max_perplexity = (n_samples - 1) / 3.0
    actual_perplexity = perplexity
    if perplexity >= max_perplexity:
        actual_perplexity = max(1.0, max_perplexity - 0.1)
        logger.warning(
            f"perplexity reduced from {perplexity} to {actual_perplexity} "
            f"(must be < (n_samples-1)/3 = {max_perplexity:.1f})"
        )

    logger.info(f"Performing t-SNE with {n_components} components, perplexity={actual_perplexity}")

    # Ensure max_iter is at least 250 (sklearn's minimum)
    actual_max_iter = max(250, max_iter)
    if actual_max_iter != max_iter:
        logger.warning(f"max_iter increased from {max_iter} to {actual_max_iter} (sklearn minimum)")

    # For large datasets, reduce the input to at most 50 principal components first (t-SNE then starts at random)
    if X.shape[0] > 5000:
        logger.info("Large dataset detected, running t-SNE on at most 50 principal components")
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

    _add_embedding_to_obs(result, tsne_coords)

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
        ValidationError: If data is not SingleCellData or n_components, n_neighbors or min_dist is out of range
    """
    if not HAS_UMAP:
        raise ImportError("UMAP requires umap-learn package")

    validation.validate_type(data, SingleCellData, "data")
    validation.validate_range(n_components, min_val=1, max_val=3, name="n_components")
    validation.validate_range(n_neighbors, min_val=5, max_val=100, name="n_neighbors")
    validation.validate_range(min_dist, min_val=0.0, max_val=1.0, name="min_dist")

    logger.info(f"Performing UMAP with {n_components} components, n_neighbors={n_neighbors}")

    # Create copy
    result = data.copy()

    if "neighbors" not in result.uns:
        # ``result`` is already this function's private copy
        result = compute_neighbors(result, n_neighbors=n_neighbors, metric=metric, inplace=True)

    # Get expression matrix
    X = result.X.toarray() if hasattr(result.X, "toarray") else result.X

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

    _add_embedding_to_obs(result, umap_coords)

    result.obsm = result.obsm if hasattr(result, "obsm") and result.obsm is not None else {}
    result.obsm["X_umap"] = X_umap

    # Store UMAP metadata
    result.uns["umap"] = {
        "n_components": n_components,
        "n_neighbors": n_neighbors,
        "min_dist": min_dist,
        "random_state": random_state,
        "metric": metric,
        "params": {
            "n_components": n_components,
            "n_neighbors": n_neighbors,
            "min_dist": min_dist,
            "random_state": random_state,
            "metric": metric,
            "n_epochs": n_epochs,
        },
    }

    logger.info("UMAP completed")
    return result


def diffusion_map_reduction(
    data: SingleCellData, n_components: int = 10, n_neighbors: int = 15, alpha: float = 1.0
) -> SingleCellData:
    """Perform diffusion map dimensionality reduction.

    The cells are the nodes of a complete graph whose edge weights are the Gaussian kernel
    ``W_ij = exp(-d_ij**2 / (2 * sigma**2))`` of their Euclidean distances, with ``sigma`` the median pairwise
    distance. The diffusion operator is the Markov matrix ``P = D**-1 W`` of the random walk on that graph
    (``D`` holds the row sums of ``W``). ``P`` is not symmetric, so it is not decomposed directly: the symmetric
    matrix ``S = D**-1/2 W D**-1/2`` has the same eigenvalues, a symmetric solver returns its orthonormal
    eigenvectors ``phi``, and the right eigenvectors of ``P`` are ``psi = D**-1/2 phi``, so that ``P psi = lambda psi``.

    Args:
        data: SingleCellData object with expression matrix
        n_components: Number of diffusion components (less than the number of cells)
        n_neighbors: Accepted and validated for compatibility, but not used: the kernel is dense over all cell
            pairs, so every cell is a neighbor of every other.
        alpha: Weight of the random-walk step in the lazy walk ``alpha * P + (1 - alpha) * I``; 1 is the plain
            random walk. The eigenvectors do not depend on it, the eigenvalues do (``alpha * lambda + 1 - alpha``).
            Must be greater than 0: with 0 the operator is the identity, whose eigenvectors are arbitrary.

    Returns:
        SingleCellData with ``obs["DC1"]`` ... ``obs["DC<n_components>"]``: the right eigenvectors ``psi_k`` of the
        Markov matrix for the ``n_components`` largest non-trivial eigenvalues, in decreasing order of eigenvalue,
        normalized so that ``sum_i D_ii psi_k(i)**2 = 1``. The trivial constant eigenvector (eigenvalue 1) is left
        out. Each sign is fixed (the entry of largest magnitude is positive). Diffusion-map coordinates at
        diffusion time ``t`` are ``eigenvalues**t * psi``; the eigenvalues are in ``uns["diffusion_map"]``, next
        to ``sigma`` and the parameters.

    Raises:
        ValidationError: If data is not SingleCellData; if n_components, n_neighbors or alpha is out of range;
            if n_components is not below the number of cells; if the expression matrix has a non-finite value;
            or if all cells are identical (the kernel width would be zero)
    """
    validation.validate_type(data, SingleCellData, "data")
    validation.validate_range(n_components, min_val=2, max_val=50, name="n_components")
    validation.validate_range(n_neighbors, min_val=5, max_val=100, name="n_neighbors")
    validation.validate_range(alpha, min_val=0.0, max_val=1.0, name="alpha")
    if alpha == 0:
        raise errors.ValidationError("alpha must be greater than 0: with alpha = 0 the operator is the identity")

    logger.info(f"Performing diffusion map with {n_components} components")

    # Create copy
    result = data.copy()

    # Get expression matrix
    X = np.asarray(data.X.toarray() if hasattr(data.X, "toarray") else data.X, dtype=float)
    n_cells = X.shape[0]
    if n_components >= n_cells:
        raise errors.ValidationError(f"n_components ({n_components}) must be less than the number of cells ({n_cells})")
    # Checked here, not left to the eigensolver: whether LAPACK raises or returns NaN for a non-finite matrix
    # depends on the build
    if not np.all(np.isfinite(X)):
        raise errors.ValidationError("The expression matrix has non-finite values (NaN or infinity)")

    # 1. Pairwise distances and the Gaussian kernel (symmetric, unit diagonal) with an adaptive width
    from scipy.spatial.distance import pdist, squareform

    distances = squareform(pdist(X, metric="euclidean"))
    sigma = float(np.median(distances))
    if sigma <= 0.0:
        raise errors.ValidationError("All cells are identical: the median pairwise distance (kernel width) is zero")
    kernel = np.exp(-(distances**2) / (2 * sigma**2))

    # 2. Symmetric conjugate of the Markov matrix P = D^-1 W: S = D^-1/2 W D^-1/2 (degrees are at least 1, the diagonal)
    degrees = kernel.sum(axis=1)
    inverse_root_degrees = 1.0 / np.sqrt(degrees)
    symmetric = kernel * inverse_root_degrees[:, np.newaxis] * inverse_root_degrees[np.newaxis, :]
    symmetric = (symmetric + symmetric.T) / 2.0  # exactly symmetric: the two roundings of the product can differ

    # 3. Lazy walk: alpha * S + (1 - alpha) * I has the eigenvectors of S and eigenvalues alpha * lambda + 1 - alpha
    if alpha != 1:
        symmetric = alpha * symmetric + (1 - alpha) * np.eye(n_cells)

    # 4. Leading eigenpairs of the symmetric matrix, then the right eigenvectors of P: psi = D^-1/2 phi
    eigenvals, phi = _leading_eigenpairs(symmetric, n_components + 1)
    right_eigenvectors = inverse_root_degrees[:, np.newaxis] * _canonical_signs(phi)

    # Skip the first, trivial eigenvector (constant, eigenvalue 1)
    diffusion_coords = right_eigenvectors[:, 1 : n_components + 1]

    # Store diffusion map results
    dm_coords = pd.DataFrame(
        diffusion_coords,
        index=result.obs.index if result.obs is not None else None,
        columns=[f"DC{i+1}" for i in range(n_components)],
    )

    _add_embedding_to_obs(result, dm_coords)

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

    The cell-by-cell Euclidean distance matrix computed from the expression
    matrix is passed to MDS as precomputed dissimilarities.

    Args:
        data: SingleCellData object with expression matrix
        n_components: Number of dimensions
        metric: Whether to use metric MDS
        random_state: Random seed for reproducibility

    Returns:
        SingleCellData with MDS coordinates added

    Raises:
        ValidationError: If data is not SingleCellData or n_components is not 2 or 3
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

    mds = MDS(
        n_components=n_components,
        metric_mds=metric,
        metric="precomputed",
        init="random",
        random_state=random_state,
        max_iter=300,
        eps=1e-6,
    )

    X_mds = mds.fit_transform(distances)

    # Store MDS results
    mds_coords = pd.DataFrame(
        X_mds,
        index=result.obs.index if result.obs is not None else None,
        columns=[f"MDS{i+1}" for i in range(n_components)],
    )

    _add_embedding_to_obs(result, mds_coords)

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


def compute_dimensionality_metrics(
    data: SingleCellData, embedding_cols: List[str], *, random_state: int | None = None
) -> Dict[str, Any]:
    """Compute quality metrics for dimensionality reduction embeddings.

    Args:
        data: SingleCellData object with embeddings
        embedding_cols: Column names of the embedding coordinates
        random_state: Seed for the cell subsample (at most 1000 cells) behind the distance
            correlation. It only matters above 1000 cells; ``None`` draws a fresh, non-reproducible
            subsample. numpy's global random state is never used.

    Returns:
        Dictionary with dimensionality reduction quality metrics: ``distance_correlation`` (Pearson correlation of
        the original and embedded pairwise distances, NaN when either set of distances is constant),
        ``embedding_dimensions``, ``n_samples``, ``n_features`` and, for at most 1000 cells,
        ``trustworthiness`` and ``continuity`` (Venna and Kaski, 2001) with their ``neighborhood_size``
        (12, or fewer neighbors for a very small dataset, because the measures need fewer than half the cells)

    Raises:
        ValidationError: If data is not SingleCellData, has no obs, or lacks the embedding columns
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
    indices = np.random.default_rng(random_state).choice(X.shape[0], size=n_samples, replace=False)

    orig_dist_flat = original_distances[indices][:, indices].flatten()
    embed_dist_flat = embedding_distances[indices][:, indices].flatten()

    # Remove self-distances; a correlation needs both sets of distances to vary (else it is undefined, not 0)
    mask = orig_dist_flat > 0
    orig_dists, embed_dists = orig_dist_flat[mask], embed_dist_flat[mask]
    if orig_dists.size > 1 and np.ptp(orig_dists) > 0 and np.ptp(embed_dists) > 0:
        distance_correlation = np.corrcoef(orig_dists, embed_dists)[0, 1]
    else:
        distance_correlation = np.nan

    # Compute embedding quality metrics
    metrics = {
        "distance_correlation": float(distance_correlation),
        "embedding_dimensions": len(embedding_cols),
        "n_samples": X.shape[0],
        "n_features": X.shape[1],
    }

    # Trustworthiness and continuity; these are computationally expensive for large datasets
    if X.shape[0] <= 1000:
        metrics.update(_compute_trustworthiness_continuity(X, embedding))

    logger.info(f"Dimensionality metrics computed: distance correlation = {distance_correlation:.3f}")
    return metrics


def _compute_trustworthiness_continuity(original: np.ndarray, embedding: np.ndarray, k: int = 12) -> Dict[str, float]:
    """Trustworthiness and continuity of an embedding (Venna and Kaski, 2001).

    Trustworthiness penalizes cells that are among the ``k`` nearest neighbors in the embedding but not in the
    original space, by how far down the original-space ranking they are; continuity penalizes the reverse (original
    neighbors that the embedding pushed away), so it is trustworthiness with the two spaces swapped. Both are 1 for
    an embedding that preserves every neighborhood. They are different numbers: counting the shared neighbors
    gives one value, which cannot be told apart from itself under two names.

    The measures need ``k < n_cells / 2``, so ``k`` is lowered to ``(n_cells - 1) // 2`` for a small dataset. With
    fewer than 3 cells there is no neighborhood and the measures are NaN.
    """
    from sklearn.manifold import trustworthiness

    n_samples = original.shape[0]
    k = min(k, (n_samples - 1) // 2)
    if k < 1:
        return {"trustworthiness": float("nan"), "continuity": float("nan"), "neighborhood_size": max(k, 0)}

    return {
        "trustworthiness": float(trustworthiness(original, embedding, n_neighbors=k)),
        "continuity": float(trustworthiness(embedding, original, n_neighbors=k)),
        "neighborhood_size": k,
    }


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
    method: str = "sklearn",
    inplace: bool = False,
) -> SingleCellData:
    """Compute neighbor graph for single-cell data.

    Like the other functions of this module, this returns a modified copy and leaves ``data``
    untouched unless ``inplace=True``.

    Args:
        data: SingleCellData object
        n_neighbors: Number of neighbors
        metric: Distance metric
        use_rep: Representation to build the graph on: ``"X"`` (the expression matrix) or the name of an array in
            ``data.obsm`` (for example ``"X_pca"`` or ``"X_umap"``). The default ``"X_pca"`` falls back to
            ``"X"`` when no PCA has been computed; that is recorded in ``uns["neighbors"]["params"]``
            (``use_rep`` is the representation that was used, ``use_rep_requested`` the one asked for).
        n_pcs: Number of leading columns of an ``obsm`` representation to use (optional)
        method: Neighbor-search backend. Currently supports ``"sklearn"``.
        inplace: Store the graph in ``data.uns["neighbors"]`` and return ``data`` itself instead of a
            copy (saves the copy of the expression matrix on large datasets).

    Returns:
        SingleCellData with neighbor graph (a copy of ``data`` unless ``inplace=True``)

    Raises:
        ImportError: If scikit-learn is unavailable
        ValueError: If ``method`` is not ``"sklearn"``
        ValidationError: If ``use_rep`` is neither ``"X"``, nor ``"X_pca"``, nor a key of ``data.obsm``
    """
    if not SKLEARN_AVAILABLE:
        raise ImportError("scikit-learn required for neighbor computation")
    if method != "sklearn":
        raise ValueError("compute_neighbors currently supports only method='sklearn'")

    from sklearn.neighbors import NearestNeighbors

    # Get the representation to use. Only the default "X_pca" may fall back to X (when no PCA exists); any other
    # name must exist, since building the graph on a different space than the one asked for would be silent.
    obsm = data.obsm if getattr(data, "obsm", None) else {}
    if use_rep != "X" and use_rep in obsm:
        used_rep = use_rep
        X = np.asarray(obsm[use_rep])
        # Limit to n_pcs if specified
        if n_pcs is not None and X.ndim == 2 and X.shape[1] > n_pcs:
            X = X[:, :n_pcs]
    elif use_rep == "X" or use_rep == "X_pca":
        used_rep = "X"
        if use_rep == "X_pca":
            logger.info("No X_pca in data.obsm: building the neighbor graph on X")
        X = data.X
        if hasattr(X, "toarray"):
            X = X.toarray()
    else:
        raise errors.ValidationError(
            f"Representation '{use_rep}' not found: use 'X', 'X_pca' or a key of data.obsm ({sorted(obsm)})"
        )

    # Compute neighbors
    nbrs = NearestNeighbors(n_neighbors=n_neighbors, metric=metric)
    nbrs.fit(X)

    distances, indices = nbrs.kneighbors(X)

    # Build connectivity matrix (sparse)
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
    result = data if inplace else data.copy()
    result.uns["neighbors"] = {
        "indices": indices,
        "distances": distances,
        "connectivities": connectivities,
        "params": {
            "n_neighbors": n_neighbors,
            "metric": metric,
            "use_rep": used_rep,
            "use_rep_requested": use_rep,
            "n_pcs": n_pcs,
            "method": method,
        },
    }

    logger.info(f"Computed neighbor graph with {n_neighbors} neighbors for {len(X)} cells")
    return result


def _checked_neighbor_arrays(data: SingleCellData, *keys: str) -> list[np.ndarray]:
    """The per-cell arrays ``keys`` of ``data.uns["neighbors"]``, checked to describe exactly the cells of ``data``.

    A graph that was computed before cells were removed (or for other data altogether) indexes the wrong
    cells, and using it would give a wrong answer without any error.

    Raises:
        ValidationError: If an array is missing, or its rows (or the cell indices in it) do not match ``data.n_obs``
    """
    graph = data.uns["neighbors"]
    missing = [key for key in keys if key not in graph]
    if missing:
        raise errors.ValidationError(
            f"uns['neighbors'] has no {' and '.join(repr(key) for key in missing)}: recompute it with compute_neighbors"
        )
    arrays = [np.asarray(graph[key]) for key in keys]
    for array in arrays:
        if array.shape[0] != data.n_obs:
            raise errors.ValidationError(
                f"uns['neighbors'] describes {array.shape[0]} cells but data has {data.n_obs}: "
                "recompute it with compute_neighbors"
            )
    if "indices" in keys:
        listed = arrays[keys.index("indices")]
        if listed.size and listed.max() >= data.n_obs:
            raise errors.ValidationError(
                "uns['neighbors']['indices'] refers to cells that data does not have: "
                "recompute it with compute_neighbors"
            )
        if "distances" in keys and arrays[keys.index("distances")].shape != listed.shape:
            raise errors.ValidationError("uns['neighbors'] 'indices' and 'distances' differ in shape")
    return arrays


def _knn_affinity_matrix(indices: np.ndarray, distances: np.ndarray) -> sparse.csr_matrix:
    """Symmetric affinity matrix of a k-nearest-neighbor graph, ``1 / (1 + distance)`` per edge, no self-loops.

    Row ``i`` of ``indices``/``distances`` lists the neighbors of cell ``i`` and their distances (a cell's own
    entry is ignored). An edge exists between two cells when either lists the other. The directed graph is made
    symmetric with the entrywise maximum, which keeps one weight per pair: adding the transpose would give two
    cells that list each other twice the weight of a pair at the same distance that only one of them lists.
    """
    n_cells, n_listed = indices.shape
    sources = np.repeat(np.arange(n_cells), n_listed)
    targets = indices.reshape(-1)
    weights = 1.0 / (1.0 + distances.reshape(-1))
    others = sources != targets
    directed = sparse.csr_matrix((weights[others], (sources[others], targets[others])), shape=(n_cells, n_cells))
    return sparse.csr_matrix(directed.maximum(directed.T))


def compute_diffusion_map(
    data: SingleCellData,
    n_components: int = 10,
    n_neighbors: int = 15,
    alpha: float = 1.0,
    *,
    inplace: bool = False,
) -> SingleCellData:
    """Compute diffusion map for single-cell data.

    Like the other functions of this module, this returns a modified copy and leaves ``data``
    untouched unless ``inplace=True``.

    The affinity between two cells is ``1 / (1 + distance)`` for every pair in which at least one cell is among
    the other's nearest neighbors in ``uns["neighbors"]``: the graph is symmetrized by the entrywise maximum, so
    two cells that list each other are joined by one edge of that weight, not by two added up. The eigenvectors
    are those of the diffusion operator below, largest eigenvalue first (the trivial first one is skipped), each
    with a fixed sign (its entry of largest magnitude is positive); they are not rescaled by the eigenvalues.

    Args:
        data: SingleCellData object
        n_components: Number of diffusion components
        n_neighbors: Number of neighbors for graph construction (only used when ``data`` has no
            ``uns["neighbors"]`` graph yet)
        alpha: Normalization of the affinity matrix ``A`` with degrees ``d``. ``1`` (the default): the symmetric
            normalization ``D**-1/2 A D**-1/2`` of the random walk's Markov matrix, whose eigenvalues lie in
            [-1, 1] and whose first (trivial) eigenvalue is 1. ``0``: the raw affinity matrix, no normalization.
            In between: ``D**-alpha A D**-alpha``.
        inplace: Write ``obsm["X_diffmap"]``, ``uns["diffmap"]`` (and the neighbor graph, when it had to
            be built) into ``data`` and return ``data`` itself instead of a copy.

    Returns:
        SingleCellData with ``obsm["X_diffmap"]`` and ``uns["diffmap"]`` (the eigenvalues and parameters), a copy of
        ``data`` unless ``inplace=True``

    Raises:
        ImportError: If scikit-learn is unavailable
        ValidationError: If ``uns["neighbors"]`` has no ``indices`` and ``distances``, or describes another number of
            cells than ``data`` has (a stale graph: recompute it with ``compute_neighbors``)
    """
    if not SKLEARN_AVAILABLE:
        raise ImportError("scikit-learn required for diffusion maps")

    result = data if inplace else data.copy()

    # First compute neighbors if not already done (``result`` is already this function's own object)
    if "neighbors" not in result.uns:
        result = compute_neighbors(result, n_neighbors, inplace=True)

    # Get neighbor information (refused if it describes other cells than the ones in data)
    n_cells = result.n_obs
    indices, distances = _checked_neighbor_arrays(result, "indices", "distances")
    adjacency = _knn_affinity_matrix(indices, distances)

    # Normalize adjacency matrix (create diffusion operator)
    degrees = np.array(adjacency.sum(axis=1)).flatten()
    degrees[degrees == 0] = 1  # Avoid division by zero

    diffusion_operator: sparse.spmatrix
    if alpha == 0:
        # Raw affinities, no normalization
        diffusion_operator = adjacency
    elif alpha == 1:
        # Symmetric normalization of the random walk's Markov matrix
        D_inv_sqrt = sparse.diags(1.0 / np.sqrt(degrees))
        diffusion_operator = D_inv_sqrt @ adjacency @ D_inv_sqrt
    else:
        # Generalized normalization
        D_alpha = sparse.diags(degrees ** (-alpha))
        diffusion_operator = D_alpha @ adjacency @ D_alpha

    # Leading eigenpairs (diffusion components); the solver is chosen by problem size, not by catching its errors
    eigenvalues, eigenvectors = _leading_eigenpairs(diffusion_operator, min(n_components + 1, n_cells))
    eigenvectors = _canonical_signs(eigenvectors)

    # Store results (skip first eigenvector, the trivial one)
    result.obsm["X_diffmap"] = eigenvectors[:, 1 : n_components + 1]
    result.uns["diffmap"] = {
        "eigenvalues": eigenvalues[1 : n_components + 1],
        "params": {"n_components": n_components, "n_neighbors": n_neighbors, "alpha": alpha},
    }

    logger.info(f"Computed diffusion map with {n_components} components")
    return result
