"""Multi-omic clustering methods for sample stratification.

Provides methods for clustering biological samples using integrated
multi-omic data.  Includes concatenation-based, late integration, and
SNF-based approaches as well as consensus clustering for determining
optimal cluster number and multi-view spectral clustering.

Algorithms:
    - multi_omic_clustering: SNF, concatenation, or late integration.
    - consensus_clustering: resampling-based consensus matrix for optimal k.
    - multi_view_spectral: spectral clustering across multiple similarity views.
    - evaluate_integration: per-omic silhouette and cross-omic ARI.
"""

from __future__ import annotations

import math
import random
from collections import Counter
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Optional scientific dependencies
try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None  # type: ignore[assignment]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _simple_kmeans(X: Any, k: int, max_iter: int = 50, seed: int = 42) -> list[int]:
    """Simple k-means clustering on a numpy array.

    Args:
        X: Data matrix (n x d).
        k: Number of clusters.
        max_iter: Maximum iterations.
        seed: Random seed for centroid initialisation.

    Returns:
        List of integer cluster labels.
    """
    n = X.shape[0]
    rng = np.random.RandomState(seed)
    indices = rng.choice(n, size=min(k, n), replace=False)
    centroids = X[indices].copy()

    labels = np.zeros(n, dtype=int)
    for _ in range(max_iter):
        # Assign each point to nearest centroid
        for i in range(n):
            dists = [float(np.sum((X[i] - centroids[c]) ** 2)) for c in range(k)]
            labels[i] = int(np.argmin(dists))
        # Recompute centroids
        new_centroids = np.zeros_like(centroids)
        for c in range(k):
            members = X[labels == c]
            if len(members) > 0:
                new_centroids[c] = members.mean(axis=0)
            else:
                new_centroids[c] = centroids[c]
        if np.allclose(centroids, new_centroids):
            break
        centroids = new_centroids

    return labels.tolist()


def _silhouette_from_data(X: Any, labels: list[int]) -> float:
    """Compute silhouette score from data matrix and labels.

    Args:
        X: Data matrix (n x d).
        labels: Cluster labels.

    Returns:
        Mean silhouette score in [-1, 1].
    """
    n = len(labels)
    unique = list(set(labels))
    if n < 2 or len(unique) < 2:
        return 0.0

    # Pairwise euclidean distances
    dists = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = float(np.sqrt(np.sum((X[i] - X[j]) ** 2)))
            dists[i, j] = d
            dists[j, i] = d

    silhouettes: list[float] = []
    for i in range(n):
        li = labels[i]
        same = [j for j in range(n) if labels[j] == li and j != i]
        if not same:
            silhouettes.append(0.0)
            continue
        a_i = float(np.mean(dists[i, same]))
        b_i = float("inf")
        for other in unique:
            if other == li:
                continue
            others = [j for j in range(n) if labels[j] == other]
            if others:
                b_i = min(b_i, float(np.mean(dists[i, others])))
        denom = max(a_i, b_i)
        silhouettes.append((b_i - a_i) / denom if denom > 0 else 0.0)

    return float(np.mean(silhouettes))


def _silhouette_from_similarity(similarity: Any, labels: list[int]) -> float:
    """Compute silhouette score from a similarity matrix and labels.

    Args:
        similarity: Symmetric similarity matrix (n x n).
        labels: Cluster labels.

    Returns:
        Mean silhouette score.
    """
    n = len(labels)
    unique = list(set(labels))
    if n < 2 or len(unique) < 2:
        return 0.0

    distance = 1.0 - similarity
    np.fill_diagonal(distance, 0.0)

    silhouettes: list[float] = []
    for i in range(n):
        li = labels[i]
        same = [j for j in range(n) if labels[j] == li and j != i]
        if not same:
            silhouettes.append(0.0)
            continue
        a_i = float(np.mean(distance[i, same]))
        b_i = float("inf")
        for other in unique:
            if other == li:
                continue
            others = [j for j in range(n) if labels[j] == other]
            if others:
                b_i = min(b_i, float(np.mean(distance[i, others])))
        denom = max(a_i, b_i)
        silhouettes.append((b_i - a_i) / denom if denom > 0 else 0.0)

    return float(np.mean(silhouettes))


def _adjusted_rand_index(labels_a: list[int], labels_b: list[int]) -> float:
    """Compute Adjusted Rand Index between two label assignments.

    Args:
        labels_a: First set of cluster labels.
        labels_b: Second set of cluster labels.

    Returns:
        ARI value in [-1, 1]. 1.0 indicates perfect agreement.
    """
    n = len(labels_a)
    if n != len(labels_b):
        raise ValueError("Label lists must have the same length")
    if n == 0:
        return 0.0

    # Build contingency table
    classes_a = sorted(set(labels_a))
    classes_b = sorted(set(labels_b))
    map_a = {c: i for i, c in enumerate(classes_a)}
    map_b = {c: i for i, c in enumerate(classes_b)}

    nij = np.zeros((len(classes_a), len(classes_b)), dtype=np.float64)
    for idx in range(n):
        nij[map_a[labels_a[idx]], map_b[labels_b[idx]]] += 1

    a_sums = nij.sum(axis=1)
    b_sums = nij.sum(axis=0)

    def _comb2(x: float) -> float:
        return x * (x - 1) / 2.0

    sum_comb_nij = sum(_comb2(nij[i, j]) for i in range(len(classes_a)) for j in range(len(classes_b)))
    sum_comb_a = sum(_comb2(a_sums[i]) for i in range(len(classes_a)))
    sum_comb_b = sum(_comb2(b_sums[j]) for j in range(len(classes_b)))
    comb_n = _comb2(n)

    expected = sum_comb_a * sum_comb_b / comb_n if comb_n > 0 else 0.0
    max_index = 0.5 * (sum_comb_a + sum_comb_b)
    denom = max_index - expected

    if abs(denom) < 1e-12:
        return 1.0 if abs(sum_comb_nij - expected) < 1e-12 else 0.0

    return (sum_comb_nij - expected) / denom


# ---------------------------------------------------------------------------
# Multi-omic clustering
# ---------------------------------------------------------------------------


def multi_omic_clustering(
    data_matrices: dict[str, Any],
    n_clusters: int,
    method: str = "snf",
) -> dict[str, Any]:
    """Cluster samples using integrated multi-omic data.

    Supports three integration strategies:
        - *concatenation*: horizontally concatenate all omic features, then
          k-means on the joint matrix.
        - *snf*: build per-omic similarity networks, fuse with SNF, then
          spectral cluster on fused network.
        - *late_integration*: cluster each omic independently, then determine
          consensus labels via majority voting.

    Args:
        data_matrices: Mapping from omic name to data matrix (numpy array
            or list-of-lists, n_samples x n_features).
        n_clusters: Number of clusters.
        method: Integration strategy, one of ``"snf"``, ``"concatenation"``,
            or ``"late_integration"``.

    Returns:
        Dictionary with keys:
            labels: Final cluster labels (list of int).
            silhouette: Silhouette score of final clustering.
            omic_contributions: Per-omic information about how each layer
                contributed to the final clustering.

    Raises:
        ValueError: If method is unsupported, n_clusters < 2, or matrices
            are incompatible.
        ImportError: If numpy is not available.
    """
    valid_methods = {"snf", "concatenation", "late_integration"}
    if method not in valid_methods:
        raise ValueError(f"method must be one of {valid_methods}, got '{method}'")
    if n_clusters < 2:
        raise ValueError(f"n_clusters must be >= 2, got {n_clusters}")
    if not data_matrices:
        raise ValueError("data_matrices must contain at least one entry")

    if not HAS_NUMPY:
        raise ImportError("numpy is required for multi-omic clustering. " "Install with: uv pip install numpy")

    logger.info("Running multi-omic clustering: method=%s, k=%d", method, n_clusters)

    # Convert and validate
    views: dict[str, Any] = {}
    n_samples: int | None = None
    for name, mat in data_matrices.items():
        arr = np.asarray(mat, dtype=np.float64)
        if arr.ndim != 2:
            raise ValueError(f"Matrix '{name}' must be 2-D")
        if n_samples is None:
            n_samples = arr.shape[0]
        elif arr.shape[0] != n_samples:
            raise ValueError(f"Sample-dimension mismatch: expected {n_samples}, " f"got {arr.shape[0]} for '{name}'")
        views[name] = arr

    assert n_samples is not None

    omic_contributions: dict[str, Any] = {}

    if method == "concatenation":
        # Standardise each view and concatenate
        parts: list[Any] = []
        for name, arr in views.items():
            mean = arr.mean(axis=0, keepdims=True)
            std = arr.std(axis=0, keepdims=True) + 1e-10
            parts.append((arr - mean) / std)
            omic_contributions[name] = {"n_features": arr.shape[1]}
        combined = np.hstack(parts)
        labels = _simple_kmeans(combined, n_clusters)
        sil = _silhouette_from_data(combined, labels)
        for name in views:
            omic_contributions[name]["feature_fraction"] = views[name].shape[1] / combined.shape[1]

    elif method == "snf":
        from metainformant.multiomics.methods.factorization import similarity_network_fusion

        # Build per-omic correlation-based similarity
        networks: list[Any] = []
        for name, arr in views.items():
            sim = np.corrcoef(arr)
            # Clip to [0, 1] for similarity
            sim = (sim + 1.0) / 2.0
            np.fill_diagonal(sim, 1.0)
            networks.append(sim)
            omic_contributions[name] = {"n_features": arr.shape[1]}

        snf_result = similarity_network_fusion(networks, k_neighbors=min(20, n_samples - 1))
        fused = np.asarray(snf_result["fused_network"])

        # Spectral clustering on fused network
        labels, sil = _spectral_cluster_similarity(fused, n_clusters)

        for i, name in enumerate(views):
            omic_contributions[name]["snf_contribution"] = float(np.mean(np.abs(np.asarray(networks[i]) - fused)))

    else:  # late_integration
        per_omic_labels: dict[str, list[int]] = {}
        for name, arr in views.items():
            mean = arr.mean(axis=0, keepdims=True)
            std = arr.std(axis=0, keepdims=True) + 1e-10
            normalised = (arr - mean) / std
            lbl = _simple_kmeans(normalised, n_clusters)
            per_omic_labels[name] = lbl
            omic_contributions[name] = {
                "n_features": arr.shape[1],
                "omic_labels": lbl,
                "omic_silhouette": _silhouette_from_data(normalised, lbl),
            }

        # Majority voting (co-association matrix)
        co_assoc = np.zeros((n_samples, n_samples), dtype=np.float64)
        for lbl in per_omic_labels.values():
            for i in range(n_samples):
                for j in range(n_samples):
                    if lbl[i] == lbl[j]:
                        co_assoc[i, j] += 1.0
        co_assoc /= len(per_omic_labels)
        labels, sil = _spectral_cluster_similarity(co_assoc, n_clusters)

    logger.info("Multi-omic clustering complete: silhouette=%.4f", sil)

    return {
        "labels": labels,
        "silhouette": sil,
        "omic_contributions": omic_contributions,
    }


def _spectral_cluster_similarity(similarity: Any, n_clusters: int) -> tuple[list[int], float]:
    """Spectral clustering from similarity matrix.

    Args:
        similarity: Symmetric similarity matrix.
        n_clusters: Number of clusters.

    Returns:
        Tuple of (labels, silhouette_score).
    """
    n = similarity.shape[0]
    if n < n_clusters:
        return list(range(n)), 0.0

    D = np.diag(similarity.sum(axis=1))
    L = D - similarity
    D_inv_sqrt = np.diag(1.0 / (np.sqrt(np.diag(D)) + 1e-10))
    L_norm = D_inv_sqrt @ L @ D_inv_sqrt

    eigenvalues, eigenvectors = np.linalg.eigh(L_norm)
    idx = np.argsort(eigenvalues)[:n_clusters]
    embedding = eigenvectors[:, idx]

    row_norms = np.linalg.norm(embedding, axis=1, keepdims=True) + 1e-10
    embedding = embedding / row_norms

    labels = _simple_kmeans(embedding, n_clusters)
    sil = _silhouette_from_similarity(similarity, labels)

    return labels, sil


# ---------------------------------------------------------------------------
# Consensus clustering
# ---------------------------------------------------------------------------


def consensus_clustering(
    data: Any,
    k_range: range | list[int] | None = None,
    n_resamples: int = 100,
    proportion: float = 0.8,
) -> dict[str, Any]:
    """Consensus clustering with resampling to determine optimal k.

    For each candidate k in *k_range*, performs *n_resamples* rounds of
    sub-sampling (keeping *proportion* of samples each time), clusters
    the subsample, and builds a consensus matrix recording the fraction
    of times each pair of samples was placed in the same cluster.

    The optimal k is selected by the Proportion of Ambiguous Clustering
    (PAC) score -- lower PAC indicates cleaner consensus.

    Args:
        data: Data matrix (n_samples x n_features) as numpy array or
            list-of-lists.
        k_range: Range of k values to evaluate.  Defaults to range(2, 7).
        n_resamples: Number of resampling rounds per k.
        proportion: Fraction of samples to include in each subsample.

    Returns:
        Dictionary with keys:
            optimal_k: The k with lowest PAC score.
            labels: Cluster labels for optimal_k (hierarchical clustering
                of the consensus matrix).
            consensus_matrix: Consensus matrix for optimal_k (n x n).
            cdf_area: Area under the CDF curve for each k.
            pac_score: PAC score for each k.

    Raises:
        ValueError: If data is empty, proportion not in (0,1], or k_range
            is empty.
        ImportError: If numpy is not available.
    """
    if not HAS_NUMPY:
        raise ImportError("numpy is required for consensus clustering. " "Install with: uv pip install numpy")

    X = np.asarray(data, dtype=np.float64)
    if X.ndim != 2:
        raise ValueError(f"data must be 2-D, got {X.ndim}-D")

    n = X.shape[0]
    if n < 3:
        raise ValueError(f"Need at least 3 samples for consensus clustering, got {n}")

    if k_range is None:
        k_range = list(range(2, min(7, n)))
    k_range = list(k_range)
    if not k_range:
        raise ValueError("k_range must be non-empty")
    if proportion <= 0 or proportion > 1.0:
        raise ValueError(f"proportion must be in (0, 1], got {proportion}")

    logger.info(
        "Running consensus clustering: k_range=%s, n_resamples=%d, proportion=%.2f",
        k_range,
        n_resamples,
        proportion,
    )

    rng = np.random.RandomState(42)
    subsample_size = max(2, int(n * proportion))

    consensus_matrices: dict[int, Any] = {}
    cdf_area: dict[int, float] = {}
    pac_score: dict[int, float] = {}

    for k in k_range:
        if k > subsample_size:
            logger.warning("k=%d > subsample_size=%d; skipping", k, subsample_size)
            continue

        # Count how many times pairs co-occur in a subsample
        co_count = np.zeros((n, n), dtype=np.float64)
        co_cluster = np.zeros((n, n), dtype=np.float64)

        for _ in range(n_resamples):
            idx = rng.choice(n, size=subsample_size, replace=False)
            sub_X = X[idx]

            # Standardise
            mean = sub_X.mean(axis=0, keepdims=True)
            std = sub_X.std(axis=0, keepdims=True) + 1e-10
            sub_X_norm = (sub_X - mean) / std

            labels_sub = _simple_kmeans(sub_X_norm, k, seed=rng.randint(0, 100000))

            for ii in range(subsample_size):
                for jj in range(ii + 1, subsample_size):
                    real_i = idx[ii]
                    real_j = idx[jj]
                    co_count[real_i, real_j] += 1
                    co_count[real_j, real_i] += 1
                    if labels_sub[ii] == labels_sub[jj]:
                        co_cluster[real_i, real_j] += 1
                        co_cluster[real_j, real_i] += 1

        # Consensus matrix
        C = np.divide(co_cluster, co_count, out=np.zeros_like(co_cluster), where=co_count > 0)
        np.fill_diagonal(C, 1.0)
        consensus_matrices[k] = C

        # CDF of consensus values
        upper_tri = C[np.triu_indices(n, k=1)]
        n_pairs = len(upper_tri)
        if n_pairs == 0:
            cdf_area[k] = 0.0
            pac_score[k] = 1.0
            continue

        sorted_vals = np.sort(upper_tri)
        cdf_vals = np.arange(1, n_pairs + 1) / n_pairs

        # Area under CDF
        area = float(np.trapz(cdf_vals, sorted_vals))
        cdf_area[k] = area

        # PAC: fraction of entries in (0.1, 0.9) -- ambiguous zone
        ambiguous = np.sum((upper_tri > 0.1) & (upper_tri < 0.9))
        pac_score[k] = float(ambiguous) / n_pairs

    if not pac_score:
        raise ValueError("No valid k values could be evaluated")

    # Select optimal k (lowest PAC)
    optimal_k = min(pac_score, key=lambda kk: pac_score[kk])
    C_opt = consensus_matrices[optimal_k]

    # Final labels from spectral clustering on consensus matrix
    labels_final, _ = _spectral_cluster_similarity(C_opt, optimal_k)

    logger.info(
        "Consensus clustering complete: optimal_k=%d, PAC=%.4f",
        optimal_k,
        pac_score[optimal_k],
    )

    return {
        "optimal_k": optimal_k,
        "labels": labels_final,
        "consensus_matrix": C_opt.tolist(),
        "cdf_area": cdf_area,
        "pac_score": pac_score,
    }


# ---------------------------------------------------------------------------
# Multi-view spectral clustering
# ---------------------------------------------------------------------------


def multi_view_spectral(
    similarity_matrices: list[Any],
    n_clusters: int,
    method: str = "average",
) -> dict[str, Any]:
    """Multi-view spectral clustering on multiple similarity matrices.

    Combines multiple view-specific similarity matrices into a single
    unified similarity, then performs spectral clustering.

    Combination strategies:
        - *average*: element-wise average of all similarity matrices.
        - *product*: element-wise product (Hadamard) then re-normalise.
        - *max*: element-wise maximum.

    Args:
        similarity_matrices: List of symmetric similarity matrices (numpy
            arrays or list-of-lists), one per view.
        n_clusters: Number of clusters.
        method: Combination method, one of ``"average"``, ``"product"``,
            or ``"max"``.

    Returns:
        Dictionary with keys:
            labels: Cluster labels (list of int).
            eigenvalues: Eigenvalues of the normalised Laplacian (list).
            eigenvectors: Eigenvectors used for embedding (list-of-lists).

    Raises:
        ValueError: If no similarity matrices provided, matrices have
            different shapes, or method is unsupported.
        ImportError: If numpy is not available.
    """
    valid_methods = {"average", "product", "max"}
    if method not in valid_methods:
        raise ValueError(f"method must be one of {valid_methods}, got '{method}'")
    if not similarity_matrices:
        raise ValueError("At least one similarity matrix is required")
    if n_clusters < 2:
        raise ValueError(f"n_clusters must be >= 2, got {n_clusters}")

    if not HAS_NUMPY:
        raise ImportError("numpy is required for multi-view spectral clustering. " "Install with: uv pip install numpy")

    logger.info(
        "Running multi-view spectral clustering: %d views, k=%d, method=%s",
        len(similarity_matrices),
        n_clusters,
        method,
    )

    arrays: list[Any] = []
    n: int | None = None
    for i, mat in enumerate(similarity_matrices):
        arr = np.asarray(mat, dtype=np.float64)
        if arr.ndim != 2 or arr.shape[0] != arr.shape[1]:
            raise ValueError(f"Matrix {i} must be square")
        if n is None:
            n = arr.shape[0]
        elif arr.shape[0] != n:
            raise ValueError(f"Matrix {i} has size {arr.shape[0]}, expected {n}")
        arrays.append(arr)

    assert n is not None

    # Combine
    if method == "average":
        combined = np.mean(arrays, axis=0)
    elif method == "product":
        combined = np.ones((n, n), dtype=np.float64)
        for arr in arrays:
            combined *= arr
        # Re-normalise to [0, 1]
        c_min = combined.min()
        c_max = combined.max()
        if c_max > c_min:
            combined = (combined - c_min) / (c_max - c_min)
    else:  # max
        combined = np.maximum.reduce(arrays)

    # Spectral clustering
    D = np.diag(combined.sum(axis=1))
    L = D - combined
    D_inv_sqrt = np.diag(1.0 / (np.sqrt(np.diag(D)) + 1e-10))
    L_norm = D_inv_sqrt @ L @ D_inv_sqrt

    eigenvalues, eigenvectors = np.linalg.eigh(L_norm)
    idx = np.argsort(eigenvalues)[:n_clusters]
    selected_eigenvalues = eigenvalues[idx]
    embedding = eigenvectors[:, idx]

    row_norms = np.linalg.norm(embedding, axis=1, keepdims=True) + 1e-10
    embedding = embedding / row_norms

    labels = _simple_kmeans(embedding, n_clusters)

    logger.info("Multi-view spectral clustering complete")

    return {
        "labels": labels,
        "eigenvalues": selected_eigenvalues.tolist(),
        "eigenvectors": embedding.tolist(),
    }


# ---------------------------------------------------------------------------
# Integration quality evaluation
# ---------------------------------------------------------------------------


def evaluate_integration(
    labels: list[int],
    omic_data: dict[str, Any],
) -> dict[str, Any]:
    """Evaluate clustering quality across multiple omic layers.

    Computes per-omic silhouette scores using the integrated labels and
    Adjusted Rand Index (ARI) between each omic-specific clustering and
    the integrated labels.  Also provides an aggregate integration metric.

    Args:
        labels: Integrated cluster labels from multi-omic clustering.
        omic_data: Mapping from omic name to data matrix (n_samples x
            n_features) as numpy array or list-of-lists.

    Returns:
        Dictionary with keys:
            silhouette_per_omic: Per-omic silhouette scores {name: float}.
            ari_per_omic: ARI between omic-specific and integrated labels
                {name: float}.
            mean_silhouette: Average silhouette across omics.
            mean_ari: Average ARI across omics.
            integration_metric: Combined quality metric (mean of normalised
                silhouette and ARI).

    Raises:
        ValueError: If labels length does not match data dimensions.
        ImportError: If numpy is not available.
    """
    if not omic_data:
        raise ValueError("omic_data must contain at least one entry")

    if not HAS_NUMPY:
        raise ImportError("numpy is required for integration evaluation. " "Install with: uv pip install numpy")

    logger.info("Evaluating integration quality across %d omics", len(omic_data))

    n_clusters = len(set(labels))
    if n_clusters < 2:
        logger.warning("Only %d cluster(s) found; metrics may not be meaningful", n_clusters)

    silhouette_per_omic: dict[str, float] = {}
    ari_per_omic: dict[str, float] = {}

    for name, mat in omic_data.items():
        arr = np.asarray(mat, dtype=np.float64)
        if arr.ndim != 2:
            raise ValueError(f"Matrix '{name}' must be 2-D")
        if arr.shape[0] != len(labels):
            raise ValueError(f"Matrix '{name}' has {arr.shape[0]} samples but labels has {len(labels)}")

        # Standardise
        mean = arr.mean(axis=0, keepdims=True)
        std = arr.std(axis=0, keepdims=True) + 1e-10
        arr_norm = (arr - mean) / std

        # Silhouette of integrated labels on this omic
        sil = _silhouette_from_data(arr_norm, labels)
        silhouette_per_omic[name] = sil

        # Omic-specific clustering, then ARI
        omic_labels = _simple_kmeans(arr_norm, n_clusters)
        ari = _adjusted_rand_index(labels, omic_labels)
        ari_per_omic[name] = ari

    mean_sil = float(np.mean(list(silhouette_per_omic.values())))
    mean_ari = float(np.mean(list(ari_per_omic.values())))

    # Integration metric: average of normalised silhouette ([0,1] mapped) and ARI
    # Silhouette is in [-1, 1], map to [0, 1]
    norm_sil = (mean_sil + 1.0) / 2.0
    # ARI is in [-1, 1], map to [0, 1]
    norm_ari = (mean_ari + 1.0) / 2.0
    integration_metric = (norm_sil + norm_ari) / 2.0

    logger.info(
        "Integration evaluation: mean_silhouette=%.4f, mean_ARI=%.4f, metric=%.4f",
        mean_sil,
        mean_ari,
        integration_metric,
    )

    return {
        "silhouette_per_omic": silhouette_per_omic,
        "ari_per_omic": ari_per_omic,
        "mean_silhouette": mean_sil,
        "mean_ari": mean_ari,
        "integration_metric": integration_metric,
    }
