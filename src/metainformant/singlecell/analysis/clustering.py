"""Single-cell clustering algorithms and cluster analysis.

This module provides functions for clustering single-cell data using various
algorithms including graph-based methods (Leiden, Louvain) and traditional
methods (K-means, hierarchical). It also includes cluster evaluation and
marker gene identification.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Any, Tuple, Union
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist

from metainformant.core import logging, errors, validation

# Optional scientific dependencies
try:
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score

    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False
    KMeans = None
    silhouette_score = None
    calinski_harabasz_score = None
    davies_bouldin_score = None

# Try to import networkx and community detection libraries
try:
    import networkx as nx

    HAS_NETWORKX = True
except ImportError:
    nx = None
    HAS_NETWORKX = False

try:
    import igraph as ig

    HAS_IGRAPH = True
except ImportError:
    ig = None
    HAS_IGRAPH = False

try:
    import leidenalg

    HAS_LEIDEN = True
except ImportError:
    leidenalg = None
    HAS_LEIDEN = False

try:
    import community as community_louvain

    HAS_LOUVAIN = True
except ImportError:
    community_louvain = None
    HAS_LOUVAIN = False

logger = logging.get_logger(__name__)

# Import our SingleCellData
from metainformant.singlecell.data.preprocessing import SingleCellData


def leiden_clustering(
    data: SingleCellData,
    resolution: float = 1.0,
    n_neighbors: int = 15,
    random_state: int | None = None,
    use_weights: bool = True,
) -> SingleCellData:
    """Perform Leiden clustering on single-cell data.

    Args:
        data: SingleCellData object with preprocessed expression matrix
        resolution: Resolution parameter for Leiden algorithm
        n_neighbors: Number of neighbors for kNN graph construction
        random_state: Random seed for reproducibility
        use_weights: Whether to use edge weights in clustering

    Returns:
        SingleCellData with cluster assignments added to obs

    Raises:
        ImportError: If required packages are not available
        TypeError: If data is not SingleCellData
    """
    validation.validate_type(data, SingleCellData, "data")

    if not HAS_NETWORKX or not HAS_IGRAPH or not HAS_LEIDEN:
        raise ImportError("Leiden clustering requires networkx, igraph, and leidenalg packages")

    logger.info(f"Performing Leiden clustering with resolution={resolution}, n_neighbors={n_neighbors}")

    # Create copy
    result = data.copy()

    # Get expression matrix
    X = data.X.toarray() if hasattr(data.X, "toarray") else data.X

    # Construct kNN graph
    n_cells = X.shape[0]
    distances = _compute_pairwise_distances(X)

    # Build adjacency matrix
    adjacency = _build_knn_adjacency(distances, n_neighbors)

    # Convert to networkx graph
    G = nx.from_scipy_sparse_matrix(adjacency) if sparse.issparse(adjacency) else nx.from_numpy_matrix(adjacency)

    # Perform Leiden clustering
    partition = leidenalg.find_partition(
        ig.Graph.from_networkx(G),
        leidenalg.ModularityVertexPartition,
        resolution_parameter=resolution,
        seed=random_state,
    )

    # Get cluster assignments
    clusters = np.array(partition.membership)

    # Add to obs
    result.obs = result.obs.copy() if result.obs is not None else pd.DataFrame()
    result.obs["leiden_cluster"] = clusters

    # Store clustering parameters
    result.uns["leiden_clustering"] = {
        "resolution": resolution,
        "n_neighbors": n_neighbors,
        "random_state": random_state,
        "n_clusters": len(np.unique(clusters)),
        "use_weights": use_weights,
    }

    logger.info(f"Leiden clustering completed: found {len(np.unique(clusters))} clusters")
    return result


def louvain_clustering(
    data: SingleCellData,
    resolution: float = 1.0,
    n_neighbors: int = 15,
    random_state: int | None = None,
    use_weights: bool = True,
) -> SingleCellData:
    """Perform Louvain clustering on single-cell data.

    Args:
        data: SingleCellData object with preprocessed expression matrix
        resolution: Resolution parameter for Louvain algorithm
        n_neighbors: Number of neighbors for kNN graph construction
        random_state: Random seed for reproducibility
        use_weights: Whether to use edge weights in clustering

    Returns:
        SingleCellData with cluster assignments added to obs

    Raises:
        ImportError: If required packages are not available
        TypeError: If data is not SingleCellData
    """
    validation.validate_type(data, SingleCellData, "data")

    if not HAS_NETWORKX or not HAS_LOUVAIN:
        raise ImportError("Louvain clustering requires networkx and python-louvain packages")

    logger.info(f"Performing Louvain clustering with resolution={resolution}, n_neighbors={n_neighbors}")

    # Create copy
    result = data.copy()

    # Get expression matrix
    X = data.X.toarray() if hasattr(data.X, "toarray") else data.X

    # Construct kNN graph
    distances = _compute_pairwise_distances(X)
    adjacency = _build_knn_adjacency(distances, n_neighbors)

    # Convert to networkx graph
    G = nx.from_scipy_sparse_matrix(adjacency) if sparse.issparse(adjacency) else nx.from_numpy_matrix(adjacency)

    # Perform Louvain clustering
    partition = community_louvain.best_partition(G, resolution=resolution, random_state=random_state)
    clusters = np.array([partition[i] for i in range(len(partition))])

    # Add to obs
    result.obs = result.obs.copy() if result.obs is not None else pd.DataFrame()
    result.obs["louvain_cluster"] = clusters

    # Store clustering parameters
    result.uns["louvain_clustering"] = {
        "resolution": resolution,
        "n_neighbors": n_neighbors,
        "random_state": random_state,
        "n_clusters": len(np.unique(clusters)),
        "use_weights": use_weights,
    }

    logger.info(f"Louvain clustering completed: found {len(np.unique(clusters))} clusters")
    return result


def kmeans_clustering(
    data: SingleCellData, n_clusters: int = 10, random_state: int | None = None, n_init: int = 10
) -> SingleCellData:
    """Perform K-means clustering on single-cell data.

    Args:
        data: SingleCellData object with preprocessed expression matrix
        n_clusters: Number of clusters
        random_state: Random seed for reproducibility
        n_init: Number of initializations for K-means

    Returns:
        SingleCellData with cluster assignments added to obs

    Raises:
        TypeError: If data is not SingleCellData
        ImportError: If scikit-learn not available
    """
    if not HAS_SKLEARN:
        raise ImportError(
            "scikit-learn is required for K-means clustering. " "Install with: uv pip install scikit-learn"
        )

    validation.validate_type(data, SingleCellData, "data")
    validation.validate_range(n_clusters, min_val=2, name="n_clusters")

    logger.info(f"Performing K-means clustering with {n_clusters} clusters")

    # Create copy
    result = data.copy()

    # Get expression matrix
    X = data.X.toarray() if hasattr(data.X, "toarray") else data.X

    # Perform K-means clustering
    kmeans = KMeans(n_clusters=n_clusters, random_state=random_state, n_init=n_init)
    clusters = kmeans.fit_predict(X)

    # Add to obs
    result.obs = result.obs.copy() if result.obs is not None else pd.DataFrame()
    result.obs["kmeans_cluster"] = clusters

    # Store clustering parameters and results
    result.uns["kmeans_clustering"] = {
        "n_clusters": n_clusters,
        "random_state": random_state,
        "n_init": n_init,
        "inertia": kmeans.inertia_,
        "n_iter": kmeans.n_iter_,
        "cluster_centers_shape": kmeans.cluster_centers_.shape,
    }

    logger.info(f"K-means clustering completed: found {n_clusters} clusters")
    return result


def hierarchical_clustering(
    data: SingleCellData, n_clusters: int = 10, linkage_method: str = "ward", metric: str = "euclidean"
) -> SingleCellData:
    """Perform hierarchical clustering on single-cell data.

    Args:
        data: SingleCellData object with preprocessed expression matrix
        n_clusters: Number of clusters
        linkage_method: Linkage method for hierarchical clustering
        metric: Distance metric

    Returns:
        SingleCellData with cluster assignments added to obs

    Raises:
        TypeError: If data is not SingleCellData
        ValueError: If linkage_method or metric is invalid
    """
    validation.validate_type(data, SingleCellData, "data")
    validation.validate_range(n_clusters, min_val=2, name="n_clusters")

    valid_linkages = ["single", "complete", "average", "weighted", "centroid", "median", "ward"]
    if linkage_method not in valid_linkages:
        raise errors.ValidationError(f"Invalid linkage method: {linkage_method}")

    logger.info(f"Performing hierarchical clustering with {n_clusters} clusters using {linkage_method} linkage")

    # Create copy
    result = data.copy()

    # Get expression matrix
    X = data.X.toarray() if hasattr(data.X, "toarray") else data.X

    # For Ward linkage, metric must be euclidean
    if linkage_method == "ward":
        metric = "euclidean"

    # Compute distance matrix
    if X.shape[0] > 1000:
        # For large datasets, use random subset for hierarchical clustering
        logger.warning(f"Large dataset ({X.shape[0]} cells), using random subset for hierarchical clustering")
        subset_indices = np.random.choice(X.shape[0], size=min(1000, X.shape[0]), replace=False)
        X_subset = X[subset_indices]

        # Compute linkage on subset
        distances = pdist(X_subset, metric=metric)
        linkage_matrix = linkage(distances, method=linkage_method)

        # Cut tree to get clusters for subset
        subset_clusters = fcluster(linkage_matrix, n_clusters, criterion="maxclust")

        # Assign clusters back to full dataset using nearest neighbors
        from sklearn.neighbors import NearestNeighbors

        nbrs = NearestNeighbors(n_neighbors=1).fit(X_subset)
        distances, indices = nbrs.kneighbors(X)
        clusters = subset_clusters[indices.flatten()]

    else:
        # Compute full hierarchical clustering
        distances = pdist(X, metric=metric)
        linkage_matrix = linkage(distances, method=linkage_method)
        clusters = fcluster(linkage_matrix, n_clusters, criterion="maxclust")

    # Add to obs
    result.obs = result.obs.copy() if result.obs is not None else pd.DataFrame()
    result.obs["hierarchical_cluster"] = clusters

    # Store clustering parameters
    result.uns["hierarchical_clustering"] = {
        "n_clusters": n_clusters,
        "linkage_method": linkage_method,
        "metric": metric,
        "linkage_matrix_shape": linkage_matrix.shape if "linkage_matrix" in locals() else None,
    }

    logger.info(f"Hierarchical clustering completed: found {n_clusters} clusters")
    return result


def find_marker_genes(data: SingleCellData, groupby: str, method: str = "t-test", n_genes: int = 100) -> pd.DataFrame:
    """Find marker genes for clusters or groups.

    Args:
        data: SingleCellData object with cluster assignments
        groupby: Column in obs containing group/cluster labels
        method: Statistical method for differential expression
        n_genes: Number of top marker genes to return per group

    Returns:
        DataFrame with marker gene statistics

    Raises:
        TypeError: If data is not SingleCellData
        ValueError: If groupby column not found
    """
    validation.validate_type(data, SingleCellData, "data")

    if data.obs is None or groupby not in data.obs.columns:
        raise errors.ValidationError(f"Group column '{groupby}' not found in data.obs")

    logger.info(f"Finding marker genes using {method} method")

    # Get expression matrix
    X = data.X.toarray() if hasattr(data.X, "toarray") else data.X
    groups = data.obs[groupby].values
    unique_groups = np.unique(groups)

    marker_results = []

    for group in unique_groups:
        # Split data by group
        group_mask = groups == group
        in_group = X[group_mask]
        out_group = X[~group_mask]

        # Calculate statistics for each gene
        for gene_idx in range(X.shape[1]):
            gene_name = data.var.index[gene_idx] if data.var is not None else f"gene_{gene_idx}"

            expr_in = in_group[:, gene_idx]
            expr_out = out_group[:, gene_idx]

            if method == "t-test":
                # Simple t-test implementation
                mean_in = np.mean(expr_in)
                mean_out = np.mean(expr_out)
                std_in = np.std(expr_in, ddof=1)
                std_out = np.std(expr_out, ddof=1)

                n_in = len(expr_in)
                n_out = len(expr_out)

                if std_in == 0 and std_out == 0:
                    t_stat = 0
                    p_val = 1.0
                elif std_in == 0 or std_out == 0:
                    # Large difference if one group has no variation
                    t_stat = 1000 if mean_in > mean_out else -1000
                    p_val = 1e-10
                else:
                    # Calculate t-statistic
                    se_diff = np.sqrt(std_in**2 / n_in + std_out**2 / n_out)
                    t_stat = (mean_in - mean_out) / se_diff if se_diff > 0 else 0

                    # Approximate p-value (simplified)
                    df = min(n_in, n_out) - 1
                    p_val = 2 * (1 - _student_t_cdf(abs(t_stat), df))

            elif method == "wilcoxon":
                # Simplified Wilcoxon rank-sum test
                from scipy.stats import ranksums

                try:
                    t_stat, p_val = ranksums(expr_in, expr_out)
                except (ValueError, TypeError, RuntimeError):
                    t_stat, p_val = 0, 1.0

            else:
                raise errors.ValidationError(f"Unsupported method: {method}")

            # Calculate fold change and other metrics
            mean_expr = np.mean(expr_in)
            mean_expr_other = np.mean(expr_out)
            fold_change = mean_expr / (mean_expr_other + 1e-10)
            log_fold_change = np.log2(fold_change) if fold_change > 0 else 0

            pct_expressed = np.mean(expr_in > 0) * 100

            marker_results.append(
                {
                    "gene": gene_name,
                    "group": group,
                    "mean_expr": mean_expr,
                    "mean_expr_other": mean_expr_other,
                    "log_fold_change": log_fold_change,
                    "pct_expressed": pct_expressed,
                    "statistic": t_stat,
                    "p_value": p_val,
                    "method": method,
                }
            )

    # Convert to DataFrame and sort
    results_df = pd.DataFrame(marker_results)
    results_df = results_df.sort_values(["group", "p_value"])

    # Keep top N genes per group
    top_markers = []
    for group in unique_groups:
        group_markers = results_df[results_df["group"] == group].head(n_genes)
        top_markers.append(group_markers)

    final_results = pd.concat(top_markers, ignore_index=True)

    logger.info(f"Found {len(final_results)} marker genes across {len(unique_groups)} groups")
    return final_results


def compute_cluster_composition(data: SingleCellData, groupby: str, cluster_col: str = "cluster") -> pd.DataFrame:
    """Compute composition of clusters by categorical variables.

    Args:
        data: SingleCellData object with cluster assignments
        groupby: Column in obs to compute composition by
        cluster_col: Column containing cluster assignments

    Returns:
        DataFrame with cluster composition statistics

    Raises:
        TypeError: If data is not SingleCellData
        ValueError: If required columns not found
    """
    validation.validate_type(data, SingleCellData, "data")

    if data.obs is None:
        raise errors.ValidationError("data.obs is required for cluster composition analysis")

    if groupby not in data.obs.columns:
        raise errors.ValidationError(f"Group column '{groupby}' not found in data.obs")

    if cluster_col not in data.obs.columns:
        raise errors.ValidationError(f"Cluster column '{cluster_col}' not found in data.obs")

    logger.info(f"Computing cluster composition by {groupby}")

    # Create contingency table
    contingency = pd.crosstab(data.obs[cluster_col], data.obs[groupby])

    # Convert to proportions
    proportions = contingency.div(contingency.sum(axis=1), axis=0)

    # Calculate statistics
    composition_stats = []

    for cluster in contingency.index:
        cluster_data = contingency.loc[cluster]
        cluster_props = proportions.loc[cluster]

        for category in contingency.columns:
            count = cluster_data[category]
            proportion = cluster_props[category]

            composition_stats.append(
                {
                    "cluster": cluster,
                    "category": category,
                    "count": count,
                    "proportion": proportion,
                    "percentage": proportion * 100,
                }
            )

    results_df = pd.DataFrame(composition_stats)

    logger.info("Cluster composition analysis completed")
    return results_df


def compute_cluster_silhouette(
    data: SingleCellData, cluster_col: str = "cluster", metric: str = "euclidean", sample_size: int | None = None
) -> Dict[str, float]:
    """Compute silhouette scores for cluster evaluation.

    Args:
        data: SingleCellData object with cluster assignments
        cluster_col: Column containing cluster assignments
        metric: Distance metric for silhouette calculation
        sample_size: Sample size for large datasets (for performance)

    Returns:
        Dictionary with silhouette scores and statistics

    Raises:
        TypeError: If data is not SingleCellData
        ValueError: If cluster column not found
        ImportError: If scikit-learn not available
    """
    if not HAS_SKLEARN:
        raise ImportError(
            "scikit-learn is required for silhouette analysis. " "Install with: uv pip install scikit-learn"
        )

    validation.validate_type(data, SingleCellData, "data")

    if data.obs is None or cluster_col not in data.obs.columns:
        raise errors.ValidationError(f"Cluster column '{cluster_col}' not found in data.obs")

    logger.info("Computing cluster silhouette scores")

    # Get data
    X = data.X.toarray() if hasattr(data.X, "toarray") else data.X
    clusters = data.obs[cluster_col].values

    # Sample data if requested (for performance)
    if sample_size is not None and sample_size < X.shape[0]:
        indices = np.random.choice(X.shape[0], size=sample_size, replace=False)
        X = X[indices]
        clusters = clusters[indices]

    # Check if we have enough clusters
    n_unique_clusters = len(np.unique(clusters))
    if n_unique_clusters < 2:
        logger.warning("Need at least 2 clusters for silhouette analysis")
        return {"error": "Need at least 2 clusters"}

    try:
        # Compute silhouette scores
        silhouette_avg = silhouette_score(X, clusters, metric=metric)
        sample_silhouette_values = silhouette_score(X, clusters, metric=metric)

        # Per-cluster silhouette scores
        per_cluster_scores = {}
        for cluster_id in np.unique(clusters):
            cluster_mask = clusters == cluster_id
            if np.sum(cluster_mask) > 1:  # Need at least 2 samples per cluster
                cluster_silhouettes = sample_silhouette_values[cluster_mask]
                per_cluster_scores[int(cluster_id)] = {
                    "mean": float(np.mean(cluster_silhouettes)),
                    "std": float(np.std(cluster_silhouettes)),
                    "min": float(np.min(cluster_silhouettes)),
                    "max": float(np.max(cluster_silhouettes)),
                }

        results = {
            "overall_silhouette_score": float(silhouette_avg),
            "n_clusters": n_unique_clusters,
            "n_samples": X.shape[0],
            "metric": metric,
            "per_cluster_scores": per_cluster_scores,
        }

        # Additional cluster quality metrics
        try:
            ch_score = calinski_harabasz_score(X, clusters)
            db_score = davies_bouldin_score(X, clusters)

            results["calinski_harabasz_score"] = float(ch_score)
            results["davies_bouldin_score"] = float(db_score)

        except Exception as e:
            logger.debug(f"Could not compute additional cluster metrics: {e}")

        logger.info(f"Silhouette analysis completed: overall score = {silhouette_avg:.3f}")
        return results

    except Exception as e:
        logger.error(f"Silhouette analysis failed: {e}")
        return {"error": str(e)}


def _compute_pairwise_distances(X: np.ndarray) -> np.ndarray:
    """Compute pairwise distances between cells."""
    # Use Euclidean distance for simplicity
    # For large datasets, consider using approximate methods
    from scipy.spatial.distance import cdist

    return cdist(X, X, metric="euclidean")


def _build_knn_adjacency(distances: np.ndarray, n_neighbors: int) -> sparse.csr_matrix:
    """Build kNN adjacency matrix from distance matrix."""
    n_cells = distances.shape[0]

    # Find k nearest neighbors for each cell
    adjacency = sparse.lil_matrix((n_cells, n_cells))

    for i in range(n_cells):
        # Sort distances and get indices of k+1 smallest (including self)
        neighbor_indices = np.argsort(distances[i])[: n_neighbors + 1]

        # Remove self from neighbors
        neighbor_indices = neighbor_indices[neighbor_indices != i][:n_neighbors]

        # Add edges with weights (inverse distance)
        for j in neighbor_indices:
            weight = 1.0 / (distances[i, j] + 1e-10)  # Avoid division by zero
            adjacency[i, j] = weight
            adjacency[j, i] = weight  # Undirected graph

    return adjacency.tocsr()


def _student_t_cdf(t: float, df: int) -> float:
    """Approximate cumulative distribution function for Student's t-distribution."""
    # Simplified approximation
    if df <= 1:
        return 0.5
    elif df == 2:
        return 0.5 + (t / np.sqrt(2 + t**2)) / (2 * np.sqrt(2))
    else:
        # Use normal approximation for large df
        return 0.5 * (1 + np.sign(t) * np.sqrt(1 - np.exp(-2 * t**2 / np.pi)))


def evaluate_clustering_performance(
    data: SingleCellData, cluster_col: str = "cluster", ground_truth_col: Optional[str] = None
) -> Dict[str, Any]:
    """Evaluate clustering performance with various metrics.

    Args:
        data: SingleCellData object with cluster assignments
        cluster_col: Column containing cluster assignments
        ground_truth_col: Optional column with ground truth labels for supervised evaluation

    Returns:
        Dictionary with clustering performance metrics

    Raises:
        TypeError: If data is not SingleCellData
    """
    validation.validate_type(data, SingleCellData, "data")

    if data.obs is None or cluster_col not in data.obs.columns:
        raise errors.ValidationError(f"Cluster column '{cluster_col}' not found in data.obs")

    logger.info("Evaluating clustering performance")

    clusters = data.obs[cluster_col].values
    n_clusters = len(np.unique(clusters))

    results = {
        "n_clusters": n_clusters,
        "cluster_sizes": [int(np.sum(clusters == c)) for c in np.unique(clusters)],
    }

    # Unsupervised metrics
    X = data.X.toarray() if hasattr(data.X, "toarray") else data.X

    try:
        silhouette = compute_cluster_silhouette(data, cluster_col)
        results["silhouette_analysis"] = silhouette
    except Exception as e:
        logger.debug(f"Silhouette analysis failed: {e}")
        results["silhouette_error"] = str(e)

    # Supervised metrics if ground truth available
    if ground_truth_col and ground_truth_col in data.obs.columns:
        from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, homogeneity_score

        true_labels = data.obs[ground_truth_col].values

        try:
            ari = adjusted_rand_score(true_labels, clusters)
            nmi = normalized_mutual_info_score(true_labels, clusters)
            homogeneity = homogeneity_score(true_labels, clusters)

            results["supervised_metrics"] = {
                "adjusted_rand_index": float(ari),
                "normalized_mutual_info": float(nmi),
                "homogeneity_score": float(homogeneity),
            }
        except Exception as e:
            logger.debug(f"Supervised metrics failed: {e}")
            results["supervised_metrics_error"] = str(e)

    logger.info("Clustering evaluation completed")
    return results
