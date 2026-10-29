"""Single-cell clustering algorithms and cluster analysis.

This module provides functions for clustering single-cell data using various
algorithms including graph-based methods (Leiden, Louvain) and traditional
methods (K-means, hierarchical). It also includes cluster evaluation and
marker gene identification.
"""

from __future__ import annotations

import warnings
from typing import Any, Dict, Literal, NamedTuple, Optional, get_args

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist
from scipy.stats import ranksums

from metainformant.core.data import validation
from metainformant.core.utils import errors, logging
from metainformant.singlecell.data.preprocessing import SingleCellData
from metainformant.singlecell.differential.expression import _welch_t_test

# Optional scientific dependencies
try:
    from sklearn.cluster import KMeans
    from sklearn.metrics import (
        calinski_harabasz_score,
        davies_bouldin_score,
        silhouette_samples,
        silhouette_score,
    )

    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False
    KMeans = None
    silhouette_score = None
    silhouette_samples = None
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

# leidenalg 0.10.2 contains invalid escape sequences in docstrings. Under the
# repository's warning-as-error policy Python can surface those as SyntaxError
# during import; keep Leiden an explicit optional capability in that case.
try:
    import leidenalg

    HAS_LEIDEN = True
except (ImportError, SyntaxError):
    leidenalg = None
    HAS_LEIDEN = False

try:
    import community as community_louvain

    HAS_LOUVAIN = True
except ImportError:
    community_louvain = None
    HAS_LOUVAIN = False

logger = logging.get_logger(__name__)

_LinkageMethod = Literal["single", "complete", "average", "weighted", "centroid", "median", "ward"]
_VALID_LINKAGE_METHODS: tuple[_LinkageMethod, ...] = get_args(_LinkageMethod)


def _as_dense_matrix(matrix: Any) -> np.ndarray:
    """Convert supported matrix-like inputs to a dense float array."""
    dense = matrix.toarray() if hasattr(matrix, "toarray") else matrix
    return np.asarray(dense, dtype=float)


def _cluster_feature_matrix(
    feature_matrix: Any,
    *,
    resolution: float = 1.0,
    random_state: int | None = None,
) -> np.ndarray:
    """Deterministically cluster a raw feature matrix for lightweight callers."""
    X = _as_dense_matrix(feature_matrix)
    if X.ndim != 2:
        raise ValueError("feature_matrix must be a 2D matrix")
    if X.shape[0] == 0:
        return np.array([], dtype=int)
    if X.shape[0] == 1:
        return np.zeros(1, dtype=int)

    n_clusters = max(2, int(round(max(resolution, 0.1) * 4)))
    n_clusters = min(n_clusters, X.shape[0], 10)

    if HAS_SKLEARN:
        kmeans = KMeans(n_clusters=n_clusters, random_state=0 if random_state is None else random_state, n_init=10)
        labels: np.ndarray = kmeans.fit_predict(X).astype(int)
        return labels

    centered = X - np.mean(X, axis=0)
    try:
        _, _, vh = np.linalg.svd(centered, full_matrices=False)
        scores = centered @ vh[0]
    except np.linalg.LinAlgError:
        scores = centered[:, 0]

    quantiles = np.linspace(0.0, 1.0, n_clusters + 1)[1:-1]
    thresholds = np.quantile(scores, quantiles)
    digitized: np.ndarray = np.digitize(scores, thresholds).astype(int)
    return digitized


def _leiden_stack_available() -> bool:
    """Whether networkx, python-igraph and leidenalg can all be imported."""
    return bool(HAS_NETWORKX and HAS_IGRAPH and HAS_LEIDEN)


class _LeidenProblem(NamedTuple):
    """Everything that is handed to leidenalg: the graph, the edge weights and the algorithm's parameters."""

    n_vertices: int
    edges: list[tuple[int, int]]
    weights: list[float] | None
    resolution: float
    seed: int | None


def _knn_graph_edges(X: np.ndarray, n_neighbors: int) -> tuple[list[tuple[int, int]], list[float]]:
    """Undirected kNN graph over the rows of ``X``: each edge once (``i < j``), with its inverse-distance weight."""
    adjacency = _build_knn_adjacency(_compute_pairwise_distances(X), n_neighbors)
    upper = sparse.triu(adjacency, k=1).tocoo()
    return list(zip(upper.row.tolist(), upper.col.tolist())), upper.data.tolist()


def _leiden_problem(
    X: np.ndarray, *, resolution: float, n_neighbors: int, random_state: int | None, use_weights: bool
) -> _LeidenProblem:
    """The kNN graph and parameters for leidenalg; the edge weights are included only when ``use_weights``."""
    edges, weights = _knn_graph_edges(X, n_neighbors)
    return _LeidenProblem(
        n_vertices=X.shape[0],
        edges=edges,
        weights=weights if use_weights else None,
        resolution=resolution,
        seed=random_state,
    )


def _leiden_labels(
    X: np.ndarray, *, resolution: float, n_neighbors: int, random_state: int | None, use_weights: bool = True
) -> np.ndarray:
    """Leiden communities of the kNN graph over the rows of ``X`` (needs the Leiden stack)."""
    problem = _leiden_problem(
        X, resolution=resolution, n_neighbors=n_neighbors, random_state=random_state, use_weights=use_weights
    )

    # Perform Leiden clustering on the weighted (or, without use_weights, unweighted) graph
    partition = leidenalg.find_partition(
        ig.Graph(n=problem.n_vertices, edges=problem.edges, directed=False),
        leidenalg.RBConfigurationVertexPartition,
        weights=problem.weights,
        resolution_parameter=problem.resolution,
        seed=problem.seed,
    )

    # Get cluster assignments
    return np.array(partition.membership)


def leiden_clustering(
    data: SingleCellData | Any | None = None,
    resolution: float = 1.0,
    n_neighbors: int = 15,
    random_state: int | None = None,
    use_weights: bool = True,
    *,
    feature_matrix: Any | None = None,
    allow_kmeans_fallback: bool = False,
) -> SingleCellData | np.ndarray:
    """Perform Leiden clustering on single-cell data.

    Builds a k-nearest-neighbour graph over the cells (or the rows of a plain matrix) and partitions it
    with the Leiden algorithm, which needs networkx, python-igraph and leidenalg. Another algorithm is
    never substituted silently: without those packages the function raises ``ImportError``, unless a
    matrix input is combined with ``allow_kmeans_fallback=True``.

    Args:
        data: SingleCellData object with preprocessed expression matrix, or a
            matrix-like object for compatibility callers
        resolution: Resolution parameter for Leiden algorithm
        n_neighbors: Number of neighbors for kNN graph construction
        random_state: Random seed for reproducibility
        use_weights: Weight each edge of the kNN graph by the inverse distance between its two cells (the
            weights are passed to leidenalg) when True; give every edge weight 1 when False. Recorded in ``uns``.
        feature_matrix: Matrix-style alias for callers that only need labels
        allow_kmeans_fallback: Matrix input only. When the Leiden stack is unavailable, cluster the rows
            with k-means (the number of clusters follows ``resolution``) instead of raising, and emit a
            ``UserWarning`` and a log warning that say so. A label array cannot carry the method, so
            rely on that warning (or on ``HAS_LEIDEN``) to know what ran. For a SingleCellData call
            :func:`kmeans_clustering` directly.

    Returns:
        SingleCellData with cluster assignments added to ``obs["leiden_cluster"]`` and the parameters,
        including ``method`` (always ``"leiden"``), in ``uns["leiden_clustering"]``; or a label array
        for matrix-style input.

    Raises:
        ImportError: If the Leiden stack is unavailable (and, for matrix input, the k-means fallback
            was not requested)
        ValueError: If a matrix input is not two-dimensional
    """
    if feature_matrix is None and isinstance(data, SingleCellData):
        cells: SingleCellData | None = data
        X = data.X.toarray() if hasattr(data.X, "toarray") else data.X
    else:
        cells = None
        X = _as_dense_matrix(feature_matrix if feature_matrix is not None else data)
        if X.ndim != 2:
            raise ValueError("feature_matrix must be a 2D matrix")

    if _leiden_stack_available():
        logger.info(f"Performing Leiden clustering with resolution={resolution}, n_neighbors={n_neighbors}")
        clusters = _leiden_labels(
            X, resolution=resolution, n_neighbors=n_neighbors, random_state=random_state, use_weights=use_weights
        )
    elif cells is None and allow_kmeans_fallback:
        message = (
            "Leiden clustering is unavailable (it needs networkx, python-igraph and leidenalg): "
            "clustered the feature matrix with k-means instead"
        )
        logger.warning(message)
        warnings.warn(message, UserWarning, stacklevel=2)
        return _cluster_feature_matrix(X, resolution=resolution, random_state=random_state)
    else:
        if cells is not None:
            hint = "Install them, or call kmeans_clustering to cluster a SingleCellData with k-means."
        else:
            hint = "Install them, or pass allow_kmeans_fallback=True to cluster the matrix with k-means instead."
        raise ImportError(f"Leiden clustering requires the networkx, python-igraph and leidenalg packages. {hint}")

    if cells is None:
        return clusters

    # Create copy
    result = cells.copy()

    # Add to obs
    result.obs = result.obs.copy() if result.obs is not None else pd.DataFrame()
    result.obs["leiden_cluster"] = clusters

    # Store clustering parameters
    result.uns["leiden_clustering"] = {
        "method": "leiden",
        "resolution": resolution,
        "n_neighbors": n_neighbors,
        "random_state": random_state,
        "n_clusters": len(np.unique(clusters)),
        "use_weights": use_weights,
    }

    logger.info(f"Leiden clustering completed: found {len(np.unique(clusters))} clusters")
    return result


def _louvain_graph(X: np.ndarray, n_neighbors: int, *, use_weights: bool) -> Any:
    """kNN graph over the rows of ``X`` as a networkx graph whose ``weight`` edge attribute python-louvain reads.

    The attribute is the inverse distance between the two cells with ``use_weights`` and 1 on every edge without.
    """
    adjacency = _build_knn_adjacency(_compute_pairwise_distances(X), n_neighbors)
    if not use_weights:
        adjacency.data[:] = 1.0
    return _graph_from_adjacency(adjacency)


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
        use_weights: Weight each edge of the kNN graph by the inverse distance between its two cells when True
            (python-louvain reads the ``weight`` edge attribute); give every edge weight 1 when False

    Returns:
        SingleCellData with cluster assignments added to obs

    Raises:
        ImportError: If required packages are not available
        ValidationError: If data is not SingleCellData
    """
    validation.validate_type(data, SingleCellData, "data")

    if not HAS_NETWORKX or not HAS_LOUVAIN:
        raise ImportError("Louvain clustering requires networkx and python-louvain packages")

    logger.info(f"Performing Louvain clustering with resolution={resolution}, n_neighbors={n_neighbors}")

    # Create copy
    result = data.copy()

    # Get expression matrix
    X = data.X.toarray() if hasattr(data.X, "toarray") else data.X

    # Construct the kNN graph as a networkx graph (unit edge weights unless use_weights)
    G = _louvain_graph(X, n_neighbors, use_weights=use_weights)

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
        ValidationError: If data is not SingleCellData or n_clusters is below 2
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


def find_markers(
    expression_matrix: Any,
    cluster_labels: Any,
    *,
    n_markers: int = 10,
    gene_names: list[str] | None = None,
) -> dict[Any, list[dict[str, float | int | str]]]:
    """Find simple marker genes from an expression matrix and cluster labels."""
    X = _as_dense_matrix(expression_matrix)
    labels = np.asarray(cluster_labels)

    if X.ndim != 2:
        raise ValueError("expression_matrix must be a 2D matrix")
    if X.shape[0] != labels.shape[0]:
        raise ValueError("Number of rows in expression_matrix must match cluster_labels")
    if n_markers <= 0:
        return {cluster.item() if hasattr(cluster, "item") else cluster: [] for cluster in np.unique(labels)}

    names = gene_names if gene_names is not None else [f"gene_{idx}" for idx in range(X.shape[1])]
    if len(names) != X.shape[1]:
        raise ValueError("gene_names length must match number of genes")

    markers: dict[Any, list[dict[str, float | int | str]]] = {}
    eps = 1e-10

    for cluster in np.unique(labels):
        cluster_mask = labels == cluster
        in_group = X[cluster_mask]
        out_group = X[~cluster_mask]
        mean_in = np.mean(in_group, axis=0)
        mean_out = np.mean(out_group, axis=0) if out_group.size else np.zeros(X.shape[1])
        log_fold_change = np.log2((mean_in + eps) / (mean_out + eps))
        scores = mean_in - mean_out
        ranked = np.argsort(scores)[::-1][: min(n_markers, X.shape[1])]

        cluster_key = cluster.item() if hasattr(cluster, "item") else cluster
        markers[cluster_key] = [
            {
                "gene": names[int(gene_idx)],
                "gene_index": int(gene_idx),
                "score": float(scores[gene_idx]),
                "mean_expr": float(mean_in[gene_idx]),
                "mean_expr_other": float(mean_out[gene_idx]),
                "log_fold_change": float(log_fold_change[gene_idx]),
            }
            for gene_idx in ranked
        ]

    return markers


def hierarchical_clustering(
    data: SingleCellData,
    n_clusters: int = 10,
    linkage_method: str = "ward",
    metric: str = "euclidean",
    *,
    random_state: int | None = None,
) -> SingleCellData:
    """Perform hierarchical clustering on single-cell data.

    Above 1000 cells the linkage is built on a random subset of 1000 cells and every other cell takes the
    cluster of its nearest subset cell.

    Args:
        data: SingleCellData object with preprocessed expression matrix
        n_clusters: Number of clusters
        linkage_method: Linkage method for hierarchical clustering
        metric: Distance metric
        random_state: Seed for the subset drawn above 1000 cells (``None`` draws a fresh,
            non-reproducible subset). numpy's global random state is never used.

    Returns:
        SingleCellData with cluster assignments added to obs

    Raises:
        ValidationError: If data is not SingleCellData, n_clusters is below 2 or linkage_method is invalid
        ValueError: If metric is not a distance metric that scipy knows (raised by scipy; ward linkage always uses
            euclidean)
    """
    validation.validate_type(data, SingleCellData, "data")
    validation.validate_range(n_clusters, min_val=2, name="n_clusters")

    if linkage_method not in _VALID_LINKAGE_METHODS:
        raise errors.ValidationError(f"Invalid linkage method: {linkage_method}")

    logger.info(f"Performing hierarchical clustering with {n_clusters} clusters using {linkage_method} linkage")

    # Create copy
    result = data.copy()

    # Get expression matrix
    X = data.X.toarray() if hasattr(data.X, "toarray") else data.X

    # For Ward linkage, metric must be euclidean
    if linkage_method == "ward":
        metric = "euclidean"

    # scipy rejects unknown metric names itself (ValueError); its stubs only accept a Literal of the known names.
    metric_name: Any = metric

    linkage_matrix = None
    # Compute distance matrix
    if X.shape[0] > 1000:
        # For large datasets, use random subset for hierarchical clustering
        logger.warning(f"Large dataset ({X.shape[0]} cells), using random subset for hierarchical clustering")
        subset_indices = np.random.default_rng(random_state).choice(
            X.shape[0], size=min(1000, X.shape[0]), replace=False
        )
        X_subset = X[subset_indices]

        # Compute linkage on subset
        distances = pdist(X_subset, metric=metric_name)
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
        distances = pdist(X, metric=metric_name)
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
        "linkage_matrix_shape": linkage_matrix.shape if linkage_matrix is not None else None,
        "random_state": random_state,
    }

    logger.info(f"Hierarchical clustering completed: found {n_clusters} clusters")
    return result


def _ranksums_outcome(a: np.ndarray, b: np.ndarray) -> tuple[float, float, str]:
    """Wilcoxon rank-sum z statistic, two-sided p-value and status of two samples.

    NaN, NaN and the reason are returned when the test is undefined.

    The test is undefined when a sample is empty, a value is not finite, or all values of both samples are tied
    (the rank sums then have no spread): no p-value is made up for those.
    """
    if a.size == 0 or b.size == 0:
        return np.nan, np.nan, "empty_group: a group has no observations"
    if not (np.all(np.isfinite(a)) and np.all(np.isfinite(b))):
        return np.nan, np.nan, "non_finite_values: the values contain NaN or infinity"
    if min(a.min(), b.min()) == max(a.max(), b.max()):
        return np.nan, np.nan, "zero_variance: every value of both groups is tied, so the rank sums have no spread"
    statistic, p_value = ranksums(a, b)
    return float(statistic), float(p_value), "ok"


def find_marker_genes(data: SingleCellData, groupby: str, method: str = "t-test", n_genes: int = 100) -> pd.DataFrame:
    """Find marker genes for clusters or groups.

    Each group is compared with all the other cells, gene by gene. ``"t-test"`` is Welch's t-test (exact Student
    t tail on the Welch-Satterthwaite degrees of freedom) and ``"wilcoxon"`` the rank-sum z test without a tie
    correction. A gene whose test is undefined for a group (constant in both, a group with fewer than two cells
    for the t-test, non-finite values) has ``p_value = NaN`` and a ``test_status`` that says why; it is ranked
    after every gene that has a p-value, and no p-value is made up for it.

    Args:
        data: SingleCellData object with cluster assignments
        groupby: Column in obs containing group/cluster labels
        method: Statistical method for differential expression ("t-test" or "wilcoxon")
        n_genes: Number of top marker genes to return per group, ranked by p-value (ties, such as p-values that
            underflow to 0, by the size of the statistic)

    Returns:
        DataFrame with marker gene statistics: gene, group, mean_expr, mean_expr_other, log_fold_change
        (``log2((mean_expr + 1e-10) / (mean_expr_other + 1e-10))``, NaN when that ratio is not positive),
        pct_expressed, statistic, p_value, method and test_status

    Raises:
        ValidationError: If data is not SingleCellData, the groupby column is not in data.obs or method is unsupported
    """
    validation.validate_type(data, SingleCellData, "data")

    if data.obs is None or groupby not in data.obs.columns:
        raise errors.ValidationError(f"Group column '{groupby}' not found in data.obs")
    if method not in ("t-test", "wilcoxon"):
        raise errors.ValidationError(f"Unsupported method: {method}")

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
                t_stat, p_val, test_status = _welch_t_test(expr_in, expr_out)
            else:
                t_stat, p_val, test_status = _ranksums_outcome(expr_in, expr_out)

            # Calculate fold change and other metrics. A gene absent from one side gets a large finite log ratio; a
            # ratio that is not positive (negative means, as in scaled data) has no logarithm: NaN, not 0
            mean_expr = np.mean(expr_in) if expr_in.size else np.nan
            mean_expr_other = np.mean(expr_out) if expr_out.size else np.nan
            ratio = (mean_expr + 1e-10) / (mean_expr_other + 1e-10)
            log_fold_change = np.log2(ratio) if ratio > 0 else np.nan

            pct_expressed = np.mean(expr_in > 0) * 100 if expr_in.size else np.nan

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
                    "test_status": test_status,
                }
            )

    # Convert to DataFrame and sort; an empty result set (e.g. zero genes)
    # would otherwise crash sort_values on missing columns
    if not marker_results:
        return pd.DataFrame(
            columns=[
                "gene",
                "group",
                "mean_expr",
                "mean_expr_other",
                "log_fold_change",
                "pct_expressed",
                "statistic",
                "p_value",
                "method",
                "test_status",
            ]
        )
    results_df = pd.DataFrame(marker_results)
    # Rank by p-value, then (p-values that underflow to 0 tie) by the size of the statistic; NaN sorts last
    results_df["_abs_statistic"] = results_df["statistic"].abs()
    results_df = results_df.sort_values(["group", "p_value", "_abs_statistic"], ascending=[True, True, False])
    results_df = results_df.drop(columns="_abs_statistic")

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
        ValidationError: If data is not SingleCellData, has no obs, or lacks the groupby or cluster column
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
    data: SingleCellData,
    cluster_col: str = "cluster",
    metric: str = "euclidean",
    sample_size: int | None = None,
    *,
    random_state: int | None = None,
) -> Dict[str, Any]:
    """Compute silhouette scores for cluster evaluation.

    Args:
        data: SingleCellData object with cluster assignments
        cluster_col: Column containing cluster assignments
        metric: Distance metric for silhouette calculation
        sample_size: Sample size for large datasets (for performance)
        random_state: Seed for the subsample drawn when ``sample_size`` is smaller than the number of
            cells (``None`` draws a fresh, non-reproducible subsample). numpy's global random state is
            never used.

    Returns:
        Dictionary with silhouette scores and statistics

    Raises:
        ValidationError: If data is not SingleCellData or the cluster column is not in data.obs
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
        indices = np.random.default_rng(random_state).choice(X.shape[0], size=sample_size, replace=False)
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
        sample_silhouette_values = silhouette_samples(X, clusters, metric=metric)

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


def _graph_from_adjacency(adjacency: Any) -> Any:
    """Build an undirected networkx graph from a weighted adjacency matrix.

    networkx 3.0 removed ``from_scipy_sparse_matrix``/``from_numpy_matrix``
    in favour of the array-based constructors; fall back to the legacy
    names so both major networkx lines are supported.
    """
    if sparse.issparse(adjacency):
        builder = getattr(nx, "from_scipy_sparse_array", None) or nx.from_scipy_sparse_matrix
    else:
        builder = getattr(nx, "from_numpy_array", None) or nx.from_numpy_matrix
    return builder(adjacency)


def _compute_pairwise_distances(X: np.ndarray) -> np.ndarray:
    """Compute pairwise distances between cells."""
    # Use Euclidean distance for simplicity
    # For large datasets, consider using approximate methods
    from scipy.spatial.distance import cdist

    pairwise: np.ndarray = cdist(X, X, metric="euclidean")
    return pairwise


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
        ValidationError: If data is not SingleCellData or the cluster column is not in data.obs
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

    try:
        silhouette = compute_cluster_silhouette(data, cluster_col)
        results["silhouette_analysis"] = silhouette
    except Exception as e:
        logger.debug(f"Silhouette analysis failed: {e}")
        results["silhouette_error"] = str(e)

    # Supervised metrics if ground truth available
    if ground_truth_col and ground_truth_col in data.obs.columns:
        from sklearn.metrics import adjusted_rand_score, homogeneity_score, normalized_mutual_info_score

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
