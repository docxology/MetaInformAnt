"""Spatially-aware clustering for spatial transcriptomics.

Implements spatial graph construction (KNN, Delaunay, radius) and community
detection algorithms (Leiden, Louvain) that incorporate spatial proximity
into transcriptomic clustering. Also provides a BayesSpace-inspired spatial
domain identification approach.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Optional dependencies
try:
    import numpy as np
    from numpy.typing import NDArray
except ImportError:
    np = None  # type: ignore[assignment]
    NDArray = None  # type: ignore[assignment,misc]

try:
    from scipy import sparse as sp_sparse
    from scipy.spatial import Delaunay, KDTree
except ImportError:
    sp_sparse = None  # type: ignore[assignment]
    Delaunay = None  # type: ignore[assignment,misc]
    KDTree = None  # type: ignore[assignment,misc]

try:
    from sklearn.cluster import KMeans
    from sklearn.decomposition import PCA
    from sklearn.neighbors import NearestNeighbors
    from sklearn.preprocessing import StandardScaler
except ImportError:
    KMeans = None  # type: ignore[assignment,misc]
    PCA = None  # type: ignore[assignment,misc]
    NearestNeighbors = None  # type: ignore[assignment,misc]
    StandardScaler = None  # type: ignore[assignment,misc]


@dataclass
class SpatialClusterResult:
    """Result of spatial clustering.

    Attributes:
        labels: Cluster label per spot/cell (integer array of length n).
        n_clusters: Number of clusters found.
        method: Clustering method used.
        modularity: Modularity score (for graph-based methods).
        spatial_graph: Adjacency matrix used for clustering.
        metadata: Additional result metadata.
    """

    labels: Any  # np.ndarray of shape (n,)
    n_clusters: int
    method: str
    modularity: float = 0.0
    spatial_graph: Any = None  # scipy sparse matrix
    metadata: dict[str, Any] = field(default_factory=dict)


def build_spatial_graph(
    coordinates: Any,
    method: Literal["knn", "delaunay", "radius"] = "knn",
    n_neighbors: int = 6,
    radius: float | None = None,
) -> Any:
    """Build a spatial neighborhood graph from coordinates.

    Args:
        coordinates: Array of shape (n, 2) with spatial coordinates.
        method: Graph construction method:
            - "knn": K-nearest neighbors graph.
            - "delaunay": Delaunay triangulation graph.
            - "radius": Fixed-radius neighborhood graph.
        n_neighbors: Number of neighbors for KNN (default 6 for hexagonal Visium grid).
        radius: Radius for radius-based graph (required if method="radius").

    Returns:
        Sparse adjacency matrix (scipy CSR) of shape (n, n).

    Raises:
        ImportError: If scipy or sklearn is not installed.
        ValueError: If method is "radius" but no radius is provided.
    """
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")
    if sp_sparse is None:
        raise ImportError("scipy is required: uv pip install scipy")

    coords = np.asarray(coordinates, dtype=np.float64)
    n = coords.shape[0]

    if method == "knn":
        if NearestNeighbors is None:
            raise ImportError("scikit-learn is required for KNN graph: uv pip install scikit-learn")

        # Build symmetric KNN graph
        nn = NearestNeighbors(n_neighbors=min(n_neighbors + 1, n), metric="euclidean")
        nn.fit(coords)
        distances, indices = nn.kneighbors(coords)

        rows: list[int] = []
        cols: list[int] = []
        vals: list[float] = []

        for i in range(n):
            for j_idx in range(1, indices.shape[1]):  # skip self (index 0)
                j = indices[i, j_idx]
                d = distances[i, j_idx]
                if d > 0:
                    weight = 1.0 / (1.0 + d)  # distance-to-similarity
                    rows.append(i)
                    cols.append(j)
                    vals.append(weight)

        adj = sp_sparse.csr_matrix((vals, (rows, cols)), shape=(n, n))
        # Make symmetric
        adj = adj + adj.T
        adj.data[:] = np.minimum(adj.data, 1.0)  # cap at 1

    elif method == "delaunay":
        if Delaunay is None:
            raise ImportError("scipy.spatial is required for Delaunay: uv pip install scipy")

        if n < 3:
            # Delaunay needs at least 3 non-collinear points
            return sp_sparse.csr_matrix((n, n))

        tri = Delaunay(coords)
        edges: set[tuple[int, int]] = set()
        for simplex in tri.simplices:
            for k in range(3):
                i, j = int(simplex[k]), int(simplex[(k + 1) % 3])
                edge = (min(i, j), max(i, j))
                edges.add(edge)

        rows_list: list[int] = []
        cols_list: list[int] = []
        for i, j in edges:
            d = float(np.linalg.norm(coords[i] - coords[j]))
            rows_list.extend([i, j])
            cols_list.extend([j, i])

        vals_list = [1.0] * len(rows_list)
        adj = sp_sparse.csr_matrix((vals_list, (rows_list, cols_list)), shape=(n, n))

    elif method == "radius":
        if radius is None:
            raise ValueError("radius must be specified for method='radius'")
        if KDTree is None:
            raise ImportError("scipy.spatial is required: uv pip install scipy")

        tree = KDTree(coords)
        pairs = tree.query_pairs(r=radius)

        rows_list = []
        cols_list = []
        for i, j in pairs:
            rows_list.extend([i, j])
            cols_list.extend([j, i])

        vals_list = [1.0] * len(rows_list)
        adj = sp_sparse.csr_matrix((vals_list, (rows_list, cols_list)), shape=(n, n))

    else:
        raise ValueError(f"Unknown graph method: {method}. Use 'knn', 'delaunay', or 'radius'.")

    n_edges = adj.nnz // 2
    logger.info(f"Built spatial graph ({method}): {n} nodes, {n_edges} edges")
    return adj


def leiden_clustering(
    adjacency_matrix: Any,
    resolution: float = 1.0,
    n_iterations: int = -1,
    seed: int = 42,
) -> tuple[Any, float]:
    """Leiden community detection on a graph.

    Implements the Leiden algorithm using a modularity-optimization approach.
    Falls back to a greedy modularity maximization if the leidenalg package
    is not available.

    Args:
        adjacency_matrix: Sparse adjacency matrix (n x n).
        resolution: Resolution parameter (higher = more clusters).
        n_iterations: Number of iterations (-1 for until convergence).
        seed: Random seed for reproducibility.

    Returns:
        Tuple of (labels, modularity) where labels is an integer array
        and modularity is the partition quality score.

    Raises:
        ImportError: If neither leidenalg nor required fallback deps are available.
    """
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")
    if sp_sparse is None:
        raise ImportError("scipy is required: uv pip install scipy")

    adj = sp_sparse.csr_matrix(adjacency_matrix)
    n = adj.shape[0]

    # Try leidenalg first (optimal implementation)
    try:
        import igraph as ig
        import leidenalg

        # Convert scipy sparse to igraph
        sources, targets = adj.nonzero()
        weights = np.array(adj[sources, targets]).flatten()
        mask = sources < targets  # upper triangle only for undirected
        edges = list(zip(sources[mask].tolist(), targets[mask].tolist()))
        edge_weights = weights[mask].tolist()

        g = ig.Graph(n=n, edges=edges, directed=False)
        g.es["weight"] = edge_weights

        partition = leidenalg.find_partition(
            g,
            leidenalg.RBConfigurationVertexPartition,
            weights="weight",
            resolution_parameter=resolution,
            n_iterations=n_iterations if n_iterations > 0 else 2,
            seed=seed,
        )

        labels = np.array(partition.membership, dtype=np.int32)
        modularity = float(partition.modularity)
        logger.info(f"Leiden clustering: {len(set(labels))} clusters, modularity={modularity:.4f}")
        return labels, modularity

    except ImportError:
        logger.info("leidenalg not available; using greedy modularity fallback")

    # Fallback: greedy modularity optimization (simplified Louvain-like)
    return _greedy_modularity_clustering(adj, resolution=resolution, seed=seed)


def louvain_clustering(
    adjacency_matrix: Any,
    resolution: float = 1.0,
    seed: int = 42,
) -> tuple[Any, float]:
    """Louvain community detection on a graph.

    Implements the Louvain method for modularity-based community detection.
    Uses community_louvain if available, otherwise falls back to a greedy approach.

    Args:
        adjacency_matrix: Sparse adjacency matrix (n x n).
        resolution: Resolution parameter (higher = more clusters).
        seed: Random seed.

    Returns:
        Tuple of (labels, modularity).
    """
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")
    if sp_sparse is None:
        raise ImportError("scipy is required: uv pip install scipy")

    adj = sp_sparse.csr_matrix(adjacency_matrix)
    n = adj.shape[0]

    # Try python-louvain (community package) first
    try:
        import community as community_louvain
        import networkx as nx

        # Convert to networkx graph
        G = nx.Graph()
        G.add_nodes_from(range(n))
        sources, targets = adj.nonzero()
        weights = np.array(adj[sources, targets]).flatten()
        for s, t, w in zip(sources.tolist(), targets.tolist(), weights.tolist()):
            if s < t:
                G.add_edge(s, t, weight=w)

        partition = community_louvain.best_partition(
            G, resolution=resolution, random_state=seed
        )
        labels = np.array([partition[i] for i in range(n)], dtype=np.int32)
        modularity = float(community_louvain.modularity(partition, G))
        logger.info(f"Louvain clustering: {len(set(labels))} clusters, modularity={modularity:.4f}")
        return labels, modularity

    except ImportError:
        logger.info("python-louvain not available; using greedy modularity fallback")

    return _greedy_modularity_clustering(adj, resolution=resolution, seed=seed)


def _greedy_modularity_clustering(
    adj: Any,
    resolution: float = 1.0,
    seed: int = 42,
    max_iter: int = 50,
) -> tuple[Any, float]:
    """Greedy modularity maximization as fallback clustering.

    Implements a simplified agglomerative modularity optimization:
    1. Start with each node in its own community.
    2. Repeatedly move nodes to the community that maximizes modularity gain.
    3. Stop when no improvement is found.

    This is a simplified version of the Louvain first phase.

    Args:
        adj: Sparse adjacency matrix.
        resolution: Resolution parameter.
        seed: Random seed.
        max_iter: Maximum iterations.

    Returns:
        Tuple of (labels, modularity).
    """
    rng = np.random.RandomState(seed)
    n = adj.shape[0]
    m = adj.sum() / 2.0  # total edge weight

    if m == 0:
        labels = np.zeros(n, dtype=np.int32)
        return labels, 0.0

    # Degree (strength) of each node
    degrees = np.array(adj.sum(axis=1)).flatten()

    # Initialize: each node in its own community
    labels = np.arange(n, dtype=np.int32)

    for iteration in range(max_iter):
        improved = False
        node_order = rng.permutation(n)

        for node in node_order:
            current_comm = labels[node]

            # Get neighbors
            row = adj.getrow(node)
            neighbor_indices = row.indices
            neighbor_weights = row.data

            # Compute delta-Q for moving node to each neighbor community
            # Sum of weights to each community
            comm_weights: dict[int, float] = {}
            for ni, nw in zip(neighbor_indices, neighbor_weights):
                nc = labels[ni]
                comm_weights[nc] = comm_weights.get(nc, 0.0) + nw

            ki = degrees[node]
            best_comm = current_comm
            best_delta_q = 0.0

            for comm, w_ic in comm_weights.items():
                if comm == current_comm:
                    continue

                # Sum of degrees in target community
                comm_mask = labels == comm
                sigma_tot = degrees[comm_mask].sum()

                # Sum of degrees in current community (excluding node)
                cur_mask = labels == current_comm
                cur_sigma = degrees[cur_mask].sum() - ki

                # Weight from node to current community
                w_current = comm_weights.get(current_comm, 0.0)

                # Delta modularity (resolution-scaled)
                delta_q = (
                    (w_ic - w_current) / m
                    - resolution * ki * (sigma_tot - cur_sigma) / (2.0 * m * m)
                )

                if delta_q > best_delta_q:
                    best_delta_q = delta_q
                    best_comm = comm

            if best_comm != current_comm:
                labels[node] = best_comm
                improved = True

        if not improved:
            break

    # Relabel communities to 0..k-1
    unique_labels = np.unique(labels)
    label_map = {old: new for new, old in enumerate(unique_labels)}
    labels = np.array([label_map[l] for l in labels], dtype=np.int32)

    # Compute modularity
    modularity = _compute_modularity(adj, labels, m)
    n_clusters = len(unique_labels)
    logger.info(
        f"Greedy modularity clustering: {n_clusters} clusters, "
        f"modularity={modularity:.4f} (resolution={resolution})"
    )
    return labels, modularity


def _compute_modularity(adj: Any, labels: Any, m: float) -> float:
    """Compute Newman-Girvan modularity Q for a partition.

    Q = (1/2m) * sum_ij [ A_ij - k_i*k_j/(2m) ] * delta(c_i, c_j)

    Args:
        adj: Sparse adjacency matrix.
        labels: Cluster label array.
        m: Total edge weight (sum of all edges).

    Returns:
        Modularity score in [-1, 1].
    """
    if m == 0:
        return 0.0

    degrees = np.array(adj.sum(axis=1)).flatten()
    Q = 0.0

    for c in np.unique(labels):
        mask = labels == c
        # Sum of edges within community
        subgraph = adj[np.ix_(mask, mask)]
        e_c = subgraph.sum() / 2.0
        # Sum of degrees in community
        a_c = degrees[mask].sum()
        Q += e_c / m - (a_c / (2.0 * m)) ** 2

    return float(Q)


def spatial_cluster(
    expression: Any,
    coordinates: Any,
    n_clusters: int | None = None,
    method: Literal["leiden", "louvain", "kmeans"] = "leiden",
    *,
    n_neighbors: int = 6,
    graph_method: Literal["knn", "delaunay", "radius"] = "knn",
    resolution: float = 1.0,
    spatial_weight: float = 0.5,
    n_pcs: int = 15,
    seed: int = 42,
) -> SpatialClusterResult:
    """Spatially-aware clustering combining expression and spatial proximity.

    Constructs a joint graph that balances transcriptomic similarity with spatial
    proximity, then applies community detection or KMeans.

    For graph-based methods (leiden, louvain):
        1. Build expression similarity graph (KNN in PCA space).
        2. Build spatial proximity graph.
        3. Combine: A_combined = (1-w)*A_expr + w*A_spatial.
        4. Run community detection on combined graph.

    For KMeans:
        1. PCA on expression.
        2. Concatenate scaled PCA coordinates with scaled spatial coordinates.
        3. Run KMeans on concatenated features.

    Args:
        expression: Expression matrix (n_spots x n_genes), dense or sparse.
        coordinates: Spatial coordinates (n_spots x 2).
        n_clusters: Number of clusters (required for KMeans, optional for graph methods).
        method: Clustering algorithm ("leiden", "louvain", "kmeans").
        n_neighbors: Number of spatial neighbors.
        graph_method: How to build the spatial graph.
        resolution: Resolution parameter for graph methods.
        spatial_weight: Weight for spatial graph in [0, 1]. 0 = expression only, 1 = spatial only.
        n_pcs: Number of PCA components for expression.
        seed: Random seed.

    Returns:
        SpatialClusterResult with labels and metadata.
    """
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")
    if sp_sparse is None:
        raise ImportError("scipy is required: uv pip install scipy")

    # Convert expression to dense if sparse
    if sp_sparse.issparse(expression):
        expr_dense = expression.toarray()
    else:
        expr_dense = np.asarray(expression, dtype=np.float64)

    coords = np.asarray(coordinates, dtype=np.float64)
    n = expr_dense.shape[0]

    if method == "kmeans":
        if KMeans is None or PCA is None or StandardScaler is None:
            raise ImportError("scikit-learn is required for KMeans: uv pip install scikit-learn")

        if n_clusters is None:
            raise ValueError("n_clusters must be specified for KMeans")

        # PCA on expression
        n_components = min(n_pcs, expr_dense.shape[1], n - 1)
        if n_components < 1:
            n_components = 1
        pca = PCA(n_components=n_components, random_state=seed)
        expr_pca = pca.fit_transform(StandardScaler().fit_transform(expr_dense))

        # Scale coordinates
        coords_scaled = StandardScaler().fit_transform(coords)

        # Concatenate with weighting
        features = np.hstack([
            expr_pca * (1.0 - spatial_weight),
            coords_scaled * spatial_weight,
        ])

        km = KMeans(n_clusters=n_clusters, random_state=seed, n_init=10)
        labels = km.fit_predict(features)

        return SpatialClusterResult(
            labels=labels.astype(np.int32),
            n_clusters=n_clusters,
            method="kmeans",
            metadata={
                "spatial_weight": spatial_weight,
                "n_pcs": n_components,
                "inertia": float(km.inertia_),
            },
        )

    # Graph-based methods (leiden, louvain)
    # 1. Build expression graph using KNN in PCA space
    if PCA is not None and StandardScaler is not None and NearestNeighbors is not None:
        n_components = min(n_pcs, expr_dense.shape[1], n - 1)
        if n_components < 1:
            n_components = 1
        pca = PCA(n_components=n_components, random_state=seed)
        expr_pca = pca.fit_transform(StandardScaler().fit_transform(expr_dense))

        nn = NearestNeighbors(n_neighbors=min(n_neighbors + 1, n), metric="euclidean")
        nn.fit(expr_pca)
        distances, indices = nn.kneighbors(expr_pca)

        rows_e: list[int] = []
        cols_e: list[int] = []
        vals_e: list[float] = []
        for i in range(n):
            for j_idx in range(1, indices.shape[1]):
                j = indices[i, j_idx]
                d = distances[i, j_idx]
                if d > 0:
                    w = 1.0 / (1.0 + d)
                    rows_e.append(i)
                    cols_e.append(j)
                    vals_e.append(w)

        expr_graph = sp_sparse.csr_matrix((vals_e, (rows_e, cols_e)), shape=(n, n))
        expr_graph = expr_graph + expr_graph.T
        # Normalize
        max_val = expr_graph.max()
        if max_val > 0:
            expr_graph = expr_graph / max_val
    else:
        raise ImportError("scikit-learn is required for expression PCA: uv pip install scikit-learn")

    # 2. Build spatial graph
    spatial_graph = build_spatial_graph(coords, method=graph_method, n_neighbors=n_neighbors)
    # Normalize
    max_val = spatial_graph.max()
    if max_val > 0:
        spatial_graph = spatial_graph / max_val

    # 3. Combine
    combined = (1.0 - spatial_weight) * expr_graph + spatial_weight * spatial_graph

    # 4. Cluster
    if method == "leiden":
        labels, modularity = leiden_clustering(combined, resolution=resolution, seed=seed)
    elif method == "louvain":
        labels, modularity = louvain_clustering(combined, resolution=resolution, seed=seed)
    else:
        raise ValueError(f"Unknown method: {method}")

    n_found = len(np.unique(labels))

    return SpatialClusterResult(
        labels=labels,
        n_clusters=n_found,
        method=method,
        modularity=modularity,
        spatial_graph=combined,
        metadata={
            "spatial_weight": spatial_weight,
            "n_pcs": n_pcs,
            "resolution": resolution,
            "graph_method": graph_method,
        },
    )


def spatial_domains(
    expression: Any,
    coordinates: Any,
    n_domains: int | None = None,
    *,
    n_pcs: int = 15,
    n_neighbors: int = 6,
    max_iterations: int = 100,
    seed: int = 42,
) -> SpatialClusterResult:
    """Identify spatial domains using a BayesSpace-inspired iterative approach.

    Algorithm:
    1. PCA on expression data.
    2. Initialize domains via KMeans on PCA + spatial features.
    3. Iteratively refine assignments by incorporating spatial neighbors:
       a. For each spot, compute mean expression profile of its spatial neighbors.
       b. Blend spot's own profile with its neighborhood mean.
       c. Re-cluster the blended profiles.
    4. Repeat until convergence or max_iterations.

    This is a simplified version of the BayesSpace spatial smoothing approach,
    using iterative neighborhood averaging instead of full MCMC.

    Args:
        expression: Expression matrix (n x genes), dense or sparse.
        coordinates: Spatial coordinates (n x 2).
        n_domains: Number of spatial domains. If None, auto-selects via gap statistic heuristic.
        n_pcs: Number of PCA components.
        n_neighbors: Spatial neighbors for smoothing.
        max_iterations: Maximum refinement iterations.
        seed: Random seed.

    Returns:
        SpatialClusterResult with domain labels.
    """
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")
    if sp_sparse is None:
        raise ImportError("scipy is required: uv pip install scipy")
    if KMeans is None or PCA is None or StandardScaler is None:
        raise ImportError("scikit-learn is required: uv pip install scikit-learn")

    # Convert to dense
    if sp_sparse.issparse(expression):
        expr_dense = expression.toarray()
    else:
        expr_dense = np.asarray(expression, dtype=np.float64)

    coords = np.asarray(coordinates, dtype=np.float64)
    n = expr_dense.shape[0]

    # PCA
    n_components = min(n_pcs, expr_dense.shape[1], n - 1)
    if n_components < 1:
        n_components = 1
    pca = PCA(n_components=n_components, random_state=seed)
    expr_pca = pca.fit_transform(StandardScaler().fit_transform(expr_dense))

    # Auto-select n_domains if not specified using elbow heuristic
    if n_domains is None:
        n_domains = _estimate_n_clusters(expr_pca, coords, max_k=15, seed=seed)
        logger.info(f"Auto-selected {n_domains} spatial domains")

    # Build spatial neighbor adjacency (binary, row-normalized)
    spatial_adj = build_spatial_graph(coords, method="knn", n_neighbors=n_neighbors)
    # Row-normalize for averaging
    row_sums = np.array(spatial_adj.sum(axis=1)).flatten()
    row_sums[row_sums == 0] = 1.0
    diag_inv = sp_sparse.diags(1.0 / row_sums)
    spatial_adj_norm = diag_inv @ spatial_adj

    # Initial clustering with spatial features
    coords_scaled = StandardScaler().fit_transform(coords)
    init_features = np.hstack([expr_pca, coords_scaled * 0.3])
    km = KMeans(n_clusters=n_domains, random_state=seed, n_init=10)
    labels = km.fit_predict(init_features)

    # Iterative spatial smoothing refinement
    smoothing_factor = 0.3  # blend ratio for neighbor averaging
    prev_labels = labels.copy()

    for iteration in range(max_iterations):
        # Smooth expression PCA by neighborhood
        neighbor_mean = spatial_adj_norm @ expr_pca
        smoothed = (1.0 - smoothing_factor) * expr_pca + smoothing_factor * neighbor_mean

        # Add spatial coordinates
        features = np.hstack([smoothed, coords_scaled * 0.3])

        # Re-cluster
        km = KMeans(n_clusters=n_domains, random_state=seed, n_init=5, init="k-means++")
        labels = km.fit_predict(features)

        # Check convergence
        changed = np.sum(labels != prev_labels)
        if changed == 0:
            logger.info(f"Spatial domains converged after {iteration + 1} iterations")
            break
        prev_labels = labels.copy()

    labels = labels.astype(np.int32)
    n_found = len(np.unique(labels))

    return SpatialClusterResult(
        labels=labels,
        n_clusters=n_found,
        method="spatial_domains",
        metadata={
            "n_pcs": n_components,
            "n_neighbors": n_neighbors,
            "iterations": iteration + 1 if 'iteration' in dir() else max_iterations,
            "smoothing_factor": smoothing_factor,
        },
    )


def _estimate_n_clusters(
    expr_pca: Any,
    coords: Any,
    max_k: int = 15,
    seed: int = 42,
) -> int:
    """Estimate optimal number of clusters using the elbow method with spatial features.

    Uses the second derivative of the within-cluster sum of squares (WCSS)
    curve to find the elbow point.

    Args:
        expr_pca: PCA-transformed expression (n x pcs).
        coords: Spatial coordinates (n x 2).
        max_k: Maximum number of clusters to try.
        seed: Random seed.

    Returns:
        Estimated optimal number of clusters (minimum 2).
    """
    coords_scaled = StandardScaler().fit_transform(coords)
    features = np.hstack([expr_pca, coords_scaled * 0.3])

    max_k = min(max_k, features.shape[0] - 1)
    if max_k < 2:
        return 2

    inertias: list[float] = []
    k_range = range(2, max_k + 1)
    for k in k_range:
        km = KMeans(n_clusters=k, random_state=seed, n_init=5)
        km.fit(features)
        inertias.append(float(km.inertia_))

    if len(inertias) < 3:
        return 2

    # Find elbow using second derivative
    inertias_arr = np.array(inertias)
    first_deriv = np.diff(inertias_arr)
    second_deriv = np.diff(first_deriv)

    # Elbow is where second derivative is maximum (most rapid change in slope)
    elbow_idx = int(np.argmax(second_deriv)) + 2  # +2 because we started at k=2 and took two diffs
    return max(2, min(elbow_idx, max_k))
