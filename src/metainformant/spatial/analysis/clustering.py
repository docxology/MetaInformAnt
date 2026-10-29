"""Spatially-aware clustering for spatial transcriptomics.

Implements spatial graph construction (KNN, Delaunay, radius) and community
detection algorithms (Leiden, Louvain) that incorporate spatial proximity
into transcriptomic clustering. Also provides a BayesSpace-inspired spatial
domain identification approach.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Literal, NamedTuple, cast

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Optional dependencies. Type checkers model them as installed; at runtime a missing
# dependency leaves the name bound to None and the functions that need it raise ImportError.
try:
    import numpy as np
except ImportError:
    if not TYPE_CHECKING:
        np = None

try:
    from scipy import sparse as sp_sparse
    from scipy.spatial import Delaunay, KDTree
except ImportError:
    if not TYPE_CHECKING:
        sp_sparse = None
        Delaunay = None
        KDTree = None

try:
    from sklearn.cluster import KMeans
    from sklearn.decomposition import PCA
    from sklearn.neighbors import NearestNeighbors
    from sklearn.preprocessing import StandardScaler
except ImportError:
    KMeans = None
    PCA = None
    NearestNeighbors = None
    StandardScaler = None


@dataclass
class SpatialClusterResult:
    """Result of spatial clustering.

    Attributes:
        labels: Cluster label per spot/cell (integer array of length n).
        n_clusters: Number of clusters found.
        method: Clustering method that actually ran: ``"leiden"``, ``"louvain"``, ``"kmeans"``,
            ``"spatial_domains"``, or ``"greedy_modularity"`` when a graph method was requested and its package was
            missing and the greedy fallback was explicitly allowed.
        modularity: Modularity score (for graph-based methods).
        spatial_graph: Adjacency matrix used for clustering.
        metadata: Additional result metadata. For the graph methods it holds ``requested_method`` (what the caller
            asked for) and ``fallback_used`` (whether ``method`` differs from it).
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
        Sparse adjacency matrix (scipy CSR) of shape (n, n). The "knn" graph is symmetric and weights each
        edge by ``1 / (1 + distance)`` whether the two spots are mutual neighbors or only one lists the
        other; "delaunay" and "radius" edges have weight 1.

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
        # Make symmetric. An edge's weight depends only on the distance between its two endpoints, so the
        # entrywise maximum keeps one weight per pair; adding the transpose would give mutual neighbors
        # twice the weight of a pair at the same distance that only one side lists.
        adj = sp_sparse.csr_matrix(adj.maximum(adj.T))

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


def _leidenalg_iterations(n_iterations: int) -> int:
    """Translate the public ``n_iterations`` into the value leidenalg expects.

    leidenalg runs until an iteration brings no improvement when given a negative number, which is how
    the public API documents ``-1``. Zero iterations would hand back the initial all-singletons partition,
    so it is rejected rather than quietly changed into some other count.
    """
    if n_iterations == 0:
        raise ValueError("n_iterations must be positive, or negative to iterate until convergence")
    return -1 if n_iterations < 0 else int(n_iterations)


class _GraphClustering(NamedTuple):
    """Outcome of a graph clustering: the labels, their modularity and the method that actually produced them."""

    labels: Any
    modularity: float
    method: str


def _greedy_fallback(
    adj: Any, *, requested: str, missing: str, resolution: float, seed: int, allow_greedy_fallback: bool
) -> _GraphClustering:
    """Run the greedy modularity optimization in place of an unavailable ``requested`` method, if that is allowed.

    Raises:
        ImportError: Unless ``allow_greedy_fallback`` is set. The greedy optimization is a different, simpler algorithm
            than Leiden or Louvain and finds different partitions, so it is never substituted silently.
    """
    if not allow_greedy_fallback:
        raise ImportError(
            f"{requested.capitalize()} clustering requires {missing}, which could not be imported. Install them, or "
            "pass allow_greedy_fallback=True to cluster with the in-package greedy modularity optimization instead "
            "(a different, simpler algorithm)."
        )
    message = (
        f"{requested.capitalize()} clustering is unavailable ({missing} could not be imported): clustered with the "
        "greedy modularity optimization instead, which finds different partitions"
    )
    logger.warning(message)
    warnings.warn(message, UserWarning, stacklevel=4)
    labels, modularity = _greedy_modularity_clustering(adj, resolution=resolution, seed=seed)
    return _GraphClustering(labels, modularity, "greedy_modularity")


def _leiden(
    adj: Any, *, resolution: float, iterations: int, seed: int, allow_greedy_fallback: bool
) -> _GraphClustering:
    """Leiden partition of the graph with ``adj`` as its weighted adjacency, or the explicitly allowed fallback."""
    try:
        import igraph as ig
        import leidenalg
    except (ImportError, SyntaxError):
        # Some leidenalg releases fail during import because of an invalid escape in an upstream docstring
        return _greedy_fallback(
            adj,
            requested="leiden",
            missing="python-igraph and leidenalg",
            resolution=resolution,
            seed=seed,
            allow_greedy_fallback=allow_greedy_fallback,
        )

    # Convert scipy sparse to igraph
    sources, targets = adj.nonzero()
    weights = np.array(adj[sources, targets]).flatten()
    mask = sources < targets  # upper triangle only for undirected
    edges = list(zip(sources[mask].tolist(), targets[mask].tolist()))
    edge_weights = weights[mask].tolist()

    g = ig.Graph(n=adj.shape[0], edges=edges, directed=False)
    g.es["weight"] = edge_weights

    partition = leidenalg.find_partition(
        g,
        leidenalg.RBConfigurationVertexPartition,
        weights="weight",
        resolution_parameter=resolution,
        n_iterations=iterations,
        seed=seed,
    )

    labels = np.array(partition.membership, dtype=np.int32)
    modularity = float(partition.modularity)
    logger.info(f"Leiden clustering: {len(set(labels))} clusters, modularity={modularity:.4f}")
    return _GraphClustering(labels, modularity, "leiden")


def _louvain(adj: Any, *, resolution: float, seed: int, allow_greedy_fallback: bool) -> _GraphClustering:
    """Louvain partition of the graph with ``adj`` as its weighted adjacency, or the explicitly allowed fallback."""
    try:
        import community as community_louvain
        import networkx as nx
    except ImportError:
        return _greedy_fallback(
            adj,
            requested="louvain",
            missing="python-louvain and networkx",
            resolution=resolution,
            seed=seed,
            allow_greedy_fallback=allow_greedy_fallback,
        )

    n = adj.shape[0]

    # Convert to networkx graph
    G = nx.Graph()
    G.add_nodes_from(range(n))
    sources, targets = adj.nonzero()
    weights = np.array(adj[sources, targets]).flatten()
    for s, t, w in zip(sources.tolist(), targets.tolist(), weights.tolist()):
        if s < t:
            G.add_edge(s, t, weight=w)

    partition = community_louvain.best_partition(G, resolution=resolution, random_state=seed)
    labels = np.array([partition[i] for i in range(n)], dtype=np.int32)
    modularity = float(community_louvain.modularity(partition, G))
    logger.info(f"Louvain clustering: {len(set(labels))} clusters, modularity={modularity:.4f}")
    return _GraphClustering(labels, modularity, "louvain")


def leiden_clustering(
    adjacency_matrix: Any,
    resolution: float = 1.0,
    n_iterations: int = -1,
    seed: int = 42,
    *,
    allow_greedy_fallback: bool = False,
) -> tuple[Any, float]:
    """Leiden community detection on a graph.

    Runs the Leiden algorithm of leidenalg on the weighted graph, which needs python-igraph and leidenalg. A
    different algorithm is never substituted silently: without those packages the function raises
    ``ImportError``, unless ``allow_greedy_fallback=True``.

    Args:
        adjacency_matrix: Sparse adjacency matrix (n x n).
        resolution: Resolution parameter (higher = more clusters).
        n_iterations: Number of leidenalg iterations; any negative value (the default, -1) iterates until
            an iteration brings no improvement. Zero is rejected. The greedy fallback ignores it: it
            sweeps the nodes until none moves, at most 50 times.
        seed: Random seed for reproducibility.
        allow_greedy_fallback: When python-igraph or leidenalg is unavailable, cluster with the in-package greedy
            modularity optimization instead of raising, and emit a ``UserWarning`` and a log warning that say so.
            The greedy optimization is a simpler algorithm and finds different partitions; the returned tuple does
            not say which ran, so rely on that warning. :func:`spatial_cluster` records the method used in
            ``SpatialClusterResult.method``. The flag does nothing when Leiden is available.

    Returns:
        Tuple of (labels, modularity) where labels is an integer array
        and modularity is the partition quality score.

    Raises:
        ImportError: If python-igraph or leidenalg is unavailable and ``allow_greedy_fallback`` is not set.
        ValueError: If ``n_iterations`` is zero.
    """
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")
    if sp_sparse is None:
        raise ImportError("scipy is required: uv pip install scipy")

    iterations = _leidenalg_iterations(n_iterations)
    adj = sp_sparse.csr_matrix(adjacency_matrix)
    outcome = _leiden(
        adj, resolution=resolution, iterations=iterations, seed=seed, allow_greedy_fallback=allow_greedy_fallback
    )
    return outcome.labels, outcome.modularity


def louvain_clustering(
    adjacency_matrix: Any,
    resolution: float = 1.0,
    seed: int = 42,
    *,
    allow_greedy_fallback: bool = False,
) -> tuple[Any, float]:
    """Louvain community detection on a graph.

    Runs the Louvain method of python-louvain on the weighted graph, which needs python-louvain and networkx. A
    different algorithm is never substituted silently: without those packages the function raises
    ``ImportError``, unless ``allow_greedy_fallback=True``.

    Args:
        adjacency_matrix: Sparse adjacency matrix (n x n).
        resolution: Resolution parameter (higher = more clusters).
        seed: Random seed.
        allow_greedy_fallback: When python-louvain or networkx is unavailable, cluster with the in-package greedy
            modularity optimization instead of raising, and emit a ``UserWarning`` and a log warning that say so.
            The greedy optimization is a simpler algorithm and finds different partitions; the returned tuple does
            not say which ran, so rely on that warning. :func:`spatial_cluster` records the method used in
            ``SpatialClusterResult.method``. The flag does nothing when Louvain is available.

    Returns:
        Tuple of (labels, modularity).

    Raises:
        ImportError: If python-louvain or networkx is unavailable and ``allow_greedy_fallback`` is not set.
    """
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")
    if sp_sparse is None:
        raise ImportError("scipy is required: uv pip install scipy")

    adj = sp_sparse.csr_matrix(adjacency_matrix)
    outcome = _louvain(adj, resolution=resolution, seed=seed, allow_greedy_fallback=allow_greedy_fallback)
    return outcome.labels, outcome.modularity


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
                delta_q = (w_ic - w_current) / m - resolution * ki * (sigma_tot - cur_sigma) / (2.0 * m * m)

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
    labels = np.array([label_map[label] for label in labels], dtype=np.int32)

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
    radius: float | None = None,
    resolution: float = 1.0,
    spatial_weight: float = 0.5,
    n_pcs: int = 15,
    seed: int = 42,
    allow_greedy_fallback: bool = False,
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
        n_neighbors: Number of spatial neighbors (used by the "knn" graph; also the expression-graph size).
        graph_method: How to build the spatial graph.
        radius: Neighborhood radius, in the units of ``coordinates``, for ``graph_method="radius"``. It is
            required with that graph method (there is no data-derived default) and ignored otherwise.
        resolution: Resolution parameter for graph methods.
        spatial_weight: Weight for spatial graph in [0, 1]. 0 = expression only, 1 = spatial only.
        n_pcs: Number of PCA components for expression.
        seed: Random seed.
        allow_greedy_fallback: Graph methods only. When the package behind "leiden" (python-igraph and leidenalg)
            or "louvain" (python-louvain and networkx) is unavailable, cluster with the in-package greedy modularity
            optimization instead of raising, with a ``UserWarning``. The result then says so: ``method`` is
            ``"greedy_modularity"``, ``metadata["requested_method"]`` the method asked for and
            ``metadata["fallback_used"]`` is True. It does nothing when the package is available.

    Returns:
        SpatialClusterResult with labels and metadata. ``method`` is the method that actually ran.

    Raises:
        ValueError: If a graph-based method is used with ``graph_method="radius"`` and no ``radius``, or
            if KMeans is requested without ``n_clusters``.
        ImportError: If the package behind a graph method is unavailable and ``allow_greedy_fallback`` is not set.
    """
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")
    if sp_sparse is None:
        raise ImportError("scipy is required: uv pip install scipy")
    if method != "kmeans" and graph_method == "radius" and radius is None:
        raise ValueError("graph_method='radius' needs the radius argument (a distance in the units of coordinates)")

    # Convert expression to dense if sparse
    if sp_sparse.issparse(expression):
        # scipy-stubs narrow issparse() to a method-less SparseABC; every sparse format provides toarray()
        expr_dense = cast("sp_sparse.csr_matrix", expression).toarray()
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
        features = np.hstack(
            [
                expr_pca * (1.0 - spatial_weight),
                coords_scaled * spatial_weight,
            ]
        )

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
        # Symmetrize by the entrywise maximum (adding the transpose would double mutual-neighbor edges)
        expr_graph = sp_sparse.csr_matrix(expr_graph.maximum(expr_graph.T))
        # Normalize
        max_val = expr_graph.max()
        if max_val > 0:
            expr_graph = expr_graph / max_val
    else:
        raise ImportError("scikit-learn is required for expression PCA: uv pip install scikit-learn")

    # 2. Build spatial graph
    spatial_graph = build_spatial_graph(coords, method=graph_method, n_neighbors=n_neighbors, radius=radius)
    # Normalize
    max_val = spatial_graph.max()
    if max_val > 0:
        spatial_graph = spatial_graph / max_val

    # 3. Combine
    combined = (1.0 - spatial_weight) * expr_graph + spatial_weight * spatial_graph

    # 4. Cluster
    combined = sp_sparse.csr_matrix(combined)
    if method == "leiden":
        outcome = _leiden(
            combined,
            resolution=resolution,
            iterations=_leidenalg_iterations(-1),
            seed=seed,
            allow_greedy_fallback=allow_greedy_fallback,
        )
    elif method == "louvain":
        outcome = _louvain(combined, resolution=resolution, seed=seed, allow_greedy_fallback=allow_greedy_fallback)
    else:
        raise ValueError(f"Unknown method: {method}")

    labels = outcome.labels
    n_found = len(np.unique(labels))

    return SpatialClusterResult(
        labels=labels,
        n_clusters=n_found,
        method=outcome.method,
        modularity=outcome.modularity,
        spatial_graph=combined,
        metadata={
            "spatial_weight": spatial_weight,
            "n_pcs": n_pcs,
            "resolution": resolution,
            "graph_method": graph_method,
            "radius": radius,
            "requested_method": method,
            "fallback_used": outcome.method != method,
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
        n_domains: Number of spatial domains. If None, it is chosen with the gap statistic (Tibshirani,
            Walther and Hastie, 2001) on the same PCA-plus-scaled-coordinates features that seed the
            clustering: between 2 and 15 domains, the smallest k whose gap is within one standard error
            of the next one's. The reference datasets are drawn with ``seed``. With only tens of spots per
            domain the statistic can prefer a few extra domains, so pass the number when it is known. How
            the number was chosen is recorded in ``metadata["n_domains_selection"]``.
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
        expr_dense = cast("sp_sparse.csr_matrix", expression).toarray()
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

    # Auto-select n_domains if not specified, using the gap statistic
    n_domains_selection = "user"
    if n_domains is None:
        n_domains = _estimate_n_clusters(expr_pca, coords, max_k=15, seed=seed)
        n_domains_selection = "gap_statistic"
        logger.info(f"Auto-selected {n_domains} spatial domains with the gap statistic")

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
    iterations_run = max_iterations

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
            iterations_run = iteration + 1
            break
        prev_labels = labels.copy()

    labels = labels.astype(np.int32)
    n_found = len(np.unique(labels))

    return SpatialClusterResult(
        labels=labels,
        n_clusters=n_found,
        method="spatial_domains",
        metadata={
            "iterations": iterations_run,
            "n_neighbors": n_neighbors,
            "smoothing_factor": smoothing_factor,
            "n_domains_selection": n_domains_selection,
        },
    )


def _log_dispersion(data: Any, k: int, seed: int) -> float:
    """log W_k: log of the within-cluster sum of squares of a k-means fit with ``k`` clusters."""
    inertia = float(KMeans(n_clusters=k, random_state=seed, n_init=5).fit(data).inertia_)
    return float(np.log(max(inertia, float(np.finfo(np.float64).tiny))))


def _gap_and_spread(observed: Any, references: list[Any], k: int, seed: int) -> tuple[float, float]:
    """Gap(k) and s_k of Tibshirani, Walther and Hastie (2001) for a given set of reference datasets.

    ``Gap(k) = mean_b(log W*_kb) - log W_k`` and ``s_k = sd_b(log W*_kb) * sqrt(1 + 1/B)``, where
    ``W_k`` is the within-cluster sum of squares of ``k`` clusters on the data, ``W*_kb`` the same on the
    b-th of ``B`` reference datasets, and sd is the (population) standard deviation over the references.
    """
    reference_logs = np.array([_log_dispersion(reference, k, seed) for reference in references])
    gap = float(reference_logs.mean() - _log_dispersion(observed, k, seed))
    spread = float(reference_logs.std() * np.sqrt(1.0 + 1.0 / len(references)))
    return gap, spread


def _first_k_by_gap_rule(gaps: list[float], spreads: list[float]) -> int | None:
    """Smallest k with ``Gap(k) >= Gap(k + 1) - s_(k+1)`` (Tibshirani et al., 2001), or None if none yet.

    ``gaps[i]`` and ``spreads[i]`` belong to ``k = i + 1``; only k up to ``len(gaps) - 1`` can be judged,
    because the rule compares k with k + 1.
    """
    for k in range(1, len(gaps)):
        if gaps[k - 1] >= gaps[k] - spreads[k]:
            return k
    return None


def _estimate_n_clusters(
    expr_pca: Any,
    coords: Any,
    max_k: int = 15,
    seed: int = 42,
    n_references: int = 10,
) -> int:
    """Choose the number of clusters with the gap statistic (Tibshirani, Walther and Hastie, 2001).

    The features are the PCA-transformed expression plus the scaled coordinates (the ones that seed the
    domain clustering). For k = 1, 2, ... the log within-cluster sum of squares of a k-means fit is
    compared with its expectation under a null of no cluster structure, estimated on ``n_references``
    reference datasets drawn uniformly over a box aligned with the principal axes of the features (their
    rotation leaves k-means unchanged, so the rotated space is clustered directly). The answer is the
    smallest k with ``Gap(k) >= Gap(k + 1) - s_(k+1)`` (R's ``clusGap`` "firstSEmax"); evaluation stops as
    soon as that holds. Like the gap statistic in general it can prefer a few extra clusters when there are
    only tens of spots per domain; pass ``n_domains`` to ``spatial_domains`` when the number is known.

    Args:
        expr_pca: PCA-transformed expression (n x pcs).
        coords: Spatial coordinates (n x 2).
        max_k: Maximum number of clusters to try.
        seed: Random seed for the reference datasets and the k-means fits.
        n_references: Number of reference datasets B.

    Returns:
        Estimated number of clusters, at least 2 (the rule may find none, k = 1, which is no domain
        structure; a domain analysis still needs two) and at most ``max_k``.
    """
    coords_scaled = StandardScaler().fit_transform(coords)
    features = np.hstack([expr_pca, coords_scaled * 0.3])

    max_k = min(max_k, features.shape[0] - 1)
    if max_k < 2:
        return 2

    # Principal-axis frame: rotation and translation do not change k-means inertia
    centred = features - features.mean(axis=0)
    _, _, axes = np.linalg.svd(centred, full_matrices=False)
    rotated = centred @ axes.T
    low, high = rotated.min(axis=0), rotated.max(axis=0)
    rng = np.random.RandomState(seed)
    references = [rng.uniform(low, high, size=rotated.shape) for _ in range(n_references)]

    gaps: list[float] = []
    spreads: list[float] = []
    for k in range(1, max_k + 1):
        gap, spread = _gap_and_spread(rotated, references, k, seed)
        gaps.append(gap)
        spreads.append(spread)
        chosen = _first_k_by_gap_rule(gaps, spreads)
        if chosen is not None:
            return max(2, chosen)
    return max_k
