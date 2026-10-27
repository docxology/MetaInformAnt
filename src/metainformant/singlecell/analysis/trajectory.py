"""Trajectory inference methods for single-cell data.

This module provides functions for inferring developmental trajectories and
pseudotemporal ordering from single-cell expression data.

Diffusion pseudotime (:func:`compute_diffusion_pseudotime`, :func:`dpt_trajectory`) implements the DPT
algorithm of Haghverdi et al. (2016). The cluster-level functions (PAGA-style connectivity,
Slingshot-style paths) and the embedding-based pseudotime are simplified approximations; see their
docstrings.
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import eigsh

from metainformant.core.data import validation
from metainformant.core.utils import errors, logging
from metainformant.singlecell.analysis.nonlinear_methods import compute_neighbors
from metainformant.singlecell.data.preprocessing import SingleCellData

# Optional scientific dependencies
try:
    from sklearn.neighbors import NearestNeighbors

    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False
    NearestNeighbors = None

logger = logging.get_logger(__name__)

# Eigenproblems up to this many cells use the dense symmetric solver (exact, and fast at this size);
# larger ones ask ARPACK for just the leading eigenpairs.
_DENSE_EIGEN_LIMIT = 1000


def _default_root_cell(X: np.ndarray) -> int:
    """Heuristic root: the cell with the largest gene-variance-weighted expression sum.

    This is only a convenience default for callers who give no ``root_cell``. DPT measures distance
    *from* the root, so the direction of the pseudotime is whatever this guess (or the caller) says.
    """
    variances = np.var(X, axis=0)
    return int(np.argmax(np.sum(X * variances, axis=1)))


def _dpt_transition_operator(affinities: sparse.csr_matrix) -> sparse.csr_matrix:
    """Symmetrised, density-normalised random-walk operator ``Z K Z`` of a connected affinity graph.

    ``K = D^-1 W D^-1`` (Coifman & Lafon, alpha = 1) removes the dependence on how densely the cells were
    sampled, and ``Z`` holds the inverse square roots of the row sums of ``K``. The operator has the same
    eigenvalues as the transition matrix ``T`` of the random walk on ``K``; its eigenvectors are
    orthonormal.
    """
    degree = np.asarray(affinities.sum(axis=1)).ravel()
    inverse_degree = sparse.diags(1.0 / degree)
    density_normalised = inverse_degree @ affinities @ inverse_degree
    inverse_root = sparse.diags(1.0 / np.sqrt(np.asarray(density_normalised.sum(axis=1)).ravel()))
    return sparse.csr_matrix(inverse_root @ density_normalised @ inverse_root)


def _leading_eigenpairs(operator: sparse.csr_matrix, k: int) -> tuple[np.ndarray, np.ndarray]:
    """The ``k`` algebraically largest eigenpairs of a symmetric matrix, largest eigenvalue first."""
    n = operator.shape[0]
    if n <= _DENSE_EIGEN_LIMIT or k >= n - 1:
        values, vectors = np.linalg.eigh(operator.toarray())
    else:
        # Fixed start vector: ARPACK's own default is drawn from internal state that changes between calls
        values, vectors = eigsh(operator, k=k, which="LA", v0=np.random.default_rng(0).standard_normal(n))
    order = np.argsort(values)[::-1][:k]
    return values[order], vectors[:, order]


def _dpt_distances(operator: sparse.csr_matrix, root: int, n_components: int) -> tuple[np.ndarray, np.ndarray]:
    """DPT distances from ``root`` to every cell, and the eigenvalues they were built from.

    With eigenvalues ``lambda_i`` and orthonormal eigenvectors ``psi_i`` of the transition operator, the
    distance between cells x and y is ``sqrt(sum_i (lambda_i / (1 - lambda_i))**2 * (psi_i(x) - psi_i(y))**2)``
    over the leading non-trivial eigenvectors: the Euclidean distance between the rows of the accumulated
    transition matrix ``M = sum_{t>=1} (T**t - steady state)`` (Haghverdi et al., 2016). The trivial
    eigenvector (eigenvalue 1, the steady state) is excluded.
    """
    n = operator.shape[0]
    values, vectors = _leading_eigenpairs(operator, min(n_components + 1, n))
    values, vectors = values[1:], vectors[:, 1:]
    # lambda / (1 - lambda); an eigenvalue indistinguishable from 1 would otherwise divide by zero
    weights = values / (1.0 - np.minimum(values, 1.0 - 1e-12))
    differences = (vectors - vectors[root]) * weights
    return np.sqrt(np.sum(differences**2, axis=1)), values


def compute_diffusion_pseudotime(
    data: SingleCellData,
    root_cell: int | None = None,
    n_components: int = 10,
    *,
    n_neighbors: int = 15,
) -> SingleCellData:
    """Compute diffusion pseudotime (DPT) for trajectory inference.

    Implements Haghverdi et al. (2016, Nature Methods 13:845). The cells are connected by a
    k-nearest-neighbour graph with Gaussian affinities: the graph in ``data.uns["neighbors"]`` when there
    is one (see :func:`~metainformant.singlecell.analysis.nonlinear_methods.compute_neighbors`),
    otherwise one built from ``data.X`` with ``n_neighbors`` neighbours and stored in the returned
    object. The affinities are density-normalised into a random-walk transition matrix ``T``. The DPT
    distance between two cells is the Euclidean distance between their rows of the accumulated transition
    matrix ``sum_{t>=1} (T**t - steady state)``, evaluated on the leading ``n_components`` non-trivial
    eigenvectors of ``T``. The pseudotime of a cell is its DPT distance from the root cell, divided by
    the largest such distance, so it spans [0, 1] with the root at 0.

    Because the distance follows the graph, pseudotime keeps increasing along a curved trajectory, where
    the straight-line (Euclidean) distance from the root would turn back.

    Args:
        data: SingleCellData object with expression matrix
        root_cell: Index of root cell (starting point of trajectory). When omitted, a heuristic picks one
            (the cell with the largest gene-variance-weighted expression sum) and the pseudotime
            direction is only as meaningful as that guess; the chosen cell is recorded in
            ``uns["dpt"]["root_cell"]``.
        n_components: Number of non-trivial diffusion components (eigenvectors) in the DPT distance
        n_neighbors: Neighbours per cell for the graph that is built when ``data.uns["neighbors"]`` is
            absent

    Returns:
        SingleCellData with ``obs["dpt_pseudotime"]`` and ``uns["dpt"]``. Cells that the neighbour graph
        does not connect to the root cell have no finite distance to it and get ``inf``
        (``uns["dpt"]["n_unreachable"]`` counts them).

    Raises:
        ValidationError: If data is not SingleCellData; if root_cell, n_components or n_neighbors is out of range; or if
            ``uns["neighbors"]`` lacks ``connectivities``, does not match the number of cells (a stale graph: recompute
            it with ``compute_neighbors``) or holds negative affinities
    """
    validation.validate_type(data, SingleCellData, "data")
    validation.validate_range(n_components, min_val=1, name="n_components")
    validation.validate_range(n_neighbors, min_val=1, name="n_neighbors")

    if root_cell is not None:
        validation.validate_range(root_cell, min_val=0, max_val=data.n_obs - 1, name="root_cell")

    logger.info(f"Computing diffusion pseudotime with {n_components} components")

    # Create copy
    result = data.copy()
    n_cells = result.n_obs

    # Affinity graph between cells
    if "neighbors" in result.uns:
        neighbors_source = "uns"
        if "connectivities" not in result.uns["neighbors"]:
            raise errors.ValidationError(
                "uns['neighbors'] has no 'connectivities' matrix: recompute it with compute_neighbors"
            )
        stored_shape = result.uns["neighbors"]["connectivities"].shape
        if stored_shape != (n_cells, n_cells):
            raise errors.ValidationError(
                f"uns['neighbors'] describes {stored_shape[0]} cells but data has {n_cells}: "
                "recompute it with compute_neighbors"
            )
    else:
        neighbors_source = "X"
        result = compute_neighbors(result, n_neighbors=min(n_neighbors, n_cells), use_rep="X", inplace=True)
    # A private copy: eliminate_zeros() edits in place and must not touch the stored graph
    affinities = sparse.csr_matrix(result.uns["neighbors"]["connectivities"], dtype=float, copy=True)
    affinities.eliminate_zeros()
    if affinities.nnz and affinities.data.min() < 0:
        raise errors.ValidationError("uns['neighbors']['connectivities'] must not contain negative affinities")

    # Root cell
    if root_cell is None:
        X = result.X.toarray() if hasattr(result.X, "toarray") else result.X
        root_cell = _default_root_cell(X)
        root_selection = "heuristic"
    else:
        root_selection = "user"

    # Distances are only defined inside the root's connected component
    _, component_labels = connected_components(affinities, directed=False)
    component = np.flatnonzero(component_labels == component_labels[root_cell])
    pseudotime = np.full(n_cells, np.inf)
    eigenvalues = np.zeros(0)
    if component.size == 1:
        pseudotime[root_cell] = 0.0
    else:
        operator = _dpt_transition_operator(sparse.csr_matrix(affinities[component][:, component]))
        distances, eigenvalues = _dpt_distances(operator, int(np.searchsorted(component, root_cell)), n_components)
        farthest = float(distances.max())
        pseudotime[component] = distances / farthest if farthest > 0 else 0.0

    n_unreachable = int(n_cells - component.size)
    if n_unreachable:
        logger.warning(
            f"{n_unreachable} cells are not connected to root cell {root_cell} in the neighbor graph; "
            "their pseudotime is inf"
        )

    # Add to obs
    result.obs = result.obs.copy() if result.obs is not None else pd.DataFrame()
    result.obs["dpt_pseudotime"] = pseudotime

    # Store metadata
    reachable = pseudotime[np.isfinite(pseudotime)]
    result.uns["dpt"] = {
        "method": "diffusion_pseudotime",
        "root_cell": root_cell,
        "root_selection": root_selection,
        "n_components": int(eigenvalues.size),
        "eigenvalues": eigenvalues.tolist(),
        "n_neighbors": result.uns["neighbors"].get("params", {}).get("n_neighbors"),
        "neighbors_source": neighbors_source,
        "n_unreachable": n_unreachable,
        "min_pseudotime": float(np.min(reachable)),
        "max_pseudotime": float(np.max(reachable)),
    }

    logger.info(f"Diffusion pseudotime computed: range [{np.min(reachable):.3f}, {np.max(reachable):.3f}]")
    return result


def dpt_trajectory(data: SingleCellData, root_cell: int | None = None) -> SingleCellData:
    """Compute Diffusion Pseudotime (DPT) trajectory.

    Runs :func:`compute_diffusion_pseudotime` with its default of 10 diffusion components and 15
    neighbours; see there for the algorithm, the root-cell default and the ``inf`` convention.

    Args:
        data: SingleCellData object with expression matrix
        root_cell: Index of root cell (starting point of trajectory)

    Returns:
        SingleCellData with ``obs["dpt_pseudotime"]`` and ``uns["dpt"]``

    Raises:
        ValidationError: If data is not SingleCellData, root_cell is out of range, or the stored neighbour graph is
            unusable (see :func:`compute_diffusion_pseudotime`)
    """
    return compute_diffusion_pseudotime(data, root_cell, n_components=10)


def paga_trajectory(data: SingleCellData, groups: str) -> SingleCellData:
    """Compute Partition-based Graph Abstraction (PAGA) trajectory.

    Args:
        data: SingleCellData object with cluster assignments
        groups: Column name containing cluster/group labels

    Returns:
        SingleCellData with PAGA connectivity information

    Raises:
        ValidationError: If data is not SingleCellData or the groups column is not in data.obs
    """
    validation.validate_type(data, SingleCellData, "data")

    if data.obs is None or groups not in data.obs.columns:
        raise errors.ValidationError(f"Groups column '{groups}' not found in data.obs")

    logger.info(f"Computing PAGA trajectory using groups from '{groups}'")

    # Create copy
    result = data.copy()

    # Get cluster labels
    cluster_labels = data.obs[groups].values
    unique_clusters = np.unique(cluster_labels)

    # Compute cluster centroids
    X = data.X.toarray() if hasattr(data.X, "toarray") else data.X
    cluster_centroids = {}

    for cluster in unique_clusters:
        mask = cluster_labels == cluster
        if np.sum(mask) > 0:
            cluster_centroids[cluster] = np.mean(X[mask], axis=0)

    # Compute connectivity between clusters
    connectivity_matrix = _compute_cluster_connectivity(cluster_labels, unique_clusters)

    # Store PAGA results
    result.uns["paga"] = {
        "groups": groups,
        "unique_groups": unique_clusters.tolist(),
        "connectivity_matrix": connectivity_matrix.tolist(),
        "n_groups": len(unique_clusters),
    }

    # Add cluster-level pseudotime (simplified)
    # In practice, PAGA would compute a trajectory through the abstracted graph
    cluster_pseudotime = _compute_cluster_pseudotime(connectivity_matrix, unique_clusters)

    # Assign pseudotime to cells based on their cluster
    cell_pseudotime = np.array([cluster_pseudotime[cluster] for cluster in cluster_labels])

    result.obs = result.obs.copy() if result.obs is not None else pd.DataFrame()
    result.obs["paga_pseudotime"] = cell_pseudotime

    logger.info(f"PAGA trajectory computed for {len(unique_clusters)} groups")
    return result


def slingshot_trajectory(data: SingleCellData, start_cluster: str, end_clusters: List[str]) -> SingleCellData:
    """Compute Slingshot-like trajectory inference.

    Args:
        data: SingleCellData object with cluster assignments
        start_cluster: Starting cluster for trajectory
        end_clusters: List of ending clusters for trajectory branches

    Returns:
        SingleCellData with trajectory pseudotime

    Raises:
        ValidationError: If data is not SingleCellData, no cluster column is found in data.obs, or start_cluster or
            end_cluster is not among the clusters
    """
    validation.validate_type(data, SingleCellData, "data")

    # Assume clusters are in a column named 'cluster' or similar
    cluster_col = None
    if data.obs is not None:
        # Find a column that looks like cluster assignments
        for col in data.obs.columns:
            if "cluster" in col.lower() or "leiden" in col.lower() or "louvain" in col.lower():
                cluster_col = col
                break

    if cluster_col is None:
        raise errors.ValidationError("No cluster column found in data.obs")

    logger.info(f"Computing Slingshot trajectory from {start_cluster} to {end_clusters}")

    # Create copy
    result = data.copy()

    cluster_labels = data.obs[cluster_col].values

    # Compute cluster centroids and distances
    X = data.X.toarray() if hasattr(data.X, "toarray") else data.X
    unique_clusters = np.unique(cluster_labels)

    # Check if start and end clusters exist
    if start_cluster not in unique_clusters:
        raise errors.ValidationError(f"Start cluster '{start_cluster}' not found in data")

    for end_cluster in end_clusters:
        if end_cluster not in unique_clusters:
            raise errors.ValidationError(f"End cluster '{end_cluster}' not found in data")

    # Compute simplified trajectory
    trajectory_info = _compute_slingshot_trajectory(X, cluster_labels, start_cluster, end_clusters)

    # Assign pseudotime to cells
    cell_pseudotime = np.zeros(data.n_obs)
    for i, cluster in enumerate(cluster_labels):
        if cluster in trajectory_info["cluster_pseudotime"]:
            cell_pseudotime[i] = trajectory_info["cluster_pseudotime"][cluster]

    result.obs = result.obs.copy() if result.obs is not None else pd.DataFrame()
    result.obs["slingshot_pseudotime"] = cell_pseudotime

    # Store trajectory information
    result.uns["slingshot"] = {
        "start_cluster": start_cluster,
        "end_clusters": end_clusters,
        "trajectory_info": trajectory_info,
    }

    logger.info(f"Slingshot trajectory computed with {len(trajectory_info['branches'])} branches")
    return result


def compute_pseudotime_from_dimensionality_reduction(
    data: SingleCellData, dim_red_cols: List[str], root_cell: int | None = None
) -> SingleCellData:
    """Compute pseudotime from dimensionality reduction coordinates.

    The pseudotime is the straight-line (Euclidean) distance from the root cell in the embedding,
    scaled to [0, 1]. A shortest path over the complete graph of pairwise Euclidean distances is the
    direct edge (triangle inequality), so this is not a geodesic along the data manifold and it turns
    back on a curved trajectory; use :func:`compute_diffusion_pseudotime` for a graph-based pseudotime.

    Args:
        data: SingleCellData object with dimensionality reduction coordinates
        dim_red_cols: Column names containing the reduced dimensions
        root_cell: Index of root cell (default: the cell closest to the origin of the embedding)

    Returns:
        SingleCellData with pseudotime computed from embedding

    Raises:
        ValidationError: If data is not SingleCellData, has no obs, or lacks the dimension columns
    """
    validation.validate_type(data, SingleCellData, "data")

    if data.obs is None:
        raise errors.ValidationError("data.obs required for pseudotime computation")

    missing_cols = [col for col in dim_red_cols if col not in data.obs.columns]
    if missing_cols:
        raise errors.ValidationError(f"Dimension columns not found: {missing_cols}")

    logger.info(f"Computing pseudotime from {len(dim_red_cols)}D embedding")

    # Create copy
    result = data.copy()

    # Get embedding coordinates
    embedding = data.obs[dim_red_cols].values

    # Compute pseudotime as the Euclidean distance from the root cell
    if root_cell is None:
        # Choose cell with minimum norm as root
        norms = np.sum(embedding**2, axis=1)
        root_cell = int(np.argmin(norms))

    # Compute pairwise distances in embedding space
    from scipy.spatial.distance import pdist, squareform

    distances = squareform(pdist(embedding, metric="euclidean"))

    # Shortest paths over the complete distance graph: each is the direct edge, i.e. the Euclidean distance
    pseudotime = _compute_shortest_path_distances(distances, root_cell)

    # Normalize pseudotime to [0, 1]
    pseudotime = (pseudotime - np.min(pseudotime)) / (np.max(pseudotime) - np.min(pseudotime))

    result.obs = result.obs.copy() if result.obs is not None else pd.DataFrame()
    result.obs["embedding_pseudotime"] = pseudotime

    # Store metadata
    result.uns["embedding_pseudotime"] = {
        "dim_red_cols": dim_red_cols,
        "root_cell": root_cell,
        "method": "euclidean_distance",
        "min_pseudotime": float(np.min(pseudotime)),
        "max_pseudotime": float(np.max(pseudotime)),
    }

    logger.info(f"Embedding-based pseudotime computed: root cell = {root_cell}")
    return result


def find_trajectory_branches(data: SingleCellData, pseudotime_col: str, min_branch_size: int = 10) -> Dict[str, Any]:
    """Identify branches in a computed trajectory.

    Args:
        data: SingleCellData object with pseudotime
        pseudotime_col: Column name containing pseudotime values
        min_branch_size: Minimum number of cells per branch

    Returns:
        Dictionary with branch information

    Raises:
        ValidationError: If data is not SingleCellData or the pseudotime column is not in data.obs
    """
    validation.validate_type(data, SingleCellData, "data")

    if data.obs is None or pseudotime_col not in data.obs.columns:
        raise errors.ValidationError(f"Pseudotime column '{pseudotime_col}' not found in data.obs")

    logger.info("Identifying trajectory branches")

    pseudotime = data.obs[pseudotime_col].values

    # Simple branch detection based on local density changes
    # Sort cells by pseudotime
    sorted_indices = np.argsort(pseudotime)
    sorted_pseudotime = pseudotime[sorted_indices]

    # Detect branches using local minima in density
    branches = _detect_trajectory_branches(sorted_pseudotime, sorted_indices, min_branch_size)

    # Assign branch labels to cells
    branch_labels = np.full(data.n_obs, -1, dtype=int)
    for branch_id, cell_indices in branches.items():
        branch_labels[cell_indices] = branch_id

    result = {
        "n_branches": len(branches),
        "branch_sizes": [len(indices) for indices in branches.values()],
        "branch_labels": branch_labels,
        "branches": branches,
    }

    logger.info(f"Found {len(branches)} trajectory branches")
    return result


def compute_trajectory_entropy(data: SingleCellData, pseudotime_col: str, window_size: int = 100) -> Dict[str, Any]:
    """Compute entropy along a trajectory to identify differentiation points.

    Args:
        data: SingleCellData object with pseudotime
        pseudotime_col: Column name containing pseudotime values
        window_size: Size of sliding window for entropy calculation

    Returns:
        Dictionary with trajectory entropy information

    Raises:
        ValidationError: If data is not SingleCellData, the pseudotime column is not in data.obs, or window_size is not
            between 1 and the number of cells
    """
    validation.validate_type(data, SingleCellData, "data")

    if data.obs is None or pseudotime_col not in data.obs.columns:
        raise errors.ValidationError(f"Pseudotime column '{pseudotime_col}' not found in data.obs")

    logger.info(f"Computing trajectory entropy with window size {window_size}")

    # Get data
    X = data.X.toarray() if hasattr(data.X, "toarray") else data.X
    pseudotime = data.obs[pseudotime_col].values

    # Sort by pseudotime
    sorted_indices = np.argsort(pseudotime)
    X_sorted = X[sorted_indices]
    pseudotime_sorted = pseudotime[sorted_indices]
    if window_size < 1 or window_size > len(X_sorted):
        raise errors.ValidationError(
            f"window_size must be between 1 and the number of cells ({len(X_sorted)}), got {window_size}"
        )

    # Compute entropy in sliding windows
    n_windows = len(X_sorted) - window_size + 1
    entropy_profile = np.zeros(n_windows)

    for i in range(n_windows):
        window_data = X_sorted[i : i + window_size]

        # Compute gene expression entropy in this window
        gene_entropies = []
        for gene_idx in range(X.shape[1]):
            gene_expr = window_data[:, gene_idx]
            # Compute expression distribution entropy
            if np.sum(gene_expr) > 0:
                # Normalize to probability distribution
                probs = gene_expr / np.sum(gene_expr)
                probs = probs[probs > 0]  # Remove zeros
                entropy = -np.sum(probs * np.log2(probs))
                gene_entropies.append(entropy)

        # Average entropy across genes
        entropy_profile[i] = np.mean(gene_entropies) if gene_entropies else 0

    # Find differentiation points (local maxima in entropy)
    differentiation_points = _find_local_maxima(entropy_profile)

    result = {
        "entropy_profile": entropy_profile,
        "pseudotime_windows": pseudotime_sorted[window_size // 2 : window_size // 2 + n_windows],
        "differentiation_points": differentiation_points,
        "window_size": window_size,
        "n_windows": n_windows,
    }

    logger.info(f"Trajectory entropy computed: {len(differentiation_points)} differentiation points found")
    return result


def _compute_shortest_path_distances(distances: np.ndarray, root: int) -> np.ndarray:
    """Compute shortest path distances from root using Dijkstra's algorithm."""
    n = distances.shape[0]
    pseudotime = np.full(n, np.inf)
    pseudotime[root] = 0
    visited = np.zeros(n, dtype=bool)

    for _ in range(n):
        # Find unvisited node with smallest distance
        unvisited = np.where(~visited)[0]
        if len(unvisited) == 0:
            break

        current = unvisited[np.argmin(pseudotime[unvisited])]

        visited[current] = True

        # Update distances to neighbors
        for neighbor in range(n):
            if not visited[neighbor]:
                new_distance = pseudotime[current] + distances[current, neighbor]
                if new_distance < pseudotime[neighbor]:
                    pseudotime[neighbor] = new_distance

    # Replace infinite distances with maximum finite distance
    max_finite = np.max(pseudotime[pseudotime < np.inf])
    pseudotime[pseudotime == np.inf] = max_finite

    return pseudotime


def _compute_cluster_connectivity(cluster_labels: np.ndarray, unique_clusters: np.ndarray) -> np.ndarray:
    """Compute connectivity matrix between clusters."""
    n_clusters = len(unique_clusters)
    connectivity = np.zeros((n_clusters, n_clusters))

    cluster_to_idx = {cluster: i for i, cluster in enumerate(unique_clusters)}

    # Count transitions between clusters
    for i in range(len(cluster_labels) - 1):
        from_cluster = cluster_labels[i]
        to_cluster = cluster_labels[i + 1]

        from_idx = cluster_to_idx[from_cluster]
        to_idx = cluster_to_idx[to_cluster]

        connectivity[from_idx, to_idx] += 1

    # Normalize by row sums
    row_sums = connectivity.sum(axis=1)
    row_sums = np.where(row_sums == 0, 1, row_sums)  # Avoid division by zero
    connectivity = connectivity / row_sums[:, np.newaxis]

    return connectivity


def _compute_cluster_pseudotime(connectivity: np.ndarray, unique_clusters: np.ndarray) -> Dict[Any, float]:
    """Compute pseudotime for clusters based on connectivity."""
    n_clusters = len(unique_clusters)

    # Simple topological sort - assume linear progression
    # In practice, this would be more sophisticated
    cluster_pseudotime = {}

    # Find starting cluster (most incoming connections)
    incoming = connectivity.sum(axis=0)
    start_idx = int(np.argmin(incoming))  # Cluster with fewest incoming connections

    # Assign pseudotime based on graph distance from start
    distances = np.full(n_clusters, np.inf)
    distances[start_idx] = 0

    # Simple BFS to compute distances
    queue = [start_idx]
    visited = set([start_idx])

    while queue:
        current = queue.pop(0)

        for neighbor in range(n_clusters):
            if connectivity[current, neighbor] > 0 and neighbor not in visited:
                distances[neighbor] = distances[current] + 1
                visited.add(neighbor)
                queue.append(neighbor)

    # Normalize distances
    max_dist = np.max(distances[distances < np.inf])
    if max_dist > 0:
        distances = distances / max_dist

    # Map back to cluster labels
    for i, cluster in enumerate(unique_clusters):
        cluster_pseudotime[cluster] = distances[i]

    return cluster_pseudotime


def _compute_slingshot_trajectory(
    X: np.ndarray, cluster_labels: np.ndarray, start_cluster: str, end_clusters: List[str]
) -> Dict[str, Any]:
    """Compute simplified Slingshot-like trajectory."""
    unique_clusters = np.unique(cluster_labels)

    # Compute cluster centroids
    centroids = {}
    for cluster in unique_clusters:
        mask = cluster_labels == cluster
        centroids[cluster] = np.mean(X[mask], axis=0)

    # Build trajectory as path from start to end clusters
    branches = []

    for end_cluster in end_clusters:
        # Simple path: start -> end
        branch = {
            "start": start_cluster,
            "end": end_cluster,
            "path": [start_cluster, end_cluster],
            "length": np.linalg.norm(centroids[end_cluster] - centroids[start_cluster]),
        }
        branches.append(branch)

    # Assign pseudotime based on distance from start
    cluster_pseudotime = {}
    start_centroid = centroids[start_cluster]

    for cluster in unique_clusters:
        distance = np.linalg.norm(centroids[cluster] - start_centroid)
        max_distance = max(np.linalg.norm(centroids[c] - start_centroid) for c in unique_clusters)
        cluster_pseudotime[cluster] = distance / max_distance if max_distance > 0 else 0

    return {
        "branches": branches,
        "cluster_pseudotime": cluster_pseudotime,
        "centroids": centroids,
    }


def _detect_trajectory_branches(
    sorted_pseudotime: np.ndarray, sorted_indices: np.ndarray, min_branch_size: int
) -> Dict[int, np.ndarray]:
    """Detect branches in pseudotime-sorted data."""
    # Simple branch detection using density changes
    # This is a highly simplified implementation

    branches = {}
    current_branch = []
    branch_id = 0

    # Group cells into branches based on pseudotime gaps
    pseudotime_diff = np.diff(sorted_pseudotime)
    gap_threshold = np.percentile(pseudotime_diff, 90)  # Top 10% of gaps

    for i, idx in enumerate(sorted_indices):
        current_branch.append(idx)

        # Check if we should start a new branch
        if i < len(sorted_indices) - 1:
            gap = sorted_pseudotime[i + 1] - sorted_pseudotime[i]
            if gap > gap_threshold and len(current_branch) >= min_branch_size:
                if current_branch:
                    branches[branch_id] = np.array(current_branch)
                    branch_id += 1
                    current_branch = []

    # Add remaining cells
    if current_branch:
        branches[branch_id] = np.array(current_branch)

    return branches


def _find_local_maxima(data: np.ndarray, min_distance: int = 10) -> List[int]:
    """Find local maxima in a 1D array."""
    maxima = []

    for i in range(min_distance, len(data) - min_distance):
        if data[i] > data[i - min_distance] and data[i] > data[i + min_distance]:
            # Check if it's a true local maximum
            is_max = True
            for j in range(-min_distance, min_distance + 1):
                if j != 0 and data[i + j] >= data[i]:
                    is_max = False
                    break
            if is_max:
                maxima.append(i)

    return maxima
