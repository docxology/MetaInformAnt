"""Trajectory inference methods for single-cell data.

This module provides functions for inferring developmental trajectories and
pseudotemporal ordering from single-cell expression data. Methods include
diffusion pseudotime, DPT, PAGA, and simplified trajectory reconstruction.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Any, Tuple, Set
import numpy as np
import pandas as pd
from scipy import sparse

from metainformant.core import logging, errors, validation

# Optional scientific dependencies
try:
    from sklearn.neighbors import NearestNeighbors
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False
    NearestNeighbors = None

logger = logging.get_logger(__name__)

# Import our SingleCellData
from metainformant.singlecell.data.preprocessing import SingleCellData


def compute_diffusion_pseudotime(data: SingleCellData,
                                root_cell: int | None = None,
                                n_components: int = 10) -> SingleCellData:
    """Compute diffusion pseudotime for trajectory inference.

    Args:
        data: SingleCellData object with expression matrix
        root_cell: Index of root cell (starting point of trajectory)
        n_components: Number of diffusion components to compute

    Returns:
        SingleCellData with pseudotime added to obs

    Raises:
        TypeError: If data is not SingleCellData
        ValueError: If root_cell is invalid
    """
    validation.validate_type(data, SingleCellData, "data")

    if root_cell is not None:
        validation.validate_range(root_cell, min_val=0, max_val=data.n_obs - 1, name="root_cell")

    logger.info(f"Computing diffusion pseudotime with {n_components} components")

    # Create copy
    result = data.copy()

    # Get expression matrix
    X = data.X.toarray() if hasattr(data.X, 'toarray') else data.X

    # Compute diffusion map (simplified version)
    pseudotime = _compute_simple_diffusion_pseudotime(X, root_cell, n_components)

    # Add to obs
    result.obs = result.obs.copy() if result.obs is not None else pd.DataFrame()
    result.obs['dpt_pseudotime'] = pseudotime

    # Store metadata
    result.uns['dpt'] = {
        'method': 'diffusion_pseudotime',
        'root_cell': root_cell,
        'n_components': n_components,
        'min_pseudotime': float(np.min(pseudotime)),
        'max_pseudotime': float(np.max(pseudotime)),
    }

    logger.info(f"Diffusion pseudotime computed: range [{np.min(pseudotime):.3f}, {np.max(pseudotime):.3f}]")
    return result


def dpt_trajectory(data: SingleCellData, root_cell: int | None = None) -> SingleCellData:
    """Compute Diffusion Pseudotime (DPT) trajectory.

    Args:
        data: SingleCellData object with expression matrix
        root_cell: Index of root cell (starting point of trajectory)

    Returns:
        SingleCellData with DPT pseudotime and ordering

    Raises:
        TypeError: If data is not SingleCellData
    """
    # For this implementation, DPT is similar to diffusion pseudotime
    # In practice, DPT involves more sophisticated eigenvector analysis
    return compute_diffusion_pseudotime(data, root_cell, n_components=10)


def paga_trajectory(data: SingleCellData, groups: str) -> SingleCellData:
    """Compute Partition-based Graph Abstraction (PAGA) trajectory.

    Args:
        data: SingleCellData object with cluster assignments
        groups: Column name containing cluster/group labels

    Returns:
        SingleCellData with PAGA connectivity information

    Raises:
        TypeError: If data is not SingleCellData
        ValueError: If groups column not found
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
    X = data.X.toarray() if hasattr(data.X, 'toarray') else data.X
    cluster_centroids = {}

    for cluster in unique_clusters:
        mask = cluster_labels == cluster
        if np.sum(mask) > 0:
            cluster_centroids[cluster] = np.mean(X[mask], axis=0)

    # Compute connectivity between clusters
    connectivity_matrix = _compute_cluster_connectivity(cluster_labels, unique_clusters)

    # Store PAGA results
    result.uns['paga'] = {
        'groups': groups,
        'unique_groups': unique_clusters.tolist(),
        'connectivity_matrix': connectivity_matrix.tolist(),
        'n_groups': len(unique_clusters),
    }

    # Add cluster-level pseudotime (simplified)
    # In practice, PAGA would compute a trajectory through the abstracted graph
    cluster_pseudotime = _compute_cluster_pseudotime(connectivity_matrix, unique_clusters)

    # Assign pseudotime to cells based on their cluster
    cell_pseudotime = np.array([cluster_pseudotime[cluster] for cluster in cluster_labels])

    result.obs = result.obs.copy() if result.obs is not None else pd.DataFrame()
    result.obs['paga_pseudotime'] = cell_pseudotime

    logger.info(f"PAGA trajectory computed for {len(unique_clusters)} groups")
    return result


def slingshot_trajectory(data: SingleCellData, start_cluster: str,
                        end_clusters: List[str]) -> SingleCellData:
    """Compute Slingshot-like trajectory inference.

    Args:
        data: SingleCellData object with cluster assignments
        start_cluster: Starting cluster for trajectory
        end_clusters: List of ending clusters for trajectory branches

    Returns:
        SingleCellData with trajectory pseudotime

    Raises:
        TypeError: If data is not SingleCellData
        ValueError: If cluster columns not found or invalid
    """
    validation.validate_type(data, SingleCellData, "data")

    # Assume clusters are in a column named 'cluster' or similar
    cluster_col = None
    if data.obs is not None:
        # Find a column that looks like cluster assignments
        for col in data.obs.columns:
            if 'cluster' in col.lower() or 'leiden' in col.lower() or 'louvain' in col.lower():
                cluster_col = col
                break

    if cluster_col is None:
        raise errors.ValidationError("No cluster column found in data.obs")

    logger.info(f"Computing Slingshot trajectory from {start_cluster} to {end_clusters}")

    # Create copy
    result = data.copy()

    cluster_labels = data.obs[cluster_col].values

    # Compute cluster centroids and distances
    X = data.X.toarray() if hasattr(data.X, 'toarray') else data.X
    unique_clusters = np.unique(cluster_labels)

    # Check if start and end clusters exist
    if start_cluster not in unique_clusters:
        raise errors.ValidationError(f"Start cluster '{start_cluster}' not found in data")

    for end_cluster in end_clusters:
        if end_cluster not in unique_clusters:
            raise errors.ValidationError(f"End cluster '{end_cluster}' not found in data")

    # Compute simplified trajectory
    trajectory_info = _compute_slingshot_trajectory(
        X, cluster_labels, start_cluster, end_clusters
    )

    # Assign pseudotime to cells
    cell_pseudotime = np.zeros(data.n_obs)
    for i, cluster in enumerate(cluster_labels):
        if cluster in trajectory_info['cluster_pseudotime']:
            cell_pseudotime[i] = trajectory_info['cluster_pseudotime'][cluster]

    result.obs = result.obs.copy() if result.obs is not None else pd.DataFrame()
    result.obs['slingshot_pseudotime'] = cell_pseudotime

    # Store trajectory information
    result.uns['slingshot'] = {
        'start_cluster': start_cluster,
        'end_clusters': end_clusters,
        'trajectory_info': trajectory_info,
    }

    logger.info(f"Slingshot trajectory computed with {len(trajectory_info['branches'])} branches")
    return result


def compute_pseudotime_from_dimensionality_reduction(data: SingleCellData,
                                                   dim_red_cols: List[str],
                                                   root_cell: int | None = None) -> SingleCellData:
    """Compute pseudotime from dimensionality reduction coordinates.

    Args:
        data: SingleCellData object with dimensionality reduction coordinates
        dim_red_cols: Column names containing the reduced dimensions
        root_cell: Index of root cell

    Returns:
        SingleCellData with pseudotime computed from embedding

    Raises:
        TypeError: If data is not SingleCellData
        ValueError: If dimension columns not found
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

    # Compute pseudotime as geodesic distance from root
    if root_cell is None:
        # Choose cell with minimum norm as root
        norms = np.sum(embedding ** 2, axis=1)
        root_cell = np.argmin(norms)

    # Compute pairwise distances in embedding space
    from scipy.spatial.distance import pdist, squareform
    distances = squareform(pdist(embedding, metric='euclidean'))

    # Compute shortest path distances from root (simplified Dijkstra)
    pseudotime = _compute_shortest_path_distances(distances, root_cell)

    # Normalize pseudotime to [0, 1]
    pseudotime = (pseudotime - np.min(pseudotime)) / (np.max(pseudotime) - np.min(pseudotime))

    result.obs['embedding_pseudotime'] = pseudotime

    # Store metadata
    result.uns['embedding_pseudotime'] = {
        'dim_red_cols': dim_red_cols,
        'root_cell': root_cell,
        'method': 'shortest_path',
        'min_pseudotime': float(np.min(pseudotime)),
        'max_pseudotime': float(np.max(pseudotime)),
    }

    logger.info(f"Embedding-based pseudotime computed: root cell = {root_cell}")
    return result


def find_trajectory_branches(data: SingleCellData, pseudotime_col: str,
                           min_branch_size: int = 10) -> Dict[str, Any]:
    """Identify branches in a computed trajectory.

    Args:
        data: SingleCellData object with pseudotime
        pseudotime_col: Column name containing pseudotime values
        min_branch_size: Minimum number of cells per branch

    Returns:
        Dictionary with branch information

    Raises:
        TypeError: If data is not SingleCellData
        ValueError: If pseudotime column not found
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
        'n_branches': len(branches),
        'branch_sizes': [len(indices) for indices in branches.values()],
        'branch_labels': branch_labels,
        'branches': branches,
    }

    logger.info(f"Found {len(branches)} trajectory branches")
    return result


def compute_trajectory_entropy(data: SingleCellData, pseudotime_col: str,
                              window_size: int = 100) -> Dict[str, Any]:
    """Compute entropy along a trajectory to identify differentiation points.

    Args:
        data: SingleCellData object with pseudotime
        pseudotime_col: Column name containing pseudotime values
        window_size: Size of sliding window for entropy calculation

    Returns:
        Dictionary with trajectory entropy information

    Raises:
        TypeError: If data is not SingleCellData
        ValueError: If pseudotime column not found
    """
    validation.validate_type(data, SingleCellData, "data")

    if data.obs is None or pseudotime_col not in data.obs.columns:
        raise errors.ValidationError(f"Pseudotime column '{pseudotime_col}' not found in data.obs")

    logger.info(f"Computing trajectory entropy with window size {window_size}")

    # Get data
    X = data.X.toarray() if hasattr(data.X, 'toarray') else data.X
    pseudotime = data.obs[pseudotime_col].values

    # Sort by pseudotime
    sorted_indices = np.argsort(pseudotime)
    X_sorted = X[sorted_indices]
    pseudotime_sorted = pseudotime[sorted_indices]

    # Compute entropy in sliding windows
    n_windows = len(X_sorted) - window_size + 1
    entropy_profile = np.zeros(n_windows)

    for i in range(n_windows):
        window_data = X_sorted[i:i + window_size]

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
        'entropy_profile': entropy_profile,
        'pseudotime_windows': pseudotime_sorted[window_size//2 : window_size//2 + n_windows],
        'differentiation_points': differentiation_points,
        'window_size': window_size,
        'n_windows': n_windows,
    }

    logger.info(f"Trajectory entropy computed: {len(differentiation_points)} differentiation points found")
    return result


def _compute_simple_diffusion_pseudotime(X: np.ndarray, root_cell: int | None,
                                       n_components: int) -> np.ndarray:
    """Compute simplified diffusion pseudotime."""
    if root_cell is None:
        # Choose cell with highest expression variance as root
        variances = np.var(X, axis=0)
        root_cell = np.argmax(np.sum(X * variances, axis=1))

    # Compute pairwise distances
    from scipy.spatial.distance import pdist, squareform
    distances = squareform(pdist(X, metric='euclidean'))

    # Compute shortest path distances from root
    pseudotime = _compute_shortest_path_distances(distances, root_cell)

    return pseudotime


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


def _compute_cluster_connectivity(cluster_labels: np.ndarray,
                                unique_clusters: np.ndarray) -> np.ndarray:
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


def _compute_cluster_pseudotime(connectivity: np.ndarray,
                              unique_clusters: np.ndarray) -> Dict[Any, float]:
    """Compute pseudotime for clusters based on connectivity."""
    n_clusters = len(unique_clusters)

    # Simple topological sort - assume linear progression
    # In practice, this would be more sophisticated
    cluster_pseudotime = {}

    # Find starting cluster (most incoming connections)
    incoming = connectivity.sum(axis=0)
    start_idx = np.argmin(incoming)  # Cluster with fewest incoming connections

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


def _compute_slingshot_trajectory(X: np.ndarray, cluster_labels: np.ndarray,
                                start_cluster: str, end_clusters: List[str]) -> Dict[str, Any]:
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
            'start': start_cluster,
            'end': end_cluster,
            'path': [start_cluster, end_cluster],
            'length': np.linalg.norm(centroids[end_cluster] - centroids[start_cluster]),
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
        'branches': branches,
        'cluster_pseudotime': cluster_pseudotime,
        'centroids': centroids,
    }


def _detect_trajectory_branches(sorted_pseudotime: np.ndarray, sorted_indices: np.ndarray,
                              min_branch_size: int) -> Dict[int, np.ndarray]:
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






