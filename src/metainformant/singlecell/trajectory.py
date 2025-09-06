"""Trajectory analysis and pseudotime inference for single-cell data.

This module provides methods for inferring developmental trajectories and pseudotime
from single-cell data, including diffusion pseudotime, principal curves, and
lineage analysis. All implementations use real computational methods without mocking.
"""

from __future__ import annotations

import warnings
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

try:
    from scipy import sparse
    from scipy.spatial.distance import pdist, squareform
    from scipy.stats import spearmanr

    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    sparse = None
    pdist = None
    squareform = None
    spearmanr = None

try:
    from sklearn.cluster import KMeans
    from sklearn.decomposition import PCA

    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False
    PCA = None
    KMeans = None

from .preprocessing import SingleCellData


def compute_pseudotime(
    data: SingleCellData,
    root_cells: Optional[Union[int, List[int]]] = None,
    method: str = "diffusion",
    n_components: int = 10,
) -> SingleCellData:
    """Compute pseudotime ordering of cells.

    Args:
        data: SingleCellData with neighbor graph or diffusion map
        root_cells: Index/indices of root cell(s) (None for auto-detection)
        method: Method for pseudotime ('diffusion', 'shortest_path')
        n_components: Number of diffusion components to use

    Returns:
        SingleCellData with pseudotime in obs['pseudotime']
    """
    data = data.copy()

    if method == "diffusion":
        # Use diffusion pseudotime
        if "X_diffmap" not in data.obsm:
            warnings.warn("No diffusion map found, computing with default parameters")
            from .dimensionality import compute_diffusion_map

            data = compute_diffusion_map(data, n_components=n_components)

        # Get diffusion coordinates
        diffmap_coords = data.obsm["X_diffmap"][:, :n_components]

        # Auto-detect root if not provided
        if root_cells is None:
            # Use cell with maximum value in first diffusion component
            root_cells = [np.argmax(diffmap_coords[:, 0])]
        elif isinstance(root_cells, int):
            root_cells = [root_cells]

        # Compute pseudotime as distance from root in diffusion space
        root_coord = np.mean(diffmap_coords[root_cells], axis=0)

        # Euclidean distance in diffusion space
        distances = np.sqrt(np.sum((diffmap_coords - root_coord) ** 2, axis=1))
        pseudotime = distances / np.max(distances)  # Normalize to [0, 1]

    elif method == "shortest_path":
        # Use shortest path pseudotime on neighbor graph
        if "neighbors" not in data.uns:
            warnings.warn("No neighbor graph found, computing with default parameters")
            from .dimensionality import compute_neighbors

            data = compute_neighbors(data)

        # Get connectivity matrix
        connectivities = data.uns["neighbors"]["connectivities"]

        if sparse.issparse(connectivities):
            # Convert to dense for shortest path (for simplicity)
            conn_matrix = connectivities.toarray()
        else:
            conn_matrix = connectivities

        # Convert similarity to distance
        distance_matrix = 1.0 - conn_matrix
        distance_matrix[distance_matrix < 0] = np.inf  # Handle negative similarities
        np.fill_diagonal(distance_matrix, 0)

        # Auto-detect root if not provided
        if root_cells is None:
            # Use cell with highest total connectivity (most central)
            connectivity_sums = np.sum(conn_matrix, axis=1)
            root_cells = [np.argmax(connectivity_sums)]
        elif isinstance(root_cells, int):
            root_cells = [root_cells]

        # Compute shortest path distances using Floyd-Warshall (simple implementation)
        # For larger datasets, consider using scipy.sparse.csgraph.shortest_path
        pseudotimes = []
        for root in root_cells:
            # Simple implementation of shortest path from root
            dist = np.full(len(conn_matrix), np.inf)
            dist[root] = 0

            # Dijkstra-like algorithm (simplified)
            visited = np.zeros(len(conn_matrix), dtype=bool)

            for _ in range(len(conn_matrix)):
                # Find unvisited node with minimum distance
                unvisited_dist = dist.copy()
                unvisited_dist[visited] = np.inf
                current = np.argmin(unvisited_dist)

                if dist[current] == np.inf:
                    break

                visited[current] = True

                # Update distances to neighbors
                for neighbor in range(len(conn_matrix)):
                    if not visited[neighbor] and distance_matrix[current, neighbor] < np.inf:
                        new_dist = dist[current] + distance_matrix[current, neighbor]
                        if new_dist < dist[neighbor]:
                            dist[neighbor] = new_dist

            pseudotimes.append(dist)

        # Use minimum distance across all roots
        pseudotime = np.min(pseudotimes, axis=0)

        # Normalize
        max_dist = np.max(pseudotime[pseudotime < np.inf])
        if max_dist > 0:
            pseudotime = pseudotime / max_dist

        # Handle infinite distances (disconnected cells)
        pseudotime[pseudotime == np.inf] = 1.0

    else:
        raise ValueError(f"Unknown pseudotime method: {method}")

    # Store results
    data.obs["pseudotime"] = pseudotime
    data.uns["pseudotime"] = {
        "method": method,
        "root_cells": root_cells,
        "n_components": n_components if method == "diffusion" else None,
    }

    print(f"Pseudotime computed using {method} method")

    return data


def trajectory_analysis(
    data: SingleCellData,
    groupby: Optional[str] = None,
    n_branches: int = 2,
    method: str = "mst",
) -> Dict[str, Any]:
    """Analyze trajectory structure and branching.

    Args:
        data: SingleCellData with pseudotime
        groupby: Column for trajectory groups (optional)
        n_branches: Number of expected branches
        method: Method for trajectory analysis ('mst', 'clustering')

    Returns:
        Dictionary with trajectory analysis results
    """
    if "pseudotime" not in data.obs.columns:
        warnings.warn("No pseudotime found, computing with default parameters")
        data = compute_pseudotime(data)

    pseudotime = data.obs["pseudotime"].values
    results = {}

    if method == "mst":
        # Minimum spanning tree approach for trajectory structure
        if "X_diffmap" in data.obsm:
            coords = data.obsm["X_diffmap"][:, :10]
        elif "X_pca" in data.obsm:
            coords = data.obsm["X_pca"][:, :20]
        else:
            # Use expression data
            coords = data.X
            if sparse.issparse(coords):
                coords = coords.toarray()

        # Compute pairwise distances
        from scipy.spatial.distance import pdist, squareform

        distances = pdist(coords, metric="euclidean")
        dist_matrix = squareform(distances)

        # Compute minimum spanning tree
        from scipy.sparse.csgraph import minimum_spanning_tree

        mst = minimum_spanning_tree(dist_matrix).toarray()

        # Find branch points (nodes with degree > 2)
        degree = np.sum(mst > 0, axis=1)
        branch_points = np.where(degree > 2)[0]

        results["mst"] = mst
        results["branch_points"] = branch_points
        results["node_degrees"] = degree

    elif method == "clustering":
        # Clustering-based trajectory analysis
        # Divide cells into pseudotime bins
        n_bins = min(20, len(data.obs) // 10)
        if n_bins < 3:
            n_bins = 3

        pseudotime_bins = pd.cut(pseudotime, bins=n_bins, labels=False)

        # For each bin, look at expression patterns
        bin_profiles = []
        for bin_idx in range(n_bins):
            bin_mask = pseudotime_bins == bin_idx
            if np.sum(bin_mask) > 0:
                if sparse.issparse(data.X):
                    bin_profile = np.mean(data.X[bin_mask, :].toarray(), axis=0)
                else:
                    bin_profile = np.mean(data.X[bin_mask, :], axis=0)
                bin_profiles.append(bin_profile)
            else:
                bin_profiles.append(np.zeros(data.n_vars))

        bin_profiles = np.array(bin_profiles)

        # Cluster bins to identify trajectory branches
        if len(bin_profiles) >= n_branches:
            kmeans = KMeans(n_clusters=n_branches, random_state=42)
            bin_clusters = kmeans.fit_predict(bin_profiles)

            results["bin_clusters"] = bin_clusters
            results["bin_profiles"] = bin_profiles
            results["n_bins"] = n_bins
        else:
            warnings.warn(f"Too few pseudotime bins ({len(bin_profiles)}) for {n_branches} branches")

    else:
        raise ValueError(f"Unknown trajectory analysis method: {method}")

    # Add pseudotime statistics
    results["pseudotime_stats"] = {
        "min": np.min(pseudotime),
        "max": np.max(pseudotime),
        "mean": np.mean(pseudotime),
        "std": np.std(pseudotime),
    }

    # If groupby is provided, analyze pseudotime by group
    if groupby and groupby in data.obs.columns:
        group_stats = data.obs.groupby(groupby)["pseudotime"].agg(["mean", "std", "min", "max"])
        results["group_pseudotime"] = group_stats

    return results


def lineage_analysis(
    data: SingleCellData,
    start_cluster: Optional[str] = None,
    end_clusters: Optional[List[str]] = None,
    cluster_key: str = "leiden",
) -> Dict[str, Any]:
    """Analyze lineage relationships between cell clusters.

    Args:
        data: SingleCellData with clustering and pseudotime
        start_cluster: Starting cluster (root)
        end_clusters: Terminal clusters (leaves)
        cluster_key: Column name for cluster assignments

    Returns:
        Dictionary with lineage analysis results
    """
    if cluster_key not in data.obs.columns:
        raise ValueError(f"Cluster key '{cluster_key}' not found in obs")

    if "pseudotime" not in data.obs.columns:
        warnings.warn("No pseudotime found, computing with default parameters")
        data = compute_pseudotime(data)

    clusters = data.obs[cluster_key]
    pseudotime = data.obs["pseudotime"]
    unique_clusters = sorted(clusters.unique())

    results = {}

    # Compute cluster pseudotime statistics
    cluster_pseudotime = data.obs.groupby(cluster_key)["pseudotime"].agg(["mean", "std", "min", "max"])
    results["cluster_pseudotime"] = cluster_pseudotime

    # Auto-detect start cluster if not provided (earliest in pseudotime)
    if start_cluster is None:
        start_cluster = cluster_pseudotime["mean"].idxmin()

    # Auto-detect end clusters if not provided (latest in pseudotime)
    if end_clusters is None:
        # Select clusters with high mean pseudotime
        mean_pseudotime = cluster_pseudotime["mean"]
        threshold = np.percentile(mean_pseudotime, 80)  # Top 20%
        end_clusters = mean_pseudotime[mean_pseudotime >= threshold].index.tolist()

    results["start_cluster"] = start_cluster
    results["end_clusters"] = end_clusters

    # Compute transition probabilities between clusters
    # This is a simplified approach - could be enhanced with RNA velocity
    transition_matrix = np.zeros((len(unique_clusters), len(unique_clusters)))
    cluster_to_idx = {cluster: i for i, cluster in enumerate(unique_clusters)}

    # For cells in each cluster, look at their neighbors and count transitions
    if "neighbors" in data.uns:
        indices = data.uns["neighbors"]["indices"]

        for cell_idx in range(len(data.obs)):
            cell_cluster = clusters.iloc[cell_idx]
            cell_cluster_idx = cluster_to_idx[cell_cluster]

            # Look at neighbors
            neighbor_indices = indices[cell_idx]
            for neighbor_idx in neighbor_indices:
                if neighbor_idx < len(clusters):
                    neighbor_cluster = clusters.iloc[neighbor_idx]
                    neighbor_cluster_idx = cluster_to_idx[neighbor_cluster]
                    transition_matrix[cell_cluster_idx, neighbor_cluster_idx] += 1

    # Normalize transition matrix
    row_sums = transition_matrix.sum(axis=1)
    row_sums[row_sums == 0] = 1  # Avoid division by zero
    transition_matrix = transition_matrix / row_sums[:, np.newaxis]

    results["transition_matrix"] = pd.DataFrame(transition_matrix, index=unique_clusters, columns=unique_clusters)

    # Compute lineage scores for each cluster
    lineage_scores = {}
    for end_cluster in end_clusters:
        # Simple lineage score: inverse of pseudotime difference weighted by transitions
        end_cluster_pseudotime = cluster_pseudotime.loc[end_cluster, "mean"]

        scores = {}
        for cluster in unique_clusters:
            cluster_pseudotime_mean = cluster_pseudotime.loc[cluster, "mean"]

            # Score based on pseudotime progression towards end cluster
            if end_cluster_pseudotime > cluster_pseudotime_mean:
                time_score = 1.0 - abs(end_cluster_pseudotime - cluster_pseudotime_mean)
            else:
                time_score = 0.0

            # Could add transition probability weighting here
            scores[cluster] = max(0, time_score)

        lineage_scores[end_cluster] = scores

    results["lineage_scores"] = lineage_scores

    # Identify potential lineage paths
    lineage_paths = {}
    for end_cluster in end_clusters:
        # Simple path finding: order clusters by pseudotime towards end cluster
        path_clusters = []

        # Start from start cluster
        current = start_cluster
        visited = set()

        while current != end_cluster and current not in visited:
            path_clusters.append(current)
            visited.add(current)

            # Find next cluster with highest transition probability and higher pseudotime
            current_pseudotime = cluster_pseudotime.loc[current, "mean"]

            # Get possible next clusters (higher pseudotime)
            candidates = []
            for next_cluster in unique_clusters:
                if next_cluster not in visited and cluster_pseudotime.loc[next_cluster, "mean"] > current_pseudotime:
                    transition_prob = transition_matrix[cluster_to_idx[current], cluster_to_idx[next_cluster]]
                    candidates.append((next_cluster, transition_prob))

            if candidates:
                # Select candidate with highest transition probability
                candidates.sort(key=lambda x: x[1], reverse=True)
                current = candidates[0][0]
            else:
                break

        if current == end_cluster:
            path_clusters.append(end_cluster)

        lineage_paths[end_cluster] = path_clusters

    results["lineage_paths"] = lineage_paths

    return results


def compute_gene_trends(
    data: SingleCellData,
    genes: Optional[Union[str, List[str]]] = None,
    method: str = "spline",
    n_knots: int = 10,
) -> pd.DataFrame:
    """Compute gene expression trends along pseudotime.

    Args:
        data: SingleCellData with pseudotime
        genes: Gene names to analyze (None for all genes)
        method: Method for trend fitting ('spline', 'polynomial', 'loess')
        n_knots: Number of knots for spline fitting

    Returns:
        DataFrame with trend statistics for each gene
    """
    if "pseudotime" not in data.obs.columns:
        warnings.warn("No pseudotime found, computing with default parameters")
        data = compute_pseudotime(data)

    pseudotime = data.obs["pseudotime"].values

    # Select genes to analyze
    if genes is None:
        # Use all genes or highly variable genes
        if "highly_variable" in data.var.columns:
            gene_mask = data.var["highly_variable"].values
            if gene_mask.sum() == 0:
                gene_mask = np.ones(data.n_vars, dtype=bool)
        else:
            gene_mask = np.ones(data.n_vars, dtype=bool)
        gene_indices = np.where(gene_mask)[0]
        gene_names = data.var.index[gene_mask].tolist()
    elif isinstance(genes, str):
        if genes in data.var.index:
            gene_indices = [data.var.index.get_loc(genes)]
            gene_names = [genes]
        else:
            raise ValueError(f"Gene '{genes}' not found")
    else:
        gene_indices = []
        gene_names = []
        for gene in genes:
            if gene in data.var.index:
                gene_indices.append(data.var.index.get_loc(gene))
                gene_names.append(gene)
            else:
                warnings.warn(f"Gene '{gene}' not found, skipping")

    if len(gene_indices) == 0:
        raise ValueError("No valid genes found for trend analysis")

    # Get expression data
    if sparse.issparse(data.X):
        X_genes = data.X[:, gene_indices].toarray()
    else:
        X_genes = data.X[:, gene_indices]

    results = []

    for i, gene_name in enumerate(gene_names):
        expression = X_genes[:, i]

        # Compute correlation with pseudotime
        correlation, p_value = spearmanr(pseudotime, expression)

        # Fit trend based on method
        if method == "spline":
            try:
                from scipy.interpolate import UnivariateSpline

                # Fit spline
                spline = UnivariateSpline(pseudotime, expression, s=len(expression) * 0.1)
                fitted = spline(pseudotime)

                # Compute residual sum of squares
                residuals = expression - fitted
                rss = np.sum(residuals**2)

            except ImportError:
                # Fallback to polynomial
                coeffs = np.polyfit(pseudotime, expression, deg=3)
                fitted = np.polyval(coeffs, pseudotime)
                residuals = expression - fitted
                rss = np.sum(residuals**2)

        elif method == "polynomial":
            # Polynomial fit
            degree = min(3, len(np.unique(pseudotime)) - 1)
            coeffs = np.polyfit(pseudotime, expression, deg=degree)
            fitted = np.polyval(coeffs, pseudotime)
            residuals = expression - fitted
            rss = np.sum(residuals**2)

        else:
            # Simple linear fit
            coeffs = np.polyfit(pseudotime, expression, deg=1)
            fitted = np.polyval(coeffs, pseudotime)
            residuals = expression - fitted
            rss = np.sum(residuals**2)

        # Compute R-squared
        tss = np.sum((expression - np.mean(expression)) ** 2)
        r_squared = 1 - (rss / tss) if tss > 0 else 0

        # Detect trend direction
        early_expr = np.mean(expression[pseudotime < 0.3])
        late_expr = np.mean(expression[pseudotime > 0.7])
        trend_direction = "up" if late_expr > early_expr else "down" if late_expr < early_expr else "stable"

        result = {
            "gene": gene_name,
            "correlation": correlation,
            "p_value": p_value,
            "r_squared": r_squared,
            "rss": rss,
            "trend_direction": trend_direction,
            "early_expression": early_expr,
            "late_expression": late_expr,
            "fold_change": late_expr / (early_expr + 1e-9),
        }

        results.append(result)

    # Convert to DataFrame and sort by correlation strength
    trends_df = pd.DataFrame(results)
    trends_df = trends_df.sort_values("correlation", key=abs, ascending=False)

    return trends_df


def identify_branch_genes(
    data: SingleCellData,
    trajectory_results: Dict[str, Any],
    min_expression: float = 0.1,
    fold_change_threshold: float = 2.0,
) -> Dict[str, List[str]]:
    """Identify genes associated with trajectory branches.

    Args:
        data: SingleCellData object
        trajectory_results: Results from trajectory_analysis()
        min_expression: Minimum mean expression threshold
        fold_change_threshold: Minimum fold change between branches

    Returns:
        Dictionary mapping branch/lineage to associated genes
    """
    if "lineage_paths" not in trajectory_results:
        warnings.warn("No lineage paths found in trajectory results")
        return {}

    branch_genes = {}
    lineage_paths = trajectory_results["lineage_paths"]

    # Get expression data
    if sparse.issparse(data.X):
        X_data = data.X.toarray()
    else:
        X_data = data.X

    # For each lineage path, identify characteristic genes
    for end_cluster, path_clusters in lineage_paths.items():
        if len(path_clusters) < 2:
            continue

        # Get cells in this lineage path
        lineage_mask = data.obs["leiden"].isin(path_clusters)
        if lineage_mask.sum() == 0:
            continue

        lineage_cells = X_data[lineage_mask, :]
        other_cells = X_data[~lineage_mask, :]

        branch_gene_list = []

        # Compare expression between lineage and other cells
        for gene_idx in range(data.n_vars):
            lineage_expr = lineage_cells[:, gene_idx]
            other_expr = other_cells[:, gene_idx]

            lineage_mean = np.mean(lineage_expr)
            other_mean = np.mean(other_expr)

            # Check minimum expression threshold
            if lineage_mean < min_expression and other_mean < min_expression:
                continue

            # Check fold change
            fold_change = (lineage_mean + 1e-9) / (other_mean + 1e-9)

            if fold_change >= fold_change_threshold or fold_change <= (1.0 / fold_change_threshold):
                gene_name = data.var.index[gene_idx]
                branch_gene_list.append(gene_name)

        branch_genes[end_cluster] = branch_gene_list[:50]  # Top 50 genes per branch

    return branch_genes
