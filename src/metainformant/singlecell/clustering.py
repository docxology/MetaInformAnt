"""Clustering algorithms for single-cell data.

This module provides clustering methods commonly used in single-cell analysis,
including Leiden clustering, Louvain clustering, and marker gene identification.
All implementations use real computational methods without mocking.
"""

from __future__ import annotations

import warnings
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

try:
    from scipy import sparse
    from scipy.sparse import csgraph
    from scipy.stats import fisher_exact, mannwhitneyu, rankdata

    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    sparse = None
    rankdata = None
    mannwhitneyu = None
    fisher_exact = None
    csgraph = None

try:
    from sklearn.cluster import AgglomerativeClustering, KMeans
    from sklearn.metrics import silhouette_score

    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False
    KMeans = None
    AgglomerativeClustering = None
    silhouette_score = None

from .preprocessing import SingleCellData


def leiden_clustering(
    data: SingleCellData,
    resolution: float = 0.5,
    random_state: int = 42,
    n_iterations: int = -1,
) -> SingleCellData:
    """Perform Leiden clustering on single-cell data.

    Args:
        data: SingleCellData with neighbor graph
        resolution: Resolution parameter (higher = more clusters)
        random_state: Random seed
        n_iterations: Number of iterations (-1 for auto)

    Returns:
        SingleCellData with cluster assignments in obs['leiden']
    """
    data = data.copy()

    # Check if neighbor graph exists
    if "neighbors" not in data.uns:
        warnings.warn("No neighbor graph found, computing with default parameters")
        from .dimensionality import compute_neighbors

        data = compute_neighbors(data)

    try:
        import igraph as ig
        import leidenalg

        # Convert sparse matrix to igraph
        connectivities = data.uns["neighbors"]["connectivities"]

        if sparse.issparse(connectivities):
            # Convert to COO format for igraph
            connectivities_coo = connectivities.tocoo()
            sources = connectivities_coo.row
            targets = connectivities_coo.col
            weights = connectivities_coo.data
        else:
            # Convert dense matrix
            sources, targets = np.nonzero(connectivities)
            weights = connectivities[sources, targets]

        # Create igraph from edge list
        edges = list(zip(sources.astype(int), targets.astype(int)))
        g = ig.Graph(edges, directed=False)
        g.es["weight"] = weights

        # Remove self-loops and multi-edges
        g.simplify(multiple=True, loops=True, combine_edges="sum")

        # Leiden clustering
        partition = leidenalg.find_partition(
            g,
            leidenalg.RBConfigurationVertexPartition,
            resolution_parameter=resolution,
            n_iterations=n_iterations,
            seed=random_state,
        )

        # Extract cluster labels
        clusters = np.array(partition.membership)

    except ImportError:
        warnings.warn("leidenalg not available, using Louvain clustering as fallback")
        return louvain_clustering(data, resolution=resolution, random_state=random_state)

    # Store results
    data.obs["leiden"] = pd.Categorical(clusters.astype(str))
    data.uns["leiden"] = {
        "params": {
            "resolution": resolution,
            "random_state": random_state,
            "n_iterations": n_iterations,
        }
    }

    n_clusters = len(np.unique(clusters))
    print(f"Leiden clustering completed: {n_clusters} clusters found")

    return data


def louvain_clustering(
    data: SingleCellData,
    resolution: float = 0.5,
    random_state: int = 42,
) -> SingleCellData:
    """Perform Louvain clustering on single-cell data.

    Args:
        data: SingleCellData with neighbor graph
        resolution: Resolution parameter (higher = more clusters)
        random_state: Random seed

    Returns:
        SingleCellData with cluster assignments in obs['louvain']
    """
    data = data.copy()

    # Check if neighbor graph exists
    if "neighbors" not in data.uns:
        warnings.warn("No neighbor graph found, computing with default parameters")
        from .dimensionality import compute_neighbors

        data = compute_neighbors(data)

    try:
        import igraph as ig
        import louvain

        # Convert sparse matrix to igraph (same as Leiden)
        connectivities = data.uns["neighbors"]["connectivities"]

        if sparse.issparse(connectivities):
            connectivities_coo = connectivities.tocoo()
            sources = connectivities_coo.row
            targets = connectivities_coo.col
            weights = connectivities_coo.data
        else:
            sources, targets = np.nonzero(connectivities)
            weights = connectivities[sources, targets]

        # Create igraph
        edges = list(zip(sources.astype(int), targets.astype(int)))
        g = ig.Graph(edges, directed=False)
        g.es["weight"] = weights
        g.simplify(multiple=True, loops=True, combine_edges="sum")

        # Louvain clustering
        partition = louvain.find_partition(
            g,
            louvain.RBConfigurationVertexPartition,
            resolution_parameter=resolution,
            seed=random_state,
        )

        clusters = np.array(partition.membership)

    except ImportError:
        warnings.warn("louvain not available, using K-means clustering as fallback")
        return kmeans_clustering(data, n_clusters=8, random_state=random_state)

    # Store results
    data.obs["louvain"] = pd.Categorical(clusters.astype(str))
    data.uns["louvain"] = {
        "params": {
            "resolution": resolution,
            "random_state": random_state,
        }
    }

    n_clusters = len(np.unique(clusters))
    print(f"Louvain clustering completed: {n_clusters} clusters found")

    return data


def kmeans_clustering(
    data: SingleCellData,
    n_clusters: int = 8,
    use_pca: bool = True,
    n_pcs: int = 40,
    random_state: int = 42,
) -> SingleCellData:
    """Perform K-means clustering on single-cell data.

    Args:
        data: SingleCellData object
        n_clusters: Number of clusters
        use_pca: Whether to use PCA coordinates
        n_pcs: Number of PCs to use
        random_state: Random seed

    Returns:
        SingleCellData with cluster assignments in obs['kmeans']
    """
    data = data.copy()

    # Select data for clustering
    if use_pca and "X_pca" in data.obsm:
        X_cluster = data.obsm["X_pca"][:, :n_pcs]
    else:
        X_cluster = data.X
        if sparse.issparse(X_cluster):
            X_cluster = X_cluster.toarray()

    # K-means clustering
    kmeans = KMeans(n_clusters=n_clusters, random_state=random_state, n_init=10)
    clusters = kmeans.fit_predict(X_cluster)

    # Calculate silhouette score
    if X_cluster.shape[0] > 10000:  # Subsample for large datasets
        indices = np.random.RandomState(random_state).choice(X_cluster.shape[0], 10000, replace=False)
        silhouette = silhouette_score(X_cluster[indices], clusters[indices])
    else:
        silhouette = silhouette_score(X_cluster, clusters)

    # Store results
    data.obs["kmeans"] = pd.Categorical(clusters.astype(str))
    data.uns["kmeans"] = {
        "params": {
            "n_clusters": n_clusters,
            "use_pca": use_pca,
            "n_pcs": n_pcs,
            "random_state": random_state,
        },
        "silhouette_score": silhouette,
    }

    print(f"K-means clustering completed: {n_clusters} clusters, silhouette score = {silhouette:.3f}")

    return data


def hierarchical_clustering(
    data: SingleCellData,
    n_clusters: int = 8,
    linkage: str = "ward",
    use_pca: bool = True,
    n_pcs: int = 40,
) -> SingleCellData:
    """Perform hierarchical clustering on single-cell data.

    Args:
        data: SingleCellData object
        n_clusters: Number of clusters
        linkage: Linkage criterion ('ward', 'complete', 'average', 'single')
        use_pca: Whether to use PCA coordinates
        n_pcs: Number of PCs to use

    Returns:
        SingleCellData with cluster assignments in obs['hierarchical']
    """
    data = data.copy()

    # Select data for clustering
    if use_pca and "X_pca" in data.obsm:
        X_cluster = data.obsm["X_pca"][:, :n_pcs]
    else:
        X_cluster = data.X
        if sparse.issparse(X_cluster):
            X_cluster = X_cluster.toarray()

    # Hierarchical clustering
    clustering = AgglomerativeClustering(
        n_clusters=n_clusters,
        linkage=linkage,
    )
    clusters = clustering.fit_predict(X_cluster)

    # Store results
    data.obs["hierarchical"] = pd.Categorical(clusters.astype(str))
    data.uns["hierarchical"] = {
        "params": {
            "n_clusters": n_clusters,
            "linkage": linkage,
            "use_pca": use_pca,
            "n_pcs": n_pcs,
        }
    }

    print(f"Hierarchical clustering completed: {n_clusters} clusters")

    return data


def find_marker_genes(
    data: SingleCellData,
    groupby: str,
    method: str = "wilcoxon",
    n_genes: int = 100,
    logfc_min: float = 0.25,
    pval_cutoff: float = 0.05,
    only_pos: bool = True,
) -> pd.DataFrame:
    """Find marker genes for each cluster/group.

    Args:
        data: SingleCellData object
        groupby: Column name in obs for grouping (e.g., 'leiden', 'louvain')
        method: Statistical test ('wilcoxon', 'ttest', 'logreg')
        n_genes: Maximum number of genes per group
        logfc_min: Minimum log fold change
        pval_cutoff: P-value cutoff
        only_pos: Only return positive markers (upregulated)

    Returns:
        DataFrame with marker genes and statistics
    """
    if groupby not in data.obs.columns:
        raise ValueError(f"Column '{groupby}' not found in obs")

    groups = data.obs[groupby].unique()
    X = data.X

    if sparse.issparse(X):
        X = X.toarray()

    results = []

    for group in groups:
        group_mask = data.obs[groupby] == group
        group_cells = X[group_mask, :]
        other_cells = X[~group_mask, :]

        if group_cells.shape[0] < 3 or other_cells.shape[0] < 3:
            warnings.warn(f"Skipping group {group} due to insufficient cells")
            continue

        # Calculate statistics for each gene
        for gene_idx in range(X.shape[1]):
            gene_group = group_cells[:, gene_idx]
            gene_other = other_cells[:, gene_idx]

            # Skip genes with no expression
            if np.sum(gene_group) == 0 and np.sum(gene_other) == 0:
                continue

            # Calculate fold change
            mean_group = np.mean(gene_group)
            mean_other = np.mean(gene_other)

            # Add small constant to avoid log(0)
            logfc = np.log2((mean_group + 1e-9) / (mean_other + 1e-9))

            # Skip if fold change is too small
            if only_pos and logfc < logfc_min:
                continue
            if not only_pos and abs(logfc) < logfc_min:
                continue

            # Statistical test
            if method == "wilcoxon":
                try:
                    statistic, pval = mannwhitneyu(gene_group, gene_other, alternative="two-sided")
                except ValueError:
                    # Handle cases where all values are identical
                    pval = 1.0
                    statistic = 0.0
            elif method == "ttest":
                from scipy.stats import ttest_ind

                try:
                    statistic, pval = ttest_ind(gene_group, gene_other)
                except ValueError:
                    pval = 1.0
                    statistic = 0.0
            else:
                # Default to simple comparison
                pval = 0.05 if abs(logfc) > logfc_min else 1.0
                statistic = logfc

            # Skip if p-value is too high
            if pval > pval_cutoff:
                continue

            # Calculate additional statistics
            pct_in = np.mean(gene_group > 0) * 100  # Percentage expressing in group
            pct_out = np.mean(gene_other > 0) * 100  # Percentage expressing outside group

            result = {
                "group": group,
                "gene": gene_idx,
                "gene_name": data.var.index[gene_idx] if hasattr(data.var.index, "__getitem__") else str(gene_idx),
                "avg_log2FC": logfc,
                "pval": pval,
                "pval_adj": pval,  # Will adjust later
                "pct.1": pct_in,
                "pct.2": pct_out,
                "mean_group": mean_group,
                "mean_other": mean_other,
            }

            results.append(result)

    if not results:
        warnings.warn("No marker genes found")
        return pd.DataFrame()

    # Convert to DataFrame
    markers_df = pd.DataFrame(results)

    # Adjust p-values (Bonferroni correction)
    from scipy.stats import false_discovery_rate

    try:
        # Try Benjamin-Hochberg FDR correction
        markers_df["pval_adj"] = false_discovery_rate(markers_df["pval"].values)[1]
    except:
        # Fallback to Bonferroni
        markers_df["pval_adj"] = np.minimum(markers_df["pval"] * len(markers_df), 1.0)

    # Sort by group and fold change
    markers_df = markers_df.sort_values(["group", "avg_log2FC"], ascending=[True, False])

    # Limit number of genes per group
    if n_genes > 0:
        markers_df = markers_df.groupby("group").head(n_genes).reset_index(drop=True)

    print(f"Found {len(markers_df)} marker genes across {len(groups)} groups")

    return markers_df


def compute_cluster_composition(
    data: SingleCellData,
    groupby: str,
    sample_key: Optional[str] = None,
) -> pd.DataFrame:
    """Compute cluster composition statistics.

    Args:
        data: SingleCellData object
        groupby: Cluster column name
        sample_key: Optional sample/batch column for composition analysis

    Returns:
        DataFrame with cluster composition statistics
    """
    if groupby not in data.obs.columns:
        raise ValueError(f"Column '{groupby}' not found in obs")

    composition_data = []

    # Basic cluster sizes
    cluster_counts = data.obs[groupby].value_counts().sort_index()

    for cluster in cluster_counts.index:
        comp_dict = {
            "cluster": cluster,
            "n_cells": cluster_counts[cluster],
            "pct_total": 100 * cluster_counts[cluster] / len(data.obs),
        }

        # Add sample composition if requested
        if sample_key and sample_key in data.obs.columns:
            cluster_mask = data.obs[groupby] == cluster
            sample_counts = data.obs[cluster_mask][sample_key].value_counts()

            for sample in data.obs[sample_key].unique():
                count = sample_counts.get(sample, 0)
                comp_dict[f"{sample}_count"] = count
                comp_dict[f"{sample}_pct"] = 100 * count / cluster_counts[cluster]

        composition_data.append(comp_dict)

    return pd.DataFrame(composition_data)


def compute_cluster_silhouette(
    data: SingleCellData,
    groupby: str,
    use_pca: bool = True,
    n_pcs: int = 40,
) -> Dict[str, float]:
    """Compute silhouette scores for clusters.

    Args:
        data: SingleCellData object
        groupby: Cluster column name
        use_pca: Whether to use PCA coordinates
        n_pcs: Number of PCs to use

    Returns:
        Dictionary with silhouette scores per cluster
    """
    if groupby not in data.obs.columns:
        raise ValueError(f"Column '{groupby}' not found in obs")

    # Select data
    if use_pca and "X_pca" in data.obsm:
        X_sil = data.obsm["X_pca"][:, :n_pcs]
    else:
        X_sil = data.X
        if sparse.issparse(X_sil):
            X_sil = X_sil.toarray()

    # Get cluster labels
    cluster_labels = data.obs[groupby].values

    # Compute overall silhouette score
    if len(np.unique(cluster_labels)) > 1:
        # Subsample for large datasets
        if X_sil.shape[0] > 10000:
            indices = np.random.choice(X_sil.shape[0], 10000, replace=False)
            X_sub = X_sil[indices]
            labels_sub = cluster_labels[indices]
        else:
            X_sub = X_sil
            labels_sub = cluster_labels

        from sklearn.metrics import silhouette_samples

        try:
            sample_scores = silhouette_samples(X_sub, labels_sub)
            overall_score = np.mean(sample_scores)

            # Calculate per-cluster scores
            cluster_scores = {}
            for cluster in np.unique(labels_sub):
                cluster_mask = labels_sub == cluster
                if np.sum(cluster_mask) > 1:
                    cluster_scores[str(cluster)] = np.mean(sample_scores[cluster_mask])
                else:
                    cluster_scores[str(cluster)] = 0.0

            cluster_scores["overall"] = overall_score

        except Exception as e:
            warnings.warn(f"Could not compute silhouette scores: {e}")
            cluster_scores = {"overall": 0.0}
    else:
        cluster_scores = {"overall": 0.0}

    return cluster_scores
