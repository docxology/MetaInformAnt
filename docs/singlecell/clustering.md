# Single-Cell Clustering

`metainformant.singlecell.analysis.clustering` groups the cells of a `SingleCellData` object and evaluates the result. The four clustering functions are fitted on the expression matrix `data.X` itself, not on stored PCA or UMAP coordinates. Each returns a modified copy and leaves its input untouched (a `SingleCellData` goes in, a `SingleCellData` comes out), writes the labels to a column of `data.obs` and records its parameters in `data.uns`:

| Function | Labels in `data.obs` | Parameters in `data.uns` | Needs |
|----------|----------------------|--------------------------|-------|
| `leiden_clustering` | `leiden_cluster` | `leiden_clustering` | `networkx`, `igraph`, `leidenalg` |
| `louvain_clustering` | `louvain_cluster` | `louvain_clustering` | `networkx`, `python-louvain` |
| `kmeans_clustering` | `kmeans_cluster` | `kmeans_clustering` | `scikit-learn` |
| `hierarchical_clustering` | `hierarchical_cluster` | `hierarchical_clustering` | `scipy` (and `scikit-learn` above 1000 cells) |

The label column name is fixed (there is no `key_added` or `use_rep` argument). To keep several clusterings side by side, copy the column out of each result, as in [Multi-resolution clustering](#multi-resolution-clustering).

## Graph-Based Clustering

### leiden_clustering()

Leiden community detection on a k-nearest-neighbor graph of the cells:

```python
from metainformant.singlecell.analysis.clustering import leiden_clustering

data = leiden_clustering(
    data,
    resolution=0.5,       # Higher = more clusters
    n_neighbors=15,       # Neighbors per cell in the kNN graph
    random_state=42,      # Seed of the Leiden optimiser
)

clusters = data.obs["leiden_cluster"]
print(data.uns["leiden_clustering"])
# {'method': 'leiden', 'resolution': 0.5, 'n_neighbors': 15, 'random_state': 42,
#  'n_clusters': ..., 'use_weights': True}
```

How it works:

- The kNN graph is built from the Euclidean distances between the rows of `data.X`. The full cells × cells distance matrix is computed, so memory grows with the square of the number of cells.
- The partition optimises RB configuration modularity (`leidenalg.RBConfigurationVertexPartition`) at the requested `resolution`; `random_state` is the optimiser seed.
- `use_weights` is only recorded in `uns["leiden_clustering"]["use_weights"]`. Edge weights are not passed to leidenalg, so the partition is computed on the unweighted graph whatever its value.
- `uns["leiden_clustering"]["method"]` is always `"leiden"`: for a `SingleCellData` the function never substitutes another algorithm.

The Leiden stack (`networkx`, `igraph` and `leidenalg`) is optional. The `singlecell` extra installs it (`leidenalg` depends on `igraph`). Without it, `leiden_clustering` on a `SingleCellData` raises an `ImportError` that names the required packages; call `kmeans_clustering` if a k-means result is what you want.

#### Matrix input and the k-means fallback

`leiden_clustering` also accepts a plain cells × features matrix (dense, sparse or array-like, passed as `data` or as `feature_matrix=`) and then returns a label array instead of a `SingleCellData`. With the Leiden stack installed this is the same Leiden clustering. Without it the call raises `ImportError`, unless `allow_kmeans_fallback=True`:

```python
from metainformant.singlecell.analysis.clustering import leiden_clustering

labels = leiden_clustering(embedding, resolution=0.5, random_state=0, allow_kmeans_fallback=True)
```

With the fallback enabled and the stack missing, the rows are clustered with k-means and the function emits a `UserWarning` and a log warning saying so. The number of clusters is `round(4 * resolution)`, at least 2 and at most 10 (and never more than the number of rows); `random_state=None` uses seed 0, and `n_neighbors` and `use_weights` have no effect. A label array cannot record which algorithm produced it, so watch for the warning or check the module flags `HAS_NETWORKX`, `HAS_IGRAPH` and `HAS_LEIDEN` of `metainformant.singlecell.analysis.clustering`. The flag is ignored for a `SingleCellData` argument, which always raises when the stack is missing. A matrix with other than two dimensions raises `ValueError`.

### louvain_clustering()

Louvain community detection on the same kNN graph:

```python
from metainformant.singlecell.analysis.clustering import louvain_clustering

data = louvain_clustering(data, resolution=0.5, n_neighbors=15, random_state=42)
clusters = data.obs["louvain_cluster"]
```

The signature is `louvain_clustering(data, resolution=1.0, n_neighbors=15, random_state=None, use_weights=True)`. `uns["louvain_clustering"]` holds `resolution`, `n_neighbors`, `random_state`, `n_clusters` and `use_weights` (recorded only; the argument is never read). It needs `networkx` and `python-louvain`; install `python-louvain` (not `louvain`), which provides the `community` module the function imports. Without them it raises `ImportError`.

## Centroid-Based Clustering

### kmeans_clustering()

K-means for a predetermined number of clusters:

```python
from metainformant.singlecell.analysis.clustering import kmeans_clustering

data = kmeans_clustering(data, n_clusters=8, random_state=42, n_init=20)
```

**Parameters:**

- `n_clusters`: number of clusters, at least 2 (default 10)
- `random_state`: seed of the initialisation (default `None`)
- `n_init`: number of random initialisations (default 10)

`uns["kmeans_clustering"]` records the parameters plus `inertia`, `n_iter` and `cluster_centers_shape`.

## Hierarchical Clustering

### hierarchical_clustering()

Agglomerative clustering with scipy's linkage:

```python
from metainformant.singlecell.analysis.clustering import hierarchical_clustering

data = hierarchical_clustering(
    data,
    n_clusters=10,
    linkage_method="ward",   # the keyword is linkage_method
    metric="euclidean",
    random_state=0,
)
```

**Linkage methods:** `"single"`, `"complete"`, `"average"`, `"weighted"`, `"centroid"`, `"median"` and `"ward"`; any other value raises `ValidationError`. Ward linkage always uses the Euclidean metric.

Above 1000 cells the linkage is built on a random subset of 1000 cells and every other cell takes the cluster of its nearest subset cell. `random_state` seeds that subset (`None`, the default, draws a different subset on each call); it only matters above 1000 cells. Labels start at 1 (scipy's `fcluster`), and `uns["hierarchical_clustering"]` records the parameters, `linkage_matrix_shape` and `random_state`.

## Marker Gene Analysis

Two functions rank the genes that characterise each cluster.

### find_marker_genes()

Marker statistics for the groups in a column of `data.obs`:

```python
from metainformant.singlecell.analysis.clustering import find_marker_genes

markers = find_marker_genes(data, groupby="leiden_cluster", method="wilcoxon", n_genes=50)
print(markers.head())
# Columns: gene, group, mean_expr, mean_expr_other, log_fold_change, pct_expressed,
#          statistic, p_value, method
```

**Methods:** `"t-test"` (the default; its p-value uses a simplified approximation of the Student t distribution) and `"wilcoxon"` (Wilcoxon rank-sum test from scipy). Any other value, such as `"ttest"` or `"logistic"`, raises `ValidationError`. Each group keeps its `n_genes` genes with the smallest p-values (default 100); the p-values are not adjusted for multiple testing. Gene names come from `data.var.index`.

### find_markers()

A lighter ranking that works on any cells × genes matrix plus one label per cell:

```python
from metainformant.singlecell.analysis.clustering import find_markers

markers = find_markers(
    data.X,
    data.obs["leiden_cluster"].to_numpy(),
    n_markers=10,
    gene_names=data.var.index.tolist(),
)
# {cluster: [{'gene', 'gene_index', 'score', 'mean_expr', 'mean_expr_other', 'log_fold_change'}, ...]}
top_genes = {cluster: [row["gene"] for row in rows] for cluster, rows in markers.items()}
```

Genes are ranked by the difference between the mean expression inside and outside the cluster. No statistical test is run, so the result has no p-values. `n_markers` and `gene_names` are keyword-only; without `gene_names` the genes are called `gene_0`, `gene_1`, ... The function raises `ValueError` when `cluster_labels` and the matrix disagree on the number of cells.

### compute_cluster_composition()

How each cluster is composed of a categorical column of `obs` such as sample or batch:

```python
from metainformant.singlecell.analysis.clustering import compute_cluster_composition

composition = compute_cluster_composition(data, groupby="sample", cluster_col="leiden_cluster")
# Columns: cluster, category, count, proportion, percentage
```

`cluster_col` defaults to `"cluster"`, which none of the clustering functions write, so pass the column name explicitly.

## Clustering Validation

### compute_cluster_silhouette()

Silhouette and related scores for an existing clustering:

```python
from metainformant.singlecell.analysis.clustering import compute_cluster_silhouette

scores = compute_cluster_silhouette(
    data,
    cluster_col="leiden_cluster",
    metric="euclidean",
    sample_size=5000,     # score a random subsample of the cells
    random_state=0,
)

print(f"Average silhouette score: {scores['overall_silhouette_score']:.3f}")
# Also available: scores['per_cluster_scores'] ({cluster: {'mean', 'std', 'min', 'max'}}),
# scores['calinski_harabasz_score'], scores['davies_bouldin_score'],
# scores['n_clusters'], scores['n_samples'], scores['metric']
```

`sample_size` only draws a subsample when it is smaller than the number of cells; `random_state` seeds that draw (`None` gives a different, non-reproducible subsample on each call). With fewer than two clusters the function returns `{"error": ...}` instead of raising.

### evaluate_clustering_performance()

```python
from metainformant.singlecell.analysis.clustering import evaluate_clustering_performance

report = evaluate_clustering_performance(data, cluster_col="leiden_cluster", ground_truth_col="cell_type")
# n_clusters, cluster_sizes, silhouette_analysis and, with a ground-truth column,
# supervised_metrics: adjusted_rand_index, normalized_mutual_info, homogeneity_score
```

## Resolution Optimization

### Finding the Optimal Resolution

```python
import matplotlib.pyplot as plt
from metainformant.singlecell.analysis.clustering import compute_cluster_silhouette, leiden_clustering

resolutions = [0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0]
n_clusters = []
silhouette_avg = []

for resolution in resolutions:
    clustered = leiden_clustering(data, resolution=resolution)   # a copy; `data` is unchanged
    n_clusters.append(clustered.uns["leiden_clustering"]["n_clusters"])
    scores = compute_cluster_silhouette(clustered, cluster_col="leiden_cluster")
    silhouette_avg.append(scores.get("overall_silhouette_score"))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
ax1.plot(resolutions, n_clusters, "bo-")
ax1.set_xlabel("Resolution")
ax1.set_ylabel("Number of Clusters")
ax2.plot(resolutions, silhouette_avg, "ro-")
ax2.set_xlabel("Resolution")
ax2.set_ylabel("Silhouette Score")
plt.tight_layout()
plt.show()
```

### Multi-resolution clustering

```python
from metainformant.singlecell.analysis.clustering import leiden_clustering

for resolution in [0.2, 0.5, 1.0]:
    clustered = leiden_clustering(data, resolution=resolution)
    data.obs[f"leiden_{resolution}"] = clustered.obs["leiden_cluster"].to_numpy()

import pandas as pd
print(pd.crosstab(data.obs["leiden_0.2"], data.obs["leiden_0.5"]))
```

### Subclustering

`SingleCellData` does not support slicing, so build the subset object explicitly and cluster it:

```python
from metainformant.singlecell.analysis.clustering import leiden_clustering
from metainformant.singlecell.data.preprocessing import SingleCellData

# 1. Extract the cells of the cluster of interest
mask = (data.obs["leiden_cluster"] == 3).to_numpy()
subset = SingleCellData(X=data.X[mask], obs=data.obs.loc[mask].copy(), var=data.var.copy())

# 2. Cluster the subset on its own
subset = leiden_clustering(subset, resolution=0.8, n_neighbors=10)

# 3. Map the sub-labels back onto the full object
data.obs["subcluster"] = "other"
data.obs.loc[mask, "subcluster"] = ["3." + str(label) for label in subset.obs["leiden_cluster"]]
```

## Cluster Annotation

### Manual Annotation

```python
# Map cluster labels to cell types
cluster_annotations = {0: "T cells", 1: "B cells", 2: "NK cells", 3: "Monocytes", 4: "Dendritic cells"}
data.obs["cell_type"] = data.obs["leiden_cluster"].map(cluster_annotations).fillna("Unknown")
```

### Marker-Based Annotation

`metainformant.singlecell.celltyping.annotation.annotate_by_markers` scores every cell against marker gene sets, and `find_marker_genes` or `find_markers` (above) give the genes to compare against a reference list. See [Cell Type Annotation](./celltyping.md).

## Clustering Visualization

Plot the clusters on an embedding with the plotting functions (see [Visualization](./visualization.md)):

```python
from metainformant.singlecell.visualization.visualization import plot_cluster_comparison, plot_umap

fig = plot_umap(data, color="leiden_cluster")            # needs the UMAP columns written by compute_umap
fig = plot_cluster_comparison(data, cluster_cols=["leiden_cluster", "kmeans_cluster"])
```

## Performance Considerations

### Large Datasets

`leiden_clustering` and `louvain_clustering` build a dense cells × cells distance matrix with `scipy.spatial.distance.cdist`, which needs `n² × 8` bytes (100 000 cells need about 80 GB). For large datasets:

- cluster a subsample or one cell population at a time (build a smaller `SingleCellData` as in [Subclustering](#subclustering)), or
- use `kmeans_clustering`, which works directly on `data.X`, or `hierarchical_clustering`, which fits a 1000-cell subset above that size.

Every clustering function densifies `data.X` and uses all of its genes, so a `SingleCellData` restricted to the genes you care about (for example the ones flagged by `select_hvgs`) is cheaper to cluster than the full matrix.

## Common Issues and Solutions

### Too Many or Too Few Clusters

```python
n = data.obs["leiden_cluster"].nunique()
if n > 20:
    print("Many clusters. Consider a lower resolution, stricter filtering or batch correction.")
if n < 5:
    print("Few clusters. Consider a higher resolution or more neighbors in the graph.")
```

### Cluster Stability

`random_state` of `leiden_clustering` seeds the optimiser (the kNN graph itself is deterministic), so the agreement of runs with different seeds measures how stable the partition is:

```python
from sklearn.metrics import adjusted_rand_score
from metainformant.singlecell.analysis.clustering import leiden_clustering

runs = [
    leiden_clustering(data, resolution=0.5, random_state=seed).obs["leiden_cluster"].to_numpy()
    for seed in range(5)
]
scores = [adjusted_rand_score(runs[i], runs[j]) for i in range(5) for j in range(i + 1, 5)]
print(f"Mean pairwise adjusted Rand index: {sum(scores) / len(scores):.3f}")
```

## Testing

Clustering is covered by `tests/singlecell/test_singlecell_clustering_trajectory.py`, `test_singlecell_clustering_diffusion_contracts.py`, `test_singlecell_leiden_fallback.py` (no silent k-means substitution) and `test_singlecell_randomness.py` (seeded subsampling):

```bash
uv run pytest tests/singlecell/test_singlecell_clustering_trajectory.py tests/singlecell/test_singlecell_leiden_fallback.py -v
```

Tests that need the Leiden or Louvain stack are skipped when it is not installed.

## Related Documentation

- [Dimensionality Reduction](./dimensionality.md): For preprocessing before clustering
- [Visualization](./visualization.md): For plotting clusters and markers
- [Trajectory Analysis](./trajectory.md): For developmental trajectory analysis
- [Integration](./integration.md): For batch correction before clustering
