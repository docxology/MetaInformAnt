# Single-Cell Clustering

The clustering module provides various algorithms for identifying cell populations and subtypes in single-cell data.

## Graph-Based Clustering

### leiden_clustering()

Leiden algorithm for community detection in cell similarity graphs:

```python
from metainformant.singlecell.clustering import leiden_clustering

# Standard Leiden clustering
data = leiden_clustering(
    data,
    resolution=0.5,        # Higher = more clusters
    use_rep='neighbors',   # Use precomputed neighbor graph
    random_state=42
)

# Access cluster assignments
clusters = data.obs['leiden']  # Cluster labels
n_clusters = len(clusters.unique())
print(f"Found {n_clusters} clusters")
```

### louvain_clustering()

Louvain algorithm (predecessor to Leiden):

```python
from metainformant.singlecell.clustering import louvain_clustering

data = louvain_clustering(
    data,
    resolution=0.5,
    use_rep='neighbors',
    random_state=42
)
```

**Parameters:**
- `resolution`: Higher values yield more clusters
- `use_rep`: Data representation to use ('neighbors', 'X_pca', 'X_umap')
- `random_state`: Random seed for reproducibility

## Centroid-Based Clustering

### kmeans_clustering()

K-means clustering for predetermined number of clusters:

```python
from metainformant.singlecell.clustering import kmeans_clustering

# K-means with known number of clusters
data = kmeans_clustering(
    data,
    n_clusters=8,          # Number of clusters
    use_rep='X_pca',       # Use PCA coordinates
    random_state=42,
    n_init=20              # Multiple initializations
)
```

**Parameters:**
- `n_clusters`: Number of clusters to find
- `use_rep`: Data representation ('X_pca' recommended)
- `n_init`: Number of random initializations
- `max_iter`: Maximum iterations

## Hierarchical Clustering

### hierarchical_clustering()

Agglomerative clustering with various linkage methods:

```python
from metainformant.singlecell.clustering import hierarchical_clustering

# Hierarchical clustering
data = hierarchical_clustering(
    data,
    n_clusters=10,         # Number of clusters
    linkage='ward',        # Linkage method
    use_rep='X_pca'        # Use PCA coordinates
)
```

**Linkage methods:**
- `'ward'`: Minimizes within-cluster variance
- `'complete'`: Maximum linkage
- `'average'`: Average linkage
- `'single'`: Single linkage

## Marker Gene Analysis

### find_markers()

Identify genes that characterize each cluster:

```python
from metainformant.singlecell.clustering import find_markers

# Find marker genes for each cluster
markers = find_markers(
    data,
    groupby='leiden',      # Column with cluster labels
    method='wilcoxon',     # Statistical test
    min_fold_change=1.5,   # Minimum fold change
    min_pct=0.1           # Minimum detection percentage
)

# Access results
print(markers.head())
# Columns: gene, cluster, avg_log2FC, pct.1, pct.2, p_val, p_val_adj
```

**Methods:**
- `'wilcoxon'`: Wilcoxon rank-sum test (default)
- `'ttest'`: Student's t-test
- `'logistic'`: Logistic regression

### cluster_composition()

Analyze cluster composition and relationships:

```python
from metainformant.singlecell.clustering import cluster_composition

# Get cluster statistics
composition = cluster_composition(data, groupby='leiden')
print(composition)
# Columns: cluster, n_cells, pct_total, top_markers
```

## Clustering Validation

### silhouette_scores()

Calculate silhouette scores to assess clustering quality:

```python
from metainformant.singlecell.clustering import silhouette_scores

# Calculate silhouette scores
scores = silhouette_scores(
    data,
    cluster_column='leiden',
    use_rep='X_pca'
)

print(f"Average silhouette score: {scores.mean():.3f}")
```

## Resolution Optimization

### Finding Optimal Resolution

```python
import matplotlib.pyplot as plt

# Test multiple resolutions
resolutions = [0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0]
n_clusters = []
silhouette_avg = []

for res in resolutions:
    # Cluster at this resolution
    data_temp = leiden_clustering(data.copy(), resolution=res)
    
    # Count clusters
    n_clust = len(data_temp.obs['leiden'].unique())
    n_clusters.append(n_clust)
    
    # Calculate silhouette score
    sil_score = silhouette_scores(data_temp, 'leiden', 'X_pca').mean()
    silhouette_avg.append(sil_score)

# Plot results
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

ax1.plot(resolutions, n_clusters, 'bo-')
ax1.set_xlabel('Resolution')
ax1.set_ylabel('Number of Clusters')
ax1.set_title('Clusters vs Resolution')

ax2.plot(resolutions, silhouette_avg, 'ro-')
ax2.set_xlabel('Resolution')
ax2.set_ylabel('Silhouette Score')
ax2.set_title('Clustering Quality vs Resolution')

plt.tight_layout()
plt.show()
```

## Advanced Clustering Workflows

### Multi-Resolution Clustering

```python
# Cluster at multiple resolutions for hierarchical analysis
resolutions = [0.2, 0.5, 1.0]
for res in resolutions:
    data = leiden_clustering(data, resolution=res, key_added=f'leiden_{res}')

# Compare clusterings
import pandas as pd
cluster_comparison = pd.crosstab(
    data.obs['leiden_0.2'], 
    data.obs['leiden_0.5']
)
print(cluster_comparison)
```

### Subclustering

```python
# Subcluster specific cell populations
# 1. Extract cells from cluster of interest
cluster_of_interest = '3'
mask = data.obs['leiden'] == cluster_of_interest
subdata = data[mask, :].copy()

# 2. Recompute neighbors and cluster
from metainformant.singlecell.dimensionality import compute_neighbors
subdata = compute_neighbors(subdata, n_neighbors=10)
subdata = leiden_clustering(subdata, resolution=0.8)

# 3. Map back to original data
data.obs['subcluster'] = 'Other'
data.obs.loc[mask, 'subcluster'] = [
    f"{cluster_of_interest}.{sub}" 
    for sub in subdata.obs['leiden']
]
```

## Integration with Dimensionality Reduction

### Clustering in Different Spaces

```python
# Compare clustering in different representations
representations = ['neighbors', 'X_pca', 'X_umap']

for rep in representations:
    if rep in data.obsm or rep in data.obsp:
        data = leiden_clustering(
            data, 
            resolution=0.5, 
            use_rep=rep,
            key_added=f'leiden_{rep}'
        )

# Compare results
from metainformant.singlecell.clustering import cluster_composition
for rep in representations:
    key = f'leiden_{rep}'
    if key in data.obs.columns:
        comp = cluster_composition(data, groupby=key)
        print(f"\n{rep}: {len(comp)} clusters")
```

## Cluster Annotation

### Manual Annotation

```python
# Create mapping from cluster numbers to cell types
cluster_annotations = {
    '0': 'T cells',
    '1': 'B cells', 
    '2': 'NK cells',
    '3': 'Monocytes',
    '4': 'Dendritic cells'
}

# Apply annotations
data.obs['cell_type'] = data.obs['leiden'].map(cluster_annotations)
data.obs['cell_type'] = data.obs['cell_type'].fillna('Unknown')
```

### Automated Annotation (Framework)

```python
def annotate_clusters(data, reference_markers):
    """
    Annotate clusters based on marker gene expression.
    
    Args:
        data: SingleCellData object
        reference_markers: Dict mapping cell types to marker genes
    """
    annotations = {}
    
    for cluster in data.obs['leiden'].unique():
        cluster_markers = find_markers(data, groupby='leiden')
        cluster_genes = cluster_markers[
            cluster_markers['cluster'] == cluster
        ]['gene'].head(10).tolist()
        
        # Score overlap with reference markers
        best_match = 'Unknown'
        best_score = 0
        
        for cell_type, markers in reference_markers.items():
            overlap = len(set(cluster_genes) & set(markers))
            if overlap > best_score:
                best_score = overlap
                best_match = cell_type
        
        annotations[cluster] = best_match
    
    return annotations

# Example usage
reference_markers = {
    'T cells': ['CD3D', 'CD3E', 'CD2'],
    'B cells': ['CD19', 'MS4A1', 'CD79A'],
    'NK cells': ['GNLY', 'NKG7', 'KLRD1'],
    'Monocytes': ['CD14', 'LYZ', 'S100A9']
}

annotations = annotate_clusters(data, reference_markers)
data.obs['predicted_cell_type'] = data.obs['leiden'].map(annotations)
```

## Clustering Visualization

### Basic Cluster Plots

```python
import matplotlib.pyplot as plt

# Plot clusters on UMAP
plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
for cluster in data.obs['leiden'].unique():
    mask = data.obs['leiden'] == cluster
    plt.scatter(
        data.obsm['X_umap'][mask, 0],
        data.obsm['X_umap'][mask, 1],
        label=f'Cluster {cluster}',
        alpha=0.7
    )
plt.xlabel('UMAP 1')
plt.ylabel('UMAP 2')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.title('Clusters on UMAP')

plt.subplot(1, 2, 2)
# Plot without legend for cleaner visualization
for i, cluster in enumerate(data.obs['leiden'].unique()):
    mask = data.obs['leiden'] == cluster
    plt.scatter(
        data.obsm['X_umap'][mask, 0],
        data.obsm['X_umap'][mask, 1],
        c=f'C{i}',
        alpha=0.7
    )
plt.xlabel('UMAP 1')
plt.ylabel('UMAP 2')
plt.title('Clusters (No Legend)')

plt.tight_layout()
plt.show()
```

## Performance Considerations

### Large Dataset Clustering

```python
# For datasets > 100k cells
# 1. Use subset for parameter optimization
subset_indices = np.random.choice(data.n_obs, size=10000, replace=False)
data_subset = data[subset_indices, :].copy()

# 2. Optimize resolution on subset
optimal_res = 0.5  # Determined from subset analysis

# 3. Apply to full dataset
data = leiden_clustering(data, resolution=optimal_res)
```

### Memory Efficient Clustering

```python
# Use sparse representations and limit memory usage
import gc

# Clear unnecessary data before clustering
if 'X_pca' in data.obsm:
    # Only keep necessary number of PCs
    data.obsm['X_pca'] = data.obsm['X_pca'][:, :30]

# Force garbage collection
gc.collect()

# Run clustering
data = leiden_clustering(data, resolution=0.5)
```

## Common Issues and Solutions

### Too Many/Few Clusters

```python
# Too many clusters
if len(data.obs['leiden'].unique()) > 20:
    print("Warning: Many clusters detected. Consider:")
    print("- Lower resolution parameter")
    print("- Better preprocessing (more stringent filtering)")
    print("- Batch correction if multiple samples")

# Too few clusters  
if len(data.obs['leiden'].unique()) < 5:
    print("Warning: Few clusters detected. Consider:")
    print("- Higher resolution parameter") 
    print("- More neighbors in graph construction")
    print("- Different clustering algorithm")
```

### Cluster Stability

```python
# Test clustering stability across random seeds
def test_clustering_stability(data, n_tests=10, resolution=0.5):
    """Test how stable clustering is across random seeds."""
    results = []
    
    for seed in range(n_tests):
        data_temp = leiden_clustering(
            data.copy(), 
            resolution=resolution, 
            random_state=seed
        )
        results.append(data_temp.obs['leiden'].values)
    
    # Calculate pairwise adjusted rand index
    from sklearn.metrics import adjusted_rand_score
    
    stability_scores = []
    for i in range(n_tests):
        for j in range(i+1, n_tests):
            ari = adjusted_rand_score(results[i], results[j])
            stability_scores.append(ari)
    
    return np.mean(stability_scores)

stability = test_clustering_stability(data)
print(f"Clustering stability (ARI): {stability:.3f}")
```

## Testing

Comprehensive tests are available in `tests/test_singlecell_clustering.py`:

```bash
# Run clustering tests
uv run pytest tests/test_singlecell_clustering.py -v
```

## Related Documentation

- [Dimensionality Reduction](./dimensionality.md): For preprocessing before clustering
- [Visualization](./visualization.md): For plotting clusters and markers
- [Trajectory Analysis](./trajectory.md): For developmental trajectory analysis
- [Integration](./integration.md): For batch correction before clustering
