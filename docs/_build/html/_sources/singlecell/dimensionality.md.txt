# Single-Cell Dimensionality Reduction

The dimensionality module provides tools for reducing the complexity of high-dimensional single-cell expression data while preserving biological signal.

## Highly Variable Genes (HVGs)

### select_hvgs()

Identify genes with high biological variability across cells:

```python
from metainformant.singlecell.dimensionality import select_hvgs

# Seurat-style HVG selection (default)
data = select_hvgs(
    data,
    n_top_genes=2000,
    method='seurat',
    min_mean=0.01,      # Minimum mean expression
    max_mean=8.0,       # Maximum mean expression
    min_disp=0.5        # Minimum dispersion
)

# Variance-based selection
data = select_hvgs(data, n_top_genes=2000, method='variance')

# Cell Ranger-style selection
data = select_hvgs(data, n_top_genes=2000, method='cellranger')
```

**Methods:**
- `'seurat'`: Fits variance-mean relationship, selects genes with high residuals
- `'variance'`: Selects genes with highest variance
- `'cellranger'`: Uses standardized variance similar to 10X Cell Ranger

**Parameters:**
- `n_top_genes`: Number of HVGs to select
- `method`: HVG selection method
- `min_mean`, `max_mean`: Mean expression thresholds
- `min_disp`: Minimum dispersion threshold

**Adds to data:**
- `data.var['highly_variable']`: Boolean mask for selected genes
- `data.var['means']`: Mean expression per gene
- `data.var['dispersions']`: Dispersion values
- `data.var['dispersions_norm']`: Normalized dispersions

## Principal Component Analysis (PCA)

### compute_pca()

Perform PCA on expression data:

```python
from metainformant.singlecell.dimensionality import compute_pca

# Standard PCA
data = compute_pca(
    data,
    n_components=50,           # Number of components
    use_highly_variable=True,  # Use only HVGs
    random_state=42
)

# Access results
print(f"PC loadings shape: {data.varm['PCs'].shape}")
print(f"Explained variance: {data.uns['pca']['variance_ratio'][:5]}")
```

**Parameters:**
- `n_components`: Number of principal components to compute
- `use_highly_variable`: Whether to use only highly variable genes
- `random_state`: Random seed for reproducibility

**Adds to data:**
- `data.obsm['X_pca']`: Cell coordinates in PC space
- `data.varm['PCs']`: Gene loadings for each PC
- `data.uns['pca']['variance']`: Explained variance per PC
- `data.uns['pca']['variance_ratio']`: Explained variance ratio per PC

### Choosing Number of Components

```python
import matplotlib.pyplot as plt

# Plot explained variance
variance_ratio = data.uns['pca']['variance_ratio']
cumsum_var = np.cumsum(variance_ratio)

plt.figure(figsize=(10, 4))
plt.subplot(1, 2, 1)
plt.plot(range(1, len(variance_ratio) + 1), variance_ratio, 'bo-')
plt.xlabel('Principal Component')
plt.ylabel('Explained Variance Ratio')
plt.title('Scree Plot')

plt.subplot(1, 2, 2)
plt.plot(range(1, len(cumsum_var) + 1), cumsum_var, 'ro-')
plt.xlabel('Principal Component')
plt.ylabel('Cumulative Explained Variance')
plt.title('Cumulative Variance')
plt.axhline(y=0.8, color='k', linestyle='--', alpha=0.5)
plt.show()
```

## Neighbor Graph Construction

### compute_neighbors()

Build k-nearest neighbor graph for downstream analysis:

```python
from metainformant.singlecell.dimensionality import compute_neighbors

# Compute neighbor graph
data = compute_neighbors(
    data,
    n_neighbors=15,        # Number of neighbors
    n_pcs=40,             # Number of PCs to use
    method='umap',        # Distance metric
    random_state=42
)
```

**Parameters:**
- `n_neighbors`: Number of nearest neighbors
- `n_pcs`: Number of principal components to use
- `method`: Algorithm ('umap' or 'gauss')
- `metric`: Distance metric ('euclidean', 'cosine', etc.)

**Adds to data:**
- `data.obsp['distances']`: Distance matrix (sparse)
- `data.obsp['connectivities']`: Connectivity matrix (sparse)
- `data.uns['neighbors']`: Neighbor graph parameters

## UMAP

### compute_umap()

Uniform Manifold Approximation and Projection for non-linear dimensionality reduction:

```python
from metainformant.singlecell.dimensionality import compute_umap

# Standard UMAP
data = compute_umap(
    data,
    n_components=2,        # Usually 2 or 3 for visualization
    n_neighbors=15,        # Local neighborhood size
    min_dist=0.5,          # Minimum distance in embedding
    spread=1.0,            # Scale of embedded points
    random_state=42
)

# Access UMAP coordinates
umap_coords = data.obsm['X_umap']  # Shape: (n_cells, n_components)
```

**Parameters:**
- `n_components`: Dimensionality of embedding (2 or 3)
- `n_neighbors`: Size of local neighborhood
- `min_dist`: Controls how tightly points are packed
- `spread`: Effective scale of embedded points
- `metric`: Distance metric ('euclidean', 'cosine', 'correlation')

**Note:** Requires the optional `umap-learn` package:
```bash
uv pip install umap-learn
```

## t-SNE

### compute_tsne()

t-distributed Stochastic Neighbor Embedding:

```python
from metainformant.singlecell.dimensionality import compute_tsne

# Standard t-SNE
data = compute_tsne(
    data,
    n_components=2,        # Usually 2 for visualization
    perplexity=30,         # Local neighborhood size
    max_iter=1000,         # Maximum iterations
    random_state=42
)

# Access t-SNE coordinates
tsne_coords = data.obsm['X_tsne']  # Shape: (n_cells, 2)
```

**Parameters:**
- `n_components`: Dimensionality of embedding (2 or 3)
- `perplexity`: Balance between local and global structure
- `max_iter`: Maximum number of optimization iterations
- `learning_rate`: Step size for gradient descent

## Diffusion Maps

### compute_diffusion_map()

Diffusion maps for trajectory analysis and pseudotime inference:

```python
from metainformant.singlecell.dimensionality import compute_diffusion_map

# Compute diffusion components
data = compute_diffusion_map(
    data,
    n_components=15,       # Number of diffusion components
    alpha=1.0,             # Normalization parameter
    random_state=42
)

# Access diffusion coordinates
dc_coords = data.obsm['X_diffusion']  # Diffusion coordinates
```

**Parameters:**
- `n_components`: Number of diffusion components to compute
- `alpha`: Determines the influence of data density (0=geometric, 1=diffusion)
- `n_neighbors`: Number of neighbors for graph construction

**Adds to data:**
- `data.obsm['X_diffusion']`: Diffusion map coordinates
- `data.uns['diffusion']`: Eigenvalues and other metadata

## Visualization Integration

### Plotting Embeddings

```python
import matplotlib.pyplot as plt

# Plot UMAP colored by cell type
if 'cell_type' in data.obs.columns:
    plt.figure(figsize=(10, 8))
    for cell_type in data.obs['cell_type'].unique():
        mask = data.obs['cell_type'] == cell_type
        plt.scatter(
            data.obsm['X_umap'][mask, 0],
            data.obsm['X_umap'][mask, 1],
            label=cell_type,
            alpha=0.7
        )
    plt.xlabel('UMAP 1')
    plt.ylabel('UMAP 2')
    plt.legend()
    plt.title('UMAP by Cell Type')
    plt.show()

# Plot t-SNE colored by total counts
plt.figure(figsize=(8, 6))
scatter = plt.scatter(
    data.obsm['X_tsne'][:, 0],
    data.obsm['X_tsne'][:, 1],
    c=data.obs['total_counts'],
    cmap='viridis',
    alpha=0.7
)
plt.colorbar(scatter, label='Total Counts')
plt.xlabel('t-SNE 1')
plt.ylabel('t-SNE 2')
plt.title('t-SNE by Total Counts')
plt.show()
```

## Parameter Guidelines

### HVG Selection
- **2000-3000 genes**: Good balance for most datasets
- **More genes (5000+)**: For very heterogeneous datasets
- **Fewer genes (1000)**: For homogeneous populations

### PCA
- **50 components**: Standard for most analyses
- **Use scree plot**: To determine optimal number
- **80% variance**: Common cutoff for cumulative variance

### UMAP Parameters
- **n_neighbors**: 
  - 15-50 for most datasets
  - Higher for smoother embeddings
  - Lower for more local structure
- **min_dist**: 
  - 0.1-0.5 for tight clusters
  - 0.5-1.0 for spread-out visualization

### t-SNE Parameters
- **perplexity**: 
  - 30-50 for most datasets
  - Scale with dataset size: ~sqrt(n_cells)
  - Multiple perplexities to check robustness

## Memory and Performance

### Large Dataset Considerations

```python
# For datasets > 100k cells
data = select_hvgs(data, n_top_genes=1000)  # Fewer HVGs
data = compute_pca(data, n_components=30)   # Fewer PCs

# Use incremental PCA for very large datasets
from sklearn.decomposition import IncrementalPCA
# (Custom implementation may be needed)
```

### GPU Acceleration

For datasets > 1M cells, consider:
- Rapids cuML for GPU-accelerated UMAP
- scanpy with GPU support
- Specialized single-cell GPU tools

## Common Workflows

### Standard Dimensionality Reduction Pipeline

```python
# Complete dimensionality reduction workflow
from metainformant.singlecell.dimensionality import *

# 1. Select highly variable genes
data = select_hvgs(data, n_top_genes=2000, method='seurat')
print(f"Selected {data.var['highly_variable'].sum()} HVGs")

# 2. Principal component analysis
data = compute_pca(data, n_components=50, use_highly_variable=True)
print(f"Explained variance (first 10 PCs): "
      f"{data.uns['pca']['variance_ratio'][:10].sum():.3f}")

# 3. Neighbor graph for clustering
data = compute_neighbors(data, n_neighbors=15, n_pcs=40)

# 4. UMAP embedding
data = compute_umap(data, n_components=2, min_dist=0.5)

# 5. t-SNE for comparison
data = compute_tsne(data, n_components=2, perplexity=30)

print("Dimensionality reduction complete!")
print(f"Available embeddings: {list(data.obsm.keys())}")
```

### Trajectory Analysis Preparation

```python
# For trajectory/pseudotime analysis
data = compute_pca(data, n_components=50)
data = compute_diffusion_map(data, n_components=15)
data = compute_neighbors(data, n_neighbors=30)  # More neighbors for trajectories
```

## Quality Assessment

### Embedding Quality Metrics

```python
# Check for overclustering in UMAP
import numpy as np
from scipy.spatial.distance import pdist, squareform

# Local neighborhood preservation
def local_preservation(X_orig, X_embed, k=15):
    """Measure how well local neighborhoods are preserved."""
    # Implementation for neighborhood preservation metric
    pass

# Global structure preservation
def global_preservation(X_orig, X_embed):
    """Measure how well global structure is preserved."""
    dist_orig = pdist(X_orig[:1000])  # Sample for speed
    dist_embed = pdist(X_embed[:1000])
    correlation = np.corrcoef(dist_orig, dist_embed)[0, 1]
    return correlation
```

## Error Handling

### Common Issues

```python
try:
    data = compute_umap(data)
except ImportError:
    print("umap-learn not installed. Install with: uv pip install umap-learn")
    # Fall back to t-SNE
    data = compute_tsne(data)

# Handle cases with too few genes
try:
    data = select_hvgs(data, n_top_genes=2000)
except ValueError as e:
    print(f"Warning: {e}")
    # Reduce number of requested HVGs
    data = select_hvgs(data, n_top_genes=min(1000, data.n_vars // 2))
```

## Integration with Clustering

The dimensionality reduction results are used by clustering algorithms:

```python
from metainformant.singlecell.clustering import leiden_clustering

# Use neighbor graph for clustering
data = leiden_clustering(data, resolution=0.5)

# Or use PCA coordinates directly
data = leiden_clustering(data, use_rep='X_pca', resolution=0.5)
```

## Testing

Dimensionality reduction tests are available in `tests/test_singlecell_dimensionality.py`:

```bash
# Run dimensionality reduction tests
uv run pytest tests/test_singlecell_dimensionality.py -v
```

## Related Documentation

- [Preprocessing](./preprocessing.md): For data preparation before dimensionality reduction
- [Clustering](./clustering.md): For using dimensionality reduction results in clustering
- [Visualization](./visualization.md): For plotting embeddings
- [Trajectory](./trajectory.md): For pseudotime analysis using diffusion maps
