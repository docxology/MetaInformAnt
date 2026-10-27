# Single-Cell Trajectory Analysis

The trajectory module provides tools for inferring developmental trajectories and pseudotime ordering of cells.

## Pseudotime Inference

### compute_diffusion_pseudotime()

Infer pseudotime using diffusion maps:

```python-snippet
from metainformant.singlecell.analysis.trajectory import compute_diffusion_pseudotime

# Compute diffusion pseudotime
data = compute_diffusion_pseudotime(
    data,
    root_cell=None,            # Index of the root cell (auto-detected if None)
    n_components=10            # Number of diffusion components
)

# Access pseudotime values
pseudotime = data.obs['dpt_pseudotime']
print(f"Pseudotime range: {pseudotime.min():.3f} - {pseudotime.max():.3f}")
```

**Parameters:**
- `root_cell`: Index of the root cell (if None, auto-detected)
- `n_components`: Number of diffusion components to use

## Trajectory Analysis

### slingshot_trajectory()

Cluster-based trajectory inference from a root cluster to one or more terminal clusters:

```python-snippet
from metainformant.singlecell.analysis.trajectory import slingshot_trajectory

# Run cluster-based trajectory analysis
data = slingshot_trajectory(
    data,
    start_cluster=0,           # Root cluster (same type as the labels in obs)
    end_clusters=[2, 4]        # Terminal clusters, one branch per end cluster
)

# Results include:
# - data.obs['slingshot_pseudotime']: pseudotime per cell
# - data.uns['slingshot']: start/end clusters, branches, cluster pseudotime and centroids
```

The cluster column is detected automatically (the first `obs` column whose name contains `cluster`, `leiden` or `louvain`). For cluster-to-cluster connectivity, see `paga_trajectory(data, groups='leiden_cluster')` in the same module.

## Lineage Analysis

### find_trajectory_branches()

Identify branches along a computed pseudotime axis:

```python-snippet
from metainformant.singlecell.analysis.trajectory import find_trajectory_branches

# Find branches in pseudotime
branches = find_trajectory_branches(
    data,
    pseudotime_col='dpt_pseudotime',
    min_branch_size=10         # Minimum number of cells per branch
)

# Add branch assignments to data (-1 = cell not assigned to a branch)
data.obs['branch'] = branches['branch_labels']
print(f"Found {branches['n_branches']} branches, sizes: {branches['branch_sizes']}")
```

## Visualization

### Trajectory Visualization

```python
import matplotlib.pyplot as plt

# Plot pseudotime on UMAP
plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
scatter = plt.scatter(
    data.obsm['X_umap'][:, 0],
    data.obsm['X_umap'][:, 1],
    c=data.obs['dpt_pseudotime'],
    cmap='viridis',
    alpha=0.7
)
plt.colorbar(scatter, label='Pseudotime')
plt.xlabel('UMAP 1')
plt.ylabel('UMAP 2')
plt.title('Pseudotime on UMAP')

# Plot gene expression along pseudotime
plt.subplot(1, 2, 2)
gene_of_interest = 'GENE1'  # Replace with actual gene
if gene_of_interest in data.var.index:
    gene_expr = data.X[:, data.var.index == gene_of_interest].flatten()
    plt.scatter(data.obs['dpt_pseudotime'], gene_expr, alpha=0.5)
    plt.xlabel('Pseudotime')
    plt.ylabel(f'{gene_of_interest} Expression')
    plt.title('Gene Expression vs Pseudotime')

plt.tight_layout()
plt.show()
```

## Integration with Other Modules

### Preprocessing for Trajectory Analysis

```python-snippet
# Recommended preprocessing for trajectory analysis
from metainformant.singlecell.data.preprocessing import (
    calculate_qc_metrics,
    filter_cells,
    log_transform,
    normalize_counts,
)
from metainformant.singlecell.analysis.pca_methods import compute_pca, select_hvgs
from metainformant.singlecell.analysis.nonlinear_methods import compute_diffusion_map, compute_neighbors
from metainformant.singlecell.analysis.trajectory import compute_diffusion_pseudotime

# 1. Standard preprocessing
data = calculate_qc_metrics(data)
data = filter_cells(data, min_genes=200, max_pct_mt=20)
data = normalize_counts(data, method='total_count')
data = log_transform(data)

# 2. Feature selection and dimensionality reduction
data = select_hvgs(data, n_top_genes=2000)
data = compute_pca(data, n_components=50)

# 3. Diffusion maps for trajectory analysis
data = compute_diffusion_map(data, n_components=15)
data = compute_neighbors(data, n_neighbors=30)  # More neighbors for trajectories

# 4. Trajectory inference
data = compute_diffusion_pseudotime(data)
```

## Advanced Trajectory Methods

### Multiple Trajectory Algorithms

The module supports multiple trajectory inference methods:

```python
from metainformant.singlecell.analysis.trajectory import (
    compute_diffusion_pseudotime,
    compute_pseudotime_from_dimensionality_reduction,
    paga_trajectory,
    slingshot_trajectory,
)

# Different trajectory methods (each adds its own pseudotime column to data.obs)
data = compute_diffusion_pseudotime(data, root_cell=0)                   # dpt_pseudotime
data = paga_trajectory(data, groups='leiden_cluster')                    # paga_pseudotime
data = slingshot_trajectory(data, start_cluster=0, end_clusters=[2])     # slingshot_pseudotime
data = compute_pseudotime_from_dimensionality_reduction(
    data, dim_red_cols=['UMAP1', 'UMAP2'], root_cell=0                   # embedding_pseudotime
)
```

## Common Trajectory Patterns

### Linear Trajectories

```python
# For simple linear differentiation
data = compute_diffusion_pseudotime(data, n_components=3)  # Fewer components for linear
```

### Branching Trajectories

```python
# For complex branching processes
data = compute_diffusion_pseudotime(data, n_components=15)  # More components for complexity
branches = find_trajectory_branches(data, pseudotime_col='dpt_pseudotime')
```

### Cyclic Trajectories

```python
# For cell cycle or other periodic processes
# Use circular embedding or specialized methods
```

## Quality Control

### Trajectory Quality Assessment

```python
def assess_trajectory_quality(data):
    """Assess quality of inferred trajectory."""
    
    # 1. Pseudotime distribution
    plt.figure(figsize=(15, 4))
    
    plt.subplot(1, 3, 1)
    plt.hist(data.obs['dpt_pseudotime'], bins=50, alpha=0.7)
    plt.xlabel('Pseudotime')
    plt.ylabel('Number of Cells')
    plt.title('Pseudotime Distribution')
    
    # 2. Pseudotime vs clusters
    plt.subplot(1, 3, 2)
    for cluster in data.obs['leiden'].unique():
        mask = data.obs['leiden'] == cluster
        plt.hist(data.obs[mask]['dpt_pseudotime'], alpha=0.5, label=f'Cluster {cluster}')
    plt.xlabel('Pseudotime')
    plt.ylabel('Density')
    plt.legend()
    plt.title('Pseudotime by Cluster')
    
    # 3. Known marker genes
    plt.subplot(1, 3, 3)
    # Plot known early vs late markers if available
    plt.xlabel('Pseudotime')
    plt.ylabel('Expression')
    plt.title('Marker Gene Expression')
    
    plt.tight_layout()
    plt.show()

assess_trajectory_quality(data)
```

## Applications

### Developmental Biology

- **Embryonic development**: Track cell fate decisions
- **Organogenesis**: Understand tissue formation
- **Stem cell differentiation**: Map lineage relationships

### Disease Studies

- **Cancer progression**: Tumor evolution trajectories  
- **Drug response**: Time-course treatment effects
- **Immune responses**: T cell activation dynamics

### Method Comparison

```python
from metainformant.singlecell.analysis.trajectory import (
    compute_diffusion_pseudotime,
    compute_pseudotime_from_dimensionality_reduction,
)

# Compare different pseudotime methods
methods = ['diffusion', 'embedding']
correlations = {}

for method in methods:
    # Compute pseudotime with different methods
    if method == 'diffusion':
        data_temp = compute_diffusion_pseudotime(data.copy(), root_cell=0)
        correlations[method] = data_temp.obs['dpt_pseudotime']
    else:
        data_temp = compute_pseudotime_from_dimensionality_reduction(
            data.copy(), dim_red_cols=['UMAP1', 'UMAP2'], root_cell=0
        )
        correlations[method] = data_temp.obs['embedding_pseudotime']

# Calculate correlation between methods
import numpy as np
correlation_matrix = np.corrcoef([
    correlations[method] for method in methods
])
print("Method correlation matrix:")
print(correlation_matrix)
```

## Best Practices

1. **Quality control**: Ensure good quality cells before trajectory analysis
2. **Sufficient coverage**: Need adequate cells along trajectory
3. **Appropriate markers**: Validate with known developmental markers
4. **Multiple methods**: Compare results from different algorithms
5. **Biological validation**: Confirm with experimental validation

## Limitations

- **Requires prior knowledge**: Need to know expected trajectory structure
- **Assumes continuity**: May not work for discrete state transitions
- **Sample coverage**: Sparse sampling can miss intermediate states
- **Technical noise**: Sensitive to technical artifacts

## Testing

Trajectory analysis functionality is tested as part of the single-cell test suite:

```bash
# Run single-cell tests (includes trajectory examples)
uv run pytest tests/singlecell/test_singlecell_*.py -v
```

## Related Documentation

- [Dimensionality Reduction](./dimensionality.md): For diffusion maps and PCA
- [Clustering](./clustering.md): For identifying cell populations  
- [Visualization](./visualization.md): For trajectory plotting
- [Preprocessing](./preprocessing.md): For data preparation
