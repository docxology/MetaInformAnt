# Single-Cell Trajectory Analysis

The trajectory module provides tools for inferring developmental trajectories and pseudotime ordering of cells.

## Pseudotime Inference

### compute_pseudotime()

Infer pseudotime using diffusion maps:

```python
from metainformant.singlecell.trajectory import compute_pseudotime

# Compute pseudotime from diffusion components
data = compute_pseudotime(
    data,
    root_cells=None,           # Auto-detect or specify root cells
    n_dcs=10,                  # Number of diffusion components
    use_rep='X_diffusion'      # Use diffusion map coordinates
)

# Access pseudotime values
pseudotime = data.obs['pseudotime']
print(f"Pseudotime range: {pseudotime.min():.3f} - {pseudotime.max():.3f}")
```

**Parameters:**
- `root_cells`: Indices of root cells (if None, auto-detected)
- `n_dcs`: Number of diffusion components to use
- `use_rep`: Representation to use for trajectory inference

## Trajectory Analysis

### trajectory_analysis()

Perform comprehensive trajectory analysis:

```python
from metainformant.singlecell.trajectory import trajectory_analysis

# Run full trajectory analysis
trajectory_results = trajectory_analysis(
    data,
    groupby='leiden',          # Cluster column for trajectory
    method='mst',              # Minimum spanning tree
    root_cluster='0'           # Starting cluster
)

# Results include:
# - Trajectory graph
# - Branch assignments
# - Milestone connections
```

## Gene Expression Trends

### compute_gene_trends()

Identify genes with significant expression changes along trajectories:

```python
from metainformant.singlecell.trajectory import compute_gene_trends

# Find genes that change along pseudotime
trends = compute_gene_trends(
    data,
    pseudotime_col='pseudotime',
    n_genes=500,               # Top variable genes
    method='gam'               # Generalized additive models
)

# Access trending genes
print(f"Found {len(trends)} genes with significant trends")
```

## Lineage Analysis

### identify_lineages()

Identify distinct developmental lineages:

```python
from metainformant.singlecell.trajectory import identify_lineages

# Find developmental lineages
lineages = identify_lineages(
    data,
    trajectory_graph,          # From trajectory_analysis
    n_lineages=3              # Expected number of lineages
)

# Add lineage assignments to data
data.obs['lineage'] = lineages
```

## Branch Point Analysis

### find_branch_genes()

Identify genes associated with trajectory branch points:

```python
from metainformant.singlecell.trajectory import find_branch_genes

# Find branch-associated genes
branch_genes = find_branch_genes(
    data,
    branch_point_cells,        # Cells at branch points
    lineage_assignments,       # Cell lineage assignments
    min_fold_change=2.0
)
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
    c=data.obs['pseudotime'],
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
    plt.scatter(data.obs['pseudotime'], gene_expr, alpha=0.5)
    plt.xlabel('Pseudotime')
    plt.ylabel(f'{gene_of_interest} Expression')
    plt.title('Gene Expression vs Pseudotime')

plt.tight_layout()
plt.show()
```

## Integration with Other Modules

### Preprocessing for Trajectory Analysis

```python
# Recommended preprocessing for trajectory analysis
from metainformant.singlecell.preprocessing import *
from metainformant.singlecell.dimensionality import *

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
data = compute_pseudotime(data)
```

## Advanced Trajectory Methods

### Multiple Trajectory Algorithms

The module supports multiple trajectory inference methods:

```python
# Different trajectory methods
methods = ['diffusion', 'mst', 'force_directed']

for method in methods:
    if method == 'diffusion':
        data = compute_pseudotime(data, method='diffusion')
    elif method == 'mst':
        results = trajectory_analysis(data, method='mst')
    elif method == 'force_directed':
        results = trajectory_analysis(data, method='force_directed')
```

## Common Trajectory Patterns

### Linear Trajectories

```python
# For simple linear differentiation
data = compute_pseudotime(data, n_dcs=3)  # Fewer components for linear
```

### Branching Trajectories

```python
# For complex branching processes
data = compute_pseudotime(data, n_dcs=15)  # More components for complexity
lineages = identify_lineages(data, n_lineages=3)
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
    plt.hist(data.obs['pseudotime'], bins=50, alpha=0.7)
    plt.xlabel('Pseudotime')
    plt.ylabel('Number of Cells')
    plt.title('Pseudotime Distribution')
    
    # 2. Pseudotime vs clusters
    plt.subplot(1, 3, 2)
    for cluster in data.obs['leiden'].unique():
        mask = data.obs['leiden'] == cluster
        plt.hist(data.obs[mask]['pseudotime'], alpha=0.5, label=f'Cluster {cluster}')
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
# Compare different pseudotime methods
methods = ['diffusion', 'force_directed']
correlations = {}

for method in methods:
    # Compute pseudotime with different methods
    data_temp = compute_pseudotime(data.copy(), method=method)
    correlations[method] = data_temp.obs['pseudotime']

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
uv run pytest tests/test_singlecell_*.py -v
```

## Related Documentation

- [Dimensionality Reduction](./dimensionality.md): For diffusion maps and PCA
- [Clustering](./clustering.md): For identifying cell populations  
- [Visualization](./visualization.md): For trajectory plotting
- [Preprocessing](./preprocessing.md): For data preparation
