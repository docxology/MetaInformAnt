# Single-Cell Visualization

The visualization module provides specialized plotting functions for single-cell RNA-seq data analysis and results.

## Quality Control Plots

### plot_qc_metrics()

Visualize quality control metrics:

```python
from metainformant.singlecell.visualization.visualization import plot_qc_metrics

# Plot QC metrics overview
fig = plot_qc_metrics(
    data,
    metrics=['total_counts', 'n_genes', 'pct_mt', 'pct_ribo'],
    plot_type='violin'     # 'violin', 'box', 'hist'
)
plt.show()
```

The figure contains histograms of total counts, genes detected and mitochondrial percentage, plus a scatter plot of total counts against genes detected.

## Dimensionality Reduction Plots

### plot_umap() / plot_tsne()

Plot cells in reduced dimensional space (coordinates are read from the `obs` columns written by `compute_umap()` / `compute_tsne()`):

```python-snippet
from metainformant.singlecell.visualization.visualization import plot_tsne, plot_umap

# UMAP colored by clusters
fig = plot_umap(
    data,
    color='leiden_cluster',    # Color by cluster
    s=30                       # Point size
)

# t-SNE colored by clusters
fig = plot_tsne(data, color='leiden_cluster')

# UMAP colored by gene expression (add the gene as an obs column first)
data.obs['CD3E'] = data.to_df()['CD3E'].values
fig = plot_umap(
    data,
    color='CD3E',
    cmap='viridis'
)
plt.show()
```

### plot_pca()

Specialized PCA visualization:

```python
from metainformant.singlecell.visualization.visualization import plot_pca

# PCA with loadings
fig = plot_pca(
    data,
    components=[1, 2],     # PC1 vs PC2
    color='leiden',
    show_loadings=True,
    n_loadings=10         # Show top 10 gene loadings
)
plt.show()
```

## Gene Expression Visualization

### plot_marker_expression()

Visualize the expression of individual genes as violin plots, dot plots or heatmaps:

```python-snippet
from metainformant.singlecell.visualization.visualization import plot_marker_expression

# Violin plots of expression among expressing cells
fig = plot_marker_expression(
    data,
    marker_genes=['CD3E', 'CD19', 'CD14'],
    method='violin'
)

# Dot plot showing expression and detection
fig = plot_marker_expression(
    data,
    marker_genes=['CD3E', 'CD19', 'CD14', 'LYZ'],
    method='dotplot'
)

# Heatmap of marker genes across (a subsample of) cells
fig = plot_marker_expression(
    data,
    marker_genes=['CD3E', 'CD19', 'CD14', 'GNLY'],
    method='heatmap',
    cmap='RdBu_r'
)
plt.show()
```

Expression is pooled over all cells; the plots are not split by cluster.

## Clustering Visualization

### plot_cluster_comparison()

Cluster views on the embedding, one panel per clustering:

```python-snippet
from metainformant.singlecell.visualization.visualization import plot_cluster_comparison

# Cluster view on the UMAP embedding
fig = plot_cluster_comparison(
    data,
    cluster_cols=['leiden_cluster'],
    embedding_cols=['UMAP1', 'UMAP2']
)

# Compare clusterings at different resolutions
fig = plot_cluster_comparison(
    data,
    cluster_cols=['leiden_0.2', 'leiden_0.5', 'leiden_1.0']
)
plt.show()
```

## Trajectory Visualization

### plot_trajectory()

Visualize developmental trajectories:

```python
from metainformant.singlecell.visualization.visualization import plot_trajectory

# Pseudotime on embedding
fig = plot_trajectory(
    data,
    basis='umap',
    color='pseudotime',
    trajectory_graph=None,   # Optional trajectory overlay
    arrow_style='->',
    cmap='viridis'
)
plt.show()
```

## Customization Options

### Styling

```python
# Set global style parameters
import matplotlib.pyplot as plt

plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 12
plt.rcParams['axes.grid'] = True

# Custom color palettes
custom_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
```

### Save Options

```python
# Save high-quality figures
fig = plot_umap(data, color='leiden_cluster')
fig.savefig(
    'output/umap_clusters.pdf', 
    dpi=300, 
    bbox_inches='tight',
    transparent=True
)
```

## Interactive Plotting

### Plotly Integration

```python-snippet
try:
    import plotly.express as px
    import plotly.graph_objects as go
    
    # Interactive UMAP plot
    def plot_interactive_umap(data, color='leiden'):
        """Create interactive UMAP plot with plotly."""
        
        df_plot = pd.DataFrame({
            'UMAP_1': data.obsm['X_umap'][:, 0],
            'UMAP_2': data.obsm['X_umap'][:, 1],
            'color': data.obs[color].astype(str)
        })
        
        fig = px.scatter(
            df_plot,
            x='UMAP_1',
            y='UMAP_2', 
            color='color',
            title=f'Interactive UMAP colored by {color}',
            hover_data=['color']
        )
        
        return fig
    
    # Display interactive plot
    interactive_fig = plot_interactive_umap(data, 'leiden')
    interactive_fig.show()
    
except ImportError:
    print("Plotly not available for interactive plots")
```

## Integration with Core Visualization

The single-cell visualization module builds upon METAINFORMANT's core visualization utilities:

```python-snippet
from metainformant.visualization.config.themes import apply_theme
from metainformant.core.io import ensure_directory

# Set up consistent styling
apply_theme("publication")

# Ensure output directory exists  
output_dir = ensure_directory("output/singlecell_plots")
```

## Batch Plotting

### Generate Multiple Plots

```python
def generate_analysis_plots(data, output_dir="output/plots"):
    """Generate standard set of analysis plots."""
    
    ensure_directory(output_dir)
    
    # 1. QC overview
    fig = plot_qc_metrics(data)
    fig.savefig(f"{output_dir}/qc_metrics.pdf", dpi=300)
    plt.close()
    
    # 2. UMAP clusters
    fig = plot_umap(data, color='leiden_cluster')  
    fig.savefig(f"{output_dir}/umap_clusters.pdf", dpi=300)
    plt.close()
    
    # 3. t-SNE clusters
    if 'X_tsne' in data.obsm:
        fig = plot_tsne(data, color='leiden_cluster')
        fig.savefig(f"{output_dir}/tsne_clusters.pdf", dpi=300)
        plt.close()
    
    # 4. PCA
    fig = plot_pca(data, components=[1, 2])
    fig.savefig(f"{output_dir}/pca.pdf", dpi=300) 
    plt.close()
    
    print(f"Analysis plots saved to {output_dir}/")

# Generate all plots
generate_analysis_plots(data)
```

## Performance Considerations

### Large Dataset Visualization

```python
# For datasets > 100k cells
# 1. Subsample for plotting
n_sample = min(50000, data.n_obs)
sample_idx = np.random.choice(data.n_obs, n_sample, replace=False)
data_subset = data[sample_idx, :].copy()

# 2. Plot subset
fig = plot_umap(data_subset, color='leiden_cluster')
```

### Memory Management

```python
# Clear figure memory
import gc

def plot_with_cleanup(plot_func, *args, **kwargs):
    """Plot function with automatic memory cleanup."""
    fig = plot_func(*args, **kwargs)
    plt.show()
    plt.close('all')
    gc.collect()
    return fig
```

## Testing

Visualization tests are integrated into the single-cell test suite:

```bash
# Run visualization tests
uv run pytest tests/singlecell/test_singlecell_visualization.py -v
```

## Related Documentation

- [Preprocessing](./preprocessing.md): For data preparation
- [Dimensionality Reduction](./dimensionality.md): For embedding generation  
- [Clustering](./clustering.md): For cluster analysis
- [Core Visualization](../visualization/plots.md): For general plotting utilities
