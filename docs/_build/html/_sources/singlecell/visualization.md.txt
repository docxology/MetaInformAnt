# Single-Cell Visualization

The visualization module provides specialized plotting functions for single-cell RNA-seq data analysis and results.

## Quality Control Plots

### plot_qc_metrics()

Visualize quality control metrics:

```python
from metainformant.singlecell.visualization import plot_qc_metrics

# Plot QC metrics overview
fig = plot_qc_metrics(
    data,
    metrics=['total_counts', 'n_genes', 'pct_mt', 'pct_ribo'],
    plot_type='violin'     # 'violin', 'box', 'hist'
)
plt.show()
```

### plot_qc_scatter()

Scatter plots of QC relationships:

```python
from metainformant.singlecell.visualization import plot_qc_scatter

# Relationship between total counts and genes detected
fig = plot_qc_scatter(
    data,
    x='total_counts',
    y='n_genes', 
    color='pct_mt',
    size=20
)
plt.show()
```

## Dimensionality Reduction Plots

### plot_embedding()

Plot cells in reduced dimensional space:

```python
from metainformant.singlecell.visualization import plot_embedding

# UMAP colored by clusters
fig = plot_embedding(
    data,
    basis='umap',          # 'umap', 'tsne', 'pca'
    color='leiden',        # Color by cluster
    legend_loc='right',
    size=30
)

# UMAP colored by gene expression
fig = plot_embedding(
    data,
    basis='umap',
    color='CD3E',          # Gene name
    cmap='viridis'
)
plt.show()
```

### plot_pca()

Specialized PCA visualization:

```python
from metainformant.singlecell.visualization import plot_pca

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

### plot_gene_expression()

Visualize individual gene expression:

```python
from metainformant.singlecell.visualization import plot_gene_expression

# Violin plots by cluster
fig = plot_gene_expression(
    data,
    genes=['CD3E', 'CD19', 'CD14'],
    groupby='leiden',
    plot_type='violin'
)

# Dot plot showing expression and detection
fig = plot_gene_expression(
    data,
    genes=['CD3E', 'CD19', 'CD14', 'LYZ'],
    groupby='leiden', 
    plot_type='dot'
)
plt.show()
```

### plot_heatmap()

Expression heatmaps:

```python
from metainformant.singlecell.visualization import plot_heatmap

# Heatmap of marker genes by cluster
fig = plot_heatmap(
    data,
    genes=['CD3E', 'CD19', 'CD14', 'GNLY'],  # Marker genes
    groupby='leiden',
    standard_scale='var',   # Scale by gene
    cmap='RdBu_r'
)
plt.show()
```

## Clustering Visualization

### plot_clusters()

Comprehensive cluster visualization:

```python
from metainformant.singlecell.visualization import plot_clusters

# Multiple cluster views
fig = plot_clusters(
    data,
    cluster_column='leiden',
    embedding='umap',
    show_cluster_labels=True,
    show_statistics=True
)
plt.show()
```

### plot_cluster_composition()

Cluster composition analysis:

```python
from metainformant.singlecell.visualization import plot_cluster_composition

# Bar plot of cluster sizes
fig = plot_cluster_composition(
    data,
    groupby='leiden',
    plot_type='bar'
)

# Pie chart of cluster proportions  
fig = plot_cluster_composition(
    data,
    groupby='leiden',
    plot_type='pie'
)
plt.show()
```

## Trajectory Visualization

### plot_trajectory()

Visualize developmental trajectories:

```python
from metainformant.singlecell.visualization import plot_trajectory

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

### plot_gene_trends()

Gene expression along trajectories:

```python
from metainformant.singlecell.visualization import plot_gene_trends

# Expression trends along pseudotime
fig = plot_gene_trends(
    data,
    genes=['GENE1', 'GENE2', 'GENE3'],
    pseudotime_col='pseudotime',
    smooth=True,
    show_confidence=True
)
plt.show()
```

## Comparative Visualization

### plot_comparison()

Compare different analyses:

```python
from metainformant.singlecell.visualization import plot_comparison

# Compare clusterings at different resolutions
fig = plot_comparison(
    data,
    groupby=['leiden_0.2', 'leiden_0.5', 'leiden_1.0'],
    basis='umap',
    ncols=3
)
plt.show()
```

### plot_split()

Split visualization by categories:

```python
from metainformant.singlecell.visualization import plot_split

# Split by sample or condition
fig = plot_split(
    data,
    split_by='sample',      # Column to split by
    color='leiden',         # Color by cluster
    basis='umap',
    ncols=2
)
plt.show()
```

## Statistical Plots

### plot_marker_genes()

Visualize marker gene results:

```python
from metainformant.singlecell.visualization import plot_marker_genes

# Volcano plot of marker genes
fig = plot_marker_genes(
    markers_df,             # From find_markers()
    cluster='0',
    plot_type='volcano',
    top_n=10               # Label top 10 genes
)

# Rank plot
fig = plot_marker_genes(
    markers_df,
    cluster='0',
    plot_type='rank',
    n_genes=20
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
fig = plot_embedding(data, basis='umap', color='leiden')
fig.savefig(
    'output/umap_clusters.pdf', 
    dpi=300, 
    bbox_inches='tight',
    transparent=True
)
```

## Interactive Plotting

### Plotly Integration

```python
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

## Multi-Panel Figures

### Complex Layouts

```python
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# Create complex multi-panel figure
fig = plt.figure(figsize=(20, 12))
gs = GridSpec(3, 4, figure=fig)

# Panel A: QC metrics
ax1 = fig.add_subplot(gs[0, :2])
plot_qc_metrics(data, ax=ax1)
ax1.set_title('A. Quality Control Metrics')

# Panel B: UMAP by clusters  
ax2 = fig.add_subplot(gs[0, 2:])
plot_embedding(data, basis='umap', color='leiden', ax=ax2)
ax2.set_title('B. UMAP by Clusters')

# Panel C: Gene expression
ax3 = fig.add_subplot(gs[1, :])
plot_gene_expression(data, genes=['CD3E', 'CD19', 'CD14'], 
                    groupby='leiden', ax=ax3)
ax3.set_title('C. Marker Gene Expression')

# Panel D: Heatmap
ax4 = fig.add_subplot(gs[2, :])
plot_heatmap(data, genes=['CD3E', 'CD19', 'CD14'], 
            groupby='leiden', ax=ax4)
ax4.set_title('D. Expression Heatmap')

plt.tight_layout()
plt.savefig('output/comprehensive_analysis.pdf', dpi=300)
plt.show()
```

## Integration with Core Visualization

The single-cell visualization module builds upon METAINFORMANT's core visualization utilities:

```python
from metainformant.visualization.plots import setup_matplotlib_style
from metainformant.core.io import ensure_directory

# Set up consistent styling
setup_matplotlib_style()

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
    fig = plot_embedding(data, basis='umap', color='leiden')  
    fig.savefig(f"{output_dir}/umap_clusters.pdf", dpi=300)
    plt.close()
    
    # 3. t-SNE clusters
    if 'X_tsne' in data.obsm:
        fig = plot_embedding(data, basis='tsne', color='leiden')
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
fig = plot_embedding(data_subset, basis='umap', color='leiden')

# 3. Use density plots for very large datasets
fig = plot_embedding(data, basis='umap', plot_type='density')
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
uv run pytest tests/test_singlecell_visualization.py -v
```

## Related Documentation

- [Preprocessing](./preprocessing.md): For data preparation
- [Dimensionality Reduction](./dimensionality.md): For embedding generation  
- [Clustering](./clustering.md): For cluster analysis
- [Core Visualization](../visualization/plots.md): For general plotting utilities
