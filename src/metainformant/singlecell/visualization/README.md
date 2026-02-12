# Single-Cell Visualization

Plotting functions for single-cell RNA-seq data including embeddings, marker expression, trajectory, and QC visualizations.

## Contents

| File | Purpose |
|------|---------|
| `visualization.py` | UMAP, t-SNE, PCA plots, marker dotplots, trajectory, and QC figures |

## Key Functions

| Function | Description |
|----------|-------------|
| `plot_umap()` | UMAP embedding colored by cluster, gene, or metadata |
| `plot_tsne()` | t-SNE embedding visualization |
| `plot_pca()` | PCA scatter plot for single-cell data |
| `plot_trajectory()` | Trajectory with pseudotime coloring |
| `plot_marker_expression()` | Dotplot, heatmap, or violin for marker genes |
| `plot_qc_metrics()` | QC violin plots for gene counts, UMIs, mito fraction |
| `plot_cluster_comparison()` | Side-by-side cluster composition comparison |

## Usage

```python
from metainformant.singlecell.visualization.visualization import plot_umap, plot_marker_expression

plot_umap(data, color="cluster", output_path="output/umap.png")
plot_marker_expression(data, marker_genes=["CD4", "CD8A"], method="dotplot")
```
