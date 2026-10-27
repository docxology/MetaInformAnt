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
| `plot_trajectory()` | Embedding colored by pseudotime; `trajectory_key` (for example `"dpt"`) is the required `uns` entry to draw |
| `plot_marker_expression()` | Dotplot, heatmap, or violin for marker genes; the heatmap shows at most `max_cells` (default 1000) randomly drawn cells, seeded by `random_state` |
| `plot_qc_metrics()` | QC histograms (counts, genes, mito %) and counts-vs-genes scatter |
| `plot_cluster_comparison()` | Side-by-side scatter of the embedding colored by each clustering column |

## Usage

```python
from metainformant.singlecell.visualization.visualization import plot_umap, plot_marker_expression

fig = plot_umap(data, color="leiden_cluster")   # an obs column, for example the labels of leiden_clustering
fig.savefig("output/umap.png")
plot_marker_expression(data, marker_genes=["CD4", "CD8A"], method="dotplot")
```
