# Multi-Omics Visualization

Plotting functions for multi-omics integration results including correlation heatmaps, layer comparisons, pathway enrichment, and interactive dashboards.

## Contents

| File | Purpose |
|------|---------|
| `visualization.py` | Multi-omics correlation, integration, network, and factor analysis plots |

## Key Functions

| Function | Description |
|----------|-------------|
| `plot_multiomics_correlation_heatmap()` | Cross-omics correlation heatmap |
| `plot_integrated_omics_heatmap()` | Heatmap of the stacked omics layers; `feature_names` has one name per stacked row (a wrong length raises `ValueError`) |
| `plot_omics_layer_comparison()` | Side-by-side comparison of omics data layers |
| `plot_multiomics_pca()` | PCA scatter colored by omics layer or sample group |
| `plot_pathway_enrichment_integration()` | Pathway enrichment across multiple omics types (sorted rows; a pathway an omic reports without a `pvalue` is blank, not p = 1) |
| `plot_multiomics_network()` | Network linking features across omics layers (seeded layout, `layout_seed=42`) |
| `plot_omics_factor_analysis()` | Factor loadings from joint NMF or MOFA |
| `plot_multiomics_upset()` | UpSet plot for feature overlap across omics |
| `create_interactive_multiomics_dashboard()` | Interactive Plotly multi-omics explorer; the network panel draws `network_data` (a NetworkX graph or `{node: [neighbours]}`); the "Enrichment" panel is reserved and empty |

## Usage

```python
from metainformant.multiomics.visualization.visualization import (
    plot_multiomics_correlation_heatmap,
    plot_multiomics_pca,
)

plot_multiomics_correlation_heatmap(corr_matrix, ["rna", "protein"], output_path="output/corr.png")
plot_multiomics_pca(integrated_features, sample_labels=sample_names, output_path="output/pca.png")
```
