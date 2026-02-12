# Visualization

Publication-quality plotting functions for spatial transcriptomics data, including spatial scatter plots, tissue image overlays, expression maps, cell type maps, and deconvolution visualizations.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports plots submodule |
| `plots.py` | All spatial visualization functions using matplotlib/seaborn |

## Key Functions

| Function | Description |
|----------|-------------|
| `plots.plot_spatial_scatter()` | Spatial scatter plot colored by gene expression or metadata |
| `plots.plot_tissue_overlay()` | Overlay expression data on H&E tissue image |
| `plots.plot_gene_expression_map()` | Spatial gene expression heatmap |
| `plots.plot_cell_type_map()` | Cell type composition map across tissue |
| `plots.plot_neighborhood_graph()` | Spatial neighborhood graph visualization |
| `plots.plot_interaction_heatmap()` | Cell-cell interaction score heatmap |
| `plots.plot_deconvolution_pie()` | Pie charts showing cell type fractions at each spot |
| `plots.plot_spatial_autocorrelation()` | LISA map showing spatial autocorrelation patterns |

## Usage

```python
from metainformant.spatial.visualization import plots

plots.plot_spatial_scatter(coordinates, values, output_path="output/scatter.png")
plots.plot_tissue_overlay(coordinates, expression, image, output_path="output/overlay.png")
plots.plot_cell_type_map(coordinates, cell_types, output_path="output/celltypes.png")
plots.plot_deconvolution_pie(coordinates, fractions, output_path="output/deconv.png")
```
