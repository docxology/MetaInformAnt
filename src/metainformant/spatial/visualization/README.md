# Visualization

Publication-quality plotting functions for spatial transcriptomics data, including spatial scatter plots, tissue image overlays, expression maps, cell type maps, and deconvolution visualizations.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports plots submodule |
| `plots.py` | All spatial visualization functions using matplotlib (seaborn is used for the interaction heatmap when installed) |

## Key Functions

| Function | Description |
|----------|-------------|
| `plots.plot_spatial_scatter()` | Spatial scatter plot colored by gene expression or metadata |
| `plots.plot_tissue_overlay()` | Overlay expression data on H&E tissue image |
| `plots.plot_gene_expression_map()` | Spatial gene expression heatmap (takes a dataset) |
| `plots.plot_cell_type_map()` | Cell type map across tissue (takes a dataset or a coordinate array) |
| `plots.plot_neighborhood_graph()` | Spatial neighborhood graph visualization |
| `plots.plot_interaction_heatmap()` | Cell-cell interaction score heatmap (takes an interaction matrix, no coordinates) |
| `plots.plot_deconvolution_pie()` | Pie charts showing cell type fractions at each spot |
| `plots.plot_spatial_autocorrelation()` | LISA map showing spatial autocorrelation patterns |

## Coordinate order

Every function that draws coordinates takes `coordinate_order`: `"row_col"` (column 0 is the vertical axis, as in the Visium loader's `[pixel_row, pixel_col]`) or `"xy"` (column 0 is x, as in the Xenium and MERFISH loaders). `plot_gene_expression_map` and `plot_cell_type_map` default to the `coordinate_order` declared by the dataset object (`SpatialDataset`: `"row_col"`; `XeniumDataset` and `MERFISHDataset`: `"xy"`; `"row_col"` when nothing is declared). The other functions take bare coordinate arrays and default to `"row_col"`, so pass `coordinate_order="xy"` for Xenium and MERFISH arrays; with the wrong order the tissue appears transposed. Any other value raises `ValueError`.

## Usage

```python
from metainformant.spatial.visualization import plots

plots.plot_spatial_scatter(coordinates, values, output_path="output/scatter.png")
plots.plot_tissue_overlay(coordinates, expression, image, output_path="output/overlay.png")
plots.plot_cell_type_map(dataset, cell_types, output_path="output/celltypes.png")   # reads dataset.coordinate_order
plots.plot_deconvolution_pie(xy_coordinates, fractions, output_path="output/deconv.png", coordinate_order="xy")
```
