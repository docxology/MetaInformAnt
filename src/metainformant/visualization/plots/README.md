# Visualization Plots

General-purpose plotting library covering basic charts, animations, multidimensional visualizations, and specialized diagram types.

## Contents

| File | Purpose |
|------|---------|
| `animations.py` | Animated time series, evolution, clustering, network, and trajectory |
| `basic.py` | Line, scatter, heatmap, bar, pie, area, and step plots |
| `general.py` | Convenience wrappers: expression heatmap, PCA, volcano, Manhattan |
| `multidim.py` | Pairwise relationships, parallel coordinates, radar, 3D scatter |
| `specialized.py` | Venn, Sankey, chord, alluvial, circular bar, UpSet plots |

## Key Functions

| Function | Description |
|----------|-------------|
| `lineplot()` | Basic line plot with optional styling |
| `scatter_plot()` | 2D scatter plot with grouping support |
| `heatmap()` | Matrix heatmap with clustering |
| `animate_time_series()` | Animated line chart over time |
| `animate_evolution()` | Animated evolutionary trajectory |
| `plot_pairwise_relationships()` | Scatterplot matrix for multi-variable data |
| `plot_parallel_coordinates()` | Parallel coordinates for high-dimensional data |
| `plot_venn_diagram()` | 2-set or 3-set Venn diagram |
| `plot_sankey_diagram()` | Flow diagram between categories |
| `plot_upset_plot()` | UpSet plot for set intersection analysis |

## Usage

```python
from metainformant.visualization.plots.basic import scatter_plot, heatmap
from metainformant.visualization.plots.specialized import plot_venn_diagram

scatter_plot(x, y, output_path="output/scatter.png")
heatmap(matrix, output_path="output/heatmap.png")
plot_venn_diagram({"A": set_a, "B": set_b}, output_path="output/venn.png")
```
