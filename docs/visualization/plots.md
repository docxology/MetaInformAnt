### Visualization: Plots

The `metainformant.visualization.plots` package provides the core plotting functions
organized into five submodules: basic charts, general biological plots, specialized
diagrams, multi-dimensional visualizations, and animations.

All plot functions return a `matplotlib.axes.Axes` object and accept an optional
`output_path` parameter to save the figure to disk.

---

## Submodule Overview

| Submodule      | Description                                            | Plot Types                    |
|----------------|--------------------------------------------------------|-------------------------------|
| `basic`        | Fundamental chart types built on matplotlib             | line, scatter, heatmap, bar, pie, area, step |
| `general`      | Biology-oriented statistical plots                      | volcano, manhattan, QQ, PCA, expression heatmap, correlation |
| `specialized`  | Advanced diagram types for systems biology              | Venn, Sankey, chord, alluvial, circular bar, UpSet |
| `multidim`     | Multi-dimensional data exploration                      | pairwise, parallel coordinates, radar, 3D scatter |
| `animations`   | Animated visualizations for temporal/dynamic data       | time series, evolution, clustering, network, trajectory |

---

## Basic Plots (`plots.basic`)

### `lineplot`

```python
import numpy as np
from metainformant.visualization import lineplot

x = np.linspace(0, 10, 100)
y = np.sin(x)
ax = lineplot(x, y, label="sine wave", color="blue")
```

| Parameter     | Type                | Default | Description                         |
|---------------|---------------------|---------|-------------------------------------|
| `x`           | `np.ndarray`        | required| X-axis data                         |
| `y`           | `np.ndarray | None` | `None`  | Y-axis data (if None, x used as y)  |
| `ax`          | `Axes | None`       | `None`  | Existing axes to plot on             |
| `output_path` | `str | Path | None` | `None`  | Path to save figure                  |

### `scatter_plot`

```python
from metainformant.visualization import scatter_plot

ax = scatter_plot(np.random.randn(100), np.random.randn(100), alpha=0.6)
```

### `heatmap`

```python
from metainformant.visualization import heatmap

data = np.random.rand(10, 10)
ax = heatmap(data, cmap="viridis")
```

Requires a 2D numpy array. Includes an automatic colorbar.

### `bar_plot`

```python
from metainformant.visualization import bar_plot

categories = np.array(["A", "B", "C", "D"])
values = np.array([10, 25, 15, 30])
ax = bar_plot(categories, values, color="teal")
```

### `pie_chart`

```python
from metainformant.visualization import pie_chart

sizes = np.array([35, 25, 20, 20])
ax = pie_chart(sizes, labels=["A", "B", "C", "D"], autopct="%1.1f%%")
```

### `area_plot`

Filled area chart using `matplotlib.fill_between`.

```python
from metainformant.visualization import area_plot

ax = area_plot(np.arange(50), np.random.rand(50).cumsum(), alpha=0.4)
```

### `step_plot`

Step function plot.

```python
from metainformant.visualization import step_plot

ax = step_plot(np.arange(20), np.random.randint(0, 10, 20))
```

---

## General Plots (`plots.general`)

### `volcano_plot`

Visualize differential expression results as log2 fold change vs. -log10 p-value.

```python
import pandas as pd
from metainformant.visualization import volcano_plot

de_results = pd.DataFrame({
    "log2FoldChange": np.random.randn(1000),
    "padj": np.random.uniform(0, 1, 1000),
})
ax = volcano_plot(de_results, log2fc_col="log2FoldChange", pval_col="padj")
```

### `manhattan_plot`

Genome-wide association study results plotted by chromosomal position.

```python
from metainformant.visualization import manhattan_plot

gwas_data = pd.DataFrame({
    "CHR": np.repeat(range(1, 23), 50),
    "BP": np.random.randint(1, 1_000_000, 1100),
    "P": np.random.uniform(0, 1, 1100),
})
ax = manhattan_plot(gwas_data, pos_col="BP", pval_col="P")
```

### `qq_plot`

Quantile-quantile plot for assessing distributional assumptions and identifying
genomic inflation.

```python
from metainformant.visualization import qq_plot

pvalues = np.random.uniform(0, 1, 500)
ax = qq_plot(pvalues)
```

### `expression_heatmap`

Clustered heatmap of gene expression data. Uses seaborn if available, falls back
to matplotlib.

```python
from metainformant.visualization import expression_heatmap

expr_data = pd.DataFrame(
    np.random.randn(20, 6),
    index=[f"Gene_{i}" for i in range(20)],
    columns=[f"Sample_{i}" for i in range(6)],
)
ax = expression_heatmap(expr_data, cmap="RdBu_r")
```

### `pca_plot`

Scatter plot of pre-computed principal components with variance annotations.

```python
from metainformant.visualization import pca_plot

pca_data = pd.DataFrame({
    "PC1": np.random.randn(50),
    "PC2": np.random.randn(50),
    "PC1_variance": [25.5] * 50,
    "PC2_variance": [18.3] * 50,
    "group": np.random.choice(["Control", "Treatment"], 50),
})
ax = pca_plot(pca_data, hue="group")
```

### `correlation_heatmap`

Correlation matrix as annotated heatmap.

```python
from metainformant.visualization import correlation_heatmap

df = pd.DataFrame(np.random.randn(100, 5), columns=["A", "B", "C", "D", "E"])
ax = correlation_heatmap(df, annot=True)
```

---

## Specialized Plots (`plots.specialized`)

| Function                       | Description                                     |
|--------------------------------|-------------------------------------------------|
| `plot_venn_diagram`            | 2-3 set Venn diagram with intersection counts   |
| `plot_sankey_diagram`          | Flow diagram showing transitions between states |
| `plot_chord_diagram`           | Circular diagram for pairwise relationships     |
| `plot_alluvial_diagram`        | Parallel category flow visualization            |
| `plot_circular_barplot`        | Radial bar chart for cyclic data                |
| `plot_network_circular_layout` | Network graph with circular node arrangement    |
| `plot_upset_plot`              | UpSet plot for multi-set intersection analysis  |

```python
from metainformant.visualization import plot_venn_diagram

sets = {
    "Gene Set A": {"BRCA1", "TP53", "EGFR", "MYC"},
    "Gene Set B": {"TP53", "KRAS", "MYC", "PIK3CA"},
}
ax = plot_venn_diagram(sets)
```

---

## Multi-Dimensional Plots (`plots.multidim`)

| Function                       | Description                                        |
|--------------------------------|----------------------------------------------------|
| `plot_pairwise_relationships`  | Scatter matrix (pairplot) for all numeric columns   |
| `plot_parallel_coordinates`    | Parallel coordinate plot for multivariate data      |
| `plot_radar_chart`             | Radar/spider chart for comparing profiles           |
| `plot_3d_scatter`              | 3D scatter plot using mpl_toolkits                  |

```python
from metainformant.visualization import plot_pairwise_relationships

df = pd.DataFrame(np.random.randn(100, 4), columns=["W", "X", "Y", "Z"])
ax = plot_pairwise_relationships(df)
```

---

## Animations (`plots.animations`)

All animation functions return a `(Figure, FuncAnimation)` tuple and can save
to GIF via the `output_path` parameter.

| Function               | Description                                      |
|------------------------|--------------------------------------------------|
| `animate_time_series`  | Progressive line plot revealing data over time    |
| `animate_evolution`    | Sequence/population change across generations     |
| `animate_clustering`   | Step-by-step clustering algorithm visualization   |
| `animate_network`      | Dynamic network graph with edge/node changes      |
| `animate_trajectory`   | Point movement path animation                     |

```python
from metainformant.visualization import animate_time_series

data = np.random.randn(100, 3).cumsum(axis=0)
fig, anim = animate_time_series(data, interval=100, output_path="output/ts.gif")
```

---

## Common Parameters

All plot functions share these optional parameters:

| Parameter     | Type                | Description                                |
|---------------|---------------------|--------------------------------------------|
| `ax`          | `Axes | None`       | Existing matplotlib axes (creates new if None) |
| `output_path` | `str | Path | None` | Save figure to this path (PNG, PDF, SVG)   |
| `figsize`     | `tuple`             | Figure size as `(width, height)` in inches |
| `**kwargs`    | `Any`               | Passed through to underlying matplotlib call |

Figures are saved at 300 DPI with tight bounding boxes. Parent directories are
created automatically.

---

## See Also

- `metainformant.visualization.config.themes` -- Theme and styling management
- `metainformant.visualization.config.palettes` -- Color palettes (Wong, IBM, etc.)
- `metainformant.visualization.dashboards` -- Composite multi-panel layouts
- `metainformant.visualization.interactive_dashboards` -- Plotly-based interactive plots
