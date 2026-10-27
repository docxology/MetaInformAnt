# Domain Integration

The generic plots live in `metainformant.visualization`, and each domain package keeps its own plotting module next to its analysis code. This page shows where the domain plots are, and the helpers in the visualization module for combining figures into one layout, exporting them, and building interactive pages.

## GWAS Integration

GWAS plots live in `metainformant.gwas.visualization`. The package re-exports `manhattan_plot` and `qq_plot` and groups the remaining plots in the `genomic`, `statistical`, `population` and `interactive` subpackages. Their arguments differ from the generic plots (they take GWAS result records rather than a DataFrame); see the [GWAS Visualization Gallery](../gwas/visualization_gallery.md) for worked examples.

```python
from metainformant.gwas.visualization import manhattan_plot, qq_plot
from metainformant.gwas.visualization.genomic.genome import circular_manhattan_plot
from metainformant.gwas.visualization.genomic.regional import regional_plot
from metainformant.gwas.visualization.population.population_pca import pca_plot
from metainformant.gwas.visualization.statistical.statistical import qq_plot_stratified
```

## Single-Cell Integration

Single-cell plots live in `metainformant.singlecell.visualization.visualization`. They take a `SingleCellData` object as their first argument: `plot_qc_metrics`, the embedding plots `plot_umap`, `plot_tsne` and `plot_pca` (with an optional `color=` column), `plot_marker_expression` for marker genes, and `plot_trajectory` and `plot_cluster_comparison`. For a QC summary built from plain arrays use `plot_singlecell_qc_metrics` in [Quality Control Plots](quality.md).

```python
from metainformant.singlecell.visualization.visualization import (
    plot_marker_expression,
    plot_pca,
    plot_qc_metrics,
    plot_tsne,
    plot_umap,
)
```

## Information Theory Integration

Information theory plots are part of the visualization module itself, in `metainformant.visualization.analysis.information` (see [Information Theory Plots](information.md) for usage):

```python
from metainformant.visualization.analysis.information import (
    plot_entropy_profile,
    plot_information_landscape,
    plot_information_network,
    plot_mutual_information_matrix,
    plot_renyi_spectra,
)
```

## Life Events Integration

Life event plots live in `metainformant.life_events.visualization`; the timeline plots are in its `timeline` module and the embedding and attention plots in its `statistical` module:

```python
from metainformant.life_events.visualization.timeline import plot_event_timeline
from metainformant.life_events.visualization.statistical import (
    plot_attention_heatmap,
    plot_event_embeddings,
)
```

## Layout Utilities

### `multi_panel(panels, *, ncols=2, figsize=None, title=None, output_path=None)`

Arrange several plots in one figure grid (`metainformant.visualization.dashboards.composite`). Each entry of `panels` is a dict with `func` (a plotting function that accepts `ax=`, such as any plot function in this package), `kwargs` (its keyword arguments) and an optional `title`. Unused grid cells are hidden, `figsize` defaults to 6 x `ncols` by 5 x rows inches, and the matplotlib `Figure` is returned. A single panel works with any `ncols`; the other cells of its row stay hidden (`ncols=1` gives a figure sized for one panel).

**Example:**
```python
import numpy as np
from metainformant.visualization.dashboards.composite import multi_panel
from metainformant.visualization.plots.basic import lineplot, scatter_plot

rng = np.random.default_rng(0)
x = np.arange(20)
fig = multi_panel(
    [
        {"func": lineplot, "kwargs": {"x": x, "y": x**2}, "title": "Line"},
        {"func": scatter_plot, "kwargs": {"x": x, "y": rng.random(20)}, "title": "Scatter"},
    ],
    ncols=2,
    title="Overview",
)
```

### `genomic_overview(data, *, figsize=(16, 12), output_path=None)` and `qc_summary(qc_data, *, figsize=(14, 10), output_path=None)`

Ready-made multi-panel dashboards in the same module. Both take a dict of data and draw a panel only for the keys that are present, returning the `Figure`. `genomic_overview` reads `expression` (a DataFrame), `pca_data` (a 2-D array), `pvalues`, `log2fc` (with `pvalues`, for a volcano plot), `gc_content` and `quality_scores`. `qc_summary` reads `read_lengths`, `gc_content`, `quality_scores`, `mapping_rates` and `duplication_rates` (sample to rate dicts) and `adapter_content` (adapter to maximum percentage).

**Example:**
```python
import numpy as np
from metainformant.visualization.dashboards.composite import qc_summary

rng = np.random.default_rng(0)
fig = qc_summary({
    "read_lengths": rng.normal(150, 5, 500),
    "gc_content": rng.normal(45, 5, 500),
    "mapping_rates": {"S1": 95.2, "S2": 91.8},
    "adapter_content": {"Illumina": 1.2, "Nextera": 0.3},
})
```

## Export Utilities

### `save_figure_deterministic(fig, path, *, dpi=300, bbox_inches=None)`

Save a matplotlib figure (`metainformant.visualization.config.conventions`) with the producer and timestamp metadata blanked, so identical figures produce identical files. The output format is taken from the suffix of `path`. The parent directory must already exist (the `output_path=` argument of the plot functions creates it for you, and most plot functions in this package save through this helper).

**Example:**
```python
from pathlib import Path
import matplotlib.pyplot as plt
from metainformant.visualization.config.conventions import save_figure_deterministic

fig, ax = plt.subplots()
ax.plot([1, 2, 3])
Path("output/figures").mkdir(parents=True, exist_ok=True)
save_figure_deterministic(fig, "output/figures/line.png", dpi=300, bbox_inches="tight")
```

## Interactive Plots

### `create_interactive_scatter(x, y, labels=None, colors=None, hover_data=None, title='Scatter Plot')`

Build a JSON-serializable dict describing an interactive scatter plot (`metainformant.visualization.interactive_dashboards.dashboards`); no plotting library is needed. `x` and `y` are equal-length lists, `labels` and `colors` give one value per point, and `hover_data` maps extra field names to lists of values. The dict has the keys `plot_data`, `layout` and `config`.

### `create_interactive_heatmap(data, row_names, col_names, cluster_rows=True, cluster_cols=True)`

Build the dict for an interactive heatmap from a 2-D list of values with one name per row and column. With clustering on, rows and columns are re-ordered by a greedy nearest-neighbour ordering. The dict has the keys `plot_data`, `row_order`, `col_order`, `dendrograms` and `color_range`.

### `export_to_html(plot_data, output_path, title='Interactive Plot')`

Write a standalone HTML page (embedded SVG and plain JavaScript, no external dependencies) from the dict returned by any `create_interactive_*` function and return its absolute path. The page draws records that have `x` and `y` keys (scatter and volcano data) with hover tooltips; for other structures it only reports the number of data points.

**Example:**
```python
from metainformant.visualization.interactive_dashboards.dashboards import (
    create_interactive_heatmap,
    create_interactive_scatter,
    export_to_html,
)

scatter = create_interactive_scatter([1, 2, 3], [4, 5, 6], labels=["a", "b", "c"], title="Demo")
html_path = export_to_html(scatter, "output/interactive/scatter.html", title="Demo")

heatmap = create_interactive_heatmap(
    [[1.0, 2.0], [3.0, 4.0]], row_names=["r1", "r2"], col_names=["c1", "c2"]
)
```

Plotly-based wrappers with a static matplotlib fallback when Plotly is not installed are in `metainformant.visualization.dashboards.interactive` (`interactive_scatter`, `interactive_heatmap`, `interactive_volcano` and `interactive_manhattan`). They return a Plotly `Figure` with Plotly and a matplotlib `Axes` without it, and each takes an `output_path`:

- **With Plotly**, the figure is written as an HTML file; a suffix other than `.html` is replaced by `.html`.
- **Without Plotly**, an `.html` or `.htm` path raises `ImportError` before any figure, file or directory is created, because matplotlib cannot write HTML. Any other path is saved as a static image.

```python
import numpy as np
import pandas as pd

from metainformant.visualization.dashboards.interactive import interactive_scatter

rng = np.random.default_rng(0)
df = pd.DataFrame({"PC1": rng.normal(size=30), "PC2": rng.normal(size=30), "cluster": rng.choice(["a", "b"], 30)})

try:
    fig = interactive_scatter(df, x="PC1", y="PC2", color="cluster", output_path="output/interactive/scatter.html")
except ImportError:  # no Plotly: the HTML path is refused, an image path gives the static fallback
    ax = interactive_scatter(df, x="PC1", y="PC2", color="cluster", output_path="output/interactive/scatter.png")
```
