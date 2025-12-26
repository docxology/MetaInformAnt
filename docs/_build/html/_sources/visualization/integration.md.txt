# Domain Integration

The visualization module integrates with domain-specific modules to provide unified access to specialized visualization functions.

## GWAS Integration

Access GWAS visualization functions through the integration module:

```python
from metainformant.visualization.gwas_integration import (
    manhattan_plot,
    circular_manhattan_plot,
    qq_plot_stratified,
    regional_plot_detailed,
    pca_plot_gwas,
    GWAS_VISUALIZATION_AVAILABLE
)
```

## Single-Cell Integration

Access single-cell visualization functions:

```python
from metainformant.visualization.singlecell_integration import (
    plot_qc_metrics,
    plot_embedding,
    plot_gene_expression,
    SINGLECELL_VISUALIZATION_AVAILABLE
)
```

## Information Theory Integration

Access information theory visualization functions:

```python
from metainformant.visualization.information_integration import (
    plot_entropy_distribution,
    plot_mutual_information_matrix,
    plot_information_profile,
    INFORMATION_VISUALIZATION_AVAILABLE
)
```

## Life Events Integration

Access life events visualization functions:

```python
from metainformant.visualization.life_events_integration import (
    plot_event_timeline,
    plot_event_embeddings,
    plot_attention_heatmap,
    LIFE_EVENTS_VISUALIZATION_AVAILABLE
)
```

## Layout Utilities

### `create_subplot_grid(n_plots, ncols=3, *, figsize=None, sharex=False, sharey=False)`

Create a subplot grid for multiple plots.

### `create_multi_panel(n_panels, layout='grid', *, figsize=None, **kwargs)`

Create a multi-panel figure layout.

### `add_shared_axis_labels(fig, xlabel=None, ylabel=None)`

Add shared axis labels to a multi-panel figure.

## Export Utilities

### `save_figure(fig, path, *, dpi=300, bbox_inches='tight', format=None, **kwargs)`

Save a figure with high-resolution settings.

### `save_figure_multiformat(fig, base_path, *, formats=('png', 'pdf', 'svg'), dpi=300, **kwargs)`

Save a figure in multiple formats.

### `batch_export_figures(figures, output_dir, *, formats=('png',), dpi=300, **kwargs)`

Batch export multiple figures.

## Interactive Plots

### `create_interactive_scatter(x, y, *, labels=None, colors=None, title='Interactive Scatter Plot', **kwargs)`

Create an interactive scatter plot using Plotly (optional dependency).

### `create_interactive_heatmap(data, *, x_labels=None, y_labels=None, title='Interactive Heatmap', **kwargs)`

Create an interactive heatmap using Plotly.

