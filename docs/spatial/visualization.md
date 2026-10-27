# Visualization

Plotting functions live in `metainformant.spatial.visualization.plots`:
`plot_spatial_scatter` (a quick one‑liner for any per‑spot values) and more
specialised `plot_*` functions. Every function except `plot_interaction_heatmap`
takes spot coordinates (or a dataset object from [io](io.md)), saves the figure to
`output_path` and returns the matplotlib `Figure`.

## Coordinate order

A coordinate array is `n × 2`, but the platforms order its columns differently: the Visium loader returns `[row, col]` (`"row_col"`: column 0 is the vertical image axis), the Xenium and MERFISH loaders return `[x, y]` (`"xy"`). Every coordinate-taking plot has a `coordinate_order` argument that puts the columns on the right axes; with the wrong order the tissue appears transposed.

| Function | Takes | `coordinate_order` |
|----------|-------|--------------------|
| `plot_spatial_scatter`, `plot_tissue_overlay`, `plot_neighborhood_graph`, `plot_deconvolution_pie`, `plot_spatial_autocorrelation` | a bare coordinate array | keyword, default `'row_col'` |
| `plot_gene_expression_map`, `plot_cell_type_map` | a dataset (or, for `plot_cell_type_map`, a bare array) | keyword, default `None`: the dataset's own `coordinate_order`, else `'row_col'` |
| `plot_interaction_heatmap` | an interaction matrix | none (no coordinates) |

`SpatialDataset` (Visium) declares `"row_col"`; `XeniumDataset` and `MERFISHDataset` declare `"xy"`. Bare arrays and objects that declare nothing are read as `"row_col"`, so pass the order explicitly for Xenium and MERFISH arrays:

```python-snippet
from metainformant.spatial.io.xenium import load_xenium
from metainformant.spatial.visualization.plots import plot_cell_type_map, plot_spatial_scatter

dataset = load_xenium('xenium_run/outs/')
plot_cell_type_map(dataset, cell_types, 'cell_types.png')       # reads dataset.coordinate_order ('xy')
plot_spatial_scatter(
    dataset.coordinates, values, 'values.png',
    coordinate_order=dataset.coordinate_order,                  # or 'xy'
)
```

A value other than `'row_col'` or `'xy'` raises `ValueError`. An explicit `coordinate_order` also overrides the order declared by a dataset.

## Scatter / spot plot

```python-snippet
from metainformant.spatial.visualization.plots import plot_spatial_scatter
plot_spatial_scatter(
    coordinates,
    expression[:, gene_names.index('EPCAM')],   # numeric values or categorical labels
    'epcam.pdf',                                # output_path
    cmap='viridis',
    point_size=10.0,
    vmin=0, vmax=10,
    title='Epithelial marker',
)
```

Numeric values get a colour bar; string labels get one colour per category.

## Tissue overlay

```python-snippet
from metainformant.spatial.visualization.plots import plot_tissue_overlay
plot_tissue_overlay(
    dataset.coordinates,
    domain_labels,                  # categorical labels (or numeric values)
    dataset.image,                  # tissue image array (H, W, C)
    'domains_overlay.png',
    scale_factor=dataset.scale_factors.get('tissue_hires_scalef', 1.0),
)
```

Overlays the spots, coloured by label, on the tissue image (`dataset` is a
`SpatialDataset`, e.g. from `load_visium()`). The coordinates are multiplied by
`scale_factor` before they are drawn. Polygon outlines around each domain are
not provided.

## Expression maps

```python-snippet
from metainformant.spatial.visualization.plots import plot_gene_expression_map
for gene in ['EPCAM', 'VIM', 'CD3D']:
    plot_gene_expression_map(dataset, gene, f'{gene}.png', cmap='Reds')
```

`dataset` is a `SpatialDataset`, `XeniumDataset` or `MERFISHDataset` (or any object with `.expression`, `.coordinates` and
`.gene_names`); the dataset's `coordinate_order` decides which column is plotted on which axis. Each call writes one figure; there is no multi‑panel layout.

## QC plots

Any per‑spot quantity can be mapped with `plot_spatial_scatter`:

- **Library size per spot:** `plot_spatial_scatter(coordinates, expression.sum(axis=1), 'library_size.png')`
- **Local spatial autocorrelation (LISA map):** `plot_spatial_autocorrelation(coordinates, lisa.local_I, 'lisa.png', cluster_labels=lisa.cluster_labels)`

## Customisation

The `plot_*` functions take keyword-only options (for example `cmap`, `title`,
`point_size`, `alpha`, `figsize`; the set differs per function, so check the signature),
save the figure to `output_path` and **return** the figure for further tweaking:

```python
fig = plot_spatial_scatter(coordinates, vimentin, 'vimentin.png', title='Vimentin')
ax = fig.axes[0]
ax.set_title('Vimentin', fontsize=14, fontweight='bold')
fig.savefig('vimentin.pdf', bbox_inches='tight')
```

Common kwargs of `plot_spatial_scatter`: `cmap`, `title`, `point_size`, `alpha`,
`figsize`, `vmin`, `vmax`, `colorbar_label`, `coordinate_order`.
