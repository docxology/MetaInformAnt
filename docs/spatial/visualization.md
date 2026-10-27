# Visualization

Plotting functions live in `metainformant.spatial.visualization.plots`:
`plot_spatial_scatter` (a quick one‑liner for any per‑spot values) and more
specialised `plot_*` functions. Every function takes spot coordinates (or a
`SpatialDataset`), saves the figure to `output_path` and returns the matplotlib
`Figure`.

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
`SpatialDataset`, e.g. from `load_visium()`). Polygon outlines around each domain are
not provided.

## Expression maps

```python-snippet
from metainformant.spatial.visualization.plots import plot_gene_expression_map
for gene in ['EPCAM', 'VIM', 'CD3D']:
    plot_gene_expression_map(dataset, gene, f'{gene}.png', cmap='Reds')
```

`dataset` is a `SpatialDataset` (or any object with `.expression`, `.coordinates` and
`.gene_names`). Each call writes one figure; there is no multi‑panel layout.

## QC plots

Any per‑spot quantity can be mapped with `plot_spatial_scatter`:

- **Library size per spot:** `plot_spatial_scatter(coordinates, expression.sum(axis=1), 'library_size.png')`
- **Local spatial autocorrelation (LISA map):** `plot_spatial_autocorrelation(coordinates, lisa.local_I, 'lisa.png', cluster_labels=lisa.cluster_labels)`

## Customisation

All `plot_*` functions take keyword arguments (`cmap`, `title`, `point_size`, `alpha`,
`figsize`, ...), save the figure to `output_path` and **return** the figure for
further tweaking:

```python
fig = plot_spatial_scatter(coordinates, vimentin, 'vimentin.png', title='Vimentin')
ax = fig.axes[0]
ax.set_title('Vimentin', fontsize=14, fontweight='bold')
fig.savefig('vimentin.pdf', bbox_inches='tight')
```

Common kwargs of `plot_spatial_scatter`: `cmap`, `title`, `point_size`, `alpha`,
`figsize`, `vmin`, `vmax`, `colorbar_label`.
