# Visualization

Two‑level API: quick one‑liners (`spatial_scatter`) and fine‑grained custom
figure assembly (`plot_*` functions).

## Scatter / spot plot

```python
from metainformant.spatial.visualization import spatial_scatter
spatial_scatter(
    adata,
    color='EPCAM',                 # gene name or obs column
    cmap='viridis',
    size=adata.obs['area'] / 1000,   # spot size proportional to area
    vmin=0, vmax=10,
    title='Epithelial marker',
    show=False,
    save='epcam.pdf',
)
```

## Domain outlines

```python
from metainformant.spatial.visualization import domain_outlines
domain_outlines(
    adata,
    domain_column='domain',
    linewidth=0.7,
    palette='tab20',
    show=True,
)
```

Overlays polygon boundaries around each cluster on the tissue image (if H&E exists,
pass `img=adata.uns['tissue_image']`).

## Expression overlay

```python
from metainformant.spatial.visualization import expression_heatmap
expression_heatmap(
    adata,
    genes=['EPCAM','VIM','CD3D'],
    ncols=3,
    standard_scale='var',
)
```

A multi‑panel figure with shared colour bar.

## QC plots

- **Library size per spot:** `plot_library_size(adata)`
- **Mitochondrial %:** `plot_mt_fraction(adata)`
- **Spatial autocorrelation heatmap:** `plot_moran_summary(adata, genes_of_interest)`

## Interactive (Plotly)

```python
from metainformant.spatial.visualization.interactive import spatial_scatter_plotly
fig = spatial_scatter_plotly(
    adata,
    color='domain',
    hover_data=['n_genes_by_counts', 'pct_counts_mt'],
)
fig.write_html('spatial_interactive.html')
```

Hover shows spot barcode and metadata.

## Customisation

All `plot_*` functions accept Matplotlib `Axes` via `ax=` and **return** the figure
and axes for further tweaking:

```python
fig, ax = spatial_scatter(adata, color='VIM', show=False)
ax.set_title('Vimentin', fontsize=14, fontweight='bold')
fig.savefig('vimentin.pdf', bbox_inches='tight')
```

Common kwargs: `show`, `save`, `dpi`, `palette`, `alpha`, `size` (scalar or array).
