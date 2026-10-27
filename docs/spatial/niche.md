# Niche Analysis

Discover and annotate micro‑environments composed of co‑occurring cell types.

## Concept & input

A niche is a spatially local neighbourhood that differs from the global
composition. Input: cell‑type proportions per spot (for example `deconv.fractions`
from `deconvolve_spots()`) and the spot coordinates. Output: niche label per spot +
niche‑specific cell type signatures.

## Algorithm

1. Smooth the cell‑type proportions of each spot over its `n_neighbors` nearest
   spatial neighbours (`spatial_weight` sets the blend: 0 = no smoothing, 1 = full).
2. Cluster the smoothed proportion vectors with K‑means into `n_niches` groups.
   Result: niche IDs.
3. Compute niche signature: mean cell‑type proportions per niche (across spots),
   plus niche size and Shannon diversity.

## Usage

```python-snippet
from metainformant.spatial.niche.identification import identify_niches

niche = identify_niches(
    proportions,               # spots × cell types, rows sum to 1
    coordinates,               # spots × 2
    cell_type_names=cell_type_names,
    n_niches=5,
    n_neighbors=15,
    spatial_weight=0.3,
    random_state=42,
)
adata.obs['niche'] = niche.niche_labels
```

`niche.compositions` (niches × cell types), `niche.niche_sizes` and
`niche.niche_diversity` describe each niche.

## Visualisation

```python-snippet
from metainformant.spatial.visualization.plots import plot_spatial_scatter
plot_spatial_scatter(coordinates, niche.niche_labels.astype(str), 'niches.png', title='Niches')
```

Draws the niche labels on the tissue coordinates (one colour per niche). The
per‑niche composition in `niche.compositions` can be plotted with matplotlib, for
example as a stacked bar chart of the top contributing cell types.

## Integration with marker analysis

Compute marker genes *within* each niche:

```python-snippet
from metainformant.singlecell.analysis.clustering import find_markers
markers = find_markers(expression, niche.niche_labels, n_markers=10, gene_names=gene_names)
```

Compare to domain markers from `spatial_cluster()` to see if niches split
domains functionally further.

## Export

```python
# Niche membership + coordinates for downstream analysis
adata.obs[['x','y','niche']].to_csv('niches.tsv', sep='	')
```
