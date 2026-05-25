# Niche Analysis

Discover and annotate micro‑environments composed of co‑occurring cell types.

## Concept & input

A niche is a spatially local neighbourhood that differs from the global
composition. Input: `adata.obsm['deconvolution']` (cell‑type proportions per
spot). Output: niche label per spot + niche‑specific cell type signatures.

## Algorithm

1. Extract neighbourhood vectors: for each spot i, average deconvolution
   proportions of neighbours within R µm.
2. Dimensionality reduction (PCA → UMAP) of neighbourhood vectors.
3. Cluster in UMAP space (HDBSCAN or Leiden). Result: niche IDs.
4. Compute niche signature: mean deconvolution profile per niche (across spots).

## Usage

```python-snippet
from metainformant.spatial.analysis.niche import discover_niches

niches = discover_niches(
    adata,
    radius=150.0,
    min_cluster_size=20,
    method='hdbscan',          # or 'leiden' after neighbour graph
)
adata.obs['niche'] = niches
```

## Visualisation

```python-snippet
from metainformant.spatial.visualization import plot_niche_summary
plot_niche_summary(adata, top_n_celltypes=5)
```

Creates a stacked bar chart: each niche → top contributing cell types.

## Niche annotation

```python-snippet
from metainformant.spatial.analysis.niche import annotate_niche
annotations = annotate_niche(adata, niche_column='niche')
# Returns dict {niche_id: 'T cell rich', 'Fibroblast core', …}
adata.obs['niche_label'] = adata.obs['niche'].map(annotations)
```

Annotation rules:
- If T_cell proportion > 0.4 and < 0.6 for B cells → 'Lymphoid aggregate'
- If Fibroblast > 0.6 and Epithelial < 0.2 → 'Stromal core'
- … (extensible dictionary in `niche.py`)

## Integration with marker analysis

Compute marker genes *within* each niche:

```python-snippet
from metainformant.spatial.analysis.clustering import find_markers
markers = find_markers(adata, groupby='niche')
```

Compare to domain markers from `spatial_leiden()` to see if niches split
domains functionally further.

## Export

```python
# Niche membership + coordinates for downstream analysis
adata.obs[['x','y','niche','niche_label']].to_csv('niches.tsv', sep='	')
```
