# Getting Started: Spatial Transcriptomics

Spatial transcriptomics captures gene expression while preserving tissue geometry.
This module supports Visium (10x), Xenium, MERFISH, Slide‑seqV2, and CosMx. Key
capabilities: tissue segmentation, cell‑type deconvolution, spatial clustering,
ligand–receptor communication, domain discovery, and batch integration.

## Installation

```bash
uv pip install metainformant[spatial]          # core + visium + xenium
uv pip install metainformant[spatial,deconvolution]  # +stereoscope + tangram
```

## One‑liner demo (Xenium)

```python
from metainformant.spatial import load_xenium
import spatial.analysis as sa

adata = load_xenium('xenium_run/outs/')
# Leiden domains at 0.8 resolution
adata.obs['domain'] = sa.clustering.spatial_leiden(adata, resolution=0.8)
# Plot
sa.visualization.spatial_scatter(adata, color='domain', show=True)
```

## Loader matrix

| Platform | Function | Directory pattern |
|----------|----------|-------------------|
| 10x Visium | `load_visium(path)` | `outs/spatial/` + `filtered_feature_bc_matrix` |
| Xenium | `load_xenium(path)` | `outs/cell_feature_matrix.h5` |
| MERFISH | `load_merfish(path)` | `Vicinity/` CSV per FOV |
| Slide‑seqV2 | `load_slide_seq(path)` | `BeadLocations.csv` + `BeadCounts.h5` |
| CosMx | `load_cosmx(path)` | `exprMat.csv` + `metadata.csv` |

All return AnnData (spots × genes) with `adata.obsm['spatial']` in µm.

## Minimal analysis pipeline

```python
adata = load_visium('sample/')
# QC filter: spots with ≥500 genes, genes in ≥3 spots
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_genes=500)

# Normalise + log1p, highly‑variable genes
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)

# PCA + neighbour graph (KNN on expression + spatial)
sa.neighbors(adata, n_neighbors=15, method='knn')
# Spatial clustering (Leiden with spatial constraint)
sa.clustering.spatial_leiden(adata, resolution=0.8, key_added='domain')
# Domain marker genes
markers = sa.find_markers(adata, groupby='domain')
```

See [ARCHITECTURE.md](ARCHITECTURE.md) for design principles.
