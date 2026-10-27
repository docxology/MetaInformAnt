# Getting Started: Spatial Transcriptomics

Spatial transcriptomics captures gene expression while preserving tissue geometry.
This module loads Visium (10x), Xenium and MERFISH data. Key capabilities: spatially
aware clustering and domain discovery, spatial autocorrelation and spatially variable
genes, cell‑type deconvolution (NNLS or NMF), neighbourhood and niche analysis,
ligand–receptor communication scoring, scRNA‑seq‑to‑spatial mapping, and plotting.

## Installation

There is no `spatial` or `deconvolution` extra. The module needs NumPy, SciPy,
scikit‑learn and matplotlib, which are core dependencies of `metainformant`:

```bash
uv pip install metainformant
```

Optional packages the module uses when they are importable:

- `leidenalg` (with `igraph`) for Leiden clustering and `python-louvain` (the `community`
  module) for Louvain; the `singlecell` extra installs both and the `networks` extra
  installs `python-louvain`. Without them `leiden_clustering` and `louvain_clustering`
  fall back to a built‑in greedy modularity routine, and the result does not say so (see
  [Clustering](clustering.md#which-implementation-runs)).
- `h5py` to read `.h5` count matrices (the `scientific` extra installs it through
  anndata, or `uv pip install h5py`).
- Pillow, which matplotlib installs, to read the Visium tissue image.

There is no GPU or deep‑learning backend.

## One‑liner demo (Xenium)

```python
from metainformant.spatial.io.xenium import load_xenium
from metainformant.spatial.analysis.clustering import spatial_cluster
from metainformant.spatial.visualization.plots import plot_spatial_scatter

ds = load_xenium('xenium_run/outs/')  # XeniumDataset: .expression (cells × genes), .coordinates (n × 2)
# Leiden domains at 0.8 resolution on a joint expression + spatial neighbour graph
result = spatial_cluster(ds.expression, ds.coordinates, resolution=0.8)
# Plot: writes a PNG and returns the matplotlib Figure (string labels give a categorical legend).
# Xenium coordinates are [x, y], so tell the plot which column is which (otherwise the tissue is transposed)
plot_spatial_scatter(
    ds.coordinates, result.labels.astype(str), 'output/spatial/domains.png',
    title='Spatial domains', coordinate_order=ds.coordinate_order,
)
```

## Loader matrix

| Platform | Function | Files read from `path` |
|----------|----------|------------------------|
| 10x Visium | `spatial.io.visium.load_visium(path)` | `filtered_feature_bc_matrix.h5` or `filtered_feature_bc_matrix/` (MEX), `spatial/tissue_positions.csv` (or `tissue_positions_list.csv`); optional `spatial/tissue_hires_image.png`, `spatial/scalefactors_json.json` |
| Xenium | `spatial.io.xenium.load_xenium(path)` | `cell_feature_matrix.h5` or `cell_feature_matrix/` (MEX) and `cells.csv.gz` (or `cells.csv`), both required; optional `transcripts.csv.gz` (read only with `load_transcripts=True`), `cell_boundaries.csv.gz` (only with `load_boundaries=True`) |
| MERFISH | `spatial.io.merfish.load_merfish(path)` | `cell_by_gene.csv` (or `.tsv`, `counts.csv`) and `cell_metadata.csv` (or `.tsv`), both required; optional `detected_transcripts.csv` (read only with `load_transcripts=True`). The delimiter follows the extension: tab for `.tsv`, comma otherwise |

There are no Slide‑seqV2 or CosMx loaders. The Xenium and MERFISH loaders need the cell
positions: a missing `cells.csv(.gz)` or `cell_metadata.csv`/`.tsv`, missing id or coordinate
columns, or a cell of the matrix without a row all raise an error (`FileNotFoundError` or
`ValueError`) instead of giving the cell a placeholder `[0, 0]` position. `load_visium`
raises `ValueError` when no barcode of the matrix has a tissue position. See [io.md](io.md).

The loaders return dataclasses, not AnnData: `SpatialDataset` (Visium), `XeniumDataset`
and `MERFISHDataset`. Each has `.expression` (rows = spots or cells, columns = genes:
SciPy CSR for Visium and Xenium, a dense `float32` array for MERFISH), `.coordinates`
(an `n × 2` float array) and `.gene_names`; Visium keeps spot barcodes in `.barcodes`,
Xenium and MERFISH keep cell IDs in `.cell_ids`. Coordinates are taken as stored in the
files, with no unit conversion: Visium gives full‑resolution image pixels as
`[row, col]` (`.scale_factors` holds the scalefactors JSON), Xenium its
`x_centroid`/`y_centroid`, MERFISH its `center_x`/`center_y`.

The column order differs by platform, and each dataset class declares it in its
`coordinate_order` constant: `"row_col"` for `SpatialDataset` (column 0 is the vertical
image axis) and `"xy"` for `XeniumDataset` and `MERFISHDataset` (column 0 is x). The
plotting functions that take a dataset (`plot_gene_expression_map`, `plot_cell_type_map`)
read it from the dataset. The ones that take a bare coordinate array default to
`"row_col"`, so pass `coordinate_order=ds.coordinate_order` (or `"xy"`) for Xenium and
MERFISH arrays; see [visualization.md](visualization.md#coordinate-order).

## Minimal analysis pipeline

The module's functions take plain arrays; the QC and normalisation below use scanpy
(`scientific` or `singlecell` extra) on an AnnData built from the loaded dataset.

```python
import anndata as ad
import pandas as pd
import scanpy as sc
from metainformant.spatial.io.visium import load_visium
from metainformant.spatial.analysis.clustering import spatial_cluster

ds = load_visium('sample/')  # SpatialDataset
adata = ad.AnnData(ds.expression, obs=pd.DataFrame(index=ds.barcodes), var=pd.DataFrame(index=ds.gene_names))
adata.obsm['spatial'] = ds.coordinates

# QC filter: spots with ≥500 genes, genes in ≥3 spots
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_genes=500)

# Normalise + log1p, highly‑variable genes
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
hvg = adata[:, adata.var['highly_variable']]

# Spatial clustering: PCA + KNN graph on expression, combined with a KNN graph on the
# coordinates (n_neighbors=6, spatial_weight=0.5 by default), then Leiden
result = spatial_cluster(hvg.X, hvg.obsm['spatial'], resolution=0.8)
adata.obs['domain'] = pd.Categorical(result.labels.astype(str))
# Domain marker genes
sc.tl.rank_genes_groups(adata, groupby='domain')
```

`spatial_cluster` converts the expression matrix to dense, so pass the highly variable
genes rather than the full matrix (see [PERFORMANCE.md](PERFORMANCE.md)).

See [ARCHITECTURE.md](ARCHITECTURE.md) for design principles.
