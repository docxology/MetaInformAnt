# Data Loading: Spatial

Unified loaders for every major platform.

## Common return value

All functions return an **AnnData** object with these guaranteed slots:

| Attribute / key | Meaning |
|------------------|---------|
| `.X` | Cell/spot × gene sparse count matrix |
| `.obs` | DataFrame per-observation metadata including `x`, `y` (µm) |
| `.obsm['spatial']` | N×2 float array of µm coordinates |
| `.var` | DataFrame per-gene metadata (gene symbols, Ensembl IDs) |
| `.uns['platform']` | String: `'visium'`, `'xenium'`, `'merfish'`, `'slideseq'`, `'cosmx'` |
| `.uns['scale']` | µm per pixel (Visium only; others 1.0) |
| `.uns['sample_id']` | Derived from folder name (configurable) |

## 10x Visium (`load_visium`)

```python
from metainformant.spatial.io.visium import load_visium
adata = load_visium(
    folder='data/visium_sample/',
    min_counts=100,
    max_counts=30_000,
    filter_genes=True,
    include_tissue=True,      # keep only spots inside tissue mask
)
```

**Expected folder layout**:

```
visium_sample/
├── filtered_feature_bc_matrix.h5
├── spatial/
│   ├── tissue_lowres_image.png
│   ├── tissue_hires_image.png
│   └── scalefactors_json.json
```

The parser reads:
- `spatial/scalefactors_json.json` → spot-to-pixel scaling factors
- `spatial/tissue_positions_list.csv` → barcode→(x,y) mapping
- HDF5 gene–barcode matrix (`filtered_feature_bc_matrix.h5`)

Coordinates are transformed from pixel-space to µm using the `tissue_hires_scalef`
and stored in `.obsm['spatial']` directly.

**Troubleshooting**:
- `ValueError: missing scalefactors` → folder is not a Space Ranger output root
- All spots have `in_tissue=False` — maybe your data uses a different mask column
  name; set `tissue_mask_col='mask'`


## Xenium (`load_xenium`)

```python
from metainformant.spatial.io.xenium import load_xenium
adata = load_xenium(
    folder='data/xenium_run/',
    quality='high',          # 'high' removes low-quality cells
    use_cyto=None,           # defaults to cytosolic counts if present
)
```

**Folder structure**:

```
xenium_run/
├── cell_feature_matrix.h5
├── cells.csv.genome
├── transcripts/
│   ├── transcripts.csv.genome
│   └── ...
```

Counts are **cell-by-gene** (not spot-based). Coordinates are in µm, already
subcellular. The loader optionally filters by QC flags present in `cells.csv`.

**Performance**: ~4.8 s for 50 k cells (HDF5 read-bound).


## MERFISH (`load_merfish`)

```python
from metainformant.spatial.io.merfish import load_merfish
adata = load_merfish(
    folder='data/merfish_aio/',
    format='vicinity_aoi',    # also 'nanowells' (deprecated)
    min_transcripts=10,
)
```

Loader parses `detected_transcripts.csv` (columns: `barcode`, `x`, `y`, `target`).
Aggregates transcripts → spots (typically cell segmentation already applied).


## Slide-seqV2 (`load_slide_seq`)

```python
from metainformant.spatial.io.slideseq import load_slide_seq
adata = load_slide_seq(
    folder='data/slideseq_v2/',
    bead_locations='bead_locations.csv',
    counts='raw_feature_bc_matrix.h5',
)
```

Coordinates stored directly in µm. Spots are beads; each bead contains 1–10 cells,
so downstream deconvolution is usually unnecessary.


## CosMx / Nanostring (`load_cosmx`)

```python
from metainformant.spatial.io.cosmx import load_cosmx
adata = load_cosmx(
    folder='data/cosmx_lung/',
    expr_matrix='exprMat.csv',
    metadata='metadata.csv',
    fov='FOV1',
)
```

`exprMat.csv` is a dense table: rows = molecules, columns = gene, cell, FOV, x, y.
Loader aggregates to cell-level counts within a single FOV. Multi-FOV data should
be concatenated after loading all FOVs.

## Registry

All loaders are registered in `io/__init__.py` → `_PLATFORM_LOADERS`. To add a new
platform:

```python
# io/mypatform.py
def load_myplatform(folder: str) -> AnnData:
    # implement parsing logic
    return adata

# then in io/__init__.py
from .myplatform import load_myplatform
_PLATFORM_LOADERS['myplatform'] = load_myplatform
```

Now `load_spatial(folder, platform='myplatform')` will auto-detect.

## Tips

- Always check `adata.obs['in_tissue'].mean()` — aim > 0.7; if lower consider raising
  `min_counts`/`min_genes` filters.
- Gene symbol harmonization: loaders keep the source gene IDs (often Ensembl). To
  convert to HGNC symbols, use `adata.var_names = map_ens_to_hgnc(adata.var_names)`.
  Built-in mapper in `metainformant.core.ids` — `convert_gene_ids()`.
