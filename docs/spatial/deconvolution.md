# Deconvolution

Estimate cell‑type proportions within each spatial spot from mixed expression.

## Backend options

| Backend | Strength | Typical runtime (5 k spots × 20 types) |
|---------|----------|----------------------------------------|
| Stereoscope | Probabilistic, accurate, uncertainty estimates | 6 min CPU / 1.5 min CUDA |
| Tangram | Supervised optimal‑transport, fast | 1.5 min CPU / 45 s CUDA |
| Cell2location | Bayesian, hierarchical | Not yet implemented |

Default: **Stereoscope**; switch via `spatial.deconvolution.backend`.

## Stereoscope workflow

```python-snippet
from metainformant.spatial.analysis.deconvolution import stereoscope
# 1. Prepare single‑cell reference (already annotated cell types)
ref_adata = sc.read_h5ad('ref_panel.h5ad')
ref_adata = ref_adata[:, :2000]     # use HVGs only

# 2. Train on reference (once)
from metainformant.spatial.analysis.deconvolution import train_stereoscope
model = train_stereoscope(ref_adata, device='cuda', max_epochs=400)

# 3. Deconvolve spatial data
spatial_adata = load_visium('sample/')
deconv = stereoscope(spatial_adata, model, device='cuda')
spatial_adata.obsm['deconvolution'] = deconv
```

`deconv` is a DataFrame spots × cell‑types, rows sum to ≈1.

## Tangram workflow (supervised)

```python-snippet
from metainformant.spatial.analysis.deconvolution import tangram
spatial_adata.obsm['deconvolution'] = tangram(
    spatial_adata,
    ref_adata,
    device='cpu',
    lambda_density=1.0,
    density_prior='uniform',
)
```

Requires a marker‑gene CSV (`gene,cell_type` rows). If not provided, uses default
canonical markers from cell‑type ontology.

## Output interpretation

```python
import seaborn as sns
proportions = spatial_adata.obsm['deconvolution']
# Visualise one cell type across tissue
spatial_scatter(spatial_adata, color=proportions['T_cell'], cmap='Reds')
```

Threshold: proportion ≥ 0.2 suggests presence; ≥ 0.7 strong enrichment.

## Quality control

- **Correlation with gold standard** (if available): `proportions.corrwith(true_props)`.
- **Spatial coherence** — smoothness of each cell‑type map; run `morans_i()` on the
  proportion vector; I > 0.3 indicates non‑random clustering (expected).
- **Training ELBO** — monotonic decrease; plateau indicates convergence.

## Common failure modes

`ImportError: faiss` — install `faiss-cpu` or `faiss-gpu`.

High memory → reduce `batch_size` or number of cell types (aggregate similar
subtypes).

Very low proportions (all near 0) → check reference gene‑space matches spatial
gene names (Ensembl vs gene symbol mismatch). Use `ref_adata.var['gene_ids'] =
ref_adata.var_names` to harmonise.

## Configuration knobs (relevant)

```yaml
spatial:
  deconvolution:
    backend: stereoscope    # or tangram
    device: cpu             # or cuda
    max_epochs: 400
    lr: 0.001
    batch_size: 256
```
