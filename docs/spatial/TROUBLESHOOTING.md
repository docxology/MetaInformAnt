# Troubleshooting: Spatial

## Common errors & fixes

| Symptom | Likely cause | Fix |
|---------|--------------|-----|
| `FileNotFoundError: spatial/scalefactors_json.json` | Visium directory incomplete | Run Space Ranger 2.0+; ensure `outs/` structure intact |
| `ValueError: platform 'xenium' unmatched` | Wrong loader function | Use `load_xenium()` not `load_visium()` |
| `KeyError: 'spatial'` (adata.obsm) | Loader didn't set coords | Verify `adata.uns['platform']` matches loader; re‑run loader |
| `MemoryError` on 100 k spots | Radius graph O(N²) | Switch to `method='knn'` or subset spots |
| `RuntimeError: CUDA out of memory` (Stereoscope) | Batch size too high | Reduce `deconvolution.batch_size` to 128 or 64 |
| `ImportError: no 'faiss' module` | Optional deps missing | `uv pip install metainformant[spatial,deconvolution]` |
| `ValueError: no valid spots after QC` | Thresholds too strict | Lower `min_genes` / `min_cells` or inspect `adata.obs['pct_counts_mt']` |

## Step‑by‑step sanity check

```python
from metainformant.spatial import load_visium
adata = load_visium('sample/')
print(adata)                     # check shape: n_obs × n_vars
print(adata.obsm['spatial'][:5]) # first 5 coordinates
print(adata.uns['platform'])     # expected 'visium'
```

If `.obsm['spatial']` missing, the loader failed; check directory structure:

```
sample/
 └─ outs/
     ├─ spatial/
     │   └─ scalefactors_json.json
     └─ filtered_feature_bc_matrix/
         ├─ barcodes.tsv.gz
         ├─ features.tsv.gz
         └─ matrix.mtx.gz
```

## Deconvolution convergence failure

**Symptom** — Stereoscope loss plateaus at high value or NaNs appear.

**Causes & fixes:**

1. **Reference single‑cell not cleaned** — filter low‑quality cells & normalize before
   training: `sc.pp.filter_cells(ref, min_genes=200); ref = ref[:, :1000]`.
2. **Mismatched gene names** — ensure both AnnData objects use identical gene ID
   space (e.g., Ensembl v105). Map with `adata.var_names = adata.var['gene_id']`.
3. **Batch effects** — integrate reference and spatial first with Harmony before
   deconvolution to align expression distributions.

## Deeper diagnostics

```python
import logging
logging.getLogger('metainformant.spatial.analysis.deconvolution').setLevel('DEBUG')
```

Console will print per‑epoch ELBO, gradient norms, and convergence flags.

## Benchmarking your hardware

```python-snippet
from metainformant.spatial.performance import benchmark
benchmark.run(platform='visium', n_spots=5000, n_genes=20000, n_jobs=1)
```

Outputs a markdown table similar to this page.

## Getting help

Include: platform (OS, Python, metainformant version), loader used, dataset size
(spots × genes), full traceback, and a 1 MB sample H5 file if possible.
