# RNA Deconvolution

Cell type deconvolution from bulk RNA-seq expression data. Estimates cell type proportions using reference expression signatures via NNLS and SVR-based methods, with tools for signature matrix construction, marker gene selection, and result validation.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Re-exports `bulk_deconvolution` module |
| `bulk_deconvolution.py` | NNLS/SVR deconvolution, signature building, marker selection, validation |

## Key Functions

| Function | Description |
|----------|-------------|
| `deconvolve_nnls()` | Non-negative least squares deconvolution |
| `deconvolve_svr()` | SVR-based deconvolution (CIBERSORT-style) |
| `build_signature_matrix()` | Build cell-type signature matrix from reference profiles |
| `select_marker_genes()` | Select informative marker genes for deconvolution |
| `validate_deconvolution()` | Assess deconvolution accuracy against known proportions |
| `batch_deconvolve()` | Deconvolve multiple bulk samples in batch |

## Usage

```python
from metainformant.rna.deconvolution import bulk_deconvolution

result = bulk_deconvolution.deconvolve_nnls(mixture, signature_matrix)
sig = bulk_deconvolution.build_signature_matrix(reference_profiles, cell_types)
```
