# Deconvolution

Estimate cell‑type proportions within each spatial spot from mixed expression.

## Backend options

| `method` | Strength | Notes |
|----------|----------|-------|
| `'nnls'` | Fast, robust; one non‑negative least‑squares fit per spot against the reference profiles | Default |
| `'nmf'` | Same fixed‑reference fit with optional L2 regularisation (`alpha`) | Requires scikit‑learn |

Default: **nnls**, selected via `deconvolve_spots(..., method='nnls')`. Probabilistic
and optimal‑transport backends (Stereoscope, Tangram, Cell2location) are not
implemented.

## Reference‑based workflow

```python-snippet
from metainformant.spatial.analysis.deconvolution import create_reference_profiles, deconvolve_spots

# 1. Build cell-type reference profiles from an annotated single-cell reference
profiles, cell_type_names, ref_genes = create_reference_profiles(
    ref_expression,              # cells × genes matrix
    ref_cell_types,              # one label per cell
    gene_names=ref_gene_names,
)

# 2. Deconvolve the spatial data (shared genes are intersected automatically)
deconv = deconvolve_spots(
    spatial_expression,          # spots × genes matrix
    profiles,
    method='nnls',
    cell_type_names=cell_type_names,
    spatial_gene_names=spatial_gene_names,
    reference_gene_names=ref_genes,
)
fractions = deconv.fractions     # spots × cell types
```

`deconv` is a `DeconvolutionResult`: `weights` (raw non‑negative weights), `fractions`
(spots × cell types, rows sum to ≈1), `cell_type_names` and per‑spot `residuals`.

## Output interpretation

```python
import pandas as pd
from metainformant.spatial.visualization.plots import plot_spatial_scatter

proportions = pd.DataFrame(deconv.fractions, columns=deconv.cell_type_names)
# Visualise one cell type across tissue
plot_spatial_scatter(coordinates, proportions['T_cell'].values, 'T_cell.png', cmap='Reds', title='T_cell')
```

Threshold: proportion ≥ 0.2 suggests presence; ≥ 0.7 strong enrichment.

## Quality control

- **Correlation with gold standard** (if available): `proportions.corrwith(true_props)`.
- **Spatial coherence** — smoothness of each cell‑type map; run `morans_i()` on the
  proportion vector; I > 0.3 indicates non‑random clustering (expected).
- **Fit residuals** — `deconv.residuals` is the per‑spot fitting error; spots with
  high residuals are poorly explained by the reference profiles.

## Common failure modes

`ValueError: No shared genes between spatial and reference data` — the gene names
passed as `spatial_gene_names` and `reference_gene_names` do not overlap (usually an
Ensembl vs gene‑symbol mismatch); harmonise the identifiers first.

Very low proportions (all near 0) → check reference gene‑space matches spatial
gene names (Ensembl vs gene symbol mismatch). Use `ref_adata.var['gene_ids'] =
ref_adata.var_names` to harmonise.

## Configuration knobs (relevant)

`deconvolve_spots()` takes its settings as arguments:

```python
deconvolve_spots(
    spatial_expression,
    profiles,
    method='nmf',   # 'nnls' or 'nmf'
    alpha=0.1,      # L2 regularisation strength (nmf only)
)
```
