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
- **Spatial coherence** — smoothness of each cell‑type map; run `morans_i(proportion_vector, weights)`
  with a `spatial_weights_matrix` of the coordinates (see [Autocorrelation](autocorrelation.md));
  a clearly positive I with a small p‑value indicates non‑random clustering (expected).
- **Fit residuals** — `deconv.residuals` is the per‑spot fitting error; spots with
  high residuals are poorly explained by the reference profiles.
- **Marker agreement** — `validate_deconvolution` (see
  [Advanced module](#advanced-module-spatialdeconvolution)) correlates each type's estimated
  proportion with the mean expression of its marker genes across spots; it takes the
  dictionary result of `spatial.deconvolution.spatial_deconvolution.deconvolve_spots`.

## Common failure modes

`ValueError: No shared genes between spatial and reference data` — the gene names
passed as `spatial_gene_names` and `reference_gene_names` do not overlap (usually an
Ensembl vs gene‑symbol mismatch); harmonise the identifiers first.

Very low proportions (all near 0) → check reference gene‑space matches spatial
gene names (Ensembl vs gene symbol mismatch). Use `ref_adata.var['gene_ids'] =
ref_adata.var_names` to harmonise.

## Advanced module: `spatial.deconvolution`

`metainformant.spatial.deconvolution.spatial_deconvolution` offers a second, dictionary-based interface with reference profiles built from marker genes, a cell-type map with mixing entropy, a validation summary and a niche finder:

```python-snippet
from metainformant.spatial.deconvolution.spatial_deconvolution import (
    build_reference_profiles,
    deconvolve_spots,
    niche_identification,
    spatial_cell_type_mapping,
    validate_deconvolution,
)

# 1. Marker-based reference profiles from a cells x genes matrix and one label per cell
reference = build_reference_profiles(sc_expression, sc_cell_types, n_markers=50)
genes = reference['gene_indices']          # selected marker genes, as columns of sc_expression

# 2. Deconvolve the spots on those genes (the spatial matrix must share the column order of sc_expression)
result = deconvolve_spots(spatial_counts[:, genes], reference['profiles'], method='nnls')
proportions = result['proportions_matrix']   # spots x types; columns follow result['cell_types']

# 3. Map, validate and group into niches ((x, y) pairs for the coordinates)
mapping = spatial_cell_type_mapping(proportions, coordinates, result['cell_types'])
report = validate_deconvolution(
    result,
    spatial_markers={'T_cell': ['CD3D', 'CD3E'], 'B_cell': ['MS4A1', 'CD79A']},
    expression=spatial_counts,         # the full spots x genes matrix the proportions came from
    gene_names=spatial_gene_names,     # one name per column of that matrix
)
print(report['marker_correlations'])   # {cell type: correlation across spots}
niches = niche_identification(proportions, coordinates, n_niches=5, seed=42)
```

- `deconvolve_spots` here takes the profile dictionary (`method` is `'nnls'` or `'regression'`) and returns a dictionary with `proportions_matrix`, `cell_types`, `confidence_scores` and `spots_summary`. Cell types are sorted alphabetically.
- `spatial_cell_type_mapping(proportions, coordinates, cell_types)` needs all three arguments and returns `spatial_map` (per spot: `x`, `y`, `dominant_type`, `proportions`, `entropy`), `dominant_type_per_spot` and `mixing_entropy`.
- `validate_deconvolution(result, spatial_markers=None, *, expression=None, gene_names=None, correlation_method='pearson')` checks that the proportions are normalised and summarises the confidence and the type diversity (`normalization_check`, `mean_confidence`, `confidence_distribution`, `type_diversity`, `overall_validity`). With `spatial_markers`, a dictionary from cell type to marker gene names, it also correlates the estimated proportion of each type with the mean expression of its markers, across spots:
    - It needs the spots × genes matrix the proportions were estimated from (`expression=`, in the same spot order, used as supplied: normalise it first if the markers should be compared on a normalised scale) and one name per column (`gene_names=`, or a dataset or DataFrame that carries them; explicit names win). Pass the full matrix, not the marker-gene subset used for the deconvolution, because markers are matched to the names exactly and a marker that is not among them is ignored.
    - `marker_correlations` maps each type to the Pearson correlation, or the Spearman rank correlation with `correlation_method='spearman'` (which needs scipy). `marker_correlation_method` records the method (`None` when no markers were checked) and `marker_genes_used` the markers that were found.
    - A type that cannot be correlated is not reported as a number: it is left out of `marker_correlations` and listed with its reason in `marker_correlations_skipped` (the type is not in the result, none of its markers is among the gene names, or its proportion or its marker expression is constant across spots).
    - `ValueError` is raised for `spatial_markers` without `expression`, for markers without gene names, for an unknown `correlation_method`, for a matrix without one row per spot and one column per name, and for a used marker that names more than one column.
- `niche_identification(proportions, coordinates, n_niches=5, *, seed=None)` clusters the proportions together with scaled coordinates using k-means and returns `niche_labels`, `niche_compositions` and `spatial_coherence`. With the default `seed=None` the k-means initialisation comes from a fresh generator seeded from operating-system entropy, so the niches can differ between calls; numpy's global random state is neither used nor advanced. Pass an integer for reproducible niches.

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
