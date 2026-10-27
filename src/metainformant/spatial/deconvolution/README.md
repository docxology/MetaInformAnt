# Deconvolution

Advanced spatial deconvolution for cell type proportion estimation, extending core analysis with reference profile building, spatial cell type mapping, validation, and tissue niche identification.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports spatial_deconvolution submodule |
| `spatial_deconvolution.py` | NNLS deconvolution, reference profiling, validation, niche detection |

## Key Functions

| Function | Description |
|----------|-------------|
| `spatial_deconvolution.deconvolve_spots()` | Deconvolve spots into cell type proportions (NNLS/regression); returns a dictionary |
| `spatial_deconvolution.build_reference_profiles()` | Build cell type expression profiles from a scRNA-seq reference (`profiles`, `gene_indices`, ...) |
| `spatial_deconvolution.spatial_cell_type_mapping()` | Map estimated cell type fractions to spatial coordinates (needs `proportions`, `coordinates` and `cell_types`) |
| `spatial_deconvolution.validate_deconvolution()` | Validate deconvolution: normalisation, confidence and type diversity checks. With `spatial_markers`, `expression` (the spots × genes matrix the proportions came from) and gene names it also returns, per cell type, the Pearson (or Spearman, `correlation_method="spearman"`) correlation across spots between the estimated proportion and the mean expression of its markers (`marker_correlations`); types that cannot be correlated are listed with a reason in `marker_correlations_skipped` |
| `spatial_deconvolution.niche_identification()` | Identify tissue niches by clustering composition vectors; `seed=None` (the default) can give different niches on every call |

## Usage

```python
from metainformant.spatial.deconvolution import spatial_deconvolution

ref = spatial_deconvolution.build_reference_profiles(scrna_data, cell_types)
genes = ref["gene_indices"]          # marker genes, as columns of scrna_data
result = spatial_deconvolution.deconvolve_spots(spatial_counts[:, genes], ref["profiles"])
proportions = result["proportions_matrix"]
mapping = spatial_deconvolution.spatial_cell_type_mapping(proportions, coordinates, result["cell_types"])
report = spatial_deconvolution.validate_deconvolution(
    result, {"T_cell": ["CD3D", "CD3E"]}, expression=spatial_counts, gene_names=spatial_gene_names
)
niches = spatial_deconvolution.niche_identification(proportions, coordinates, n_niches=5, seed=42)
```
