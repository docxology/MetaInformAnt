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
| `spatial_deconvolution.deconvolve_spots()` | Deconvolve spots into cell type proportions (NNLS/regression) |
| `spatial_deconvolution.build_reference_profiles()` | Build cell type expression profiles from scRNA-seq reference |
| `spatial_deconvolution.spatial_cell_type_mapping()` | Map estimated cell type fractions to spatial coordinates |
| `spatial_deconvolution.validate_deconvolution()` | Validate deconvolution via marker gene agreement |
| `spatial_deconvolution.niche_identification()` | Identify tissue niches by clustering composition vectors |

## Usage

```python
from metainformant.spatial.deconvolution import spatial_deconvolution

ref = spatial_deconvolution.build_reference_profiles(scrna_data, cell_types)
result = spatial_deconvolution.deconvolve_spots(spatial_counts, ref)
mapping = spatial_deconvolution.spatial_cell_type_mapping(result, coordinates)
niches = spatial_deconvolution.niche_identification(result, coordinates, n_niches=5)
```
