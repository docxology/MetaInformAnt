# Multi-Sample Integration

Batch-correct and spatially align multiple tissue sections.

## Current API surface

`metainformant.spatial.integration` currently provides single-cell-to-spatial
mapping utilities in `scrna_mapping`:

```python
from metainformant.spatial.integration.scrna_mapping import (
    correlation_mapping,
    anchor_based_transfer,
    map_scrna_to_spatial,
    impute_spatial_genes,
)
```

## Label transfer from a reference

Project cell-type labels from a scRNA-seq reference onto a spatial section:

```python
transfer = anchor_based_transfer(
    ref_adata=scrna_ref,
    query_adata=spatial_section,
    label_key="cell_type",
)
spatial_section.obs["predicted_cell_type"] = transfer["labels"]
```

## Gene imputation into spatial coordinates

```python
imputed = impute_spatial_genes(
    scrna_ref=scrna_ref,
    spatial=spatial_section,
    genes=["Apoe", "Gfap"],
)
```

## Not yet implemented

General multi-sample batch correction (`spatial_batch_correct` with
BBKNN/Harmony/Scanorama), inter-section coordinate registration, and Visium
downsampling are specified but not implemented in this module. Do not use
them in workflows until the corresponding functions land; batch correction
for single-cell data lives in `metainformant.singlecell`.
