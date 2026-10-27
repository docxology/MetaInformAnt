# Integration

Spatial-scRNA-seq integration providing methods for mapping single-cell RNA-seq data to spatial transcriptomics, including label transfer, correlation-based mapping, and gene imputation.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports scrna_mapping submodule |
| `scrna_mapping.py` | Label transfer, correlation mapping, anchor-based transfer, gene imputation |

## Key Functions

| Function | Description |
|----------|-------------|
| `scrna_mapping.correlation_mapping()` | Map scRNA-seq cell types to spatial spots via correlation |
| `scrna_mapping.anchor_based_transfer()` | Transfer labels using mutual nearest neighbor anchors (`anchors` is a required positional argument; pass `None` to compute them) |
| `scrna_mapping.map_scrna_to_spatial()` | End-to-end scRNA-seq to spatial mapping pipeline |
| `scrna_mapping.impute_spatial_genes()` | Impute unmeasured genes in spatial data using scRNA-seq (needs `spatial_genes` and `scrna_genes`) |

## Usage

```python
from metainformant.spatial.integration import scrna_mapping

result = scrna_mapping.map_scrna_to_spatial(
    scrna_expression, spatial_expression, scrna_labels=scrna_cell_types
)
imputed = scrna_mapping.impute_spatial_genes(
    spatial_expression, scrna_expression, genes_to_impute,
    spatial_genes=spatial_gene_names, scrna_genes=scrna_gene_names,
)
```
