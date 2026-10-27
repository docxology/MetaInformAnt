# Multi-Sample Integration

Map single-cell RNA-seq references onto spatial sections: label transfer and gene
imputation. Multi-sample batch correction and section registration are not
implemented yet (see the last section).

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

The functions take plain matrices (spots × genes and cells × genes, dense or SciPy sparse) and
optional lists of gene names; they do not read AnnData objects. They return `MappingResult`
(`predicted_labels`, `prediction_scores`, `label_probabilities`, `cell_type_names`, `method`,
`metadata`) or `ImputationResult` (`imputed_expression`, `gene_names`, `confidence`, `method`).

## Label transfer from a reference

Project cell-type labels from a scRNA-seq reference onto a spatial section:

```python
from metainformant.spatial.integration.scrna_mapping import anchor_based_transfer

transfer = anchor_based_transfer(
    spatial_expression,             # spots x genes
    scrna_expression,               # cells x genes
    None,                           # anchors: None finds mutual-nearest-neighbour anchors itself
    scrna_labels=scrna_cell_types,  # one label per cell (required)
    spatial_genes=spatial_gene_names,
    scrna_genes=scrna_gene_names,   # the shared genes are intersected when both name lists are given
    n_pcs=30,
    k_anchor=5,
    seed=42,
)
predicted = transfer.predicted_labels            # one label per spot
probabilities = transfer.label_probabilities     # spots x types; columns follow transfer.cell_type_names
```

`anchors` is a required positional argument: pass `None` to compute the anchors, or an array of
`(spatial_index, scrna_index)` pairs. `scrna_labels` is required for the label transfer.

`correlation_mapping(spatial_expression, scrna_profiles, *, cell_type_names=None,
correlation_method='pearson')` assigns each spot the cell type whose mean profile (`scrna_profiles`
has one row per type) it correlates with most (`'pearson'` or `'spearman'`).
`map_scrna_to_spatial(scrna_data, spatial_data, method='correlation', *, scrna_labels=...)` builds
those profiles from labelled cells and runs `correlation_mapping` (or, with `method='anchor'`,
`anchor_based_transfer`) in one call.

## Gene imputation into spatial coordinates

```python
from metainformant.spatial.integration.scrna_mapping import impute_spatial_genes

imputed = impute_spatial_genes(
    spatial_expression,             # spots x spatial genes
    scrna_expression,               # cells x scRNA-seq genes
    ['Apoe', 'Gfap'],               # genes to impute; they must be among scrna_gene_names
    spatial_genes=spatial_gene_names,
    scrna_genes=scrna_gene_names,   # both name lists are required
    n_neighbors=10,
    seed=42,
)
imputed.gene_names                  # the requested genes that were found
imputed.imputed_expression          # spots x imputed genes
```

Each spot receives the inverse-distance-weighted mean expression of its `n_neighbors` nearest
scRNA-seq cells in a joint PCA space built from the genes the two datasets share.

## Not yet implemented

General multi-sample batch correction (`spatial_batch_correct` with
BBKNN/Harmony/Scanorama), inter-section coordinate registration, and Visium
downsampling are specified but not implemented in this module. Do not use
them in workflows until the corresponding functions land; batch correction
for single-cell data lives in `metainformant.singlecell`.
