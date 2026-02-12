# Cell Type Annotation

The celltyping module assigns cell type identities to single cells using marker
gene panels, transfers labels from annotated reference datasets, detects novel
cell populations, and analyzes cell type composition across experimental groups.

## Concepts

### Marker-Based Annotation

`annotate_by_markers` scores every cell against a dictionary of known marker gene
sets. Three scoring strategies are available:

- **Overlap**: Enrichment of expressed marker genes relative to the background
  expression rate across all genes.
- **Correlation**: Pearson correlation between the cell's expression at marker
  positions and an ideal all-ones profile.
- **Scoring**: Mean expression of marker genes minus a random background,
  highlighting marker-specific enrichment.

Each cell receives the label of the highest-scoring type. Cells where the best
score falls below `threshold` are labeled `"unassigned"`. Ambiguous cells (top
two scores within 10% of each other) are flagged separately.

### Label Transfer

`transfer_labels` projects reference and query expression matrices into a shared
PCA space and applies k-nearest-neighbor classification. Prediction confidence
equals the fraction of neighbors that agree on the assigned label.

### Novel Type Detection

`find_novel_types` identifies cells that do not confidently match any known
type. When marker genes are provided, cells are re-scored and those below the
confidence threshold are flagged. Novel cells are grouped with a simple k-means
clustering for summary reporting.

### Composition Analysis

`cell_type_composition` computes the proportion of each cell type per group or
sample, with Wilson score confidence intervals for the binomial proportions.

## Function Reference

### `annotate_by_markers(expression_matrix, marker_genes, gene_names, method, threshold) -> dict`

Assign cell types from marker gene sets. Returns `cell_labels`,
`confidence_scores`, `ambiguous_cells`, and `marker_stats`.

| Parameter          | Type                   | Default      | Description                        |
|--------------------|------------------------|--------------|------------------------------------|
| `expression_matrix`| Any                    | --           | Cells x genes matrix               |
| `marker_genes`     | `dict[str, list[str]]` | --           | Cell type to marker gene lists     |
| `gene_names`       | `list[str] \| None`    | `None`       | Gene names for columns             |
| `method`           | `str`                  | `"overlap"`  | `"overlap"`, `"correlation"`, or `"scoring"` |
| `threshold`        | `float`                | `0.0`        | Min score for confident assignment |

### `score_cell_type(expression, gene_names, marker_set, n_background, seed) -> float`

Score a single cell against one marker gene set. Returns marker mean minus
background mean.

| Parameter      | Type         | Default | Description                            |
|----------------|--------------|---------|----------------------------------------|
| `expression`   | `list[float]`| --      | Expression values for one cell         |
| `gene_names`   | `list[str]`  | --      | Gene names matching expression values  |
| `marker_set`   | `list[str]`  | --      | Marker genes for the target type       |
| `n_background` | `int`        | `50`    | Number of background genes to sample   |
| `seed`         | `int \| None`| `None`  | Random seed for reproducibility        |

### `transfer_labels(reference_data, reference_labels, query_data, n_neighbors, n_components) -> dict`

Transfer labels via kNN in PCA space. Returns `predicted_labels`,
`prediction_scores`, and `mapping_quality`.

| Parameter         | Type         | Default | Description                          |
|-------------------|--------------|---------|--------------------------------------|
| `reference_data`  | Any          | --      | Reference expression matrix          |
| `reference_labels`| `list[str]`  | --      | Labels for reference cells           |
| `query_data`      | Any          | --      | Query expression matrix (same genes) |
| `n_neighbors`     | `int`        | `30`    | Neighbors for kNN classifier         |
| `n_components`    | `int`        | `50`    | PCA components before kNN            |

### `find_novel_types(expression_matrix, known_labels, marker_genes, gene_names, threshold) -> dict`

Flag cells that do not match known types. Returns `novel_cell_indices`,
`novel_cell_scores`, `n_novel`, `fraction_novel`, and `cluster_summary`.

| Parameter          | Type                        | Default | Description                    |
|--------------------|-----------------------------|---------|--------------------------------|
| `expression_matrix`| Any                         | --      | Cells x genes matrix           |
| `known_labels`     | `list[str]`                 | --      | Current labels per cell        |
| `marker_genes`     | `dict[str, list[str]] \| None` | `None`  | Marker sets for re-scoring  |
| `gene_names`       | `list[str] \| None`         | `None`  | Gene names for columns         |
| `threshold`        | `float`                     | `0.5`   | Score cutoff for novelty       |

### `cell_type_composition(labels, groups, confidence_level) -> dict`

Compute proportions per group with Wilson confidence intervals. Returns
`overall`, `per_group`, `confidence_intervals`, `n_cells`, `n_types`, `n_groups`.

| Parameter          | Type              | Default | Description                       |
|--------------------|-------------------|---------|-----------------------------------|
| `labels`           | `list[str]`       | --      | Cell type label per cell          |
| `groups`           | `list[str] \| None`| `None` | Group/sample label per cell       |
| `confidence_level` | `float`           | `0.95`  | Confidence level for intervals    |

## Code Examples

```python
from metainformant.singlecell.celltyping.annotation import (
    annotate_by_markers, transfer_labels, find_novel_types,
    cell_type_composition, score_cell_type,
)

# Marker-based annotation
matrix = [[2.5, 0.1, 0.0], [0.0, 3.1, 0.2], [0.3, 0.0, 4.0]]
markers = {"T_cell": ["CD3D"], "B_cell": ["CD19"], "NK": ["NKG7"]}
genes = ["CD3D", "CD19", "NKG7"]
result = annotate_by_markers(matrix, markers, gene_names=genes, method="overlap")
print(result["cell_labels"])  # ["T_cell", "B_cell", "NK"]

# Transfer labels from reference to query
transfer = transfer_labels(
    [[2.5, 0.1], [0.0, 3.1]], ["T_cell", "B_cell"],
    [[2.0, 0.3]], n_neighbors=2, n_components=2,
)
print(transfer["predicted_labels"])

# Detect novel types and analyze composition
novel = find_novel_types(matrix, result["cell_labels"], markers, genes, threshold=0.5)
comp = cell_type_composition(
    ["T_cell", "B_cell", "T_cell", "NK"], groups=["ctrl", "ctrl", "treat", "treat"],
)
print(comp["per_group"])
```

## Configuration

All parameters are function arguments. No external config files required.

## Import Path

```python
from metainformant.singlecell.celltyping.annotation import annotate_by_markers
from metainformant.singlecell.celltyping.annotation import transfer_labels
from metainformant.singlecell.celltyping.annotation import find_novel_types
from metainformant.singlecell.celltyping.annotation import cell_type_composition
from metainformant.singlecell.celltyping.annotation import score_cell_type
```
