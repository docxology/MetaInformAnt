# Single-Cell Cell Typing

Marker-based cell type annotation and classification for single-cell data. Supports overlap scoring, correlation-based assignment, label transfer from reference datasets via kNN, novel cell type discovery, and cell type composition analysis.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Re-exports `annotation` module |
| `annotation.py` | Marker-based annotation, label transfer, novel type detection, composition |

## Key Functions

| Function | Description |
|----------|-------------|
| `annotate_by_markers()` | Assign cell types using marker gene sets and scoring methods |
| `score_cell_type()` | Score individual cells against a specific cell type signature |
| `transfer_labels()` | Transfer cell type labels from reference via kNN in PCA space |
| `find_novel_types()` | Identify cells that do not match any known cell type |
| `cell_type_composition()` | Compute cell type proportions across samples or groups |

## Usage

```python
from metainformant.singlecell.celltyping import annotation

labels = annotation.annotate_by_markers(expression_matrix, marker_genes, gene_names)
transferred = annotation.transfer_labels(reference, query, ref_labels)
composition = annotation.cell_type_composition(labels, sample_ids)
```
