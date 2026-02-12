# Single-Cell Data

Data loading, preprocessing, QC filtering, normalization, and batch integration for single-cell RNA-seq datasets.

## Contents

| File | Purpose |
|------|---------|
| `integration.py` | Batch correction: BBKNN, Harmony, Scanorama, MNN, ComBat |
| `preprocessing.py` | SingleCellData class, QC metrics, filtering, normalization, HVG selection |

## Key Classes and Functions

| Symbol | Description |
|--------|-------------|
| `SingleCellData` | Core data container for count matrix, metadata, and embeddings |
| `load_count_matrix()` | Load from h5ad, CSV, or MTX format |
| `calculate_qc_metrics()` | Compute mitochondrial %, gene counts, total UMIs |
| `filter_cells()` | Remove low-quality cells by QC thresholds |
| `filter_genes()` | Remove lowly expressed genes |
| `normalize_counts()` | Library size normalization |
| `identify_highly_variable_genes()` | Select HVGs for downstream analysis |
| `harmony_integration()` | Harmony batch correction on PCA space |
| `bbknn_integration()` | Batch-balanced kNN graph construction |
| `scanorama_integration()` | Scanorama batch integration across datasets |

## Usage

```python
from metainformant.singlecell.data.preprocessing import load_count_matrix, filter_cells
from metainformant.singlecell.data.integration import harmony_integration

data = load_count_matrix("data/counts.h5ad")
data = filter_cells(data, min_genes=200, max_mito_pct=20.0)
data = harmony_integration(data, batch_key="sample")
```
