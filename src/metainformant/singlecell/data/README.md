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
| `normalize_counts()` | Library size normalization (also `median`, `size_factors` methods) |
| `log_transform()` | Log1p transformation with configurable base |
| `scale_data()` | Z-score standardization with optional clipping |
| `identify_highly_variable_genes()` | Select HVGs for downstream analysis |
| `harmony_integration()` | Harmony batch correction on PCA space |
| `bbknn_integration()` | Batch-balanced kNN graph construction |
| `scanorama_integration()` | Scanorama batch integration across datasets |
| `combat_integration()` | ComBat-style batch correction on a `batch_key` |
| `mnn_integration()` | Mutual-nearest-neighbor integration of multiple datasets |
| `integrate_multiple_batches()` | Dispatcher for all integration methods |
| `evaluate_integration_quality()` | Batch mixing and distance-ratio metrics |

## Usage

```python
from metainformant.singlecell.data.preprocessing import load_count_matrix, filter_cells
from metainformant.singlecell.data.integration import harmony_integration

data = load_count_matrix("data/counts.h5ad")
data = filter_cells(data, min_genes=200, max_pct_mt=20.0)
data = harmony_integration(data, batch_key="sample")
```
