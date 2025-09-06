# Single-Cell Preprocessing

The preprocessing module provides the foundational tools for loading, cleaning, and preparing single-cell RNA-seq data for analysis.

## SingleCellData Class

The `SingleCellData` class is the central data structure for single-cell analysis, providing an AnnData-like interface:

```python
class SingleCellData:
    """
    Container for single-cell expression data and metadata.
    
    Attributes:
        X: Expression matrix (cells × genes), can be sparse or dense
        obs: Cell metadata (DataFrame with cell_id as index)
        var: Gene metadata (DataFrame with gene_id as index)
        obsm: Multi-dimensional cell annotations (dict of arrays)
        varm: Multi-dimensional gene annotations (dict of arrays)
        uns: Unstructured annotations (dict)
        n_obs: Number of cells
        n_vars: Number of genes
    """
```

### Initialization

```python
from metainformant.singlecell.preprocessing import SingleCellData
import numpy as np
import pandas as pd

# Create from arrays
X = np.random.poisson(2, (100, 1000))  # 100 cells, 1000 genes
obs = pd.DataFrame({'cell_type': ['TypeA'] * 50 + ['TypeB'] * 50})
var = pd.DataFrame({'gene_name': [f'Gene_{i}' for i in range(1000)]})

data = SingleCellData(X=X, obs=obs, var=var)
```

## Data Loading

### load_count_matrix()

Load count matrices from various file formats:

```python
from metainformant.singlecell.preprocessing import load_count_matrix

# From CSV (genes as rows, cells as columns)
data = load_count_matrix("counts.csv", format="csv")

# From TSV with transposition (cells as rows, genes as columns)
data = load_count_matrix("counts.tsv", format="tsv", transpose=True)

# From Matrix Market format
data = load_count_matrix("matrix.mtx", format="mtx")

# From HDF5 (requires h5py)
data = load_count_matrix("counts.h5", format="h5", h5_key="matrix")
```

**Parameters:**
- `file_path`: Path to the count matrix file
- `format`: File format ("csv", "tsv", "mtx", "h5")
- `transpose`: Whether to transpose the matrix (default: False)
- `h5_key`: HDF5 dataset key (for H5 format)
- `sep`: Delimiter for CSV/TSV files

**Returns:** `SingleCellData` object with loaded expression matrix

## Quality Control

### calculate_qc_metrics()

Calculate standard quality control metrics for cells:

```python
from metainformant.singlecell.preprocessing import calculate_qc_metrics

# Calculate QC metrics
data = calculate_qc_metrics(data)

# Access calculated metrics
print(data.obs.columns)  # ['total_counts', 'n_genes', 'pct_mt', 'pct_ribo']
```

**Calculated metrics:**
- `total_counts`: Total UMI/read counts per cell
- `n_genes`: Number of detected genes per cell
- `pct_mt`: Percentage of mitochondrial gene expression
- `pct_ribo`: Percentage of ribosomal gene expression

**Mitochondrial genes:** Identified by gene names starting with "MT-", "mt-", or "Mt-"
**Ribosomal genes:** Identified by gene names starting with "RPS", "RPL", "rps", or "rpl"

## Filtering

### filter_cells()

Filter cells based on QC metrics:

```python
from metainformant.singlecell.preprocessing import filter_cells

# Filter cells with basic thresholds
data_filtered = filter_cells(
    data,
    min_genes=200,      # Minimum genes detected
    max_genes=5000,     # Maximum genes detected
    min_counts=1000,    # Minimum total counts
    max_counts=50000,   # Maximum total counts
    max_pct_mt=20,      # Maximum mitochondrial percentage
    max_pct_ribo=None   # Maximum ribosomal percentage (optional)
)
```

**Parameters:**
- `min_genes`: Minimum number of genes per cell
- `max_genes`: Maximum number of genes per cell
- `min_counts`: Minimum total counts per cell
- `max_counts`: Maximum total counts per cell
- `max_pct_mt`: Maximum mitochondrial gene percentage
- `max_pct_ribo`: Maximum ribosomal gene percentage

### filter_genes()

Filter genes based on detection across cells:

```python
from metainformant.singlecell.preprocessing import filter_genes

# Filter rarely expressed genes
data_filtered = filter_genes(
    data,
    min_cells=3,        # Gene must be expressed in at least 3 cells
    min_counts=None     # Optional minimum total counts per gene
)
```

## Normalization

### normalize_counts()

Normalize count data to account for differences in sequencing depth:

```python
from metainformant.singlecell.preprocessing import normalize_counts

# Total count normalization (most common)
data_norm = normalize_counts(
    data,
    target_sum=1e4,           # Scale to 10,000 counts per cell
    method='total_count'      # Normalization method
)

# Median normalization
data_norm = normalize_counts(data, method='median')

# Custom scaling factors
scaling_factors = data.obs['total_counts'].median() / data.obs['total_counts']
data_norm = normalize_counts(data, scaling_factors=scaling_factors)
```

**Methods:**
- `'total_count'`: Scale each cell to `target_sum` total counts
- `'median'`: Scale each cell to the median total count across all cells
- Custom: Provide scaling factors directly

### log_transform()

Apply log transformation to normalized counts:

```python
from metainformant.singlecell.preprocessing import log_transform

# Log2 transformation with pseudocount
data_log = log_transform(data, base=2, pseudocount=1)

# Natural log transformation
data_log = log_transform(data, base='e')

# Log10 transformation
data_log = log_transform(data, base=10)
```

**Parameters:**
- `base`: Logarithm base (2, 10, 'e')
- `pseudocount`: Value added before log transformation

## Scaling

### scale_data()

Z-score standardization of expression data:

```python
from metainformant.singlecell.preprocessing import scale_data

# Standard z-score scaling (mean=0, std=1)
data_scaled = scale_data(data, center=True, scale=True)

# Center only (mean=0)
data_centered = scale_data(data, center=True, scale=False)

# Scale to specific genes only (e.g., highly variable genes)
hvg_indices = data.var['highly_variable'].values
data_scaled = scale_data(data, gene_subset=hvg_indices)
```

**Parameters:**
- `center`: Whether to center data (subtract mean)
- `scale`: Whether to scale data (divide by standard deviation)
- `gene_subset`: Boolean array to select specific genes
- `max_value`: Clip scaled values to this maximum (default: 10)

## Complete Preprocessing Pipeline

Here's a typical preprocessing workflow:

```python
from metainformant.singlecell.preprocessing import (
    load_count_matrix, calculate_qc_metrics, filter_cells, 
    filter_genes, normalize_counts, log_transform, scale_data
)

# 1. Load data
data = load_count_matrix("raw_counts.csv", transpose=True)
print(f"Loaded: {data.n_obs} cells × {data.n_vars} genes")

# 2. Calculate QC metrics
data = calculate_qc_metrics(data)

# 3. Filter cells and genes
print("Before filtering:")
print(f"  Cells: {data.n_obs}")
print(f"  Genes: {data.n_vars}")

data = filter_cells(data, min_genes=200, max_genes=5000, max_pct_mt=20)
data = filter_genes(data, min_cells=3)

print("After filtering:")
print(f"  Cells: {data.n_obs}")
print(f"  Genes: {data.n_vars}")

# 4. Normalize and transform
data = normalize_counts(data, target_sum=1e4, method='total_count')
data = log_transform(data, base=2)

# 5. Scale data (typically done after HVG selection)
data = scale_data(data)

print("Preprocessing complete!")
```

## Data Format Compatibility

### Input Formats
- **CSV/TSV files**: Standard text files with headers
- **Matrix Market**: `.mtx` files (sparse format)
- **HDF5**: Binary format for large datasets
- **Gzip compressed**: Automatic detection and decompression

### Expected Structure
- **Genes as rows, cells as columns** (default)
- **Cells as rows, genes as columns** (use `transpose=True`)
- **Gene names**: Should be in first column or as row indices
- **Cell names**: Should be in first row or as column headers

## Memory Management

### Sparse Matrix Support
The module automatically handles sparse matrices to reduce memory usage:

```python
import scipy.sparse as sparse

# Create sparse matrix
X_sparse = sparse.csr_matrix(X_dense)
data = SingleCellData(X=X_sparse, obs=obs, var=var)

# Check if matrix is sparse
print(f"Matrix is sparse: {sparse.issparse(data.X)}")
```

### Memory-Efficient Operations
- Filtering operations preserve sparsity
- Normalization can work in-place to save memory
- Large datasets automatically use sparse representations

## Error Handling

Common issues and solutions:

### File Loading Issues
```python
try:
    data = load_count_matrix("counts.csv")
except FileNotFoundError:
    print("Count matrix file not found")
except ValueError as e:
    print(f"Invalid file format: {e}")
```

### Filtering Edge Cases
```python
# Check for empty results after filtering
data_filtered = filter_cells(data, min_genes=200)
if data_filtered.n_obs == 0:
    print("Warning: No cells passed filtering criteria")
```

## Integration with Core Utilities

The preprocessing module integrates with METAINFORMANT's core utilities:

```python
from metainformant.core.io import ensure_directory

# Use core I/O utilities for output paths
output_dir = ensure_directory("output/singlecell_analysis")
```

## Performance Tips

1. **Use sparse matrices** for large, sparse datasets
2. **Filter early** to reduce computational burden
3. **Batch processing** for very large datasets
4. **Memory monitoring** for datasets > 50k cells

## Testing

Preprocessing tests are available in `tests/test_singlecell_preprocessing.py`:

```bash
# Run preprocessing tests
uv run pytest tests/test_singlecell_preprocessing.py -v
```

## Related Functions

- [Dimensionality Reduction](./dimensionality.md): For HVG selection and PCA
- [Quality Control](../quality/fastq.md): For raw sequencing QC
- [Core I/O](../core/io.md): For file handling utilities
