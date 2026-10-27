# Single-Cell Preprocessing

The preprocessing module (`metainformant.singlecell.data.preprocessing`) provides the foundational tools for loading, cleaning, and preparing single-cell RNA-seq data for analysis.

## SingleCellData Class

The `SingleCellData` class is the central data structure for single-cell analysis, providing an AnnData-like interface:

```python
class SingleCellData:
    """
    Container for single-cell expression data and metadata.

    Constructor: SingleCellData(X, obs=None, var=None, uns=None, obsm=None,
                                varm=None, obsp=None, varp=None, layers=None)

    Attributes:
        X: Expression matrix (cells × genes), can be sparse or dense
        obs: Cell metadata (DataFrame, one row per cell; empty frame by default)
        var: Gene metadata (DataFrame, one row per gene; empty frame by default)
        uns: Unstructured annotations (dict)
        obsm: Multi-dimensional cell annotations (dict of arrays)
        varm: Multi-dimensional gene annotations (dict of arrays)
        obsp: Pairwise cell annotations (dict)
        varp: Pairwise gene annotations (dict)
        layers: Alternative matrix representations (dict)
        n_obs: Number of cells
        n_vars: Number of genes
        shape: (n_obs, n_vars)

    Methods: copy(), to_df()
    """
```

### Initialization

```python-snippet
from metainformant.singlecell.data.preprocessing import SingleCellData
import numpy as np
import pandas as pd

# Create from arrays
X = np.random.poisson(2, (100, 1000))  # 100 cells, 1000 genes
obs = pd.DataFrame({'cell_type': ['TypeA'] * 50 + ['TypeB'] * 50})
var = pd.DataFrame({'gene_name': [f'Gene_{i}' for i in range(1000)]})

data = SingleCellData(X=X, obs=obs, var=var)
```

`obs` and `var` must have one row per cell and per gene respectively, otherwise the constructor raises `ValueError`.

`data.copy()` returns an independent object: the matrix, `obs`, `var`, `obsm`, `varm`, `obsp`, `varp` and `layers` are copied, and `uns` is copied recursively, so editing a nested entry such as `copy.uns["qc_summary"]` never changes `data`. The preprocessing functions below return new objects in the same way and leave their input untouched.

## Data Loading

### load_count_matrix()

Load count matrices from various file formats:

```python-snippet
from metainformant.singlecell.data.preprocessing import load_count_matrix

# From CSV with cells as rows and genes as columns (first column = cell names)
data = load_count_matrix("counts.csv", format="csv")

# From TSV with genes as rows and cells as columns (transposed on load)
data = load_count_matrix("counts.tsv", format="tsv", transpose=True)

# From Matrix Market format (genes × cells), with gene and barcode name files
data = load_count_matrix(
    "matrix.mtx", format="mtx", genes_file="genes.txt", barcodes_file="barcodes.txt"
)

# From AnnData H5AD (the default format; requires anndata)
data = load_count_matrix("counts.h5ad")
```

**Parameters:**
- `filepath`: Path to the count matrix file
- `format`: File format ("h5ad", "csv", "tsv", "mtx"); the default is "h5ad", so pass `format="csv"` for text files
- `transpose`: For CSV/TSV only. Set to `True` when the file has genes as rows and cells as columns, so that the result has cells as rows (default: False)
- `genes_file`, `barcodes_file`: For Matrix Market only, both required. Text files with one gene name and one cell barcode per line (the first comma-separated column is read)
- Any other keyword arguments are passed to the reader (`pandas.read_csv` for CSV/TSV, `anndata.read_h5ad` for H5AD) and are ignored for Matrix Market. The CSV/TSV delimiter is fixed by `format`, so there is no `sep` argument

**Returns:** `SingleCellData` object with loaded expression matrix (cells × genes). For CSV, TSV and Matrix Market input `X` is a scipy CSR sparse matrix when scipy is installed.

**Errors:** A missing file or an unsupported format raises `ValidationError`, a parse failure raises `IOError`, and H5AD without `anndata` installed raises `ConfigError`. These are the exception classes of `metainformant.core.utils.errors` (subclasses of `METAINFORMANTError`), not the Python built-ins.

## Quality Control

### calculate_qc_metrics()

Calculate standard quality control metrics for cells:

```python-snippet
from metainformant.singlecell.data.preprocessing import calculate_qc_metrics

# Calculate QC metrics
data = calculate_qc_metrics(data)

# Access calculated metrics
print(data.obs.columns)  # ['total_counts', 'n_genes', 'pct_mt', 'pct_ribo']
```

The function returns a copy of the data; the input is not modified.

**Calculated metrics (columns of `data.obs`):**
- `total_counts`: Total UMI/read counts per cell
- `n_genes`: Number of detected genes per cell
- `pct_mt`: Percentage of mitochondrial gene expression
- `pct_ribo`: Percentage of ribosomal gene expression

It also adds per-gene columns to `data.var` (`n_cells`, `total_counts`, `mean_expression`, `pct_dropout`) and a summary dictionary in `data.uns["qc_summary"]`.

**Mitochondrial genes:** Identified by gene names (the `data.var` index) starting with "MT-", "MT." or "MT_", case-insensitively (so "mt-" and "Mt-" also match)
**Ribosomal genes:** Identified by gene names starting with "RPS", "RPL", "MRPS" or "MRPL", case-insensitively

## Filtering

### filter_cells()

Filter cells based on QC metrics:

```python-snippet
from metainformant.singlecell.data.preprocessing import filter_cells

# Filter cells with basic thresholds
data_filtered = filter_cells(
    data,
    min_genes=200,      # Minimum genes detected
    max_genes=5000,     # Maximum genes detected
    min_counts=1000,    # Minimum total counts
    max_counts=50000,   # Maximum total counts
    max_pct_mt=20,      # Maximum mitochondrial percentage
)
```

**Parameters:**
- `min_genes`: Minimum number of genes per cell
- `max_genes`: Maximum number of genes per cell
- `min_counts`: Minimum total counts per cell
- `max_counts`: Maximum total counts per cell
- `max_pct_mt`: Maximum mitochondrial gene percentage (`max_mito_percent` is an alias)

There is no ribosomal-percentage filter; apply one yourself on the `pct_ribo` column if needed. When `data.obs` has no `total_counts` column (the only column it checks), `filter_cells` runs `calculate_qc_metrics` first and logs a warning.

**What the filtered object holds:**

- The matrix, `obs`, `obsm` and `layers` keep only the retained cells, and `obsp` keeps their rows and columns. `var`, `varm` and `varp` are carried over.
- The per-gene QC columns already present in `data.var` (`n_cells`, `total_counts`, `mean_expression`, `pct_dropout`) are recomputed on the retained cells, so a later `filter_genes` reads current values (when no cell is left, `mean_expression` and `pct_dropout` are `NaN`). Other `var` annotations are carried over unchanged.
- `uns` is copied recursively, so the input object is never modified. The cell-indexed `uns["neighbors"]` is not updated: it still describes the cells before filtering. Recompute the graph with `compute_neighbors` afterwards (diffusion pseudotime raises an error on a graph with the wrong number of cells).
- `uns["qc_summary"]`, when present, gains `cells_after_filtering` and `cells_filtered`.

### filter_genes()

Filter genes based on detection across cells:

```python-snippet
from metainformant.singlecell.data.preprocessing import filter_genes

# Filter rarely expressed genes
data_filtered = filter_genes(
    data,
    min_cells=3,        # Gene must be expressed in at least 3 cells
    min_counts=None     # Optional minimum total counts per gene
)
```

**Parameters:** `min_cells`, `max_cells`, `min_counts` and `max_counts` (all optional). They are compared with the per-gene `n_cells` and `total_counts` columns of `data.var` written by `calculate_qc_metrics`; when `data.var` has no `n_cells` column, `filter_genes` runs `calculate_qc_metrics` first and logs a warning. `filter_cells` recomputes those columns on the retained cells, so after `filter_cells` the gene thresholds refer to the remaining cells without another `calculate_qc_metrics` call.

The result keeps only the retained genes in the matrix, `var`, `varm` and `layers`; `uns` is copied recursively, and `uns["qc_summary"]`, when present, gains `genes_after_filtering` and `genes_filtered`.

## Normalization

### normalize_counts()

Normalize count data to account for differences in sequencing depth:

```python-snippet
from metainformant.singlecell.data.preprocessing import normalize_counts

# Total count normalization (most common)
data_norm = normalize_counts(
    data,
    target_sum=1e4,           # Scale to 10,000 counts per cell
    method='total_count'      # Normalization method
)

# Median normalization
data_norm = normalize_counts(data, method='median')

# Median-of-ratios size factors
data_norm = normalize_counts(data, method='size_factors')
```

**Methods:**
- `'total_count'` (alias `'total'`): Scale each cell to `target_sum` total counts. If `target_sum` is omitted it defaults to the median library size of the cells with non-zero counts
- `'median'`: Scale each cell to the median total count of the cells with non-zero counts (`target_sum` is ignored)
- `'size_factors'`: Divide each cell by a simplified median-of-ratios size factor, then multiply by `target_sum` (default: the median of the positive size factors)

An unknown method raises `ValueError`; passing precomputed scaling factors is not supported. The result holds a dense array in `X` even for sparse input, and `data_norm.uns["normalization"]` records `method`, `target_sum` and `size_factors`.

**The size factors are divisors, the same way for every method.** `uns["normalization"]["size_factors"]` holds one value per cell: the library size (the row sum of the counts) for `'total_count'`/`'total'` and `'median'`, and the median-of-ratios size factor for `'size_factors'`. Whichever method was used,

```
X_normalized == X / size_factors[:, None] * target_sum
```

with `target_sum` read from `uns["normalization"]["target_sum"]`. Cells whose size factor is zero (for example empty cells) are set to all zeros. The array is stored as a plain list of floats, one per cell; use `numpy.asarray` to divide by it.

### log_transform()

Apply log transformation to normalized counts:

```python-snippet
from metainformant.singlecell.data.preprocessing import log_transform

# Log2 transformation: log2(x + 1)
data_log = log_transform(data, base=2)

# Natural log transformation (the default base)
data_log = log_transform(data)

# Log10 transformation
data_log = log_transform(data, base=10)
```

**Parameters:**
- `base`: Logarithm base, a positive number other than 1 (default: `numpy.e`). The string `'e'` is not accepted

The transformation is always `log(x + 1)`: the pseudocount of 1 is fixed and cannot be changed. The result holds a dense array, and `data_log.uns["transformation"]` records the `type`, `base` and `pseudocount`.

## Scaling

### scale_data()

Z-score standardization of expression data:

```python-snippet
from metainformant.singlecell.data.preprocessing import scale_data

# Standard z-score scaling (mean=0, std=1)
data_scaled = scale_data(data)

# Divide by the standard deviation without subtracting the mean
data_uncentered = scale_data(data, zero_center=False)

# Z-score scaling, clipping values to [-10, 10]
data_clipped = scale_data(data, max_value=10)
```

**Parameters:**
- `zero_center`: Whether to center each gene (subtract its mean; default: True)
- `max_value`: Clip scaled values to `[-max_value, max_value]` (default: `None`, no clipping)

Every gene is always divided by its standard deviation (genes with zero variance are left unscaled); there is no option to skip scaling and no gene-subset argument, so all genes are scaled. The result holds a dense array, and `data_scaled.uns["scaling"]` records `zero_centered`, `max_value`, `gene_means` and `gene_stds`.

## Complete Preprocessing Pipeline

Here's a typical preprocessing workflow:

```python-snippet
from metainformant.singlecell.data.preprocessing import (
    load_count_matrix, calculate_qc_metrics, filter_cells,
    filter_genes, normalize_counts, log_transform, scale_data
)

# 1. Load data (a CSV with genes as rows and cells as columns, hence transpose=True)
data = load_count_matrix("raw_counts.csv", format="csv", transpose=True)
print(f"Loaded: {data.n_obs} cells × {data.n_vars} genes")

# 2. Calculate QC metrics
data = calculate_qc_metrics(data)

# 3. Filter cells and genes
print("Before filtering:")
print(f"  Cells: {data.n_obs}")
print(f"  Genes: {data.n_vars}")

data = filter_cells(data, min_genes=200, max_genes=5000, max_pct_mt=20)  # also refreshes the per-gene QC columns
data = filter_genes(data, min_cells=3)

print("After filtering:")
print(f"  Cells: {data.n_obs}")
print(f"  Genes: {data.n_vars}")

# 4. Normalize and transform
data = normalize_counts(data, target_sum=1e4, method='total_count')
data = log_transform(data, base=2)

# 5. Scale data (optional: compute_pca standardizes the data itself by default)
data = scale_data(data)

print("Preprocessing complete!")
```

## Data Format Compatibility

### Input Formats
- **CSV/TSV files**: Standard text files with a header row and row names in the first column
- **Matrix Market**: `.mtx` files (sparse format), with the gene and barcode names in separate text files
- **H5AD**: AnnData HDF5 files (requires `anndata`)
- **Gzip compressed**: CSV/TSV and Matrix Market files with a `.gz` extension are decompressed by pandas and scipy

### Expected Structure
- **Cells as rows, genes as columns** (default for CSV/TSV)
- **Genes as rows, cells as columns** (CSV/TSV with `transpose=True`; Matrix Market files are always read as genes × cells)
- **Row names**: The first column becomes the index (cell names, or gene names with `transpose=True`)
- **Column names**: Taken from the header row (gene names, or cell names with `transpose=True`)

## Memory Management

### Sparse Matrix Support
The module accepts sparse matrices, and CSV/TSV/Matrix Market input is loaded as a CSR matrix when scipy is installed:

```python
import scipy.sparse as sparse

# Create sparse matrix
X_sparse = sparse.csr_matrix(X_dense)
data = SingleCellData(X=X_sparse, obs=obs, var=var)

# Check if matrix is sparse
print(f"Matrix is sparse: {sparse.issparse(data.X)}")
```

### Memory-Efficient Operations
- `filter_cells` and `filter_genes` preserve sparsity
- `calculate_qc_metrics`, `normalize_counts`, `log_transform`, `scale_data` and the downstream HVG and PCA functions convert to a dense array internally, so plan for a dense cells × genes matrix; `normalize_counts`, `log_transform` and `scale_data` return dense results and never modify the input in place
- CSV, TSV and Matrix Market files are loaded as CSR matrices; H5AD keeps whatever matrix type is stored in the file

## Error Handling

Common issues and solutions:

### File Loading Issues
```python
from metainformant.core.utils.errors import METAINFORMANTError
from metainformant.singlecell.data.preprocessing import load_count_matrix

try:
    data = load_count_matrix("counts.csv", format="csv")
except METAINFORMANTError as e:
    # ValidationError (missing file, unsupported format), IOError (parse failure)
    # or ConfigError (H5AD without anndata)
    print(f"Could not load count matrix: {e}")
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

1. **Use sparse matrices** for large, sparse datasets (they help while loading and filtering; later steps densify)
2. **Filter early** to reduce computational burden
3. **Batch processing** for very large datasets
4. **Memory monitoring** for datasets > 50k cells

## Testing

Preprocessing tests are available in `tests/singlecell/test_singlecell_preprocessing.py`:

```bash
# Run preprocessing tests
uv run pytest tests/singlecell/test_singlecell_preprocessing.py -v
```

## Related Functions

- [Dimensionality Reduction](./dimensionality.md): For HVG selection and PCA
- [Quality Control](../quality/fastq.md): For raw sequencing QC
- [Core I/O](../core/io.md): For file handling utilities
