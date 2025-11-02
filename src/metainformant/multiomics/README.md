# Multi-Omics Integration Module

The `multiomics` module provides tools for integrating and analyzing data from multiple omic layers including genomics, transcriptomics, proteomics, metabolomics, and epigenomics.

## Overview

This module enables systems-level biological analysis by combining data from different molecular layers to understand complex biological processes. It handles sample alignment, data harmonization, and provides joint analysis methods for multi-omic datasets.

## Key Components

### Data Integration (`integration.py`)
Framework for combining heterogeneous biological datasets with automatic sample alignment.

**Key Features:**
- Automatic sample alignment across omics layers
- Support for multiple file formats (CSV, TSV, Excel)
- Sample and feature ID mapping
- Metadata integration

**Usage:**
```python
from metainformant.multiomics import MultiOmicsData, integrate_omics_data
import pandas as pd

# Create multi-omics data object
genomics = pd.DataFrame(np.random.randn(10, 100), index=[f"S{i}" for i in range(10)])
transcriptomics = pd.DataFrame(np.random.randn(10, 200), index=[f"S{i}" for i in range(10)])

omics_data = MultiOmicsData(genomics=genomics, transcriptomics=transcriptomics)
print(f"Aligned samples: {omics_data.n_samples}")
print(f"Available layers: {omics_data.layer_names}")

# Or integrate from file paths (paths must point to CSV/TSV files)
data_dict = {
    "genomics": "path/to/genomics.csv",
    "transcriptomics": "path/to/expression.tsv"
}
integrated = integrate_omics_data(data_dict)
```

### Joint Analysis Methods

#### Joint Principal Component Analysis (`joint_pca`)
Performs PCA on concatenated features from multiple omics layers.

```python
from metainformant.multiomics import joint_pca

# Joint PCA with optional layer weights
embeddings, loadings, variance = joint_pca(
    omics_data,
    n_components=50,
    layer_weights={"genomics": 1.0, "transcriptomics": 2.0},
    standardize=True
)

print(f"Joint embeddings shape: {embeddings.shape}")
print(f"Explained variance: {variance[:5]}")
```

#### Joint Non-negative Matrix Factorization (`joint_nmf`)
Factorizes each omics layer with shared sample factors.

```python
from metainformant.multiomics import joint_nmf

# Joint NMF with regularization
W, H = joint_nmf(
    omics_data,
    n_components=20,
    max_iter=200,
    regularization=0.01,
    random_state=42
)

# W: shared sample factors (n_samples x n_components)
# H: layer-specific feature factors
print(f"Sample factors shape: {W.shape}")
print(f"Genomics feature factors: {H['genomics'].shape}")
```

#### Canonical Correlation Analysis (`canonical_correlation`)
Finds maximally correlated patterns between two omics layers.

```python
from metainformant.multiomics import canonical_correlation

# CCA between genomics and transcriptomics
X_c, Y_c, X_w, Y_w, correlations = canonical_correlation(
    omics_data,
    layer_pair=("genomics", "transcriptomics"),
    n_components=10,
    regularization=0.01
)

print(f"Canonical correlations: {correlations}")
print(f"First canonical variate correlation: {correlations[0]}")
```

**Key Features:**
- Multi-omic dimensionality reduction
- Cross-omic pattern discovery
- Integrative biomarker identification
- Sample stratification across omics

## Integration with Other Modules

### With DNA, RNA, and Protein Modules
```python
from metainformant.multiomics import integrate_omics_data, MultiOmicsData
import pandas as pd

# Prepare omics data as DataFrames
genomics_df = pd.DataFrame(...)  # DNA/genomic data
transcriptomics_df = pd.DataFrame(...)  # RNA expression data
proteomics_df = pd.DataFrame(...)  # Protein data

# Integrate using MultiOmicsData or integrate_omics_data
omics_data = MultiOmicsData(
    genomics=genomics_df,
    transcriptomics=transcriptomics_df,
    proteomics=proteomics_df
)
```

## Performance Features

- Scalable integration algorithms
- Memory-efficient data structures
- Parallel processing for large datasets

## Testing

Comprehensive tests cover:
- Integration algorithm correctness
- Multi-omic data compatibility
- Performance with large datasets

## Dependencies

- Statistical packages for integrative analysis
- Optional: specialized multi-omic tools

This module enables comprehensive systems biology analysis through multi-omic data integration.
