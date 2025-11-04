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
from metainformant.dna import variants
from metainformant.rna import workflow
from metainformant.protein import parse_fasta
import pandas as pd

# Prepare omics data from domain modules
# Genomics: variant data from DNA module
genomic_variants = variants.parse_vcf("variants.vcf")
genomics_df = pd.DataFrame(genomic_variants)

# Transcriptomics: expression data from RNA module
expression_data = workflow.extract_expression("expression.tsv")
transcriptomics_df = pd.DataFrame(expression_data)

# Proteomics: protein abundance from protein module
proteins = parse_fasta(Path("proteome.fasta"))
proteomics_df = pd.DataFrame(proteins)

# Integrate all omics layers
omics_data = MultiOmicsData(
    genomics=genomics_df,
    transcriptomics=transcriptomics_df,
    proteomics=proteomics_df
)
```

### With Single-Cell Module
```python
from metainformant.multiomics import joint_pca, canonical_correlation
from metainformant.singlecell import compute_pca, load_count_matrix

# Integrate single-cell data with bulk omics
sc_counts = load_count_matrix("singlecell_counts.h5ad")
bulk_expression = pd.read_csv("bulk_expression.csv", index_col=0)

# Joint analysis across omics layers
omics_data = MultiOmicsData(
    singlecell=sc_counts,
    transcriptomics=bulk_expression
)

# Joint dimensionality reduction
embeddings, _, _ = joint_pca(omics_data, n_components=50)

# Canonical correlation between single-cell and bulk
X_c, Y_c, _, _, correlations = canonical_correlation(
    omics_data,
    layer_pair=("singlecell", "transcriptomics"),
    n_components=10
)
```

### With Epigenome Module
```python
from metainformant.multiomics import integrate_omics_data
from metainformant.epigenome import load_cpg_table, compute_beta_values

# Integrate epigenomics with transcriptomics
methylation_df = load_cpg_table("methylation.tsv")
beta_values = compute_beta_values(methylation_df)

# Create multi-omics dataset
omics_data = MultiOmicsData(
    epigenomics=beta_values,
    transcriptomics=expression_df
)

# Joint analysis of methylation and expression
# Find correlations between methylation and gene expression
```

### With Visualization Module
```python
from metainformant.multiomics import joint_pca
from metainformant.visualization import scatter_plot, heatmap

# Visualize joint PCA results
embeddings, loadings, variance = joint_pca(omics_data, n_components=2)

# Scatter plot of first two components
ax = scatter_plot(embeddings[:, 0], embeddings[:, 1],
                  xlabel="PC1", ylabel="PC2",
                  title="Joint PCA: Multi-Omics Integration")

# Heatmap of loadings across omics layers
ax = heatmap(loadings, title="Joint PCA Loadings")
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
