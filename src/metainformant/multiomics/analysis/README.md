# Multi-Omics Analysis

Integration methods for combining genomics, transcriptomics, proteomics, and epigenomics data layers using joint factorization and correlation approaches.

## Contents

| File | Purpose |
|------|---------|
| `integration.py` | Multi-omics data container, joint PCA, NMF, CCA, and data converters |

## Key Classes and Functions

| Symbol | Description |
|--------|-------------|
| `MultiOmicsData` | Container for multiple omics layers with sample alignment |
| `integrate_omics_data()` | Combine multiple omics DataFrames into unified representation |
| `joint_pca()` | Joint PCA across concatenated omics layers |
| `joint_nmf()` | Non-negative matrix factorization across omics layers |
| `canonical_correlation()` | Canonical correlation analysis between two omics types |
| `from_dna_variants()` | Convert VCF data to integration-ready format |
| `from_rna_expression()` | Convert expression matrix to integration-ready format |
| `from_protein_abundance()` | Convert protein quantification to integration-ready format |
| `from_epigenome_data()` | Convert methylation or ChIP data to integration-ready format |
| `from_metabolomics()` | Convert metabolomics measurements (log/z-score normalization) |
| `compute_multiomics_similarity()` | Sample-by-sample similarity from concatenated layers |
| `find_multiomics_modules()` | Joint-NMF-based co-regulated feature modules |
| `descriptive_summary` | Per-block statistics and cross-omics Spearman correlations |

## Usage

```python
from metainformant.multiomics.analysis.integration import (
    MultiOmicsData,
    integrate_omics_data,
    joint_pca,
)

multi = MultiOmicsData(data={"rna": rna_df, "protein": prot_df})
integrated = integrate_omics_data({"rna": rna_df, "protein": prot_df})
embeddings, loadings, variance = joint_pca(integrated, n_components=10)
```
