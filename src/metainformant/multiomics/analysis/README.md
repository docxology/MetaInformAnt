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
| `canonical_correlation()` | Ridge-regularized canonical correlation analysis between two omics layers: closed form via thin SVDs of the standardized layers (the same numerical core as `methods.factorization.canonical_correlation`), `regularization` defaults to 0.0 (classical CCA) and a layer of numerical rank `n_samples - 1` or more raises `ValueError` asking for `regularization > 0`, missing values are replaced by their own feature's mean, and the returned weights apply to the standardized features and reproduce the canonical variates |
| `from_dna_variants()` | Convert VCF data to integration-ready format |
| `from_rna_expression()` | Convert expression matrix to integration-ready format (`normalize`: `log1p` when every observed value is non-negative, zeros and gaps included, then a z-score per gene) |
| `from_protein_abundance()` | Convert protein quantification to integration-ready format |
| `from_epigenome_data()` | Returns a copy of methylation or ChIP data, unchanged for every `data_type` |
| `from_metabolomics()` | Convert metabolomics measurements (natural log when every observed value is positive, then a z-score per metabolite) |
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
