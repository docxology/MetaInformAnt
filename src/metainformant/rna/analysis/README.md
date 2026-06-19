# Analysis

RNA-seq expression analysis including normalization, differential expression, quality control, batch effect detection, cross-species comparison, and protein integration.

## Contents

| File | Purpose |
|------|---------|
| `expression_core.py` | Count normalization (CPM, TPM, RPKM, quantile, median-ratio) and filtering |
| `expression_analysis.py` | Differential expression (DESeq2-like, t-test, Wilcoxon), PCA, and volcano data |
| `qc_metrics.py` | Sample/gene QC metrics, outlier detection, library complexity, saturation curves |
| `qc_filtering.py` | Batch effect detection, GC bias, length bias, and QC report generation |
| `cross_species.py` | Ortholog mapping, expression conservation, divergence matrices, cross-species PCA |
| `protein_integration.py` | Translation efficiency, protein abundance prediction, ribosome profiling |
| `validation.py` | Pipeline validation: per-sample status checks and end-to-end reports |

## Key Functions

| Function | Description |
|----------|-------------|
| `normalize_counts()` | Normalize raw counts by CPM, TPM, RPKM, quantile, or median-ratio |
| `filter_low_expression()` | Remove genes below a minimum count threshold |
| `differential_expression()` | DE analysis with method selection and p-value adjustment |
| `pca_analysis()` | Principal component analysis on expression matrices |
| `compute_sample_metrics()` | Per-sample statistics: total counts, detected genes, complexity |
| `detect_outlier_samples()` | Flag outlier samples using median deviation |
| `detect_batch_effects()` | Identify confounding batch structure in expression data |
| `build_ortholog_map()` | Load one-to-one ortholog mappings between species |
| `compute_expression_conservation()` | Spearman correlation of orthologs across species |
| `calculate_translation_efficiency()` | Estimate translation efficiency from RNA and protein data using `ratio` or `correlation` |
| `predict_protein_abundance_from_rna()` | Predict protein abundance from RNA with the implemented `linear` method |
| `validate_all_samples()` | Check pipeline completion status for every sample |

## Usage

```python
from metainformant.rna.analysis.expression_core import normalize_counts, filter_low_expression
from metainformant.rna.analysis.expression_analysis import differential_expression, pca_analysis
from metainformant.rna.analysis.qc_metrics import compute_sample_metrics
from metainformant.rna.analysis.protein_integration import calculate_translation_efficiency

normalized = normalize_counts(raw_counts, method="tpm", gene_lengths=lengths)
de_results = differential_expression(counts, groups, method="deseq2_like")
pca_result = pca_analysis(normalized, n_components=3)
rna_expression = normalized.T  # samples x genes for protein integration
protein_abundance = rna_expression.copy()
translation_efficiency = calculate_translation_efficiency(rna_expression, protein_abundance, method="ratio")
```

## Validation Notes

- Count/QC matrices must be numeric, finite, and non-negative where raw counts are expected.
- Batch labels must include every expression sample; labels for extra samples are ignored.
- GC content values for matched genes must be in `[0, 1]`; matched gene lengths must be positive.
- RNA-protein integration filters NaN measurements deterministically and rejects unsupported methods instead of returning silent empty results.
