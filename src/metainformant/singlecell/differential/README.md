# Single-Cell Differential Expression

Statistical testing for differential gene expression between cell groups. Supports Wilcoxon rank-sum, Welch's t-test, and pseudobulk aggregation for multi-sample designs, plus fold change computation, volcano plot preparation, and gene set scoring.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Re-exports `expression` module |
| `expression.py` | DE testing, pseudobulk, fold change, volcano data, gene set scoring |

## Key Functions

| Function | Description |
|----------|-------------|
| `differential_expression()` | Test for DE genes between two cell groups (Wilcoxon or t-test) |
| `pseudobulk_de()` | Aggregate cells by sample and perform sample-level DE testing |
| `compute_log_fold_change()` | Compute log2 fold change between groups |
| `volcano_data()` | Prepare data for volcano plot visualization |
| `gene_set_scoring()` | Score gene set activity per cell |

## Usage

```python
from metainformant.singlecell.differential import expression

results = expression.differential_expression(matrix, groups, method="wilcoxon")
pseudobulk = expression.pseudobulk_de(matrix, groups, sample_ids)
volcano = expression.volcano_data(results)
```
