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

A gene whose test cannot be computed (a group with no or one value, constant values, non-finite values) has `p_value = NaN`, `adjusted_p = NaN` and a `test_status` that says why; it takes no part in the Benjamini-Hochberg correction and sorts last. `log2fc` is group 1 over group 2, where group 1 carries the smaller group label.

## Usage

```python
from metainformant.singlecell.differential import expression

results = expression.differential_expression(matrix, groups, gene_names, method="wilcoxon")
pseudobulk = expression.pseudobulk_de(matrix, cell_labels, sample_ids, groups)
volcano = expression.volcano_data(results)
```
