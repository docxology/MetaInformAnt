# Differential Expression Analysis

The differential expression module identifies genes whose expression levels
differ significantly between groups of single cells. It supports Wilcoxon
rank-sum tests, Welch's t-tests, pseudobulk aggregation for multi-sample
designs, volcano plot data preparation, and gene set activity scoring.

## Concepts

### Statistical Testing

Two cell-level tests are available. The **Wilcoxon rank-sum** test is a
non-parametric test that compares expression ranks between two groups. It is
robust to non-normality and outliers. The **Welch's t-test** assumes approximate
normality (reasonable for log-normalized data) and handles unequal variances.
Both tests are applied gene-by-gene and corrected for multiple testing with the
Benjamini-Hochberg procedure.

### Pseudobulk Aggregation

For multi-donor or multi-sample experiments, individual cells within a sample
are not statistically independent. `pseudobulk_de` sums single-cell expression
per sample, then runs sample-level Welch's t-tests. This produces more
conservative and statistically valid results for hierarchical experimental
designs.

### Volcano Plot Data

`volcano_data` takes a DE result table and classifies each gene as `"up"`,
`"down"`, or `"ns"` (not significant) based on fold-change and adjusted p-value
thresholds. It returns the vectors needed to plot a standard volcano plot.

### Gene Set Scoring

`gene_set_scoring` computes per-cell activity scores for curated gene sets. The
`"mean"` method subtracts the mean expression of a random background gene set
from the mean expression of the gene set genes, correcting for overall
expression level. The `"sum"` method totals expression without background
correction.

## Function Reference

### `differential_expression(expression_matrix, groups, gene_names, method, min_cells, min_log2fc) -> list[dict]`

Run gene-by-gene statistical testing between two cell groups. Returns a sorted
list of result dictionaries containing `gene`, `log2fc`, `p_value`,
`adjusted_p`, `pct_group1`, `pct_group2`, `mean_group1`, and `mean_group2`.

| Parameter          | Type         | Default      | Description                           |
|--------------------|--------------|--------------|---------------------------------------|
| `expression_matrix`| Any          | --           | Cells x genes matrix                  |
| `groups`           | `list[int]`  | --           | Group label per cell (exactly 2 unique values) |
| `gene_names`       | `list[str]`  | --           | Gene names for columns                |
| `method`           | `str`        | `"wilcoxon"` | `"wilcoxon"` or `"t_test"`            |
| `min_cells`        | `int`        | `3`          | Min expressing cells per group        |
| `min_log2fc`       | `float`      | `0.0`        | Min absolute fold change to report    |

### `pseudobulk_de(expression_matrix, cell_labels, sample_labels, groups, gene_names, min_cells_per_sample) -> list[dict]`

Aggregate cells by sample, then perform sample-level DE. Returns the same result
format as `differential_expression`.

| Parameter              | Type         | Default | Description                        |
|------------------------|--------------|---------|------------------------------------|
| `expression_matrix`    | Any          | --      | Cells x genes matrix               |
| `cell_labels`          | `list[str]`  | --      | Cell type label per cell           |
| `sample_labels`        | `list[str]`  | --      | Sample ID per cell                 |
| `groups`               | `list[int]`  | --      | Group assignment per cell          |
| `gene_names`           | `list[str]`  | `None`  | Gene names (auto-generated if None)|
| `min_cells_per_sample` | `int`        | `5`     | Min cells to include a sample      |

### `compute_log_fold_change(mean_a, mean_b, pseudocount) -> float`

Compute `log2((mean_a + pseudocount) / (mean_b + pseudocount))`.

### `volcano_data(de_results, fc_threshold, p_threshold) -> dict`

Classify DE results into up/down/ns and return plotting vectors (`genes`,
`log2fc`, `neg_log10_p`, `classification`, `n_up`, `n_down`, `n_ns`).

### `gene_set_scoring(expression_matrix, gene_sets, gene_names, method, n_background, seed) -> dict`

Score gene set activity per cell. Returns `scores` (dict of per-cell score
lists), `n_cells`, `n_gene_sets`, and `gene_set_sizes`.

| Parameter          | Type                       | Default  | Description                     |
|--------------------|----------------------------|----------|---------------------------------|
| `expression_matrix`| Any                        | --       | Cells x genes matrix            |
| `gene_sets`        | `dict[str, list[str]]`     | --       | Named gene sets                 |
| `gene_names`       | `list[str]`                | --       | Gene names for columns          |
| `method`           | `str`                      | `"mean"` | `"mean"` or `"sum"`             |
| `n_background`     | `int`                      | `50`     | Background genes for mean method|
| `seed`             | `int \| None`              | `None`   | Random seed for reproducibility |

## Code Examples

```python
from metainformant.singlecell.differential.expression import (
    differential_expression,
    pseudobulk_de,
    volcano_data,
    gene_set_scoring,
)

# Basic two-group DE with Wilcoxon test
matrix = [[1.2, 0.0, 3.4], [0.5, 2.1, 0.0], [4.0, 0.1, 1.5], [0.2, 3.3, 0.8]]
groups = [0, 0, 1, 1]
genes = ["TP53", "BRCA1", "MYC"]

results = differential_expression(matrix, groups, genes, method="wilcoxon")
for hit in results:
    if hit["adjusted_p"] < 0.05:
        print(f"{hit['gene']}: log2FC={hit['log2fc']:.2f}, padj={hit['adjusted_p']:.4f}")

# Pseudobulk DE for multi-sample experiment
pb_results = pseudobulk_de(
    expression_matrix=matrix,
    cell_labels=["T_cell"] * 4,
    sample_labels=["S1", "S1", "S2", "S2"],
    groups=[0, 0, 1, 1],
    gene_names=genes,
)

# Prepare volcano plot data
volcano = volcano_data(results, fc_threshold=1.0, p_threshold=0.05)
print(f"Up: {volcano['n_up']}, Down: {volcano['n_down']}, NS: {volcano['n_ns']}")

# Score pathway activity per cell
scores = gene_set_scoring(
    matrix,
    gene_sets={"apoptosis": ["TP53", "MYC"], "repair": ["BRCA1"]},
    gene_names=genes,
)
print(scores["scores"]["apoptosis"])  # per-cell scores
```

## Configuration

All parameters are function arguments. No external config files. Optional
dependencies (`numpy`, `scipy`) are auto-detected; pure-Python fallbacks apply.

## Import Path

```python
from metainformant.singlecell.differential.expression import differential_expression
from metainformant.singlecell.differential.expression import pseudobulk_de
from metainformant.singlecell.differential.expression import volcano_data
from metainformant.singlecell.differential.expression import gene_set_scoring
from metainformant.singlecell.differential.expression import compute_log_fold_change
```
