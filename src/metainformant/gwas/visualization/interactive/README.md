# GWAS Interactive and Composite Visualization

Plotly interactive figures, matplotlib composite panels, fine-mapping plots and phenotype plots for GWAS results. Plot functions return a dict with a `status` key and create the missing parent directories of their output path.

## Contents

| File | Purpose |
|------|---------|
| `interactive.py` | `interactive_manhattan()`, `interactive_pca()`, `interactive_volcano()`: Plotly HTML figures, with a styled HTML table when Plotly is unavailable |
| `finemapping.py` | `compute_credible_set()`, `credible_set_plot()`, `conditional_analysis_plot()`, `pip_vs_ld_plot()` |
| `composite.py` | `gwas_summary_panel()`, `population_structure_panel()`, `top_hit_detail_panel()` |
| `phenotype.py` | `phenotype_distribution()`, `phenotype_correlation_matrix()`, `genotype_phenotype_boxplot()`, `top_hits_genotype_phenotype()`, `phenotype_pca_correlation()` |
| `suite.py` | `generate_all_plots()`: the full plot set from association results |

## Interactive Manhattan Plot

`interactive_manhattan(assoc_results, output_file=None, gene_annotations=None, significance_threshold=5e-8, title="Interactive Manhattan Plot")` takes a list of result dicts. Each result needs a chromosome, a position and a positive p-value. For each value the first listed key that is present (not `None`) is used:

| Value | Accepted keys |
|-------|---------------|
| Chromosome | `chromosome`, `chrom`, `chr`, `CHROM`, `CHR` |
| Position | `position`, `pos`, `POS` |
| p-value | `p_value`, `pvalue`, `P` |

`variant_id` (or `rsid`) and `beta` (or `effect_size`) are optional and feed the hover text. A result without a chromosome, a position or a finite positive p-value is skipped and counted in `n_skipped`; it is never placed on a default chromosome. A float chromosome such as `7.0` is labelled `"7"`.

The returned dict holds `status` (`"success"`, or `"skipped"` when no result is usable), `output_path`, `n_variants` (results plotted), `n_skipped` and `interactive` (`False` when Plotly is unavailable and an HTML table of the most significant results is written instead). `interactive_volcano()` reads the optional chromosome from the same keys and returns `status`, `output_path`, `n_variants` and `interactive`. `metainformant.gwas.interactive_manhattan` is this function; the DataFrame-based function of the same name in `metainformant.visualization.dashboards.interactive` takes `chrom_col`, `pos_col` and `pvalue_col` instead.

## Credible Sets

`compute_credible_set(assoc_results, credible_level=0.95)` delegates to `metainformant.gwas.finemapping.credible_sets.compute_credible_set`, which uses Wakefield approximate Bayes factors. When every result has `beta` and a positive `se`, the Z-scores are `beta / se` and the Bayes factors use `V = se**2`. Otherwise the Z-scores are recovered from the two-sided `p_value` (or `P`) and the z-score variance model is used. It returns `status`, `pips`, `credible_set_indices` (highest PIP first), `credible_set_size` and `cumulative_probability`, or `{"status": "failed", "reason": ...}` when the input is empty or the computation fails. `credible_set_plot()` uses it.

`metainformant.gwas.compute_credible_set` is this wrapper. The Z-score function with `standard_errors=` and `sample_size=` arguments is `metainformant.gwas.finemapping.credible_sets.compute_credible_set`.

## Usage

```python
from metainformant.gwas.visualization.interactive import finemapping, interactive

results = [
    {"CHR": "1", "POS": 1000, "P": 3e-9, "beta": 0.4, "se": 0.07},
    {"CHR": "1", "POS": 2000, "P": 0.2, "beta": 0.05, "se": 0.04},
    {"CHR": "2", "POS": 500, "P": 1e-4, "beta": 0.1, "se": 0.03},
]

summary = interactive.interactive_manhattan(results, "output/gwas/manhattan.html")
print(summary["n_variants"], summary["n_skipped"], summary["interactive"])

credible = finemapping.compute_credible_set(results, credible_level=0.95)
print(credible["credible_set_indices"])
```
