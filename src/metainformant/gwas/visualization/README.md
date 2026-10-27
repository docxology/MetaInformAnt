# GWAS Visualization

Comprehensive GWAS plotting organized into subpackages: core plots (Manhattan, QQ), genomic views (regional, LD, genome-wide), population structure (PCA, geography), statistical summaries, and interactive dashboards.

## Contents

| File/Directory | Purpose |
|----------------|---------|
| `general.py` | Core plots: Manhattan, QQ, regional, PCA, kinship heatmap, effect sizes |
| `config.py` | Plot configuration, color palettes, style defaults |
| `utils.py` | Shared helpers: the key spellings of chromosome, position and p-value (`get_chromosome`, `get_position`, `standardize_result_columns`, `CHROMOSOME_KEYS`, `POSITION_KEYS`, `P_VALUE_KEYS`), chromosome labels and natural ordering, cumulative chromosome offsets, and `safe_log10_p` |
| `genomic/` | Genome-wide views, LD heatmaps, regional association, variant annotations |
| `population/` | Population PCA, admixture, geographic distribution plots |
| `statistical/` | Effect size forest plots, method comparison, statistical summary plots |
| `interactive/` | Plotly interactive Manhattan, fine-mapping, composite dashboards |

## Key Functions (general.py)

| Function | Description |
|----------|-------------|
| `manhattan_plot()` | Genome-wide -log10(p) Manhattan plot from one record, a list of records or a DataFrame, reading the chromosome, position and p-value under the spellings in `utils.py`; results without a p-value are not drawn |
| `qq_plot()` | Quantile-quantile plot of observed vs expected p-values |
| `regional_plot()` | Locus-zoom regional association plot |
| `pca_plot()` | Population structure PCA scatter with cluster coloring |
| `kinship_heatmap()` | Kinship/relatedness matrix heatmap |
| `effect_size_plot()` | Beta coefficient forest plot with confidence intervals |
| `generate_all_plots()` | Manhattan and Q-Q plots from association results (records, a DataFrame, or a `.json`, `.tsv`, `.txt` or `.csv` file) plus PCA and kinship plots from JSON files; returns a `{plot name: path}` mapping. The larger suite is `interactive/suite.py` |
| `missingness_plot()` | Sample/variant missingness rate visualization |
| `functional_enrichment_plot()` | Bar chart of the genome-wide significant variants per chromosome (`gff_path` is accepted but not read). The category-based plot, which uses annotations you supply, is `statistical/effects.py::functional_enrichment_plot` |

## Usage

```python
from metainformant.gwas.visualization.general import manhattan_plot, qq_plot

manhattan_plot(results, output_path="output/gwas/manhattan.png")
qq_plot(p_values, output_path="output/gwas/qq.png")
```

The plot functions in `general.py` create the missing parent directories of `output_path`, and `generate_all_plots()` creates its `output_dir`.
