# GWAS Visualization

Comprehensive GWAS plotting organized into subpackages: core plots (Manhattan, QQ), genomic views (regional, LD, genome-wide), population structure (PCA, geography), statistical summaries, and interactive dashboards.

## Contents

| File/Directory | Purpose |
|----------------|---------|
| `general.py` | Core plots: Manhattan, QQ, regional, PCA, kinship heatmap, effect sizes |
| `config.py` | Plot configuration, color palettes, style defaults |
| `genomic/` | Genome-wide views, LD heatmaps, regional association, variant annotations |
| `population/` | Population PCA, admixture, geographic distribution plots |
| `statistical/` | Effect size forest plots, method comparison, statistical summary plots |
| `interactive/` | Plotly interactive Manhattan, fine-mapping, composite dashboards |

## Key Functions (general.py)

| Function | Description |
|----------|-------------|
| `manhattan_plot()` | Chromosome-wide -log10(p) Manhattan plot |
| `qq_plot()` | Quantile-quantile plot of observed vs expected p-values |
| `regional_plot()` | Locus-zoom regional association plot |
| `pca_plot()` | Population structure PCA scatter with cluster coloring |
| `kinship_heatmap()` | Kinship/relatedness matrix heatmap |
| `effect_size_plot()` | Beta coefficient forest plot with confidence intervals |
| `generate_all_plots()` | Generate full GWAS plot suite from results |
| `missingness_plot()` | Sample/variant missingness rate visualization |
| `functional_enrichment_plot()` | Enrichment of hits in functional categories |

## Usage

```python
from metainformant.gwas.visualization.general import manhattan_plot, qq_plot

manhattan_plot(results, output_path="output/gwas/manhattan.png")
qq_plot(p_values, output_path="output/gwas/qq.png")
```
