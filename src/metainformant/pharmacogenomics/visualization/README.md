# Visualization

Publication-quality plots for pharmacogenomics data, including metabolizer status distributions, allele frequencies, activity score distributions, drug response heatmaps, and ACMG criteria summaries.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports plots submodule |
| `plots.py` | All pharmacogenomics visualization functions using matplotlib/seaborn |

## Key Functions

| Function | Description |
|----------|-------------|
| `plots.plot_metabolizer_status()` | Bar chart of metabolizer phenotype distributions |
| `plots.plot_allele_frequencies()` | Allele frequency bar chart for a pharmacogene |
| `plots.plot_activity_score_distribution()` | Histogram of activity scores across samples |
| `plots.plot_drug_response_heatmap()` | Heatmap of drug responses by metabolizer phenotype |
| `plots.plot_population_comparison()` | Cross-population allele frequency comparison |
| `plots.plot_acmg_criteria()` | Summary plot of ACMG evidence criteria |

## Usage

```python
from metainformant.pharmacogenomics.visualization import plots

plots.plot_metabolizer_status(phenotype_counts, output_path="output/metabolizers.png")
plots.plot_allele_frequencies(allele_data, gene="CYP2D6", output_path="output/alleles.png")
plots.plot_drug_response_heatmap(response_data, output_path="output/response.png")
```
