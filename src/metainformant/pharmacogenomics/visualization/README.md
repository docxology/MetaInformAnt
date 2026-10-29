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
| `plots.plot_allele_frequencies()` | Allele frequency chart for a pharmacogene: `plot_type="bar"` or `"pie"` (anything else raises `ValueError`) |
| `plots.plot_activity_score_distribution()` | Histogram of activity scores across samples, with CPIC phenotype boundaries for genes that have them (none for CYP1A2, which CPIC defines no phenotype for) |
| `plots.plot_drug_response_heatmap()` | Drug-gene matrix of interaction severity or phenotype for `drugs` x `genes` |
| `plots.plot_population_comparison()` | Cross-population allele frequency comparison |
| `plots.plot_acmg_criteria()` | Summary plot of ACMG evidence criteria |

## Usage

```python
from metainformant.pharmacogenomics.visualization import plots

plots.plot_metabolizer_status(phenotype_counts, output_path="output/metabolizers.png")
plots.plot_allele_frequencies(allele_data, gene="CYP2D6", output_path="output/alleles.png")
plots.plot_drug_response_heatmap(
    drugs=["codeine", "warfarin"],
    genes=["CYP2D6", "CYP2C9"],
    phenotypes={("codeine", "CYP2D6"): "Major", ("warfarin", "CYP2C9"): "Moderate"},
    output_path="output/response.png",
)
```
