# Visualization

Specialized visualization functions for metagenomic data, including Krona-style taxonomy charts, stacked bar composition plots, rarefaction curves, ordination plots, diversity boxplots, and abundance heatmaps.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports plots submodule |
| `plots.py` | All metagenomics visualization functions using matplotlib/seaborn |

## Key Functions

| Function | Description |
|----------|-------------|
| `plots.plot_krona_chart()` | Krona-style hierarchical taxonomy sunburst chart |
| `plots.plot_stacked_bar()` | Stacked bar plot of taxonomic composition across samples |
| `plots.plot_rarefaction_curves()` | Rarefaction curves for species richness estimation |
| `plots.plot_ordination()` | PCoA/NMDS ordination scatter plot |
| `plots.plot_alpha_diversity()` | Alpha diversity boxplots by group |
| `plots.plot_heatmap()` | Abundance heatmap with hierarchical clustering |

## Usage

```python
from metainformant.metagenomics.visualization import plots

plots.plot_krona_chart(taxonomy_data, output_path="output/krona.png")
plots.plot_stacked_bar(abundance_table, output_path="output/barplot.png")
plots.plot_rarefaction_curves(rarefaction_data, output_path="output/rarefaction.png")
plots.plot_ordination(coords, groups, output_path="output/pcoa.png")
```
