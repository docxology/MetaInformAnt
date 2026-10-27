# Visualization

Specialized visualization functions for metagenomic data, including Krona-style taxonomy charts, stacked bar composition plots, rarefaction curves, ordination plots, diversity boxplots, and abundance heatmaps.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports plots submodule |
| `plots.py` | All metagenomics visualization functions using matplotlib |

## Key Functions

| Function | Description |
|----------|-------------|
| `plots.plot_krona_chart()` | Krona-style hierarchical taxonomy sunburst chart |
| `plots.plot_stacked_bar()` | Stacked bar plot of taxonomic composition across samples |
| `plots.plot_rarefaction_curves()` | Rarefaction curves for species richness estimation |
| `plots.plot_ordination()` | PCoA or NMDS (`method`, case-insensitive; `seed` for NMDS) ordination scatter plot of a distance matrix |
| `plots.plot_alpha_diversity()` | Alpha diversity boxplots by group |
| `plots.plot_heatmap()` | Abundance heatmap with rows and columns ordered by greedy nearest-neighbour similarity |

## Usage

```python
from metainformant.metagenomics.visualization import plots

plots.plot_krona_chart(taxonomy_data, output_path="output/krona.png")
plots.plot_stacked_bar(abundance_table, output_path="output/barplot.png")
plots.plot_rarefaction_curves(otu_table, step=20, output_path="output/rarefaction.png", seed=42)
plots.plot_ordination(distance_matrix, groups=sample_groups, output_path="output/pcoa.png")
plots.plot_ordination(distance_matrix, method="NMDS", groups=sample_groups, seed=42, output_path="output/nmds.png")
```
