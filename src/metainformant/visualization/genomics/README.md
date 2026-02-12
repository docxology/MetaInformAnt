# Visualization Genomics

Genomics-specific visualization including Manhattan plots, volcano plots, expression heatmaps, network diagrams, and phylogenetic trees.

## Contents

| File | Purpose |
|------|---------|
| `expression.py` | Expression heatmaps, enrichment barplots, differential expression plots |
| `genomics.py` | Manhattan plots, volcano plots, regional plots, chromosome ideograms |
| `networks.py` | Network layouts: basic, circular, hierarchical, force-directed, community |
| `trees.py` | Phylogenetic tree plots: rectangular, circular, unrooted, annotated |

## Key Functions

| Function | Description |
|----------|-------------|
| `plot_expression_heatmap()` | Clustered gene expression heatmap |
| `plot_differential_expression()` | Volcano or MA plot for DE results |
| `manhattan_plot()` | Genome-wide association Manhattan plot |
| `volcano_plot()` | -log10(p) vs fold-change scatter |
| `regional_plot()` | Locus zoom regional association plot |
| `chromosome_ideogram()` | Chromosome banding pattern diagram |
| `plot_network_basic()` | Spring-layout network visualization |
| `plot_phylo_tree()` | Rectangular phylogenetic tree |
| `circular_tree_plot()` | Circular/radial phylogenetic tree layout |

## Usage

```python
from metainformant.visualization.genomics.genomics import manhattan_plot, volcano_plot
from metainformant.visualization.genomics.trees import plot_phylo_tree

manhattan_plot(gwas_results, output_path="output/manhattan.png")
volcano_plot(de_results, output_path="output/volcano.png")
plot_phylo_tree(tree, output_path="output/tree.png")
```
