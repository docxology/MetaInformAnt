# Visualization Genomics

Genomics-specific visualization including Manhattan plots, volcano plots, expression heatmaps, network diagrams, and phylogenetic trees.

## Contents

| File | Purpose |
|------|---------|
| `expression.py` | Expression heatmaps, enrichment barplots, differential expression plots |
| `genomics.py` | Manhattan plots, volcano plots, regional plots, chromosome ideograms |
| `networks.py` | Network layouts: basic, circular, hierarchical (columns by distance from the nearest root), force-directed, community |
| `trees.py` | Phylogenetic tree plots: rectangular, circular, unrooted, annotated |

## Key Functions

| Function | Description |
|----------|-------------|
| `plot_expression_heatmap()` | Clustered gene expression heatmap |
| `plot_differential_expression()` | Volcano or MA plot for DE results |
| `manhattan_plot()` | Genome-wide association Manhattan plot |
| `volcano_plot()` | -log10(p) vs fold-change scatter |
| `regional_plot()` | Locus zoom regional association plot, with an optional recombination-rate track from a `recombination_rate` DataFrame (`BP`, `RATE`, optional `CHR`) |
| `chromosome_ideogram()` | One bar per chromosome with its relative length (GRCh38 by default, or any `chromosome_lengths`); no cytogenetic bands |
| `plot_network_basic()` | Spring-layout network visualization |
| `plot_phylo_tree()` | Rectangular phylogram (branch lengths to scale, bootstrap support beside the nodes) of a `dna.phylogeny` tree dictionary, a Newick string, a `Bio.Phylo` tree or a NetworkX `DiGraph` |
| `circular_tree_plot()` | Circular/radial phylogenetic tree layout (leaves on the unit circle, branch lengths not to scale) |
| `unrooted_tree_plot()` | Unrooted tree (a two-child root is removed; seeded force-directed layout, topology only) |
| `tree_comparison_plot()` | Two trees side by side, on a new figure or inside a given `ax` |
| `tree_annotation_plot()` | A tree with per-node `color` and `label` annotations |

## Usage

```python
from metainformant.visualization.genomics.genomics import manhattan_plot, volcano_plot
from metainformant.visualization.genomics.trees import plot_phylo_tree

manhattan_plot(gwas_results, output_path="output/manhattan.png")
volcano_plot(de_results, output_path="output/volcano.png")
plot_phylo_tree(tree, output_path="output/tree.png")
```
