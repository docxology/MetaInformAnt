# Phylogenetic Trees

Phylogenetic tree visualization functions including basic tree plots, circular layouts, unrooted trees, tree comparisons, and tree annotations.

## Functions

### `plot_phylo_tree(tree, *, ax=None)`

Plot a Biopython Phylo tree to matplotlib Axes.

**Example:**
```python
from pathlib import Path
from metainformant.dna import sequences, phylogeny
from metainformant.visualization import plot_phylo_tree

seqs = sequences.read_fasta(str(Path("tests/data/dna/toy.fasta")))
tree = phylogeny.neighbor_joining_tree(seqs)
ax = plot_phylo_tree(tree)
```

### `circular_tree_plot(tree, *, ax=None)`

Plot a phylogenetic tree in circular layout.

### `unrooted_tree_plot(tree, *, ax=None)`

Plot an unrooted phylogenetic tree.

### `tree_comparison_plot(trees, labels=None, *, ncols=2, figsize=None)`

Plot multiple trees side by side for comparison.

### `tree_annotation_plot(tree, annotations=None, *, ax=None, **kwargs)`

Plot a phylogenetic tree with annotations.

## Cross-links

- Build trees in [DNA: Phylogeny](../dna/phylogeny.md) and render here.
