### Visualization: Trees

Function: `plot_phylo_tree`

```python
from pathlib import Path
from metainformant.dna import sequences, phylogeny
from metainformant.visualization import trees

seqs = sequences.read_fasta(str(Path("tests/data/dna/toy.fasta")))
tree = phylogeny.neighbor_joining_tree(seqs)
ax = trees.plot_phylo_tree(tree)

# Cross-link: build trees in [DNA: Phylogeny](../dna/phylogeny.md) and render here.
```
