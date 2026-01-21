### DNA: Phylogeny

Functions: `neighbor_joining_tree`, `upgma_tree`, `to_newick`, `to_ascii`, `bootstrap_support`, `basic_tree_stats`.

```mermaid
flowchart LR
  AalignedSequences[Aligned sequences] --> BdistanceMatrix[Distance Matrix]
  B --> C[NJ]
  B --> D[UPGMA]
  C --> E[Tree]
  D --> E
  E --> F[Newick/ASCII/Stats]
```

Example

```python
from pathlib import Path
from metainformant.dna import sequences, phylogeny

seqs = sequences.read_fasta(str(Path("tests/data/dna/toy.fasta")))
tree = phylogeny.neighbor_joining_tree(seqs)
newick = phylogeny.to_newick(tree)

# Optional: ASCII rendering and basic stats
ascii_tree = phylogeny.to_ascii(tree)
stats = phylogeny.basic_tree_stats(tree)
```
