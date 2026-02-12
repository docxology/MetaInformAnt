# Phylogeny

Phylogenetic tree construction (Neighbor-Joining, UPGMA, k-mer based) and tree analysis including Newick I/O, bootstrap support, and topology comparison.

## Contents

| File | Purpose |
|------|---------|
| `tree.py` | Re-exports all public symbols from construction and analysis submodules |
| `tree_construction.py` | NJ, UPGMA, and k-mer tree building from distance matrices |
| `tree_analysis.py` | Newick parsing, bootstrap, topology metrics, and tree manipulation |

## Key Functions

| Function | Description |
|----------|-------------|
| `neighbor_joining_tree()` | Build a tree using the Neighbor-Joining algorithm |
| `upgma_tree()` | Build an ultrametric tree using UPGMA clustering |
| `nj_tree_from_kmer()` | Alignment-free tree from k-mer distance profiles |
| `to_newick()` | Serialize a tree to Newick format string |
| `from_newick()` | Parse a Newick string into a tree data structure |
| `bootstrap_support()` | Compute bootstrap support values for tree branches |
| `to_ascii()` | Render a tree as an ASCII art diagram |
| `basic_tree_stats()` | Number of leaves, internal nodes, and total branch length |
| `robinson_foulds_distance()` | Topological distance between two trees |
| `is_monophyletic()` | Test whether a set of taxa forms a monophyletic group |
| `prune_tree()` | Remove specified taxa from a tree |
| `tree_diameter()` | Longest path between any two leaves |

## Usage

```python
from metainformant.dna.phylogeny.tree import neighbor_joining_tree, to_newick, bootstrap_support

tree = neighbor_joining_tree(distance_matrix)
newick_str = to_newick(tree)
support = bootstrap_support(sequences, n_bootstraps=100)
```
