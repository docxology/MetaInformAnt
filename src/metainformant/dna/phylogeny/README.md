# Phylogeny

Phylogenetic tree construction (Neighbor-Joining, UPGMA, k-mer based) and tree analysis including Newick I/O, bootstrap support, and topology comparison.

## Contents

| File | Purpose |
|------|---------|
| `tree.py` | Re-exports all public symbols from construction and analysis submodules |
| `tree_construction.py` | NJ, UPGMA, and k-mer tree building from `{name: sequence}` dicts (the distance matrices are computed internally) |
| `tree_analysis.py` | Newick parsing, bootstrap, topology metrics, and tree manipulation |

## Key Functions

| Function | Description |
|----------|-------------|
| `neighbor_joining_tree()` | Build a tree from a `{name: sequence}` dict using the Neighbor-Joining algorithm |
| `upgma_tree()` | Build an ultrametric tree from a `{name: sequence}` dict using UPGMA clustering |
| `nj_tree_from_kmer()` | Alignment-free tree from k-mer distance profiles |
| `to_newick()` | Serialize a tree to Newick format string |
| `from_newick()` | Parse a Newick string into a tree data structure |
| `bootstrap_support()` | Resample alignment sites, rebuild a tree per replicate (`method` `"nj"` or `"upgma"`), and return a new tree with a `bootstrap` support value (0-100) on each internal node; the input tree is not modified |
| `to_ascii()` | Render a tree as an ASCII art diagram |
| `basic_tree_stats()` | Number of leaves, internal nodes, and total branch length |
| `robinson_foulds_distance()` | Topological distance between two trees |
| `is_monophyletic()` | Test whether a set of taxa forms a monophyletic group |
| `prune_tree()` | Remove specified taxa from a tree |
| `tree_diameter()` | Longest path between any two leaves |

## Usage

```python
from metainformant.dna.phylogeny.tree import neighbor_joining_tree, to_newick, bootstrap_support

sequences = {  # aligned sequences of equal length
    "seq1": "ACGTACGTACGT",
    "seq2": "ACGTACGTACGA",
    "seq3": "ACGAACGAACGT",
    "seq4": "TCGAACGAACGT",
}

tree = neighbor_joining_tree(sequences)
newick_str = to_newick(tree)

# A new tree is returned and `tree` is left unchanged
supported = bootstrap_support(tree, sequences, n_replicates=100)
support = {node: children["bootstrap"] for node, children in supported.items() if isinstance(children, dict)}  # internal node -> 0-100
```
