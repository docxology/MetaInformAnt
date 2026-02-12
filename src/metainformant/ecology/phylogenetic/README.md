# Phylogenetic Ecology

Phylogenetic diversity metrics and community structure analysis using pure Python implementations with nested-dict tree representations.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Package exports for `diversity` module |
| `diversity.py` | Phylogenetic diversity, UniFrac, NRI/NTI, phylogenetic signal, tree construction |

## Key Functions

| Function | Description |
|----------|-------------|
| `faiths_pd()` | Compute Faith's Phylogenetic Diversity for a set of taxa |
| `compute_unifrac()` | Weighted and unweighted UniFrac distance between communities |
| `phylogenetic_beta_diversity()` | Pairwise phylogenetic beta diversity matrix |
| `nri_nti()` | Net Relatedness Index and Nearest Taxon Index for community structure |
| `phylogenetic_signal()` | Blomberg's K and Pagel's lambda for trait-phylogeny association |
| `build_simple_tree()` | UPGMA/NJ tree construction from distance matrices |

## Usage

```python
from metainformant.ecology.phylogenetic import diversity

pd_value = diversity.faiths_pd(tree, taxa_present=["sp_A", "sp_B"])
unifrac = diversity.compute_unifrac(tree, community_a, community_b)
nri, nti = diversity.nri_nti(tree, community_taxa)
signal = diversity.phylogenetic_signal(tree, trait_values)
tree = diversity.build_simple_tree(distance_matrix, labels, method="upgma")
```
