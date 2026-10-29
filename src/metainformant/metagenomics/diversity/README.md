# Diversity

Community diversity metrics for metagenomic data, providing alpha diversity (Shannon, Simpson, Chao1), beta diversity (Bray-Curtis, Jaccard, Aitchison), rarefaction analysis, PERMANOVA, and ordination (PCoA, NMDS).

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports metrics submodule |
| `metrics.py` | Alpha/beta diversity, rarefaction, PERMANOVA, ordination methods |

## Key Functions

| Function | Description |
|----------|-------------|
| `metrics.alpha_diversity()` | Compute within-sample diversity (Shannon, Simpson, Chao1, ACE, Fisher) |
| `metrics.beta_diversity()` | Compute between-sample distance matrices (Bray-Curtis, Jaccard, Aitchison) |
| `metrics.rarefaction_curve()` | Generate rarefaction curves for species richness estimation |
| `metrics.rarefy()` | Subsample abundance data to uniform depth |
| `metrics.permanova()` | PERMANOVA test for community composition differences between groups (pseudo-F is `inf` for perfectly separated groups) |
| `metrics.ordination()` | PCoA (exact eigendecomposition; all eigenvalues, negative ones included, are reported) or Kruskal's non-metric MDS (isotonic regression, SMACOF, several seeded starts; reports stress-1) |

## Usage

```python
from metainformant.metagenomics.diversity import metrics

alpha = metrics.alpha_diversity(abundances, metric="shannon")
beta = metrics.beta_diversity(samples, metric="bray_curtis")
curve = metrics.rarefaction_curve(abundances, n_iterations=20, seed=42)
perm = metrics.permanova(distance_matrix, groups, n_permutations=999, seed=42)
```

`chao1`, `ace` and `fisher_alpha` need integer counts and raise `ValueError` for proportions; the other metrics accept counts or proportions. `rarefaction_curve`, `rarefy` and `permanova` are unseeded by default (`seed=None`); `ordination` seeds its NMDS random starts with `seed=42`. PCoA `variance_explained` is relative to the sum of the positive eigenvalues; a non-Euclidean distance matrix has negative eigenvalues, listed in `all_eigenvalues` and summarised by `negative_eigenvalue_fraction`. NMDS has no eigenvalues or variance shares and reports `stress` instead. `ordination` and `permanova` reject matrices that are not square, symmetric, finite, non-negative and zero on the diagonal.
