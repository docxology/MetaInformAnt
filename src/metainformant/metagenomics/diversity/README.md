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
| `metrics.permanova()` | PERMANOVA test for community composition differences between groups |
| `metrics.ordination()` | Dimensionality reduction (PCoA, NMDS) for visualization |

## Usage

```python
from metainformant.metagenomics.diversity import metrics

alpha = metrics.alpha_diversity(abundances, metric="shannon")
beta = metrics.beta_diversity(samples, metric="bray_curtis")
curve = metrics.rarefaction_curve(abundances, steps=20)
perm = metrics.permanova(distance_matrix, groups, permutations=999)
```
