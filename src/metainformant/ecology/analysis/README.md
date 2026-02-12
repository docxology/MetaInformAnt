# Ecology Analysis

Community ecology analysis including diversity indices, functional traits, species-abundance distributions, ordination, and indicator species analysis.

## Contents

| File | Purpose |
|------|---------|
| `community.py` | Alpha/beta/gamma diversity, rarefaction, nestedness, similarity matrices |
| `functional.py` | Functional diversity: richness, evenness, divergence, Rao's Q, CWM |
| `indicators.py` | IndVal, ANOSIM, PERMANOVA, SIMPER, community clustering |
| `macroecology.py` | Species-abundance distributions, species-area relationships, distance decay |
| `ordination.py` | PCoA, NMDS, CCA, distance matrices, Procrustes analysis |

## Key Functions

| Function | Description |
|----------|-------------|
| `calculate_diversity()` | Shannon, Simpson, and other diversity indices |
| `beta_diversity()` | Bray-Curtis, Jaccard, and other pairwise dissimilarities |
| `rarefaction_curve()` | Expected species richness at subsampled depths |
| `functional_richness()` | Convex hull volume in trait space |
| `functional_dispersion()` | Abundance-weighted mean distance to centroid |
| `indval()` | Indicator value analysis for species-habitat associations |
| `permanova()` | Permutational multivariate ANOVA on distance matrices |
| `fit_logseries()` | Fit log-series species-abundance distribution |
| `pcoa()` | Principal coordinates analysis on distance matrices |
| `nmds()` | Non-metric multidimensional scaling |

## Usage

```python
from metainformant.ecology.analysis.community import calculate_diversity, beta_diversity
from metainformant.ecology.analysis.ordination import pcoa, distance_matrix

indices = calculate_diversity([10, 20, 30, 5], method="shannon")
dist = distance_matrix(communities, method="bray_curtis")
coords = pcoa(dist, n_components=2)
```
