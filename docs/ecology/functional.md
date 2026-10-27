# Functional Ecology

Trait-based diversity metrics for analyzing community structure through species functional traits. Implements Functional Richness (FRic), Functional Evenness (FEve), Functional Divergence (FDiv), Functional Dispersion (FDis), Rao's Quadratic Entropy, Community Weighted Means, and related metrics.

All functions live in `metainformant.ecology.analysis.functional`. The trait matrix is species-by-traits (`trait_matrix[i]` is the list of trait values for species `i`), and `abundances[i]` is the abundance of that same species. Species with zero abundance are ignored.

## Key Concepts

**Functional traits** are measurable characteristics of organisms (body size, leaf area, metabolic rate) that influence performance and fitness. Functional diversity describes the variety and distribution of these traits in a community.

**FRic** (Functional Richness) measures the volume of trait space occupied by the community. Computed as range (1-D), convex hull area (2-D), or bounding-box volume (3-D+).

**FEve** (Functional Evenness) quantifies regularity of abundance distribution along the minimum spanning tree in trait space. Values range from 0 (uneven) to 1 (perfectly even).

**FDiv** (Functional Divergence) measures how abundances distribute relative to the trait-space centroid. High FDiv indicates dominant species have extreme trait values.

**FDis** (Functional Dispersion) is the abundance-weighted mean distance of species to the weighted centroid (Laliberte and Legendre, 2010).

**Rao's Q** is the expected trait dissimilarity between two randomly selected individuals: Q = sum(d_ij * p_i * p_j).

## Function Reference

### `functional_richness(trait_matrix, abundances=None) -> float`

Volume of occupied trait space: the range for one trait, the convex hull area for two, and the product of the per-trait ranges (a bounding box that skips traits with zero range) for three or more. Species with zero abundance are excluded when `abundances` is given, and fewer than two remaining species give 0.0.

### `functional_evenness(trait_matrix, abundances) -> float`

Regularity of abundance along the MST in trait space (Villeger et al. 2008). Returns 0-1, and 0.0 when fewer than two species are present.

### `functional_divergence(trait_matrix, abundances) -> float`

Abundance-weighted deviation from the trait centroid (Villeger et al. 2008). Returns 0-1.

### `functional_dispersion(trait_matrix, abundances) -> float`

Abundance-weighted mean distance to centroid. The value is 0 or more and has the scale of the traits.

### `raos_quadratic_entropy(trait_matrix, abundances) -> float`

Rao's Q = sum(d_ij * p_i * p_j) over all ordered species pairs, using Euclidean trait distances (traits are not standardized).

### `community_weighted_mean(trait_matrix, abundances) -> List[float]`

CWM_t = sum(p_i * trait_it) for each trait. Returns one value per trait.

### `functional_redundancy(trait_matrix, abundances) -> float`

Simpson diversity minus Rao's Q. High values indicate many species share similar traits. Q uses unscaled Euclidean trait distances, so the result is negative when species are far apart on the trait scale (it cannot be negative when every pairwise trait distance is at most 1).

### `functional_beta_diversity(community1_traits, community1_abundances, community2_traits, community2_abundances) -> Dict[str, float]`

Decompose functional dissimilarity into `total`, `turnover`, and `nestedness`. `total` is Rao's Q of the pooled community (each community weighted equally) minus the mean of the two communities' Q, floored at 0. `nestedness` is the difference between the two Q values, capped at `total`, and `turnover` is `total - nestedness`.

### `trait_distance_matrix(trait_matrix, method="euclidean") -> List[List[float]]`

Pairwise trait distance matrix (symmetric, species by species). Methods: `euclidean`, `gower` (each trait is divided by its range and the result averaged over traits with a nonzero range, so values lie in 0-1); any other method raises `ValueError`.

### `functional_diversity_suite(trait_matrix, abundances) -> Dict`

All metrics in one call: `fric`, `feve`, `fdiv`, `fdis`, `raos_q`, `cwm`, `redundancy`.

## Usage Examples

```python
from metainformant.ecology.analysis.functional import (
    functional_richness, functional_evenness, functional_divergence,
    raos_quadratic_entropy, community_weighted_mean,
    functional_diversity_suite, trait_distance_matrix,
    functional_beta_diversity,
)

# Species traits: body size and metabolic rate
traits = [[2.0, 0.5], [4.0, 1.2], [6.0, 0.8], [3.0, 1.5], [5.0, 1.0]]
abundances = [20, 15, 10, 8, 5]

# Individual metrics
fric = functional_richness(traits, abundances)         # 1.85 (convex hull area)
feve = functional_evenness(traits, abundances)         # 0.91
fdiv = functional_divergence(traits, abundances)       # 0.70
rao_q = raos_quadratic_entropy(traits, abundances)     # 1.70
cwm = community_weighted_mean(traits, abundances)      # [3.60, 0.91]

# All metrics at once
suite = functional_diversity_suite(traits, abundances)
print(suite)  # keys: fric, feve, fdiv, fdis, raos_q, cwm, redundancy

# Trait distance matrix
dm = trait_distance_matrix(traits, method="gower")

# Functional dissimilarity between two communities
other_traits = [[8.0, 2.0], [9.0, 2.5], [10.0, 3.0]]
beta = functional_beta_diversity(traits, abundances, other_traits, [5, 5, 5])
print(beta["total"], beta["turnover"], beta["nestedness"])
```

## Configuration

Settings are passed as function arguments; the functional metrics read no environment variables.

All algorithms are pure Python with no external dependencies beyond the standard library.

## Related Modules

- `metainformant.ecology.analysis.community` -- taxonomic diversity metrics
- `metainformant.ecology.analysis.ordination` -- ordination using trait distance matrices
- `metainformant.ecology.visualization` -- plotting community results (abundance distributions, rarefaction, ordination, distance heatmaps)
