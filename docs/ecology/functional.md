# Functional Ecology

Trait-based diversity metrics for analyzing community structure through species functional traits. Implements Functional Richness (FRic), Functional Evenness (FEve), Functional Divergence (FDiv), Functional Dispersion (FDis), Rao's Quadratic Entropy, Community Weighted Means, and related metrics.

## Key Concepts

**Functional traits** are measurable characteristics of organisms (body size, leaf area, metabolic rate) that influence performance and fitness. Functional diversity describes the variety and distribution of these traits in a community.

**FRic** (Functional Richness) measures the volume of trait space occupied by the community. Computed as range (1-D), convex hull area (2-D), or bounding-box volume (3-D+).

**FEve** (Functional Evenness) quantifies regularity of abundance distribution along the minimum spanning tree in trait space. Values range from 0 (uneven) to 1 (perfectly even).

**FDiv** (Functional Divergence) measures how abundances distribute relative to the trait-space centroid. High FDiv indicates dominant species have extreme trait values.

**FDis** (Functional Dispersion) is the abundance-weighted mean distance of species to the weighted centroid (Laliberte and Legendre, 2010).

**Rao's Q** is the expected trait dissimilarity between two randomly selected individuals: Q = sum(d_ij * p_i * p_j).

## Function Reference

### `functional_richness(trait_matrix, abundances=None) -> float`

Volume of occupied trait space. Filters by abundance if provided.

### `functional_evenness(trait_matrix, abundances) -> float`

Regularity of abundance along the MST in trait space. Returns 0-1.

### `functional_divergence(trait_matrix, abundances) -> float`

Abundance-weighted deviation from the trait centroid. Returns 0-1.

### `functional_dispersion(trait_matrix, abundances) -> float`

Abundance-weighted mean distance to centroid.

### `raos_quadratic_entropy(trait_matrix, abundances) -> float`

Rao's Q = sum(d_ij * p_i * p_j) using Euclidean trait distances.

### `community_weighted_mean(trait_matrix, abundances) -> List[float]`

CWM_t = sum(p_i * trait_it) for each trait. Returns one value per trait.

### `functional_redundancy(trait_matrix, abundances) -> float`

Simpson diversity minus Rao's Q. High values indicate many species share similar traits.

### `functional_beta_diversity(c1_traits, c1_ab, c2_traits, c2_ab) -> Dict`

Decompose functional dissimilarity into `total`, `turnover`, and `nestedness`.

### `trait_distance_matrix(trait_matrix, method="euclidean") -> List[List[float]]`

Pairwise trait distance matrix. Methods: `euclidean`, `gower`.

### `functional_diversity_suite(trait_matrix, abundances) -> Dict`

All metrics in one call: `fric`, `feve`, `fdiv`, `fdis`, `raos_q`, `cwm`, `redundancy`.

## Usage Examples

```python
from metainformant.ecology import (
    functional_richness, functional_evenness, functional_divergence,
    raos_quadratic_entropy, community_weighted_mean,
    functional_diversity_suite, trait_distance_matrix,
)

# Species traits: body size and metabolic rate
traits = [[2.0, 0.5], [4.0, 1.2], [6.0, 0.8], [3.0, 1.5], [5.0, 1.0]]
abundances = [20, 15, 10, 8, 5]

# Individual metrics
fric = functional_richness(traits, abundances)
feve = functional_evenness(traits, abundances)
fdiv = functional_divergence(traits, abundances)
rao_q = raos_quadratic_entropy(traits, abundances)
cwm = community_weighted_mean(traits, abundances)

# All metrics at once
suite = functional_diversity_suite(traits, abundances)
print(suite)

# Trait distance matrix
dm = trait_distance_matrix(traits, method="gower")
```

## Configuration

Environment variable prefix: `ECO_`

All algorithms are pure Python with no external dependencies beyond the standard library.

## Related Modules

- `metainformant.ecology.community` -- taxonomic diversity metrics
- `metainformant.ecology.ordination` -- ordination using trait distance matrices
- `metainformant.ecology.visualization` -- plotting functional diversity results
