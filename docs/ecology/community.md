# Community Ecology

Community ecology analysis including alpha/beta diversity, rarefaction, species accumulation, and comprehensive biodiversity assessment.

## Key Concepts

**Alpha diversity** measures species diversity within a single community. Common indices include Shannon entropy (information-theoretic), Simpson (probability of interspecific encounter), and species richness (simple count).

**Beta diversity** measures compositional differences between communities using dissimilarity metrics such as Bray-Curtis (abundance-weighted), Jaccard (presence/absence), and Sorensen.

**Rarefaction** estimates expected richness at standardized sampling effort, enabling fair comparisons across communities with unequal sampling.

**Chao1** is a nonparametric richness estimator that uses singleton/doubleton counts to predict unobserved species.

## Function Reference

### `shannon_diversity(abundances: List[float]) -> float`

Calculate Shannon diversity index H' = -sum(p_i * ln(p_i)).

### `simpson_diversity(abundances: List[float]) -> float`

Calculate Simpson diversity index (1 - D), where D = sum(p_i^2).

### `species_richness(community_data) -> int | List[int] | Dict[str, int]`

Count species with nonzero abundance. Accepts flat lists, nested lists, or dicts.

### `species_richness_simple(abundances: List[float]) -> int`

Count of non-zero abundances in a single community vector.

### `pielou_evenness(abundances: List[float]) -> float`

Pielou's evenness J = H / ln(S). Returns value between 0 and 1.

### `chao1_estimator(abundances: List[float]) -> float`

Estimate total richness using Chao1: S_obs + f1^2 / (2*f2).

### `beta_diversity(community1, community2, method="bray_curtis") -> float`

Calculate beta diversity between two communities. Methods: `bray_curtis`, `jaccard`, `sorensen`. Returns 0 (identical) to 1 (completely different).

### `alpha_beta_gamma_diversity(communities: List[List[float]]) -> Dict[str, float]`

Partition diversity into alpha (within), beta (between), and gamma (total) components using Shannon index.

### `rarefaction_curve(abundances, max_samples=None) -> List[Tuple[int, float]]`

Generate rarefaction curve as (sample_size, expected_richness) tuples using hypergeometric expectations.

### `species_accumulation_curve(sampling_effort: List[int]) -> List[Tuple[int, float]]`

Smooth cumulative species counts into an accumulation curve.

### `community_similarity_matrix(communities, method="bray_curtis") -> List[List[float]]`

Compute pairwise similarity matrix (1 - dissimilarity) for all community pairs.

### `community_metrics(abundances: List[float]) -> Dict[str, float]`

Comprehensive metrics in one call: `shannon`, `simpson`, `richness`, `pielou`, `chao1`.

### `calculate_biodiversity_indices(community_data, indices=None) -> Dict[str, List[float]]`

Batch-calculate multiple diversity indices across communities.

### `nestedness_temperature_calculator(presence_absence_matrix) -> float`

Nestedness temperature (0 = perfectly nested, 100 = random).

### `species_area_relationship(species_counts, area_sizes) -> Dict[str, Any]`

Fit power-law species-area model S = c * A^z via log-linear regression.

### `generate_ecology_report(community_data, sample_names=None, output_path=None) -> str`

Generate a formatted text report with diversity indices, alpha-beta-gamma decomposition, and similarity statistics.

## Usage Examples

```python
from metainformant.ecology import (
    shannon_diversity, simpson_diversity, beta_diversity,
    community_metrics, rarefaction_curve, alpha_beta_gamma_diversity,
)

# Single community metrics
abundances = [45, 23, 12, 8, 5, 3, 2, 1, 1]
h = shannon_diversity(abundances)        # 1.72
d = simpson_diversity(abundances)        # 0.78
metrics = community_metrics(abundances)  # all five indices

# Beta diversity between two sites
site_a = [10, 5, 0, 3]
site_b = [0, 6, 8, 2]
bc = beta_diversity(site_a, site_b, method="bray_curtis")

# Rarefaction
curve = rarefaction_curve(abundances, max_samples=50)

# Alpha-beta-gamma partitioning
communities = [site_a, site_b, [7, 7, 7, 7]]
abg = alpha_beta_gamma_diversity(communities)
print(abg["alpha"], abg["beta"], abg["gamma"])
```

## Configuration

Environment variable prefix: `ECO_`

No external dependencies required for core community functions.

## Related Modules

- `metainformant.ecology.ordination` -- ordination of community distance matrices
- `metainformant.ecology.indicators` -- indicator species and hypothesis testing
- `metainformant.ecology.functional` -- trait-based functional diversity
- `metainformant.information` -- information-theoretic entropy measures
