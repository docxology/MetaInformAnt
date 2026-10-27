# Community Ecology

Community ecology analysis including alpha/beta diversity, rarefaction, species accumulation, and comprehensive biodiversity assessment.

All functions live in `metainformant.ecology.analysis.community`; `calculate_diversity`, `species_richness` and `beta_diversity` are also importable from `metainformant.ecology`. A community is a plain list of abundances with one entry per species, and communities that are compared must list their species in the same order.

## Key Concepts

**Alpha diversity** measures species diversity within a single community. Common indices include Shannon entropy (information-theoretic), Simpson (probability of interspecific encounter), and species richness (simple count).

**Beta diversity** measures compositional differences between communities using dissimilarity metrics such as Bray-Curtis (abundance-weighted), Jaccard (presence/absence), and Sorensen.

**Rarefaction** estimates expected richness at standardized sampling effort, enabling fair comparisons across communities with unequal sampling.

**Chao1** is a nonparametric richness estimator that uses singleton/doubleton counts to predict unobserved species.

## Function Reference

### `shannon_diversity(abundances: List[float]) -> float`

Calculate Shannon diversity index H' = -sum(p_i * ln(p_i)) over species with positive abundance (natural logarithm). Returns 0.0 for an empty community.

### `simpson_diversity(abundances: List[float]) -> float`

Calculate Simpson diversity index (1 - D), where D = sum(p_i^2).

### `calculate_diversity(species_matrix, method="shannon") -> List[float] | Dict[str, List[float]]`

Apply one index to many communities. `species_matrix` is a list of abundance lists or a dict of name -> abundance list; `method` is `shannon`, `simpson`, `invsimpson` or `richness`, and anything else raises `ValueError`. A list input gives a list with one value per community; a dict input gives `{name: [value]}`.

### `species_richness(community_data) -> int | List[int] | Dict[str, int]`

Count species with nonzero abundance. A flat list gives an `int`, a list of lists gives a list of counts, and a dict of lists gives a dict of counts.

### `species_richness_simple(abundances: List[float]) -> int`

Count of non-zero abundances in a single community vector.

### `pielou_evenness(abundances: List[float]) -> float`

Pielou's evenness J = H / ln(S). Returns a value between 0 and 1: 1.0 for a community with a single species and 0.0 for an empty one.

### `chao1_estimator(abundances: List[float]) -> float`

Estimate total richness using Chao1: S_obs + f1^2 / (2*f2), where f1 and f2 count the species with abundance exactly 1 and exactly 2 (so abundances must be counts). Without doubletons it uses the bias-corrected form S_obs + f1*(f1-1)/2.

### `beta_diversity(community1, community2, method="bray_curtis") -> float`

Calculate beta diversity (dissimilarity) between two communities given as abundance vectors over the same species; the shorter vector is padded with zeros. Methods: `bray_curtis` (abundance-weighted), `jaccard` and `sorensen` (both presence/absence); any other method raises `ValueError`. Returns 0 (identical) to 1 (no species in common).

### `alpha_beta_gamma_diversity(communities: List[List[float]]) -> Dict[str, float]`

Partition diversity using the Shannon index: `alpha` is the mean Shannon H' of the communities, `gamma` is the Shannon H' of the pooled abundances (summed per species position), and `beta` is `gamma / alpha` (0.0 when alpha is 0). An empty list gives zeros for all three.

### `rarefaction_curve(abundances, max_samples=None) -> List[Tuple[int, float]]`

Generate rarefaction curve as (sample_size, expected_richness) tuples for sample sizes 1 to `max_samples`, using hypergeometric expectations. Abundances are truncated to whole individuals, and `max_samples` defaults to the total number of individuals (at most 1000) and is capped at it.

### `species_accumulation_curve(sampling_effort: List[int]) -> List[Tuple[int, float]]`

Pair each cumulative species count with its 1-based sampling effort as (effort, richness) tuples. The counts are returned as floats without smoothing; empty input gives `[(0, 0.0)]`.

### `community_similarity_matrix(communities, method="bray_curtis") -> List[List[float]]`

Compute pairwise similarity matrix (1 - dissimilarity) for all community pairs, with 1.0 on the diagonal. It accepts the same methods as `beta_diversity`.

### `community_metrics(abundances: List[float]) -> Dict[str, float]`

Comprehensive metrics in one call: `shannon`, `simpson`, `richness`, `pielou`, `chao1`.

### `calculate_biodiversity_indices(community_data, indices=None) -> Dict[str, List[float]]`

Batch-calculate multiple diversity indices across communities (a list of abundance lists). `indices` may name `shannon`, `simpson`, `invsimpson` and `richness` (all four by default); each result is a list with one value per community.

### `nestedness_temperature_calculator(presence_absence_matrix) -> float`

Nestedness temperature for a binary sites-by-species matrix: the percentage of site-pair and species comparisons in which a species occurs at the poorer site but not at the richer one. 0 means perfectly nested and larger values mean less nested; fewer than two sites or species gives 100.0.

### `species_area_relationship(species_counts, area_sizes) -> Dict[str, Any]`

Fit power-law species-area model S = c * A^z via log-linear regression. Counts and areas must be positive and there must be at least two of each. Returns `c`, `z`, `r_squared`, `intercept`, `slope`, `model` and `equation`.

### `generate_ecology_report(community_data, sample_names=None, output_path=None) -> str`

Generate a formatted text report with diversity indices, alpha-beta-gamma decomposition, and similarity statistics. It is written to `output_path` when one is given. `sample_names` is accepted but does not appear in the report text.

## Usage Examples

```python
from metainformant.ecology.analysis.community import (
    shannon_diversity, simpson_diversity, beta_diversity,
    community_metrics, rarefaction_curve, alpha_beta_gamma_diversity,
)

# Single community metrics
abundances = [45, 23, 12, 8, 5, 3, 2, 1, 1]
h = shannon_diversity(abundances)        # 1.58
d = simpson_diversity(abundances)        # 0.72
metrics = community_metrics(abundances)  # shannon, simpson, richness, pielou, chao1

# Beta diversity between two sites
site_a = [10, 5, 0, 3]
site_b = [0, 6, 8, 2]
bc = beta_diversity(site_a, site_b, method="bray_curtis")  # 0.59

# Rarefaction
curve = rarefaction_curve(abundances, max_samples=50)  # 50 points; curve[-1] is about (50, 7.6)

# Alpha-beta-gamma partitioning
communities = [site_a, site_b, [7, 7, 7, 7]]
abg = alpha_beta_gamma_diversity(communities)
print(abg["alpha"], abg["beta"], abg["gamma"])  # about 1.11 1.23 1.38
```

## Configuration

Settings are passed as function arguments; the community functions read no environment variables.

No external dependencies required for core community functions.

## Related Modules

- `metainformant.ecology.analysis.ordination` -- ordination of community distance matrices
- `metainformant.ecology.analysis.indicators` -- indicator species and hypothesis testing
- `metainformant.ecology.analysis.functional` -- trait-based functional diversity
- `metainformant.information` -- information-theoretic entropy measures
