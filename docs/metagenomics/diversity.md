# Metagenomics Diversity Metrics

Alpha diversity (within-sample), beta diversity (between-sample), rarefaction analysis, PERMANOVA hypothesis testing, and ordination methods for metagenomic community data.

## Key Concepts

**Alpha diversity** measures species richness and evenness within a single community. Eight metrics are supported: Shannon entropy, Simpson index, inverse Simpson, Chao1 richness estimator, ACE (Abundance-based Coverage Estimator), observed species count, Fisher's alpha, and Pielou evenness.

**Beta diversity** quantifies compositional differences between communities using distance metrics. Supports Bray-Curtis (abundance-weighted), Jaccard (presence/absence), and Aitchison (CLR-transformed Euclidean) distances.

**Rarefaction** subsamples communities to standardized depths to assess whether sequencing effort is sufficient to capture full diversity. Saturation is detected when the curve plateaus.

**PERMANOVA** tests whether community composition differs between groups using a permutation-based pseudo-F statistic on the distance matrix.

**Ordination** reduces high-dimensional distance matrices to 2-3 dimensions for visualization via PCoA (eigendecomposition) or NMDS (stress minimization).

## Function Reference

### `alpha_diversity(abundances, metric="shannon") -> Dict`

Compute within-sample diversity.

**Metrics:**
| Metric | Description |
|--------|-------------|
| `shannon` | H' = -sum(p_i * ln(p_i)) |
| `simpson` | D = sum(p_i^2) (dominance) |
| `invsimpson` | 1/D (effective species) |
| `chao1` | S_obs + f1^2/(2*f2) |
| `ace` | Abundance-based Coverage Estimator |
| `observed` | Count of non-zero taxa |
| `fisher_alpha` | Fisher's log-series alpha |
| `pielou_evenness` | H'/ln(S) |

**Returns** dict with `value`, `metric`, `n_species`, `total_count`.

### `beta_diversity(samples, metric="bray_curtis") -> Dict`

Compute pairwise distance matrix between samples.

**Metrics:**
| Metric | Description |
|--------|-------------|
| `bray_curtis` | sum|a-b| / sum(a+b) |
| `jaccard` | 1 - |A and B| / |A or B| |
| `aitchison` | Euclidean on CLR-transformed data |

**Returns** dict with `distance_matrix`, `metric`, `n_samples`.

### `rarefaction_curve(abundances, depths=None, n_iterations=10, seed=None) -> Dict`

Generate a rarefaction curve by subsampling at increasing depths.

**Returns** dict with `depths`, `mean_species`, `std_species`, and `is_saturated` (True if last three points are within 5% of each other).

### `rarefy(abundances, depth, seed=None) -> List[int]`

Subsample a community to a specific depth without replacement. Returns rarefied count vector summing to `depth`.

### `permanova(distance_matrix, groups, n_permutations=999, seed=None) -> Dict`

Permutational MANOVA for testing group differences in community composition.

**Returns** dict with `pseudo_f`, `p_value`, `r_squared`, `n_permutations`.

### `ordination(distance_matrix, method="pcoa", n_components=2) -> Dict`

Dimensionality reduction for visualization.

**Methods:**
- `pcoa`: Principal Coordinates Analysis via double-centering eigendecomposition.
- `nmds`: Non-metric MDS via gradient-based stress minimization.

**Returns** dict with `coordinates` (n_samples x n_components), `eigenvalues`, `variance_explained`.

## Usage Examples

```python
from metainformant.metagenomics import (
    alpha_diversity, beta_diversity,
    rarefaction_curve, rarefy,
    permanova, ordination,
)

# Alpha diversity
counts = [100, 50, 30, 20, 10, 5, 3, 1, 1]
shannon = alpha_diversity(counts, metric="shannon")
print(f"Shannon H': {shannon['value']:.3f}")

chao1 = alpha_diversity(counts, metric="chao1")
print(f"Chao1 estimated richness: {chao1['value']:.1f}")

# Beta diversity
samples = [
    [100, 50, 30, 0, 0],
    [0, 60, 40, 20, 10],
    [80, 40, 25, 5, 0],
]
bc = beta_diversity(samples, metric="bray_curtis")
print(f"Distance matrix: {bc['distance_matrix']}")

# Rarefaction
rare = rarefaction_curve(counts, n_iterations=20, seed=42)
print(f"Saturated: {rare['is_saturated']}")

# Rarefy to even depth
rarefied = rarefy(counts, depth=100, seed=42)

# PERMANOVA
groups = ["forest", "forest", "grassland"]
result = permanova(bc["distance_matrix"], groups, n_permutations=999)
print(f"pseudo-F={result['pseudo_f']:.3f}, p={result['p_value']:.4f}")

# PCoA ordination
ord_result = ordination(bc["distance_matrix"], method="pcoa")
coords = ord_result["coordinates"]
var_exp = ord_result["variance_explained"]
```

## Configuration

Environment variable prefix: `META_`

Requires `numpy` when available for PCoA eigendecomposition; pure Python power iteration fallback otherwise.

## Related Modules

- `metainformant.metagenomics.amplicon` -- OTU/ASV data input
- `metainformant.metagenomics.shotgun` -- profiled community input
- `metainformant.metagenomics.comparative` -- differential abundance testing
- `metainformant.ecology.community` -- equivalent metrics for macroecological data
