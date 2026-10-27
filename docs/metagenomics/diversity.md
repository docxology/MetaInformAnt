# Metagenomics Diversity Metrics

Alpha diversity (within-sample), beta diversity (between-sample), rarefaction analysis, PERMANOVA hypothesis testing, and ordination methods for metagenomic community data.

## Key Concepts

**Alpha diversity** measures species richness and evenness within a single community. Eight metrics are supported: Shannon entropy, Simpson index, inverse Simpson, Chao1 richness estimator, ACE (Abundance-based Coverage Estimator), observed species count, Fisher's alpha, and Pielou evenness.

**Beta diversity** quantifies compositional differences between communities using distance metrics. Supports Bray-Curtis (abundance-weighted), Jaccard (presence/absence), and Aitchison (CLR-transformed Euclidean, pseudocount 0.5) distances.

**Rarefaction** subsamples communities to standardized depths to assess whether sequencing effort is sufficient to capture full diversity. Saturation is detected when the curve plateaus.

**PERMANOVA** tests whether community composition differs between groups using a permutation-based pseudo-F statistic on the distance matrix.

**Ordination** reduces high-dimensional distance matrices to 2-3 dimensions for visualization via PCoA (eigendecomposition) or NMDS (stress minimization).

## Function Reference

All functions live in `metainformant.metagenomics.diversity.metrics`; the `metainformant.metagenomics` package itself exposes only submodules. Abundance vectors and distance matrices are plain Python lists (`alpha_diversity` also takes a NumPy array or a pandas Series).

### `alpha_diversity(abundances, metric="shannon") -> Dict`

Compute within-sample diversity of one sample. Zero entries are ignored. `abundances` is a list, a NumPy array or a pandas Series; an empty vector or an unknown `metric` raises `ValueError`.

**Metrics:**
| Metric | Description |
|--------|-------------|
| `shannon` | H' = -sum(p_i * ln(p_i)) |
| `simpson` | D = sum(p_i^2) (dominance) |
| `invsimpson` | 1/D (effective species) |
| `chao1` | S_obs + f1^2/(2*f2); S_obs + f1*(f1-1)/2 when f2 = 0 |
| `ace` | Abundance-based Coverage Estimator (rare taxa: abundance <= 10) |
| `observed` | Count of non-zero taxa |
| `fisher_alpha` | Fisher's log-series alpha: the root of `S = alpha * ln(1 + N / alpha)` for `S` taxa and `N` individuals, found by bracketed bisection; `inf` when every taxon is a singleton (`N = S`), where no finite root exists |
| `pielou_evenness` | H'/ln(S) (1.0 when S <= 1) |

`shannon`, `simpson`, `invsimpson`, `pielou_evenness` and `observed` accept counts or proportions. `chao1`, `ace` and `fisher_alpha` are defined on integer counts (singletons, doubletons, a log-series), so they raise `ValueError` for any abundance that is not an integer (within 1e-9; floats such as `3.0` are fine). Relative abundances (proportions) therefore cannot be used with them; supply raw or rarefied counts.

**Returns** dict with `value`, `metric`, `n_species` (non-zero taxa), `total_count` (sum of the abundances).

### `beta_diversity(samples, metric="bray_curtis") -> Dict`

Compute pairwise distance matrix between samples. `samples` is a list of equal-length abundance vectors (at least two); otherwise `ValueError`.

**Metrics:**
| Metric | Description |
|--------|-------------|
| `bray_curtis` | sum(abs(a - b)) / sum(a + b) |
| `jaccard` | 1 - (taxa present in both) / (taxa present in either) |
| `aitchison` | Euclidean on CLR-transformed data |

**Returns** dict with `distance_matrix` (list of lists), `metric`, `n_samples`.

### `rarefaction_curve(abundances, depths=None, n_iterations=10, seed=None) -> Dict`

Generate a rarefaction curve by subsampling (`n_iterations` random draws per depth). The draws are unseeded by default (`seed=None`), so the means and standard deviations vary between runs; pass an integer `seed` for reproducible results. If `depths` is None, 20 evenly spaced depths up to the total count are used (fewer when the total is below 20, because duplicate depths are merged); depths above the total are dropped. Raises `ValueError` for an empty vector, negative counts, or a zero total.

**Returns** dict with `depths`, `mean_species`, `std_species`, and `is_saturated` (True if last three points are within 5% of each other).

### `rarefy(abundances, depth, seed=None) -> List[int]`

Subsample a community to a specific depth without replacement. Returns a rarefied count vector of the same length as the input, summing to `depth`. The draw is unseeded by default (`seed=None`); pass an integer `seed` for a reproducible subsample. Raises `ValueError` if `depth` is negative or exceeds the total count.

### `permanova(distance_matrix, groups, n_permutations=999, seed=None) -> Dict`

Permutational MANOVA for testing group differences in community composition. `groups` holds one label per sample (at least two distinct labels). `p_value` is (1 + number of permutations whose pseudo-F is at least the observed value) / (`n_permutations` + 1). The permutations are unseeded by default (`seed=None`), so the p-value varies between runs; pass an integer `seed` for a reproducible p-value.

**Returns** dict with `pseudo_f`, `p_value`, `r_squared`, `n_permutations`.

### `ordination(distance_matrix, method="pcoa", n_components=2, seed=42) -> Dict`

Dimensionality reduction for visualization.

**Methods:**
- `pcoa`: Principal Coordinates Analysis via double-centering eigendecomposition. At most `n_samples - 1` components are returned.
- `nmds`: gradient-descent stress minimization against the raw distances, starting from random coordinates drawn with `seed` (default 42, so repeated calls give identical results; `None` draws fresh entropy). It has no rank-order (isotonic) step, so it behaves like metric MDS rather than classical non-metric MDS. PCoA is deterministic and ignores `seed`.

**Returns** dict with `coordinates` (n_samples lists of n_components floats), `eigenvalues`, `variance_explained`. For `pcoa`, `variance_explained` is normalized over the returned components only (it sums to 1 unless no returned eigenvalue is positive), not over all positive eigenvalues, so it overstates each axis's share of the total variance. For `nmds`, `eigenvalues` holds the last stress values of the optimisation and `variance_explained` is `1 - final stress` repeated for each component.

## Usage Examples

```python
from metainformant.metagenomics.diversity.metrics import (
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
result = permanova(bc["distance_matrix"], groups, n_permutations=999, seed=42)
print(f"pseudo-F={result['pseudo_f']:.3f}, p={result['p_value']:.4f}")

# PCoA ordination
ord_result = ordination(bc["distance_matrix"], method="pcoa")
coords = ord_result["coordinates"]
var_exp = ord_result["variance_explained"]
```

## Configuration

This module reads no environment variables; the `META_` prefix used elsewhere in the package does not apply here.

PCoA uses `numpy.linalg.eigh` when numpy is importable and falls back to pure-Python power iteration otherwise.

## Related Modules

- `metainformant.metagenomics.amplicon` -- OTU/ASV data input
- `metainformant.metagenomics.shotgun` -- profiled community input
- `metainformant.metagenomics.comparative` -- differential abundance testing
- `metainformant.ecology.analysis.community` -- equivalent metrics for macroecological data
