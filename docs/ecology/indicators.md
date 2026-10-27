# Indicator Species and Multivariate Tests

Classical ecological statistical methods for analyzing community composition across groups, including indicator species analysis (IndVal), permutational hypothesis testing (ANOSIM, PERMANOVA), hierarchical clustering, SIMPER decomposition, and tests for homogeneity of dispersions (PERMDISP).

All functions live in `metainformant.ecology.analysis.indicators`. Abundance matrices are sites by species (`abundance_matrix[i][j]` is the abundance of species `j` at site `i`), distance matrices are square lists of lists (build one with `distance_matrix` from `metainformant.ecology.analysis.ordination`), and `group_labels` holds one label per site or row. `indval`, `permanova`, `simper` and `multivariate_dispersion` need at least two distinct groups and raise `ValueError` otherwise (`anosim` returns R = 0 for a single group).

## Key Concepts

**IndVal** (Dufrene and Legendre, 1997) combines specificity (relative abundance concentration) and fidelity (frequency of occurrence) to identify indicator species for habitat or treatment groups.

**ANOSIM** (Clarke, 1993) tests whether between-group dissimilarities are larger than within-group dissimilarities using a rank-based R statistic with permutation testing.

**PERMANOVA** (Anderson, 2001) partitions the total sum of squared distances into within-group and between-group components, testing significance with a pseudo-F statistic.

**SIMPER** (Clarke, 1993) decomposes Bray-Curtis dissimilarity between groups to identify which species contribute most to observed differences.

**PERMDISP** (Anderson, 2006) tests whether group dispersions (distances to centroid) differ, complementing PERMANOVA by checking the homogeneity assumption.

## Function Reference

The permutation tests (`indval`, `anosim`, `permanova`, `multivariate_dispersion`) take `n_permutations` and `seed` as keyword-only arguments. Each p-value is `(k + 1) / (n_permutations + 1)`, where `k` counts the permutations that gave a statistic at least as large as the observed one, so the smallest possible p-value is `1 / (n_permutations + 1)`, and with few sites there are few distinct label arrangements, which keeps p-values from getting small. The fixed default `seed=42` makes results reproducible.

### `indval(abundance_matrix, group_labels, *, n_permutations=999, seed=42) -> List[Dict]`

Indicator Value analysis. Returns one dict per species with `species_idx`, `indval` (0-100), `specificity`, `fidelity`, `best_group`, and `p_value`.

### `anosim(distance_matrix, group_labels, *, n_permutations=999, seed=42) -> Dict`

Analysis of Similarities. Returns `r_statistic` (-1 to 1, higher = stronger separation), `p_value`, and `n_permutations`.

### `permanova(distance_matrix, group_labels, *, n_permutations=999, seed=42) -> Dict`

Permutational MANOVA. Returns `f_statistic`, `p_value`, `r_squared`, and `n_permutations`.

### `simper(abundance_matrix, group_labels) -> List[Dict]`

Similarity Percentages. Returns one dict per species sorted by `contribution_pct`, each with `species_idx`, `avg_dissimilarity`, `cumulative_pct`, `sd`, and `avg_abundance_by_group`.

### `cluster_communities(distance_matrix, method="upgma", n_clusters=None) -> Dict`

Hierarchical agglomerative clustering. Methods: `upgma`, `single`, `complete`. Returns `dendrogram` (one dict per merge with `cluster_a`, `cluster_b`, `distance` and `size`), `cluster_labels` and `cophenetic_correlation`. `cluster_labels` is only meaningful when `n_clusters` is given (1 to the number of sites); otherwise every site gets its own label.

### `multivariate_dispersion(distance_matrix, group_labels, *, n_permutations=999, seed=42) -> Dict`

PERMDISP test. Returns `f_statistic`, `p_value`, `group_dispersions` (group label -> mean distance to the group centroid), and `n_permutations`.

## Usage Examples

```python
from metainformant.ecology.analysis.indicators import (
    indval, anosim, permanova, simper,
    cluster_communities, multivariate_dispersion,
)
from metainformant.ecology.analysis.ordination import distance_matrix

# Sample data: sites x species (five forest sites, then five grassland sites)
abundance = [
    [10, 0, 3, 1], [8, 1, 2, 0], [12, 0, 4, 2], [9, 2, 1, 1], [11, 1, 3, 0],
    [0, 7, 1, 1], [1, 9, 0, 2], [0, 8, 2, 0], [2, 10, 1, 1], [0, 6, 0, 2],
]
labels = ["forest"] * 5 + ["grassland"] * 5

# Indicator species
results = indval(abundance, labels, n_permutations=999)
for r in results:
    if r["p_value"] < 0.05:
        print(f"Species {r['species_idx']}: IndVal={r['indval']}, group={r['best_group']}")
# Species 0: IndVal=94.3396, group=forest
# Species 1: IndVal=90.9091, group=grassland

# Compute distance matrix for hypothesis tests
dm = distance_matrix(abundance, method="bray_curtis")

# ANOSIM
anosim_result = anosim(dm, labels)
print(f"R = {anosim_result['r_statistic']}, p = {anosim_result['p_value']}")  # R = 1.0, p = 0.01

# PERMANOVA
perm_result = permanova(dm, labels)
print(f"F = {perm_result['f_statistic']}, R2 = {perm_result['r_squared']}")  # F = 57.632306, R2 = 0.878109

# SIMPER decomposition
simper_result = simper(abundance, labels)
for sp in simper_result[:3]:
    print(f"Species {sp['species_idx']}: {sp['contribution_pct']:.1f}%")  # species 0, 1, 2

# Clustering
clusters = cluster_communities(dm, method="upgma", n_clusters=2)
print(f"Labels: {clusters['cluster_labels']}")  # [0, 0, 0, 0, 0, 1, 1, 1, 1, 1]

# Dispersion test
disp = multivariate_dispersion(dm, labels)
print(f"Dispersions: {disp['group_dispersions']}")  # one mean distance to centroid per group
```

## Configuration

Settings are passed as function arguments; the indicator and test functions read no environment variables.

All algorithms are pure Python with no external dependencies beyond the standard library.

## Related Modules

- `metainformant.ecology.analysis.community` -- diversity metrics used by these tests
- `metainformant.ecology.analysis.ordination` -- ordination of distance matrices
- `metainformant.metagenomics.diversity` -- diversity metrics and PERMANOVA for microbial data
