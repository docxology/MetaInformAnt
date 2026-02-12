# Indicator Species and Multivariate Tests

Classical ecological statistical methods for analyzing community composition across groups, including indicator species analysis (IndVal), permutational hypothesis testing (ANOSIM, PERMANOVA), hierarchical clustering, SIMPER decomposition, and tests for homogeneity of dispersions (PERMDISP).

## Key Concepts

**IndVal** (Dufrene and Legendre, 1997) combines specificity (relative abundance concentration) and fidelity (frequency of occurrence) to identify indicator species for habitat or treatment groups.

**ANOSIM** (Clarke, 1993) tests whether between-group dissimilarities are larger than within-group dissimilarities using a rank-based R statistic with permutation testing.

**PERMANOVA** (Anderson, 2001) partitions the total sum of squared distances into within-group and between-group components, testing significance with a pseudo-F statistic.

**SIMPER** (Clarke, 1993) decomposes Bray-Curtis dissimilarity between groups to identify which species contribute most to observed differences.

**PERMDISP** (Anderson, 2006) tests whether group dispersions (distances to centroid) differ, complementing PERMANOVA by checking the homogeneity assumption.

## Function Reference

### `indval(abundance_matrix, group_labels, n_permutations=999, seed=42) -> List[Dict]`

Indicator Value analysis. Returns per-species dicts with `indval` (0-100), `specificity`, `fidelity`, `best_group`, and `p_value`.

### `anosim(distance_matrix, group_labels, n_permutations=999, seed=42) -> Dict`

Analysis of Similarities. Returns `r_statistic` (-1 to 1, higher = stronger separation), `p_value`, and `n_permutations`.

### `permanova(distance_matrix, group_labels, n_permutations=999, seed=42) -> Dict`

Permutational MANOVA. Returns `f_statistic`, `p_value`, `r_squared`, and `n_permutations`.

### `simper(abundance_matrix, group_labels) -> List[Dict]`

Similarity Percentages. Returns species sorted by `contribution_pct` with `avg_dissimilarity`, `cumulative_pct`, `sd`, and `avg_abundance_by_group`.

### `cluster_communities(distance_matrix, method="upgma", n_clusters=None) -> Dict`

Hierarchical agglomerative clustering. Methods: `upgma`, `single`, `complete`. Returns `dendrogram`, `cluster_labels`, and `cophenetic_correlation`.

### `multivariate_dispersion(distance_matrix, group_labels, n_permutations=999, seed=42) -> Dict`

PERMDISP test. Returns `f_statistic`, `p_value`, `group_dispersions`, and `n_permutations`.

## Usage Examples

```python
from metainformant.ecology import (
    indval, anosim, permanova, simper,
    cluster_communities, multivariate_dispersion,
    distance_matrix,
)

# Sample data: sites x species
abundance = [[10, 0, 3], [8, 1, 2], [0, 7, 1], [1, 9, 0]]
labels = ["forest", "forest", "grassland", "grassland"]

# Indicator species
results = indval(abundance, labels, n_permutations=999)
for r in results:
    if r["p_value"] < 0.05:
        print(f"Species {r['species_idx']}: IndVal={r['indval']}, group={r['best_group']}")

# Compute distance matrix for hypothesis tests
dm = distance_matrix(abundance, method="bray_curtis")

# ANOSIM
anosim_result = anosim(dm, labels)
print(f"R = {anosim_result['r_statistic']}, p = {anosim_result['p_value']}")

# PERMANOVA
perm_result = permanova(dm, labels)
print(f"F = {perm_result['f_statistic']}, R2 = {perm_result['r_squared']}")

# SIMPER decomposition
simper_result = simper(abundance, labels)
for sp in simper_result[:3]:
    print(f"Species {sp['species_idx']}: {sp['contribution_pct']:.1f}%")

# Clustering
clusters = cluster_communities(dm, method="upgma", n_clusters=2)
print(f"Labels: {clusters['cluster_labels']}")

# Dispersion test
disp = multivariate_dispersion(dm, labels)
print(f"Dispersions: {disp['group_dispersions']}")
```

## Configuration

Environment variable prefix: `ECO_`

All algorithms are pure Python with no external dependencies beyond the standard library.

## Related Modules

- `metainformant.ecology.community` -- diversity metrics used by these tests
- `metainformant.ecology.ordination` -- ordination of distance matrices
- `metainformant.metagenomics.diversity` -- similar tests for microbial data
