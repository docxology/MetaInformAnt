# Ordination Methods

Multivariate ordination techniques for visualizing and analyzing patterns in ecological community data. All algorithms are implemented from first principles using numpy -- no scipy or scikit-bio dependencies required.

## Key Concepts

**PCoA (Principal Coordinates Analysis)** embeds a distance matrix into Euclidean space via Gower's double-centering and eigendecomposition. Preserves pairwise distances as faithfully as possible in low dimensions.

**NMDS (Non-metric Multidimensional Scaling)** finds a configuration whose rank-order of pairwise distances best matches the original. Uses Kruskal's stress minimization with isotonic regression.

**CCA (Canonical Correspondence Analysis)** is a constrained ordination that relates community composition to measured environmental variables via chi-square-weighted eigenanalysis.

**Procrustes rotation** optimally superimposes two ordination configurations by translation, uniform scaling, and rotation to minimize the sum of squared differences.

## Function Reference

### `distance_matrix(communities, method="bray_curtis") -> List[List[float]]`

Compute a symmetric pairwise distance matrix. Methods: `bray_curtis`, `jaccard`, `euclidean`, `manhattan`, `canberra`.

### `pcoa(distance_matrix, n_components=2) -> Dict[str, Any]`

Principal Coordinates Analysis via Gower's double-centering.

Returns `coordinates` (n x n_components), `eigenvalues`, and `variance_explained`.

### `nmds(distance_matrix, n_components=2, max_iter=300, n_init=4, tol=1e-7, seed=None) -> Dict[str, Any]`

Non-metric MDS using Kruskal's algorithm with multiple random initializations.

Returns `coordinates`, `stress` (0 = perfect), and `n_iter`.

### `cca(species_matrix, environmental_matrix, n_components=None) -> Dict[str, Any]`

Canonical Correspondence Analysis relating species to environmental variables.

Returns `site_scores`, `species_scores`, `eigenvalues`, and `variance_explained`.

### `procrustes(coords1, coords2) -> Dict[str, Any]`

Procrustes rotation to align two ordination configurations.

Returns `transformed_coords`, `m2` (Procrustes statistic), and `correlation`.

## Usage Examples

```python
from metainformant.ecology import distance_matrix, pcoa, nmds, cca, procrustes

# Compute Bray-Curtis distance matrix
communities = [[10, 0, 5], [0, 8, 2], [3, 3, 3], [1, 9, 0]]
dm = distance_matrix(communities, method="bray_curtis")

# PCoA ordination
result = pcoa(dm, n_components=2)
coords = result["coordinates"]       # 4 x 2 coordinates
var_exp = result["variance_explained"]  # per-axis variance

# NMDS ordination
nmds_result = nmds(dm, n_components=2, seed=42)
print(f"Stress: {nmds_result['stress']:.4f}")

# CCA with environmental data
species = [[10, 0, 5], [0, 8, 2], [3, 3, 3]]
env = [[1.0, 2.0], [3.0, 1.0], [2.0, 3.0]]
cca_result = cca(species, env)

# Compare two ordinations with Procrustes
proc = procrustes(result["coordinates"][:3], nmds_result["coordinates"][:3])
print(f"Procrustes correlation: {proc['correlation']:.4f}")
```

## Stress Interpretation (NMDS)

| Stress Value | Interpretation |
|-------------|----------------|
| < 0.05      | Excellent       |
| 0.05 - 0.10 | Good           |
| 0.10 - 0.20 | Fair           |
| > 0.20      | Poor           |

## Configuration

Environment variable prefix: `ECO_`

Requires `numpy` for matrix operations.

## Related Modules

- `metainformant.ecology.community` -- diversity metrics and distance calculations
- `metainformant.ecology.indicators` -- statistical tests on ordination groups
- `metainformant.ecology.visualization` -- plotting ordination results
