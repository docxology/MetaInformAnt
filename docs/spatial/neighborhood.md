# Neighborhood Analysis

Quantify spatial relationships between cell types or gene expression patterns.

The functions in `metainformant.spatial.analysis.neighborhood` take plain arrays: one cell-type label per spot or cell (`cell_types`), the matching coordinates (`n × 2`) and, where needed, an expression matrix or a spatial graph.

## Cell-type colocalisation (if deconvolution done)

After running `deconvolve_spots()`, `deconv.fractions` holds spot × cell-type
proportions (see [Deconvolution](deconvolution.md)); the dominant type per spot is one
label per spot:

```python
import numpy as np
cell_types = np.array(deconv.cell_type_names)[deconv.fractions.argmax(axis=1)]
```

### Ripley's K (spread of a single type)

```python
import matplotlib.pyplot as plt
import numpy as np
from metainformant.spatial.analysis.neighborhood import ripley_k

points = coordinates[cell_types == 'Epithelial']                         # n x 2 points of one type
area = float(np.ptp(points[:, 0]) * np.ptp(points[:, 1]))                # area of the bounding box of `points`
res = ripley_k(
    points,
    radii=[50, 100, 150, 200],   # in the units of the coordinates
    area=area,
    n_simulations=99,            # CSR simulations for the envelope
    seed=42,
    edge_correction='isotropic', # the default; 'none' counts every pair with weight 1
)
plt.plot(res.radii, res.l_values, label='Epithelial L(r)')
plt.fill_between(res.radii, res.csr_envelope_lower, res.csr_envelope_upper, alpha=0.2, label='CSR envelope')
plt.legend()
```

`points`, `radii` and `area` are required (`area` is not computed for you). Under complete spatial randomness (CSR) K(r) = π r². The result's `l_values` hold Besag's L(r) = sqrt(K(r)/π) − r, which is 0 under CSR: positive values mean clustering, negative values regularity.

**Edge correction.** `edge_correction='isotropic'` (the default) applies Ripley's isotropic correction at every number of points: each ordered pair (i, j) at distance d_ij ≤ r counts with weight 1 / (the share of the circle of radius d_ij around point i that lies inside the study window), the exact share for a rectangle, and `K(r) = area / (n (n − 1)) × the sum of the weights`.

- The window is always the bounding box of `points`, so pass the area of that box as `area`; a different `area` only rescales K (the intensity), it does not change the window.
- `edge_correction='none'` gives every pair weight 1, which is biased low near the window boundary. Any other value raises `ValueError`.
- Keep the radii well below half the shorter side of the window (a quarter is a common rule): the weights grow quickly beyond that. Memory grows with the number of point pairs within the largest radius.
- `ValueError` is also raised for empty or negative `radii`, and for the isotropic correction on points with no extent in x or in y (there is no window to correct for; use `edge_correction='none'` for such a pattern).

`csr_envelope_lower` and `csr_envelope_upper` are the 2.5th and 97.5th percentiles of L(r) over `n_simulations` CSR point sets drawn uniformly in the bounding box of `points` (the same number of points each, with the same edge correction), and `RipleyKResult` also holds `k_values`, `n_points` and `area`.

### Interaction index (pairwise)

```python
from metainformant.spatial.analysis.clustering import build_spatial_graph
from metainformant.spatial.analysis.neighborhood import compute_interaction_matrix

graph = build_spatial_graph(coordinates, method='radius', radius=100.0)   # same units as coordinates
II = compute_interaction_matrix(
    cell_types,                      # one label per spot / cell
    graph,
    normalize_by_type_frequency=True,
)
# II.interaction_matrix[i, j] > 1 → more edges than expected under random mixing; < 1 → fewer
# II.cell_type_names gives the row / column order; II.method is 'frequency_normalized' or 'raw_count'
```

Edge counts between two cell types, observed / expected under independence (from the
cell‑type frequencies). The counts are sums of the graph's edge weights, so a radius or
Delaunay graph (weight 1) gives plain edge counts, while a `'knn'` graph weighs each edge
`1/(1+distance)`. Significance by Monte‑Carlo randomisation of the labels is
provided by `neighborhood_enrichment()` below.

## Ligand–receptor communication scoring

Ligand expression in a spot, receptor expression in its neighbours within a radius.

```python
from metainformant.spatial.analysis.neighborhood import ligand_receptor_spatial

lr = ligand_receptor_spatial(
    expression,
    [('EGF', 'EGFR')],
    coordinates,
    gene_names=gene_names,
    radius=200.0,
)
print(lr['scores'])  # {('EGF', 'EGFR'): mean spatial interaction score}
```

`lr['per_spot_scores']` holds the per‑spot values and `lr['n_pairs_tested']` the number of pairs whose genes are both in `gene_names` (others are skipped). Without `radius` the `n_neighbors` (default 6) nearest spots are used. The score is not restricted to particular cell types and has no permutation test; for a label‑permutation test of cell‑type neighbourhoods use `neighborhood_enrichment()` below. See [Communication](communication.md) for the cell-type level scoring in `spatial.communication`.

## Nearest‑neighbour enrichment

```python
import pandas as pd
from metainformant.spatial.analysis.neighborhood import neighborhood_enrichment

nne = neighborhood_enrichment(
    cell_types,
    coordinates,
    n_neighbors=5,
    n_permutations=1000,
    seed=42,
)
# (n_types × n_types) matrices ordered like nne.cell_type_names:
# enrichment_matrix (z-scores), count_matrix, expected_matrix, p_values
print(pd.DataFrame(nne.enrichment_matrix, index=nne.cell_type_names, columns=nne.cell_type_names))
```

A positive z‑score means type‑A spots have type‑B neighbours more often than expected
when the cell‑type labels are permuted: z = (observed − mean of the permutations) / their standard deviation. `p_values` is the fraction of permutations whose count is at least the observed one. Pass `radius=` to define neighbours by distance instead of the `n_neighbors` nearest spots.

## Niche detection

`niche_detection(cell_types, coordinates, n_niches=5, *, n_neighbors=10, seed=42)` describes every cell by the cell-type composition of the cell and its `n_neighbors` nearest neighbours, and clusters those composition vectors with k-means:

```python
from metainformant.spatial.analysis.neighborhood import niche_detection

niches = niche_detection(cell_types, coordinates, n_niches=5, n_neighbors=10, seed=42)
# niches.niche_labels, niches.niche_compositions (n_niches × n_types), niches.cell_type_names
```

To cluster cell-type proportions instead of labels, see [Niche Analysis](niche.md).

## Workflow end‑to‑end

```python
import numpy as np

# 1. Deconvolution
from metainformant.spatial.analysis.deconvolution import deconvolve_spots
deconv = deconvolve_spots(expression, profiles, method='nnls', cell_type_names=cell_type_names)
cell_types = np.array(deconv.cell_type_names)[deconv.fractions.argmax(axis=1)]  # dominant type per spot

# 2. Colocalisation
from metainformant.spatial.analysis.clustering import build_spatial_graph
from metainformant.spatial.analysis.neighborhood import compute_interaction_matrix
from metainformant.spatial.visualization.plots import plot_interaction_heatmap
graph = build_spatial_graph(coordinates, method='radius', radius=100.0)
im = compute_interaction_matrix(cell_types, graph)
plot_interaction_heatmap(im.interaction_matrix, 'interaction_heatmap.png', cell_type_names=im.cell_type_names)

# 3. LR communication (top pairs)
from metainformant.spatial.analysis.neighborhood import ligand_receptor_spatial
from metainformant.spatial.communication.cell_communication import default_lr_database
lr_pairs = [(p['ligand'], p['receptor']) for p in default_lr_database()['pairs']]
lr = ligand_receptor_spatial(expression, lr_pairs, coordinates, gene_names=gene_names, radius=150.0)
top_pairs = sorted(lr['scores'].items(), key=lambda kv: kv[1], reverse=True)[:10]
print(top_pairs)
```

## Cost

All of these functions are single‑threaded Python loops; there is no `n_jobs` option.

- `neighborhood_enrichment` repeats a loop over every spot and its neighbours `n_permutations` times (default 1000); lower it for exploration.
- `ripley_k` evaluates K(r) once for the data and once per simulation (`n_simulations`, default 99); the cost per evaluation grows with the number of point pairs within the largest radius.
- `ligand_receptor_spatial` loops over the spots once per ligand–receptor pair.

See [Performance](PERFORMANCE.md) for measured timings of the graph and autocorrelation functions.

## References

- Ripley, B. D. (1976). The second-order analysis of stationary point processes. *Journal of Applied Probability* 13(2), 255–266.
