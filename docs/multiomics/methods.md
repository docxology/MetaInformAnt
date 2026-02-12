# Multi-Omic Integration Methods

Matrix factorization, network fusion, canonical correlation, and clustering
methods for integrating multiple omic data layers. Includes joint NMF, MOFA,
CP tensor decomposition, Similarity Network Fusion, regularised CCA,
multi-omic clustering, consensus clustering, and multi-view spectral
clustering.

## Key Concepts

**Shared-sample assumption:** All integration methods require that omic
matrices share the same sample (row) dimension. Feature dimensions may differ.

**Joint NMF** factorises non-negative matrices V_i into a shared basis W and
per-omic coefficients H_i using multiplicative update rules.

**MOFA** (Multi-Omics Factor Analysis) uses a Bayesian factor model with ARD
(Automatic Relevance Determination) priors to learn sparse shared factors
across views via EM.

**SNF** (Similarity Network Fusion) iteratively diffuses per-omic patient
similarity networks through local neighbourhood kernels to produce a single
fused network.

**CCA** (Canonical Correlation Analysis) finds linear combinations of two
omic matrices that are maximally correlated, with ridge regularisation for
high-dimensional data.

## Factorization Functions

### joint_nmf

```python
def joint_nmf(
    data_matrices: dict[str, Any],
    k: int = 10,
    max_iter: int = 200,
    tol: float = 1e-4,
) -> dict[str, Any]
```

Returns `W` (shared factor matrix), `H_dict` (per-omic coefficients),
`reconstruction_error`, `n_iter`, `converged`. Requires non-negative input.

### mofa_simple

```python
def mofa_simple(
    data_matrices: dict[str, Any],
    k: int = 10,
    max_iter: int = 100,
    tol: float = 1e-3,
) -> dict[str, Any]
```

Returns `factors` (Z), `weights_per_view`, `variance_explained` per factor
per view, and `active_factors` (indices with ARD precision below threshold).
Requires numpy.

### tensor_decomposition

```python
def tensor_decomposition(
    tensor: list[list[list[float]]] | Any,
    rank: int = 5,
    method: str = "cp",
    max_iter: int = 100,
) -> dict[str, Any]
```

CP/CANDECOMP-PARAFAC decomposition of a 3-D tensor via alternating least
squares. Returns `factors` (list of three matrices), `fit` (variance
explained), and `core_consistency` (0--100%). Requires numpy.

### similarity_network_fusion

```python
def similarity_network_fusion(
    networks: list[list[list[float]]] | list[Any],
    k_neighbors: int = 20,
    n_iter: int = 20,
    alpha: float = 0.5,
) -> dict[str, Any]
```

Fuses >= 2 similarity networks via iterative diffusion. Returns
`fused_network`, `cluster_labels` (from spectral clustering), and
`silhouette_score`. Requires numpy.

### canonical_correlation

```python
def canonical_correlation(
    X: Any, Y: Any,
    n_components: int = 2,
    regularization: float = 0.1,
) -> dict[str, Any]
```

SVD-based regularised CCA for two omic matrices. Returns `x_scores`,
`y_scores`, `correlations`, `x_loadings`, `y_loadings`. Requires numpy.

## Clustering Functions

### multi_omic_clustering

```python
def multi_omic_clustering(
    data_matrices: dict[str, Any],
    n_clusters: int,
    method: str = "snf",
) -> dict[str, Any]
```

Integrates and clusters using `"snf"`, `"concatenation"`, or
`"late_integration"`. Returns `labels`, `silhouette`, and per-omic
`omic_contributions`.

### consensus_clustering

```python
def consensus_clustering(
    data: Any,
    k_range: range | list[int] | None = None,
    n_resamples: int = 100,
    proportion: float = 0.8,
) -> dict[str, Any]
```

Determines optimal k via resampling consensus matrices and PAC (Proportion of
Ambiguous Clustering) scores. Returns `optimal_k`, `labels`,
`consensus_matrix`, `cdf_area`, `pac_score`.

### multi_view_spectral

```python
def multi_view_spectral(
    similarity_matrices: list[Any],
    n_clusters: int,
    method: str = "average",
) -> dict[str, Any]
```

Combines similarity matrices via `"average"`, `"product"`, or `"max"`, then
performs spectral clustering. Returns `labels`, `eigenvalues`, `eigenvectors`.

### evaluate_integration

```python
def evaluate_integration(
    labels: list[int],
    omic_data: dict[str, Any],
) -> dict[str, Any]
```

Evaluates clustering quality across omic layers using per-omic silhouette
scores, Adjusted Rand Index (ARI), and a combined integration metric.

## Usage Example

```python
import numpy as np
from metainformant.multiomics.methods.factorization import (
    joint_nmf, similarity_network_fusion, canonical_correlation,
)
from metainformant.multiomics.methods.clustering import (
    multi_omic_clustering, consensus_clustering,
)

# Joint NMF
data = {"rna": np.random.rand(50, 100), "protein": np.random.rand(50, 80)}
result = joint_nmf(data, k=5, max_iter=100)
print(f"Converged: {result['converged']}, Error: {result['reconstruction_error']:.4f}")

# Multi-omic clustering
clust = multi_omic_clustering(data, n_clusters=3, method="concatenation")
print(f"Silhouette: {clust['silhouette']:.4f}")

# Consensus clustering
cc = consensus_clustering(np.random.randn(30, 10), k_range=range(2, 6))
print(f"Optimal k: {cc['optimal_k']}")
```

## Related Modules

- `metainformant.multiomics.pathways` -- multi-omic pathway enrichment
- `metainformant.multiomics.survival` -- survival analysis with multi-omic data
- `metainformant.multiomics.analysis` -- integration analysis utilities
