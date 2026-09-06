# Multiomics Methods

Matrix factorization, clustering, and network fusion methods for multi-omic integration.

## Components

### factorization.py

Matrix factorization and fusion algorithms (numpy or pure-Python fallback):

| Function | Description |
|----------|-------------|
| `joint_nmf()` | Joint non-negative matrix factorization across omic layers (shared W, per-layer H) |
| `mofa_simple()` | Simplified MOFA+ Bayesian factor model with ARD priors |
| `tensor_decomposition()` | CP (CANDECOMP/PARAFAC) tensor decomposition via ALS |
| `similarity_network_fusion()` | Similarity Network Fusion of patient similarity networks |
| `canonical_correlation()` | Regularized (ridge) CCA for two omic matrices via SVD |

### clustering.py

Multi-omic clustering methods:

| Function | Description |
|----------|-------------|
| `multi_omic_clustering()` | Sample clustering via `"snf"`, `"concatenation"`, or `"late_integration"` |
| `consensus_clustering()` | Resampling-based consensus matrix and PAC-score selection of optimal k |
| `multi_view_spectral()` | Spectral clustering on averaged/product/max-combined similarity views |
| `evaluate_integration()` | Per-omic silhouette, ARI vs. integrated labels, aggregate metric |

## Usage

```python
from metainformant.multiomics.methods import clustering, factorization

# Joint NMF (non-negative matrices sharing the sample dimension)
result = factorization.joint_nmf(
    {"rna": rna_matrix, "protein": protein_matrix},
    k=10,
)
W, H_dict = result["W"], result["H_dict"]

# SNF-based multi-omic clustering
clusters = clustering.multi_omic_clustering(
    {"rna": rna_matrix, "protein": protein_matrix},
    n_clusters=4,
    method="snf",
)

# Consensus clustering to select the optimal k
consensus = clustering.consensus_clustering(rna_matrix, k_range=[2, 3, 4, 5])
```

## Related

- [analysis/integration.py](../analysis/integration.py) - Core integration
- [pathways/](../pathways/) - Pathway-level analysis
