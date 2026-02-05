# Multiomics Methods

Matrix factorization, clustering, and network fusion methods for multi-omic integration.

## 📦 Components

### factorization.py

Matrix factorization methods for discovering latent factors:

| Function | Description |
|----------|-------------|
| `joint_nmf()` | Non-negative matrix factorization across omics |
| `mofa_factors()` | Multi-Omics Factor Analysis |
| `tensor_decomposition()` | CP/Tucker tensor factorization |

### clustering.py

Multi-omic clustering methods:

| Function | Description |
|----------|-------------|
| `spectral_clustering()` | Graph-based clustering |
| `snf_clustering()` | Similarity Network Fusion |
| `consensus_clustering()` | Ensemble approach |

## 🚀 Usage

```python
from metainformant.multiomics.methods import factorization, clustering

# Joint NMF
factors = factorization.joint_nmf(
    {"rna": rna_matrix, "protein": protein_matrix},
    n_components=20
)

# SNF clustering
clusters = clustering.snf_clustering(
    {"rna": rna_data, "dna": variant_data},
    n_clusters=4
)
```

## 🔗 Related

- [analysis/integration.py](../analysis/integration.py) - Core integration
- [pathways/](../pathways/) - Pathway-level analysis
