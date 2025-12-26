# Single-Cell Genomics

METAINFORMANT's single-cell genomics module provides comprehensive tools for analyzing single-cell RNA-sequencing (scRNA-seq) data. The module follows standard scRNA-seq analysis workflows while maintaining compatibility with popular formats and tools.

## Overview

The single-cell module is designed around a lightweight `SingleCellData` class that provides an AnnData-like interface for storing and manipulating single-cell data. It includes:

- **Data Structure**: `SingleCellData` class for managing expression matrices, cell/gene metadata
- **Preprocessing**: Loading, QC metrics, filtering, normalization, scaling
- **Dimensionality Reduction**: HVG selection, PCA, UMAP, t-SNE, diffusion maps
- **Clustering**: Leiden, Louvain, K-means, hierarchical clustering
- **Trajectory Analysis**: Pseudotime, lineage inference, gene trends
- **Visualization**: QC plots, dimensionality reduction plots, expression heatmaps
- **Integration**: Batch correction and dataset merging

## Architecture

```mermaid
flowchart TB
  subgraph Data[Data Management]
    A[SingleCellData]
    B[load_count_matrix]
  end
  
  subgraph Preprocessing
    C[calculate_qc_metrics]
    D[filter_cells]
    E[filter_genes]
    F[normalize_counts]
    G[log_transform]
    H[scale_data]
  end
  
  subgraph Dimensionality[Dimensionality Reduction]
    I[select_hvgs]
    J[compute_pca]
    K[compute_neighbors]
    L[compute_umap]
    M[compute_tsne]
    N[compute_diffusion_map]
  end
  
  subgraph Clustering
    O[leiden_clustering]
    P[louvain_clustering]
    Q[kmeans_clustering]
    R[hierarchical_clustering]
    S[find_markers]
  end
  
  subgraph Analysis[Advanced Analysis]
    T[compute_pseudotime]
    U[trajectory_analysis]
    V[gene_trends]
    W[batch_correction]
  end
  
  A --> C
  C --> D
  D --> E
  E --> F
  F --> G
  G --> H
  H --> I
  I --> J
  J --> K
  K --> L
  K --> M
  K --> N
  K --> O
  O --> S
  J --> T
  T --> U
  U --> V
  A --> W
```

## Key Features

### Data Structure
- **AnnData-compatible**: Similar interface to scanpy/AnnData
- **Sparse matrix support**: Efficient memory usage for large datasets
- **Flexible metadata**: Store cell and gene annotations
- **Extensible**: Easy to add custom analysis results

### Preprocessing Pipeline
- **Quality control**: Calculate standard QC metrics (total counts, n_genes, mitochondrial/ribosomal percentages)
- **Flexible filtering**: Filter cells and genes based on various criteria
- **Multiple normalization methods**: Total count, median, custom scaling factors
- **Log transformation**: With pseudocount handling
- **Data scaling**: Z-score standardization with optional centering

### Dimensionality Reduction
- **HVG selection**: Multiple methods (Seurat-style, variance-based, Cell Ranger)
- **PCA**: Principal component analysis with explained variance
- **UMAP**: Uniform Manifold Approximation and Projection (requires umap-learn)
- **t-SNE**: t-distributed Stochastic Neighbor Embedding
- **Diffusion maps**: For trajectory analysis
- **Neighbor graphs**: K-nearest neighbors for downstream analysis

### Clustering & Analysis
- **Graph-based clustering**: Leiden and Louvain algorithms
- **Centroid-based clustering**: K-means with multiple initialization strategies
- **Hierarchical clustering**: With various linkage methods
- **Marker gene identification**: Find genes that define clusters
- **Trajectory inference**: Pseudotime and lineage analysis

## Quick Start

```python
from metainformant.singlecell.preprocessing import load_count_matrix, calculate_qc_metrics
from metainformant.singlecell.dimensionality import select_hvgs, compute_pca, compute_umap
from metainformant.singlecell.clustering import leiden_clustering

# Load data
data = load_count_matrix("counts.csv", transpose=True)

# Quality control and filtering
data = calculate_qc_metrics(data)
data = filter_cells(data, min_genes=200, max_genes=5000, max_pct_mt=20)
data = filter_genes(data, min_cells=3)

# Normalization and scaling
data = normalize_counts(data, target_sum=1e4, method='total_count')
data = log_transform(data, base=2)
data = scale_data(data)

# Dimensionality reduction
data = select_hvgs(data, n_top_genes=2000)
data = compute_pca(data, n_components=50)
data = compute_umap(data, n_components=2)

# Clustering
data = leiden_clustering(data, resolution=0.5)
```

## File Format Support

- **CSV/TSV**: Standard comma or tab-separated files
- **Matrix Market**: Sparse matrix format (.mtx files)
- **HDF5**: Binary format for large datasets (requires h5py)
- **Compressed formats**: Automatic gzip handling

## Integration with Other Tools

The module is designed to complement existing single-cell analysis tools:
- **Input compatibility**: Can read outputs from Cell Ranger, STARsolo, alevin
- **Export capabilities**: Results can be converted to formats compatible with Seurat, scanpy
- **Visualization**: Integrates with METAINFORMANT's visualization module

## Module Structure

- **[preprocessing.py](./preprocessing.md)**: Data loading, QC, normalization, scaling
- **[dimensionality.py](./dimensionality.md)**: HVG selection, PCA, UMAP, t-SNE, diffusion maps
- **[clustering.py](./clustering.md)**: Various clustering algorithms and marker identification
- **[trajectory.py](./trajectory.md)**: Pseudotime and trajectory inference
- **[visualization.py](./visualization.md)**: Plotting functions for scRNA-seq data
- **[integration.py](./integration.md)**: Batch correction and dataset integration

## Dependencies

### Required
- `numpy`: Numerical computations
- `pandas`: Data manipulation
- `scipy`: Statistical functions and sparse matrices
- `scikit-learn`: Machine learning algorithms

### Optional
- `umap-learn`: For UMAP dimensionality reduction
- `h5py`: For HDF5 file support
- `matplotlib`: For visualization
- `seaborn`: For advanced plotting

## Performance Considerations

- **Memory efficiency**: Uses sparse matrices when appropriate
- **Scalability**: Tested with datasets up to 100k+ cells
- **Parallel processing**: Some operations utilize multiple cores
- **Caching**: Intermediate results can be cached for large analyses

## Best Practices

1. **Quality control first**: Always examine QC metrics before filtering
2. **Appropriate filtering**: Be conservative with filtering thresholds initially
3. **Normalization strategy**: Consider your downstream analysis when choosing normalization
4. **HVG selection**: Number of HVGs affects downstream clustering resolution
5. **Parameter tuning**: Clustering and UMAP parameters may need optimization per dataset

## Examples and Tutorials

For detailed examples and analysis workflows, see:
- `tests/test_singlecell_preprocessing.py`: Example preprocessing workflows
- `tests/test_singlecell_dimensionality.py`: Dimensionality reduction examples
- `output/`: Generated example outputs from test runs

## Related Documentation

- [Quality Control](../quality/index.md): For raw sequencing data QC
- [DNA Sequences](../dna/sequences.md): For sequence-based analysis
- [Visualization](../visualization/index.md): For general plotting utilities
- [Core I/O](../core/io.md): For file handling utilities
