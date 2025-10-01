# Single-Cell Genomics Documentation

This directory contains comprehensive documentation for METAINFORMANT's single-cell RNA sequencing analysis capabilities.

## Overview

The single-cell domain provides tools for preprocessing, dimensionality reduction, clustering, trajectory analysis, and visualization of single-cell transcriptomic data.

## Documentation Files

### Core Single-Cell Analysis
- **`index.md`**: Single-cell domain overview and module index
- **`preprocessing.md`**: Data loading, QC, normalization, and scaling
- **`dimensionality.md`**: HVG selection, PCA, UMAP, t-SNE, diffusion maps
- **`clustering.md`**: Various clustering algorithms and marker identification
- **`trajectory.md`**: Pseudotime and trajectory inference
- **`visualization.md`**: Plotting functions for scRNA-seq data
- **`integration.md`**: Batch correction and dataset integration

## Related Source Code

- See `src/metainformant/singlecell/` for implementation details
- See `tests/test_singlecell_*.py` for comprehensive test coverage
- See `src/metainformant/singlecell/README.md` for module-specific documentation

## Usage Examples

The single-cell domain supports complete scRNA-seq workflows:

```python
from metainformant.singlecell.preprocessing import load_count_matrix, calculate_qc_metrics
from metainformant.singlecell.dimensionality import select_hvgs, compute_pca, compute_umap
from metainformant.singlecell.clustering import leiden_clustering

# Complete single-cell analysis pipeline
data = load_count_matrix("counts.csv", transpose=True)
data = calculate_qc_metrics(data)
data = select_hvgs(data, n_top_genes=2000)
data = compute_pca(data, n_components=50)
data = compute_umap(data, n_components=2)
data = leiden_clustering(data, resolution=0.5)
```

## Integration

Single-cell analysis integrates with:
- **RNA workflows** for transcriptomic context
- **Statistical methods** for differential expression
- **Visualization tools** for publication-quality figures
- **Machine learning** for cell type classification

## Testing

Comprehensive tests ensure analysis reliability:
- Data format compatibility testing
- Algorithm implementation verification
- Performance benchmarking with large datasets
- Integration testing across preprocessing steps

## Contributing

When adding new single-cell functionality:
1. Update preprocessing and analysis documentation
2. Add comprehensive algorithm tests
3. Ensure compatibility with scanpy/Seurat workflows
4. Update visualization and integration guides

This documentation provides complete coverage of METAINFORMANT's single-cell genomics capabilities.
