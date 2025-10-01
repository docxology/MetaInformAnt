# Single Cell Analysis Module

The `singlecell` module provides comprehensive tools for single-cell transcriptomic analysis, including preprocessing, dimensionality reduction, clustering, and trajectory inference.

## Overview

This module handles the complete single-cell RNA sequencing analysis pipeline from raw count matrices to biological interpretation.

## Key Components

### Preprocessing (`preprocessing.py`)
Quality control, normalization, and filtering of single-cell data.

**Key Features:**
- Quality metric calculation and filtering
- Normalization methods (CPM, TPM, etc.)
- Batch effect detection and correction
- Mitochondrial gene percentage filtering

**Usage:**
```python
from metainformant.singlecell import preprocessing

# Load and preprocess data
counts = preprocessing.load_counts("counts.h5ad")
filtered = preprocessing.filter_cells(counts, min_genes=200, max_mito=0.1)
normalized = preprocessing.normalize_counts(filtered, method="cpm")
```

### Dimensionality Reduction (`dimensionality.py`)
PCA, t-SNE, UMAP, and other dimensionality reduction techniques.

**Key Features:**
- Principal component analysis (PCA)
- t-Distributed Stochastic Neighbor Embedding (t-SNE)
- Uniform Manifold Approximation and Projection (UMAP)
- Non-negative matrix factorization (NMF)

**Usage:**
```python
from metainformant.singlecell import dimensionality

# Reduce dimensions
pca_result = dimensionality.run_pca(normalized, n_components=50)
umap_result = dimensionality.run_umap(pca_result, n_neighbors=15, min_dist=0.1)
```

### Clustering (`clustering.py`)
Cell type identification and cluster analysis.

**Key Features:**
- K-means and hierarchical clustering
- Graph-based clustering (Louvain, Leiden)
- Cluster marker identification
- Cluster stability assessment

**Usage:**
```python
from metainformant.singlecell import clustering

# Cluster cells
clusters = clustering.louvain_clustering(umap_result, resolution=1.0)
markers = clustering.find_markers(normalized, clusters)
```

### Trajectory Analysis (`trajectory.py`)
Pseudotime and developmental trajectory inference.

**Key Features:**
- Pseudotime ordering algorithms
- Trajectory branching detection
- Developmental stage identification
- Trajectory comparison across conditions

**Usage:**
```python
from metainformant.singlecell import trajectory

# Infer trajectory
trajectory_result = trajectory.infer_pseudotime(umap_result, start_cell="root")
branching = trajectory.detect_branching(trajectory_result)
```

### Integration (`integration.py`)
Multi-sample batch correction and integration.

**Key Features:**
- Canonical correlation analysis (CCA)
- Mutual nearest neighbors (MNN) correction
- Harmony batch correction
- Data integration quality assessment

**Usage:**
```python
from metainformant.singlecell import integration

# Integrate multiple samples
integrated = integration.harmony_integration([sample1, sample2, sample3])
corrected = integration.mnn_correction(integrated)
```

### Visualization (`visualization.py`)
Specialized single-cell data visualization.

**Key Features:**
- UMAP and t-SNE scatter plots
- Feature expression overlays
- Trajectory visualizations
- Cluster annotation plots

**Usage:**
```python
from metainformant.singlecell import visualization

# Create visualizations
umap_plot = visualization.plot_umap(umap_result, color_by=clusters)
trajectory_plot = visualization.plot_trajectory(trajectory_result)
```

## Integration with Other Modules

### With RNA Module
```python
from metainformant.rna import workflow
from metainformant.singlecell import preprocessing

# Process bulk RNA data for single-cell comparison
bulk_expression = workflow.extract_expression_patterns(rna_data)
single_cell = preprocessing.normalize_counts(single_cell_counts)
comparison = singlecell.compare_with_bulk(single_cell, bulk_expression)
```

### With Machine Learning Module
```python
from metainformant.singlecell import dimensionality, clustering
from metainformant.ml import classification

# Use clustering results for supervised learning
features = dimensionality.run_pca(normalized, n_components=20)
clusters = clustering.louvain_clustering(features)
classifier = classification.train_classifier(features, clusters)
```

## Performance Features

- Memory-efficient processing of large single-cell datasets
- Parallel computation for time-intensive operations
- Streaming processing for very large datasets
- GPU acceleration support where applicable

## Testing

Comprehensive tests cover:
- Algorithm implementation correctness
- Integration with real single-cell datasets
- Performance benchmarking
- Edge case handling

## Dependencies

- Scanpy for core single-cell algorithms
- Optional: scVelo for trajectory analysis, scikit-learn for additional ML methods

This module provides a complete toolkit for single-cell transcriptomic analysis.
