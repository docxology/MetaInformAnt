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
from metainformant.singlecell import (
    load_count_matrix,
    calculate_qc_metrics,
    filter_cells,
    normalize_counts
)

# Load and preprocess data
data = load_count_matrix("counts.mtx", format="mtx")
data = calculate_qc_metrics(data)
data = filter_cells(data, min_genes=200, max_pct_mt=10.0)
data = normalize_counts(data)
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
from metainformant.singlecell import (
    select_hvgs,
    compute_pca,
    compute_umap,
    compute_neighbors
)

# Select highly variable genes
data = select_hvgs(data, n_top_genes=2000)

# Reduce dimensions
data = compute_pca(data, n_components=50)
data = compute_neighbors(data, n_neighbors=15)
data = compute_umap(data, min_dist=0.1)
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
from metainformant.singlecell import (
    leiden_clustering,
    louvain_clustering,
    find_marker_genes
)

# Cluster cells
data = leiden_clustering(data, resolution=1.0)
# Or use Louvain
data = louvain_clustering(data, resolution=1.0)

# Find marker genes
markers = find_marker_genes(data, cluster_key="cluster")
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
from metainformant.singlecell import (
    compute_pseudotime,
    trajectory_analysis,
    lineage_analysis
)

# Compute pseudotime
data = compute_pseudotime(data, root_cells=0, method="diffusion")

# Trajectory analysis
traj_result = trajectory_analysis(data)

# Lineage analysis
lineages = lineage_analysis(data)
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
from metainformant.singlecell import (
    integrate_datasets,
    harmony_integration,
    batch_correction
)

# Integrate multiple samples
integrated = integrate_datasets([sample1_data, sample2_data, sample3_data])

# Harmony batch correction
corrected = harmony_integration(integrated, batch_key="batch")

# Alternative: Combat batch correction
corrected = batch_correction(integrated, batch_key="batch", method="combat")
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
from metainformant.singlecell import (
    plot_qc_metrics,
    plot_dimensionality_reduction,
    plot_gene_expression,
    plot_clusters
)

# Create visualizations
fig = plot_qc_metrics(data)
fig = plot_dimensionality_reduction(data, method="umap", color_by="cluster")
fig = plot_gene_expression(data, gene="GENE1")
fig = plot_clusters(data, cluster_key="cluster")
```

## Integration with Other Modules

### With RNA Module
```python
from metainformant.rna import workflow
from metainformant.singlecell import preprocessing

# Process single-cell data
from metainformant.singlecell import normalize_counts, load_count_matrix

single_cell_data = load_count_matrix("counts.mtx")
single_cell_data = normalize_counts(single_cell_data)

# Compare with bulk RNA data from RNA workflow
# See RNA module for bulk expression data extraction
```

### With Machine Learning Module
```python
from metainformant.singlecell import leiden_clustering, compute_pca
from metainformant.ml import BiologicalClassifier

# Use clustering results for supervised learning
data = leiden_clustering(data, resolution=1.0)
data = compute_pca(data, n_components=20)

# Extract cluster labels for classification
clusters = data.metadata.get("cluster", [])
# Train classifier with cluster labels
classifier = BiologicalClassifier(algorithm="random_forest")
# classifier.fit(data.counts, clusters)
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
