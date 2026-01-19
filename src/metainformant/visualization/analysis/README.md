# Statistical Analysis Visualization Module

The `metainformant.visualization.analysis` module provides tools for visualizing statistical and quality metrics.

## Overview

This module covers common statistical plots like histograms, box plots, and correlation heatmaps, as well as specialized visualizations for dimensional reduction (PCA, UMAP) and data quality.

## Key Functions

- **`histogram`**: Visualizing distribution of single variables.
- **`box_plot`**: Comparing distributions across categories.
- **`correlation_heatmap`**: Visualizing relationships between multiple variables.
- **`plot_pca`**: 2D/3D visualization of dimensionality reduction results.
- **`plot_quality_metrics`**: Summarizing FASTQ or assembly quality results.

## Usage Example

```python
from metainformant.visualization.analysis import plot_pca
import numpy as np

# Sample data
data = np.random.rand(100, 2)
labels = [0] * 50 + [1] * 50

# Plot PCA
ax = plot_pca(data, labels=labels)
ax.set_title("PCA Result")
```
