# Visualization: Overview

The visualization module provides comprehensive plotting and animation utilities for creating publication-quality figures across all METAINFORMANT domains.

## Module Organization

The visualization package is organized into category-specific modules:

### Core Plotting Modules
- **[Basic Plots](./basic.md)**: Line, scatter, bar, pie, area, step plots, and heatmaps
- **[Statistical Plots](./statistical.md)**: Histograms, box plots, violin plots, Q-Q plots, density plots, ROC curves
- **[Genomics Plots](./genomics.md)**: Manhattan plots, volcano plots, regional plots, chromosome ideograms
- **[Expression Plots](./expression.md)**: Expression heatmaps, enrichment plots, differential expression
- **[Dimensionality Reduction](./dimred.md)**: PCA, UMAP, t-SNE plots with loadings and scree plots
- **[Network Plots](./networks.md)**: Network graphs, circular layouts, hierarchical layouts
- **[Time Series](./timeseries.md)**: Time series plots, autocorrelation, seasonal decomposition, forecasts
- **[Multi-dimensional](./multidim.md)**: Pair plots, parallel coordinates, radar charts, 3D scatter
- **[Quality Control](./quality.md)**: QC metrics, quality scores, adapter content, sequence length
- **[Information Theory](./information.md)**: Entropy plots, mutual information, information profiles

### Specialized Modules
- **[Trees](./trees.md)**: Phylogenetic tree visualization with multiple layouts
- **[Animations](./animations.md)**: Time series, evolution, clustering, network, and trajectory animations

### Utility Modules
- **[Styling](./styling.md)**: Publication-quality styles, color palettes, font management
- **[Layout](./integration.md)**: Multi-panel figure creation and layout management
- **[Export](./integration.md)**: High-resolution export in multiple formats
- **[Interactive](./integration.md)**: Plotly integration for interactive plots

### Domain Integration
- **[GWAS Integration](./integration.md)**: Unified access to GWAS visualization functions
- **[Single-Cell Integration](./integration.md)**: Single-cell RNA-seq visualization functions
- **[Information Theory Integration](./integration.md)**: Information theory visualization functions
- **[Life Events Integration](./integration.md)**: Life course event visualization functions

## Quick Start

```python
from metainformant.visualization import (
    lineplot, scatter_plot, histogram, heatmap,
    manhattan_plot, volcano_plot, pca_plot,
    plot_phylo_tree, animate_time_series
)

# Basic plots
ax = lineplot(None, [1, 4, 2, 8, 5])
ax = scatter_plot([1, 2, 3], [4, 5, 6])
ax = histogram([1, 2, 3, 4, 5], bins=10)

# Genomic plots
ax = manhattan_plot(data, 'position', 'neg_log10_p', 'chromosome')
ax = volcano_plot(data, 'log2fc', 'neg_log10_p')

# Dimensionality reduction
ax = pca_plot(data, pc_x=1, pc_y=2, hue='group')

# Trees
from Bio import Phylo
tree = Phylo.read("tree.nwk", "newick")
ax = plot_phylo_tree(tree)

# Animations
fig, anim = animate_time_series([1, 2, 3, 2, 4, 5])
```

## Documentation

- **[Basic Plots](./basic.md)**: Simple plotting functions
- **[Statistical Plots](./statistical.md)**: Statistical visualization
- **[Genomics Plots](./genomics.md)**: Genomic visualization
- **[Expression Plots](./expression.md)**: Expression analysis
- **[Dimensionality Reduction](./dimred.md)**: PCA, UMAP, t-SNE
- **[Network Plots](./networks.md)**: Network visualization
- **[Time Series](./timeseries.md)**: Time series analysis
- **[Multi-dimensional](./multidim.md)**: Multi-dimensional plots
- **[Quality Control](./quality.md)**: Quality control plots
- **[Information Theory](./information.md)**: Information theory plots
- **[Trees](./trees.md)**: Phylogenetic trees
- **[Animations](./animations.md)**: Dynamic animations
- **[Styling](./styling.md)**: Styling and customization
- **[Integration](./integration.md)**: Domain integration patterns
- **[Examples](./examples.md)**: Comprehensive usage examples
- **[Gallery](./gallery.md)**: Visualization gallery
