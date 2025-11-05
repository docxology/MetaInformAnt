# Visualization Documentation

This directory contains comprehensive documentation for METAINFORMANT's plotting, animation, and visualization capabilities.

## Overview

The visualization domain provides unified plotting and animation utilities for creating publication-quality figures and interactive visualizations across all METAINFORMANT domains.

## Documentation Files

### Core Documentation
- **`index.md`**: Visualization domain overview and module index
- **`README.md`**: This file - overview and navigation

### Plotting Modules
- **`basic.md`**: Basic plotting functions (line, scatter, bar, pie, area, heatmap)
- **`statistical.md`**: Statistical visualization (histograms, box plots, Q-Q plots, ROC curves)
- **`genomics.md`**: Genomic visualization (Manhattan plots, volcano plots, regional plots)
- **`expression.md`**: Expression analysis plots (heatmaps, enrichment, differential expression)
- **`dimred.md`**: Dimensionality reduction (PCA, UMAP, t-SNE)
- **`networks.md`**: Network visualization (graphs, layouts, communities)
- **`timeseries.md`**: Time series visualization (plots, autocorrelation, decomposition)
- **`multidim.md`**: Multi-dimensional plots (pair plots, parallel coordinates, radar charts)
- **`quality.md`**: Quality control plots (QC metrics, quality scores, adapter content)
- **`information.md`**: Information theory plots (entropy, mutual information, profiles)

### Specialized Modules
- **`trees.md`**: Phylogenetic tree visualization
- **`animations.md`**: Time-series and dynamic data animation

### Utilities and Integration
- **`styling.md`**: Styling and customization guide
- **`integration.md`**: Domain integration patterns
- **`examples.md`**: Comprehensive usage examples
- **`gallery.md`**: Visualization gallery with code

## Related Source Code

- See `src/metainformant/visualization/` for implementation details
- See `tests/test_visualization_*.py` for comprehensive test coverage
- See `src/metainformant/visualization/README.md` for module-specific documentation

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

# Genomic plots
ax = manhattan_plot(data, 'position', 'neg_log10_p', 'chromosome')

# Dimensionality reduction
ax = pca_plot(data, pc_x=1, pc_y=2, hue='group')

# Trees
from Bio import Phylo
tree = Phylo.read("tree.nwk", "newick")
ax = plot_phylo_tree(tree)

# Animations
fig, anim = animate_time_series([1, 2, 3, 2, 4, 5])
```

## Integration

Visualization tools integrate with:
- **DNA/RNA analysis** for sequence and expression plots
- **GWAS** for genome-wide association visualization
- **Single-cell analysis** for clustering and trajectory plots
- **Information theory** for entropy and mutual information visualization
- **Life events** for event sequence visualization
- **Mathematical models** for theoretical visualization
- **Statistical methods** for data exploration

## Testing

Comprehensive tests ensure visualization reliability:
- Plot rendering correctness
- Data accuracy in visualizations
- Animation functionality
- Export format validation
- Performance benchmarking

## Contributing

When adding new visualization functionality:
1. Update plotting and rendering documentation
2. Add comprehensive visual tests
3. Ensure compatibility with multiple backends
4. Update integration examples

This documentation provides complete coverage of METAINFORMANT's visualization capabilities.
