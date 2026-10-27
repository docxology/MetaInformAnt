# Visualization: Overview

The visualization module provides comprehensive plotting and animation utilities for creating publication-quality figures across all METAINFORMANT domains.

## Module Organization

The visualization package is organized into category-specific subpackages. The package namespace `metainformant.visualization` itself exposes only these subpackages (plus the `WONG` and `OKABE_ITO` color lists), so import plot functions from their modules:

- `metainformant.visualization.plots`: `basic`, `general`, `multidim`, `specialized` and `animations`
- `metainformant.visualization.analysis`: `statistical`, `dimred`, `timeseries`, `information` and the quality-control modules
- `metainformant.visualization.genomics`: `genomics`, `expression`, `networks` and `trees`
- `metainformant.visualization.config`: themes, palettes and figure conventions
- `metainformant.visualization.dashboards` and `metainformant.visualization.interactive_dashboards`: composite figures and interactive plots

### Core Plotting Modules
- **[Basic Plots](./basic.md)**: Line, scatter, bar, pie, area, step plots, and heatmaps
- **[Statistical Plots](./statistical.md)**: Histograms, box plots, violin plots, Q-Q plots, density plots, ROC curves
- **[Genomics Plots](./genomics.md)**: Manhattan plots, volcano plots, regional plots, circular Manhattan plots, coverage and variant plots
- **[Expression Plots](./expression.md)**: Expression heatmaps, enrichment plots, differential expression
- **[Dimensionality Reduction](./dimred.md)**: PCA, UMAP, t-SNE plots with loadings and biplots
- **[Network Plots](./networks.md)**: Network graphs, circular layouts, hierarchical layouts
- **[Time Series](./timeseries.md)**: Time series plots, autocorrelation, seasonal decomposition, forecasts
- **[Multi-dimensional](./multidim.md)**: Pair plots, parallel coordinates, radar charts, 3D scatter
- **[Quality Control](./quality.md)**: QC metrics, quality scores, adapter content, read length distributions
- **[Information Theory](./information.md)**: Entropy profiles, mutual information, Rényi spectra, information landscapes and networks

### Specialized Modules
- **[Trees](./trees.md)**: Phylogenetic tree visualization with multiple layouts
- **[Animations](./animations.md)**: Time series, evolution, clustering, network, and trajectory animations

### Utility Modules
- **[Styling](./styling.md)**: Publication-quality styles, color palettes, font management
- **[Layout](./integration.md)**: Multi-panel figure creation and layout management
- **[Export](./integration.md)**: Deterministic high-resolution figure export
- **[Interactive](./integration.md)**: Interactive plot data structures with standalone HTML export, and Plotly wrappers

### Domain Integration
- **[GWAS Integration](./integration.md)**: Where the GWAS visualization functions live
- **[Single-Cell Integration](./integration.md)**: Single-cell RNA-seq visualization functions
- **[Information Theory Integration](./integration.md)**: Information theory visualization functions
- **[Life Events Integration](./integration.md)**: Life course event visualization functions

## Quick Start

```python
import io

import numpy as np
import pandas as pd
from Bio import Phylo

from metainformant.visualization.analysis.statistical import histogram
from metainformant.visualization.genomics.trees import plot_phylo_tree
from metainformant.visualization.plots import manhattan_plot, pca_plot, scatter_plot, volcano_plot
from metainformant.visualization.plots.animations import animate_time_series
from metainformant.visualization.plots.basic import lineplot

rng = np.random.default_rng(0)

# Basic plots (numeric inputs are NumPy arrays)
ax = lineplot(np.array([1, 4, 2, 8, 5]))
ax = scatter_plot(np.array([1, 2, 3]), np.array([4, 5, 6]))
ax = histogram(rng.normal(size=200), bins=10)

# Genomic plots (p-value columns hold raw p-values; -log10 is applied for you)
gwas = pd.DataFrame({
    'chromosome': np.repeat(['chr1', 'chr2'], 50),
    'position': np.tile(np.arange(50) * 1000, 2),
    'pvalue': rng.uniform(1e-6, 1, 100),
})
ax = manhattan_plot(gwas, 'position', 'pvalue', chr_col='chromosome')

de = pd.DataFrame({'log2fc': rng.normal(0, 1.5, 100), 'pvalue': rng.uniform(1e-6, 1, 100)})
ax = volcano_plot(de, 'log2fc', 'pvalue')

# Dimensionality reduction (pre-computed components and their variance percentages)
pcs = pd.DataFrame({
    'PC1': rng.normal(size=60),
    'PC2': rng.normal(size=60),
    'PC1_variance': 41.0,
    'PC2_variance': 23.5,
    'group': np.repeat(['A', 'B'], 30),
})
ax = pca_plot(pcs, hue='group')

# Trees (a Biopython tree; here parsed from a Newick string)
tree = Phylo.read(io.StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
ax = plot_phylo_tree(tree)

# Animations (returns the figure and the animation; output_path saves a GIF)
fig, anim = animate_time_series([1, 2, 3, 2, 4, 5], output_path="output/animation.gif")
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
