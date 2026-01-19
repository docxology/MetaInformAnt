# Genomics Visualization Module

The `metainformant.visualization.genomics` module provides specialized plotting functions for genomic data.

## Overview

This module includes visualizations for GWAS results, phylogenetic trees, gene expression heatmaps, and chromosome ideograms.

## Key Functions

- **`manhattan_plot`**: Standard visualization for association studies.
- **`volcano_plot`**: Visualizing differential expression or association results.
- **`plot_phylo_tree`**: Drawing phylogenetic trees from various data structures.
- **`chromosome_ideogram`**: Visualizing features across the genome.

## Usage Example

```python
from metainformant.visualization.genomics import manhattan_plot
import pandas as pd

# Load results
results = pd.read_csv("gwas_results.tsv", sep="\t")

# Plot
ax = manhattan_plot(results)
ax.figure.savefig("manhattan.png")
```
