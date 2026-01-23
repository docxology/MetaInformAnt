# Agent Directives: docs/visualization

## Role
Documentation for METAINFORMANT's comprehensive visualization system.

## Contents
- `index.md` - Visualization module overview
- `gallery.md` - Plot gallery with examples
- Domain-specific visualization docs:
  - `genomics.md` - Manhattan, QQ, regional plots
  - `statistical.md` - Statistical analysis plots
  - `networks.md` - Network visualization
  - `expression.md` - Expression heatmaps
  - `trees.md` - Phylogenetic trees
  - `dimred.md` - Dimensionality reduction plots
  - `timeseries.md` - Time series plots
  - `animations.md` - Animated visualizations

## Key Patterns
- All functions return matplotlib Figure objects (never None)
- Output to `output/visualization/` by default
- Consistent styling across plot types
- Support for both interactive and publication-quality figures
