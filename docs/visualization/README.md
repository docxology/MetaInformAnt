# Visualization Documentation

This directory contains comprehensive documentation for METAINFORMANT's plotting, animation, and visualization capabilities.

## Overview

The visualization domain provides unified plotting and animation utilities for creating publication-quality figures and interactive visualizations across all METAINFORMANT domains.

## Documentation Files

### Core Visualization Tools
- **`index.md`**: Visualization domain overview and module index
- **`plots.md`**: Statistical plotting and data visualization
- **`trees.md`**: Phylogenetic tree visualization and rendering
- **`animations.md`**: Time-series and dynamic data animation

## Related Source Code

- See `src/metainformant/visualization/` for implementation details
- See `tests/test_visualization_*.py` for comprehensive test coverage
- See `src/metainformant/visualization/README.md` for module-specific documentation

## Usage Examples

The visualization domain supports diverse plotting needs:

```python
from metainformant.visualization import plots, trees

# Statistical plots
ax = plots.lineplot([1, 4, 2, 8, 5], label="data")
ax.set_xlabel("Time")
ax.set_ylabel("Value")

# Phylogenetic trees
tree = # ... load or create tree
ax = trees.plot_phylo_tree(tree)

# Animations
from metainformant.visualization import animations
fig, anim = animations.animate_time_series([0, 1, 3, 2, 5])
```

## Integration

Visualization tools integrate with:
- **DNA/RNA analysis** for sequence and expression plots
- **Mathematical models** for theoretical visualization
- **Single-cell analysis** for clustering and trajectory plots
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
