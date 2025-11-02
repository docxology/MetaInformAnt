# Visualization Module

The `visualization` module provides unified plotting and animation utilities for creating publication-quality figures and interactive visualizations across all METAINFORMANT domains.

## Overview

This module offers a cohesive API for creating various types of plots, from statistical charts to phylogenetic trees and animations. It integrates multiple plotting backends while providing a consistent interface.

## Submodules

### Statistical Plots (`plots.py`)
Core plotting functionality for statistical data visualization.

**Key Features:**
- Line plots and scatter plots
- Heatmaps and correlation matrices
- Histograms and density plots
- Box plots and violin plots
- Publication-quality formatting

**Usage:**
```python
from metainformant.visualization import lineplot, heatmap
from metainformant.visualization.plots import scatter_plot, histogram

# Line plot
data = [1, 4, 2, 8, 5, 7]
ax = lineplot(None, data)
ax.set_xlabel("Time")
ax.set_ylabel("Value")
ax.set_title("Time Series")

# Heatmap
import numpy as np
matrix = np.random.random((10, 10))
ax = heatmap(matrix)

# Histogram
from metainformant.visualization.plots import histogram
ax = histogram(data, bins=20)
```

### Phylogenetic Trees (`trees.py`)
Specialized visualization for phylogenetic trees and evolutionary relationships.

**Key Features:**
- Newick tree parsing and rendering
- Tree annotation and labeling
- Branch length visualization
- Multiple tree comparison
- Interactive tree exploration

**Usage:**
```python
from metainformant.visualization import plot_phylo_tree
from Bio import Phylo

# Load tree from Newick file (using BioPython)
tree = Phylo.read("tree.nwk", "newick")

# Visualize tree
ax = plot_phylo_tree(tree)

# Save figure using matplotlib
ax.figure.savefig("tree.png", dpi=300)
```

### Animations (`animations.py`)
Time-series and dynamic data animation.

**Key Features:**
- Time-series animation
- Dynamic plot updates
- GIF and video export
- Real-time data visualization
- Animation controls and playback

**Usage:**
```python
from metainformant.visualization import animate_time_series

# Animate time series data
data = [1, 2, 3, 2, 4, 5, 3, 6]
fig, anim = animate_time_series(data, interval_ms=500)

# Save animation using matplotlib's animation writer
from matplotlib.animation import PillowWriter
writer = PillowWriter(fps=2)
anim.save("timeseries.gif", writer=writer)
```

## Integration Features

### Backend Support
- **Matplotlib**: Primary backend for static plots
- **Seaborn**: Enhanced statistical visualizations
- **Plotly**: Interactive web-based visualizations (optional)
- **PyQtGraph**: Real-time plotting (optional)

### Consistent API
All plotting functions follow a consistent pattern:
```python
ax = plot_function(data, **kwargs)
ax.set_xlabel("X Label")
ax.set_ylabel("Y Label")
ax.set_title("Plot Title")
```

## Advanced Features

### Publication Quality
- High-resolution export (300+ DPI)
- Consistent color schemes and styling
- LaTeX text rendering support
- Vector format export (SVG, PDF, EPS)

### Customization
- Extensive styling options
- Custom color palettes
- Font and layout control
- Legend and annotation management

## Usage Examples

### Basic Statistical Visualization
```python
import numpy as np
from metainformant.visualization import lineplot
from metainformant.visualization.plots import scatter_plot

# Generate sample data
x = np.linspace(0, 10, 100)
y = np.sin(x) + np.random.normal(0, 0.1, 100)

# Create publication-quality plot
ax = scatter_plot(x, y, alpha=0.6, s=20)
ax = lineplot(x, np.sin(x), ax=ax, color='red', style='-')
ax.set_xlabel("X Values")
ax.set_ylabel("Y Values")
ax.set_title("Sample Data with Trend Line")
ax.grid(True, alpha=0.3)
ax.figure.savefig("sample_plot.png", dpi=300)
```

### Phylogenetic Tree Visualization
```python
from metainformant.visualization import plot_phylo_tree
from Bio import Phylo
import matplotlib.pyplot as plt

# Load tree from file
tree = Phylo.read("phylogenetic_tree.nwk", "newick")

# Create visualization
fig, ax = plt.subplots(figsize=(12, 8))
ax = plot_phylo_tree(tree, ax=ax)
ax.set_title("Phylogenetic Tree")

# Save high-resolution figure
fig.savefig("tree_visualization.png", dpi=300)
plt.close(fig)
```

### Animation Creation
```python
from metainformant.visualization import animate_time_series
import numpy as np
from matplotlib.animation import PillowWriter

# Generate time series data
time_points = np.linspace(0, 4*np.pi, 50)
data = np.sin(time_points)

# Create animation
fig, anim = animate_time_series(data, interval_ms=200)

# Customize animation
ax = fig.gca()
ax.set_xlabel("Time")
ax.set_ylabel("Value")
ax.set_title("Dynamic Time Series")
ax.grid(True, alpha=0.3)

# Export animation
writer = PillowWriter(fps=5)
anim.save("dynamic_series.gif", writer=writer)
```

## Performance Considerations

- **Memory Management**: Figures are properly closed after saving
- **Vector Graphics**: SVG/PDF export for scalable graphics
- **Batch Processing**: Efficient handling of multiple plots
- **Large Dataset Optimization**: Downsampling for large datasets

## Integration with Other Modules

### With DNA Module
```python
from metainformant.dna import phylogeny
from metainformant.visualization import plot_phylo_tree

# Visualize DNA-based phylogenetic analysis
tree = phylogeny.neighbor_joining_tree(dna_sequences)
ax = plot_phylo_tree(tree)
ax.set_title("DNA-based Phylogeny")
```

### With Math Module
```python
from metainformant.math import price_equation
from metainformant.visualization import lineplot

# Visualize mathematical models
fitness = [0.2, 0.4, 0.1]
trait_parent = [1.0, 1.2, 0.9]
trait_offspring = [1.25, 1.35, 0.95]
cov, trans, total = price_equation(fitness, trait_parent, trait_offspring)

# Plot trait evolution
ax = lineplot(None, trait_offspring)
ax.set_title("Trait Evolution")
```

## Testing

Comprehensive tests ensure:
- Plot rendering correctness
- Data accuracy in visualizations
- Animation functionality
- Export format validation
- Performance benchmarks

## Dependencies

- **Required**: matplotlib, numpy
- **Optional**: seaborn (enhanced styling), plotly (interactive plots)
- **Animation**: pillow (GIF support), ffmpeg (video export)

This module provides a complete visualization toolkit for biological data analysis and publication-quality figure generation.
