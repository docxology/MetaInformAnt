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
from metainformant.visualization import plots

# Line plot
data = [1, 4, 2, 8, 5, 7]
ax = plots.lineplot(data, xlabel="Time", ylabel="Value", title="Time Series")

# Heatmap
import numpy as np
matrix = np.random.random((10, 10))
ax = plots.heatmap(matrix, xlabel="X", ylabel="Y", title="Correlation Matrix")

# Histogram
ax = plots.histogram(data, bins=20, xlabel="Value", ylabel="Frequency")
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
from metainformant.visualization import trees

# Load and visualize tree
tree = trees.load_newick("tree.nwk")
ax = trees.plot_tree(tree, show_branch_lengths=True)

# Compare multiple trees
trees.compare_trees([tree1, tree2, tree3])
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
from metainformant.visualization import animations

# Animate time series data
data = [[1, 2], [2, 3], [3, 2], [4, 4]]
fig, anim = animations.animate_time_series(data, interval_ms=500)

# Save animation
animations.save_animation(anim, "timeseries.gif", fps=2)
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
from metainformant.visualization import plots

# Generate sample data
x = np.linspace(0, 10, 100)
y = np.sin(x) + np.random.normal(0, 0.1, 100)

# Create publication-quality plot
ax = plots.scatterplot(x, y, alpha=0.6, s=20)
ax = plots.lineplot(x, np.sin(x), ax=ax, color='red', linewidth=2)
ax.set_xlabel("X Values")
ax.set_ylabel("Y Values")
ax.set_title("Sample Data with Trend Line")
ax.grid(True, alpha=0.3)
plots.save_figure(ax.figure, "sample_plot.png", dpi=300)
```

### Phylogenetic Tree Visualization
```python
from metainformant.visualization import trees

# Load tree from file
tree = trees.load_newick("phylogenetic_tree.nwk")

# Create visualization with annotations
ax = trees.plot_tree(
    tree,
    show_branch_lengths=True,
    leaf_labels=True,
    internal_labels=False,
    figsize=(12, 8)
)

# Customize appearance
ax.set_title("Phylogenetic Tree")
trees.style_tree(ax, color_scheme="viridis")

# Save high-resolution figure
trees.save_tree_plot(ax.figure, "tree_visualization.png", dpi=300)
```

### Animation Creation
```python
from metainformant.visualization import animations
import numpy as np

# Generate time series data
time_points = np.linspace(0, 4*np.pi, 50)
data = [np.sin(time_points + i*0.5) for i in range(20)]

# Create animation
fig, anim = animations.animate_multiple_series(
    data,
    interval_ms=200,
    colors=animations.get_color_palette("plasma", 20)
)

# Customize animation
ax = fig.gca()
ax.set_xlabel("Time")
ax.set_ylabel("Value")
ax.set_title("Dynamic Time Series")
ax.grid(True, alpha=0.3)

# Export animation
animations.save_animation(anim, "dynamic_series.gif", fps=5)
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
from metainformant.visualization import trees

# Visualize DNA-based phylogenetic analysis
tree = phylogeny.neighbor_joining_tree(dna_sequences)
trees.plot_tree(tree, title="DNA-based Phylogeny")
```

### With Math Module
```python
from metainformant.math import price_equation
from metainformant.visualization import lineplot

# Visualize mathematical models
results = price_equation.simulate_evolution(...)
lineplot(results["trait_values"], title="Trait Evolution")
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
