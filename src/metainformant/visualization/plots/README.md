# Plotting and Visualizations Module

The `metainformant.visualization.plots` module contains the fundamental and specialized plotting engines.

## Overview

This module provides high-level APIs for common plots (lines, scatters, heatmaps) as well as complex biological visualizations (Venn diagrams, Sankey diagrams, animations).

## Key Components

### Basic Plots
General purpose functions for standard data visualization.
- `lineplot`, `scatter_plot`, `heatmap`, `bar_plot`.

### Specialized Plots
Complex visualizations tailored for biological data.
- `plot_venn_diagram`, `plot_sankey_diagram`, `plot_chord_diagram`.

### Animations
Dynamic visualizations of time-series or evolutionary data.
- `animate_time_series`, `animate_evolution`, `animate_network`.

## Usage Example

```python
from metainformant.visualization.plots.animations import animate_evolution
import numpy as np

# Simulate evolution data
generations = [np.random.normal(0, 1, 100) for _ in range(50)]

# Create animation
fig, anim = animate_evolution(generations)
anim.save("evolution.mp4")
```
