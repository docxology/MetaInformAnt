# Basic Plots

Basic plotting functions for simple visualizations including line plots, scatter plots, bar charts, pie charts, area plots, and heatmaps.

## Functions

### `lineplot(x, y, *, label=None, ax=None, style='-', color=None)`

Create a simple line plot.

**Parameters:**
- `x`: X-axis values (if None, uses indices)
- `y`: Y-axis values
- `label`: Label for legend
- `ax`: Matplotlib axes (creates new if None)
- `style`: Line style (e.g., '-', '--', '-.')
- `color`: Line color

**Example:**
```python
from metainformant.visualization import lineplot

ax = lineplot(None, [1, 4, 2, 8, 5], label="Data")
ax.set_xlabel("Index")
ax.set_ylabel("Value")
```

### `scatter_plot(x, y, *, ax=None, color=None, size=20, alpha=0.7, xlabel=None, ylabel=None, title=None)`

Create a scatter plot.

**Example:**
```python
from metainformant.visualization import scatter_plot

ax = scatter_plot([1, 2, 3], [4, 5, 6], xlabel="X", ylabel="Y", title="Scatter")
```

### `bar_plot(x, y, *, ax=None, color=None, alpha=0.7, horizontal=False, xlabel=None, ylabel=None, title=None)`

Create a bar plot.

**Example:**
```python
from metainformant.visualization import bar_plot

ax = bar_plot(["A", "B", "C"], [10, 20, 15])
```

### `pie_chart(sizes, labels=None, *, ax=None, colors=None, autopct='%1.1f%%', title=None)`

Create a pie chart.

**Example:**
```python
from metainformant.visualization import pie_chart

ax = pie_chart([30, 25, 45], ["A", "B", "C"])
```

### `area_plot(x, y, *, ax=None, color=None, alpha=0.5, xlabel=None, ylabel=None, title=None)`

Create an area plot (filled line plot).

**Example:**
```python
from metainformant.visualization import area_plot

ax = area_plot([1, 2, 3, 4], [1, 4, 2, 3])
```

### `step_plot(x, y, *, ax=None, where='pre', label=None, color=None, xlabel=None, ylabel=None, title=None)`

Create a step plot.

**Example:**
```python
from metainformant.visualization import step_plot

ax = step_plot([1, 2, 3, 4], [1, 4, 2, 3])
```

### `heatmap(data, *, cmap='viridis', cbar=True, ax=None, annot=False)`

Create a heatmap from 2D data or DataFrame.

**Example:**
```python
from metainformant.visualization import heatmap
import numpy as np

data = np.random.random((10, 10))
ax = heatmap(data, annot=True)
```

