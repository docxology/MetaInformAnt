# Multi-dimensional Plots

Multi-dimensional visualization functions including pair plots, parallel coordinates, radar charts, scatter plot matrices, and 3D plots.

## Functions

### `pairplot_dataframe(df, *, hue=None)`

Pairplot for a tidy DataFrame, returns seaborn PairGrid or matplotlib figure.

**Example:**
```python
from metainformant.visualization import pairplot_dataframe
import pandas as pd
import numpy as np

df = pd.DataFrame(np.random.random((100, 3)), columns=['A', 'B', 'C'])
grid = pairplot_dataframe(df)
```

### `parallel_coordinates_plot(data, class_column=None, *, ax=None, title='Parallel Coordinates Plot', **kwargs)`

Create a parallel coordinates plot.

### `radar_chart(categories, values, *, ax=None, title='Radar Chart', **kwargs)`

Create a radar chart (spider chart).

### `splom_plot(data, *, hue=None, figsize=None, **kwargs)`

Create a scatter plot matrix (SPLOM).

### `scatter_3d(x, y, z, *, color=None, ax=None, title='3D Scatter Plot', **kwargs)`

Create a 3D scatter plot.

