# Multi-dimensional Plots

Multi-dimensional visualization functions including pair plots (scatter plot matrices), parallel coordinates, radar charts, and 3D scatter plots. They live in `metainformant.visualization.plots.multidim`, take a `pandas.DataFrame` (only its numeric columns are plotted), and accept `output_path=` to save the figure at 300 dpi.

## Functions

### `plot_pairwise_relationships(data, *, ax=None, output_path=None, **kwargs)`

Pairplot (scatter plot matrix) of the numeric columns of a DataFrame. Requires seaborn; extra
keyword arguments are passed to `seaborn.pairplot`, and the first axes of the grid is returned.
`seaborn.pairplot` creates its own figure, so the `ax` argument is ignored.

**Example:**
```python
from metainformant.visualization.plots.multidim import plot_pairwise_relationships
import pandas as pd
import numpy as np

df = pd.DataFrame(np.random.random((100, 3)), columns=['A', 'B', 'C'])
ax = plot_pairwise_relationships(df)
```

### `plot_parallel_coordinates(data, *, ax=None, output_path=None, **kwargs)`

Create a parallel coordinates plot: one semi-transparent line per row across the numeric columns, each column min-max scaled to 0-1 (at least two numeric columns are required). If the DataFrame has a text or categorical column (the first one, or the column named by `color='<column>'`), the rows are plotted group by group and each line is labelled with its group name in the legend; line colors cycle per row and do not encode the group. `figsize` and `color` are the only keyword arguments read; any others are ignored.

**Example:**
```python
from metainformant.visualization.plots.multidim import plot_parallel_coordinates
import pandas as pd
import numpy as np

rng = np.random.default_rng(0)
df = pd.DataFrame(rng.random((20, 4)), columns=['A', 'B', 'C', 'D'])
ax = plot_parallel_coordinates(df)
```

### `plot_radar_chart(data, *, ax=None, output_path=None, **kwargs)`

Create a radar chart (spider chart). Each numeric column of `data` is one axis of the chart (at least three are required) and each row is one polygon, labelled by its index in the legend. When `ax` is `None` the function creates polar axes. `figsize` is the only keyword argument it reads; any others are ignored.

**Example:**
```python
from metainformant.visualization.plots.multidim import plot_radar_chart
import pandas as pd

data = pd.DataFrame(
    {
        'Trait_A': [0.8, 0.6, 0.7],
        'Trait_B': [0.6, 0.8, 0.7],
        'Trait_C': [0.9, 0.5, 0.7],
        'Trait_D': [0.4, 0.6, 0.7],
    },
    index=['Sample_1', 'Sample_2', 'Sample_3'],
)
ax = plot_radar_chart(data)
```

### `plot_3d_scatter(data, *, x_col=None, y_col=None, z_col=None, ax=None, output_path=None, **kwargs)`

Create a 3D scatter plot of three numeric columns (the first three numeric columns unless `x_col`, `y_col` and `z_col` name them; at least three numeric columns are required). Pass `c='<column>'` to color the points by another column of `data`. Remaining keyword arguments (`s`, `alpha`, `cmap`, ...) are passed to the 3D `scatter`. The function creates 3D axes when `ax` is `None`.

**Example:**
```python
from metainformant.visualization.plots.multidim import plot_3d_scatter
import pandas as pd
import numpy as np

rng = np.random.default_rng(0)
df = pd.DataFrame(rng.random((50, 3)), columns=['A', 'B', 'C'])
ax = plot_3d_scatter(df, x_col='A', y_col='B', z_col='C', c='A', s=20, alpha=0.8)
```
