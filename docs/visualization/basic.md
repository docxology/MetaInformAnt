# Basic Plots

Basic plotting functions for simple visualizations including line plots, scatter plots, bar charts, pie charts, area plots, and heatmaps. They live in `metainformant.visualization.plots.basic` and share one calling convention:

- Numeric inputs must be NumPy arrays. A plain list (or `None` where an array is required) raises `metainformant.core.utils.errors.ValidationError`; arrays of different lengths raise `ValueError`.
- `ax=None` creates a new figure. Use `figsize=(width, height)` to size it (`figsize` is only valid when the function creates the figure, not together with `ax=`); pass `ax=` to draw on existing axes instead.
- `output_path=` saves the figure at 300 dpi and creates missing parent directories.
- Every other keyword argument is forwarded to the matplotlib call named under each function, so matplotlib's own argument names apply (`color`, `linestyle`, `alpha`, `label`, ...).
- The return value is the `matplotlib.axes.Axes`.

## Functions

### `lineplot(x, y=None, *, ax=None, output_path=None, **kwargs)`

Create a line plot with `Axes.plot`. `x` and `y` must be NumPy arrays of equal length. With `y=None`, `x` is plotted against its indices.

**Example:**
```python
import numpy as np
from metainformant.visualization.plots.basic import lineplot

ax = lineplot(np.array([1, 4, 2, 8, 5]), label="Data", linestyle="--", color="tab:blue")
ax.set_xlabel("Index")
ax.set_ylabel("Value")
ax.legend()
```

### `scatter_plot(x, y, *, ax=None, output_path=None, **kwargs)`

Create a scatter plot. `x` and `y` must be NumPy arrays of equal length; extra keyword
arguments (for example `color`, `s`, `alpha`) are passed to `matplotlib.axes.Axes.scatter`.

**Example:**
```python
import numpy as np
from metainformant.visualization.plots.basic import scatter_plot

ax = scatter_plot(np.array([1, 2, 3]), np.array([4, 5, 6]), alpha=0.7)
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_title("Scatter")
```

### `bar_plot(x, height, *, ax=None, output_path=None, **kwargs)`

Create a bar plot with `Axes.bar`. `x` is a list or NumPy array of bar positions or category labels (a tuple raises `ValueError`); `height` must be a NumPy array of the same length. Extra keyword arguments (`color`, `alpha`, `width`, `bottom`, `edgecolor`, ...) go to `Axes.bar`; there is no `horizontal` option.

**Example:**
```python
import numpy as np
from metainformant.visualization.plots.basic import bar_plot

ax = bar_plot(["A", "B", "C"], np.array([10, 20, 15]), color="tab:green", alpha=0.7)
ax.set_ylabel("Count")
```

### `pie_chart(sizes, labels=None, *, ax=None, output_path=None, **kwargs)`

Create a pie chart with `Axes.pie`, drawn with an equal aspect ratio. `sizes` must be a NumPy array; `labels`, if given, must be a `list` with one entry per slice. Extra keyword arguments (`autopct`, `colors`, `startangle`, ...) go to `Axes.pie`. There is no `title` option; call `ax.set_title`.

**Example:**
```python
import numpy as np
from metainformant.visualization.plots.basic import pie_chart

ax = pie_chart(np.array([30, 25, 45]), ["A", "B", "C"], autopct="%1.1f%%")
ax.set_title("Composition")
```

### `area_plot(x, y, *, ax=None, output_path=None, **kwargs)`

Create an area plot (a filled line plot) with `Axes.fill_between(x, y)`. `x` and `y` must be NumPy arrays of equal length. Extra keyword arguments (`color`, `alpha`, `label`, ...) go to `Axes.fill_between`; there are no `xlabel`, `ylabel` or `title` options, so use the axes methods.

**Example:**
```python
import numpy as np
from metainformant.visualization.plots.basic import area_plot

ax = area_plot(np.array([1, 2, 3, 4]), np.array([1, 4, 2, 3]), color="tab:orange", alpha=0.5)
ax.set_xlabel("x")
ax.set_ylabel("y")
```

### `step_plot(x, y, *, ax=None, output_path=None, **kwargs)`

Create a step plot with `Axes.step`. `x` and `y` must be NumPy arrays of equal length. Extra keyword arguments (`where='pre'`, `'mid'` or `'post'`, `label`, `color`, ...) go to `Axes.step`.

**Example:**
```python
import numpy as np
from metainformant.visualization.plots.basic import step_plot

ax = step_plot(np.array([1, 2, 3, 4]), np.array([1, 4, 2, 3]), where="post")
```

### `heatmap(data, *, cmap='viridis', ax=None, output_path=None, **kwargs)`

Create a heatmap of a 2-D NumPy array with `Axes.imshow` and add a colorbar. `data` must be a 2-D `numpy.ndarray`: a DataFrame raises `ValidationError` and an array that is not 2-D raises `ValueError`. Extra keyword arguments (`vmin`, `vmax`, `aspect`, `interpolation`, ...) go to `Axes.imshow`. There are no `annot` or `cbar` options; for an annotated correlation matrix use `correlation_heatmap` from [Statistical Plots](./statistical.md).

**Example:**
```python
import numpy as np
from metainformant.visualization.plots.basic import heatmap

data = np.random.default_rng(0).random((10, 10))
ax = heatmap(data, cmap="magma", vmin=0, vmax=1)
```
