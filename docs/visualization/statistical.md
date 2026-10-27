# Statistical Plots

Statistical plotting functions for data analysis including histograms, box plots, violin plots, Q-Q plots, correlation heatmaps, density plots, and statistical diagnostic plots. They live in `metainformant.visualization.analysis.statistical` and share one calling convention:

- Inputs are type-checked: a NumPy array (`numpy.ndarray`) for numeric data, a **list** of NumPy arrays for grouped data (box, violin and ridge plots), and a `pandas.DataFrame` for the correlation heatmap. A wrong type raises `metainformant.core.utils.errors.ValidationError`; empty or mismatched inputs raise `ValueError`.
- `ax=None` creates a new figure (`figsize=(width, height)` sizes it); pass `ax=` to draw on existing axes. `output_path=` saves the figure at 300 dpi.
- Every other keyword argument is forwarded to the matplotlib or seaborn call named under each function.
- The return value is the `matplotlib.axes.Axes`.

## Functions

### `histogram(data, *, bins=30, ax=None, output_path=None, **kwargs)`

Create a histogram with `Axes.hist`. `data` must be a non-empty NumPy array. Extra keyword arguments (`density`, `alpha`, `color`, `range`, ...) go to `Axes.hist`.

**Example:**
```python
from metainformant.visualization.analysis.statistical import histogram
import numpy as np

data = np.random.default_rng(0).normal(0, 1, 1000)
ax = histogram(data, bins=30, density=True, color="tab:blue", alpha=0.7)
```

### `box_plot(data, *, ax=None, output_path=None, **kwargs)`

Create a box plot with `Axes.boxplot`. `data` is a list of non-empty NumPy arrays, one box per array, drawn at x positions 1, 2, ... Extra keyword arguments (`notch`, `showfliers`, `patch_artist`, ...) go to `Axes.boxplot`. Matplotlib renamed the box-label argument (`labels` became `tick_labels` in 3.9, and the old name no longer works in current releases), so label the boxes with `ax.set_xticklabels`, which works on every supported matplotlib.

**Example:**
```python
from metainformant.visualization.analysis.statistical import box_plot
import numpy as np

rng = np.random.default_rng(0)
data = [rng.normal(0, 1, 100) for _ in range(3)]
ax = box_plot(data)
ax.set_xticklabels(["A", "B", "C"])
```

### `violin_plot(data, *, ax=None, output_path=None, **kwargs)`

Create a violin plot. `data` is a list of non-empty NumPy arrays, one violin per array. The arrays are stacked into a long-form table with groups numbered 1, 2, ... and drawn with `seaborn.violinplot`, so extra keyword arguments (`inner`, `cut`, `palette`, ...) are seaborn arguments. The categorical x axis sits at positions 0, 1, ... with tick labels "1", "2", ...; relabel it with `ax.set_xticks`. If seaborn cannot be imported the function logs a warning and draws a box plot instead.

**Example:**
```python
from metainformant.visualization.analysis.statistical import violin_plot
import numpy as np

rng = np.random.default_rng(0)
data = [rng.normal(0, 1, 100) for _ in range(3)]
ax = violin_plot(data, inner="quartile")
ax.set_xticks(range(3), labels=["A", "B", "C"])
```

### `qq_plot(data, *, distribution='norm', ax=None, output_path=None, **kwargs)`

Create a Q-Q plot of a sample against a theoretical distribution (`'norm'` or `'uniform'`).
`data` must be a NumPy array; use `distribution='uniform'` to check the distribution of p-values.

**Example:**
```python
from metainformant.visualization.analysis.statistical import qq_plot
import numpy as np

pvals = np.random.uniform(0, 1, 1000)
ax = qq_plot(pvals, distribution="uniform")
```

### `correlation_heatmap(data, *, ax=None, output_path=None, **kwargs)`

Create a heatmap of the Pearson correlation matrix of the numeric columns of a DataFrame. Extra
keyword arguments are passed to `seaborn.heatmap` (`annot` defaults to `True`, `cmap` to `'coolwarm'`).

**Example:**
```python
from metainformant.visualization.analysis.statistical import correlation_heatmap
import pandas as pd
import numpy as np

df = pd.DataFrame(np.random.random((10, 5)))
ax = correlation_heatmap(df)
```

### `density_plot(data, *, ax=None, output_path=None, **kwargs)`

Create a kernel density estimate plot of a non-empty NumPy array with `seaborn.kdeplot` (a normalized histogram if seaborn cannot be imported). Extra keyword arguments (`fill`, `alpha`, `color`, `bw_adjust`, ...) are seaborn arguments. Seaborn labels the y axis "Density"; no title or x label is set, so use the axes methods.

**Example:**
```python
from metainformant.visualization.analysis.statistical import density_plot
import numpy as np

data = np.random.default_rng(0).normal(0, 1, 1000)
ax = density_plot(data, fill=True, alpha=0.5)
```

### `ridge_plot(data, *, ax=None, output_path=None, **kwargs)`

Overlay one density curve per array of `data` (a list of non-empty NumPy arrays) on a single axes, labelled "Group 1", "Group 2", ... in the legend. The curves are filled and overlaid, not vertically offset, and the function fixes `fill`, `alpha` and `label` itself; other keyword arguments go to `seaborn.kdeplot` for every group.

**Example:**
```python
from metainformant.visualization.analysis.statistical import ridge_plot
import numpy as np

rng = np.random.default_rng(0)
data = [rng.normal(i, 1, 100) for i in range(3)]
ax = ridge_plot(data)
```

### `roc_curve(y_true, y_scores, *, ax=None, output_path=None, **kwargs)`

Plot a ROC curve with the AUC in the legend and the chance diagonal. `y_true` (binary labels) and `y_scores` must be NumPy arrays of equal length; scikit-learn computes the curve. Extra keyword arguments go to `Axes.plot` for the ROC line. The title is fixed to "ROC Curve" (there is no `title` option).

**Example:**
```python
from metainformant.visualization.analysis.statistical import roc_curve
import numpy as np

y_true = np.array([0, 1, 0, 1, 1])
y_scores = np.array([0.1, 0.9, 0.2, 0.8, 0.7])
ax = roc_curve(y_true, y_scores)
```

### `precision_recall_curve(y_true, y_scores, *, ax=None, output_path=None, **kwargs)`

Plot a precision-recall curve with the average precision in the legend. Inputs are the same as for `roc_curve` (NumPy arrays of equal length; scikit-learn computes the curve), extra keyword arguments go to `Axes.plot`, and the title is fixed to "Precision-Recall Curve".

**Example:**
```python
from metainformant.visualization.analysis.statistical import precision_recall_curve
import numpy as np

y_true = np.array([0, 1, 0, 1, 1])
y_scores = np.array([0.1, 0.9, 0.2, 0.8, 0.7])
ax = precision_recall_curve(y_true, y_scores)
```

### `residual_plot(y_true, y_pred, *, ax=None, output_path=None, **kwargs)`

Create a residual plot for regression diagnostics: residuals (`y_true - y_pred`) against the predicted values, with a dashed zero line. `y_true` and `y_pred` must be NumPy arrays of equal length. Extra keyword arguments go to `Axes.scatter`; the title is fixed to "Residual Plot".

**Example:**
```python
from metainformant.visualization.analysis.statistical import residual_plot
import numpy as np

y_true = np.array([1, 2, 3, 4, 5])
y_pred = np.array([1.1, 1.9, 3.2, 3.8, 5.1])
ax = residual_plot(y_true, y_pred)
```

### `leverage_plot(X, y, *, ax=None, output_path=None, **kwargs)`

Fit an ordinary least squares model with scikit-learn and plot standardized residuals against a leverage proxy. `X` (samples x features) and `y` must be NumPy arrays with the same number of rows. The leverage value is each row's share of the total squared magnitude of `X` (`np.sum(X**2, axis=1) / np.sum(X**2)`), a simple approximation rather than the hat-matrix diagonal, so treat it as a rough screen for unusual rows. Extra keyword arguments go to `Axes.scatter`; the title is fixed to "Leverage Plot".

**Example:**
```python
from metainformant.visualization.analysis.statistical import leverage_plot
import numpy as np

rng = np.random.default_rng(0)
X = rng.random((100, 2))
y = rng.random(100)
ax = leverage_plot(X, y)
```
