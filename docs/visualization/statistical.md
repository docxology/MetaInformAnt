# Statistical Plots

Statistical plotting functions for data analysis including histograms, box plots, violin plots, Q-Q plots, correlation heatmaps, density plots, and statistical diagnostic plots.

## Functions

### `histogram(data, *, bins=30, ax=None, density=False, alpha=0.7, color='blue', xlabel=None, ylabel=None, title=None)`

Create a histogram.

**Example:**
```python
from metainformant.visualization import histogram
import numpy as np

data = np.random.normal(0, 1, 1000)
ax = histogram(data, bins=30)
```

### `box_plot(data, *, ax=None, positions=None, labels=None, xlabel=None, ylabel=None, title=None)`

Create a box plot.

**Example:**
```python
from metainformant.visualization import box_plot
import numpy as np

data = [np.random.normal(0, 1, 100) for _ in range(3)]
ax = box_plot(data, labels=["A", "B", "C"])
```

### `violin_plot(data, *, ax=None, positions=None, labels=None, xlabel=None, ylabel=None, title=None)`

Create a violin plot.

**Example:**
```python
from metainformant.visualization import violin_plot
import numpy as np

data = [np.random.normal(0, 1, 100) for _ in range(3)]
ax = violin_plot(data, labels=["A", "B", "C"])
```

### `qq_plot(p_values, *, ax=None, title='Q-Q Plot', **kwargs)`

Create a Q-Q plot for p-value distribution analysis.

**Example:**
```python
from metainformant.visualization import qq_plot
import numpy as np

pvals = np.random.uniform(0, 1, 1000)
ax = qq_plot(pvals.tolist())
```

### `correlation_heatmap(data, *, method='pearson', cmap='RdBu_r', annot=True, ax=None, **kwargs)`

Create a correlation heatmap.

**Example:**
```python
from metainformant.visualization import correlation_heatmap
import pandas as pd
import numpy as np

df = pd.DataFrame(np.random.random((10, 5)))
ax = correlation_heatmap(df)
```

### `density_plot(data, *, ax=None, fill=True, alpha=0.5, color=None, xlabel=None, ylabel=None, title=None)`

Create a density plot (kernel density estimation).

**Example:**
```python
from metainformant.visualization import density_plot
import numpy as np

data = np.random.normal(0, 1, 1000)
ax = density_plot(data)
```

### `ridge_plot(data, labels=None, *, ax=None, overlap=0.5, **kwargs)`

Create a ridge plot (overlapping density plots).

**Example:**
```python
from metainformant.visualization import ridge_plot
import numpy as np

data = [np.random.normal(i, 1, 100) for i in range(3)]
ax = ridge_plot(data, labels=["A", "B", "C"])
```

### `roc_curve(y_true, y_scores, *, ax=None, title='ROC Curve', **kwargs)`

Create a ROC curve plot.

**Example:**
```python
from metainformant.visualization import roc_curve

y_true = [0, 1, 0, 1, 1]
y_scores = [0.1, 0.9, 0.2, 0.8, 0.7]
ax = roc_curve(y_true, y_scores)
```

### `precision_recall_curve(y_true, y_scores, *, ax=None, title='Precision-Recall Curve', **kwargs)`

Create a precision-recall curve plot.

**Example:**
```python
from metainformant.visualization import precision_recall_curve

y_true = [0, 1, 0, 1, 1]
y_scores = [0.1, 0.9, 0.2, 0.8, 0.7]
ax = precision_recall_curve(y_true, y_scores)
```

### `residual_plot(y_true, y_pred, *, ax=None, title='Residual Plot', **kwargs)`

Create a residual plot for regression diagnostics.

**Example:**
```python
from metainformant.visualization import residual_plot
import numpy as np

y_true = np.array([1, 2, 3, 4, 5])
y_pred = np.array([1.1, 1.9, 3.2, 3.8, 5.1])
ax = residual_plot(y_true, y_pred)
```

### `leverage_plot(X, y, *, ax=None, title='Leverage Plot', **kwargs)`

Create a leverage plot for regression diagnostics.

**Example:**
```python
from metainformant.visualization import leverage_plot
import numpy as np

X = np.random.random((100, 2))
y = np.random.random(100)
ax = leverage_plot(X, y)
```

