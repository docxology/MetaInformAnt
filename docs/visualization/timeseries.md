# Time Series Plots

Time series visualization functions including time series plots, autocorrelation plots, seasonal decomposition plots, forecast plots, and trend analysis. They live in `metainformant.visualization.analysis.timeseries`, return the `matplotlib.axes.Axes`, and accept `output_path=` to save the figure at 300 dpi. Inputs are type-checked: a `pandas.DataFrame` for `plot_time_series` and `pandas.Series` for the others, with the time axis in the index. The autocorrelation and decomposition plots need `statsmodels`, which is a core dependency.

## Functions

### `plot_time_series(data, *, ax=None, output_path=None, **kwargs)`

Create a time series plot. `data` is a non-empty DataFrame indexed by time; each column is drawn as a
separate line (labelled with the column name, with a legend when there is more than one column). Extra keyword
arguments are passed to `Axes.plot`.

**Example:**
```python
from metainformant.visualization.analysis.timeseries import plot_time_series
import numpy as np
import pandas as pd

time = np.arange(0, 10, 0.1)
data = pd.DataFrame({'signal': np.sin(time)}, index=time)
ax = plot_time_series(data)
```

### `plot_autocorrelation(data, *, ax=None, output_path=None, **kwargs)`

Plot the autocorrelation function (ACF) of a Series with `statsmodels.graphics.tsaplots.plot_acf`. Extra keyword arguments such as `lags=36` or `alpha=0.05` are passed to `plot_acf`.

**Example:**
```python
from metainformant.visualization.analysis.timeseries import plot_autocorrelation
import numpy as np
import pandas as pd

rng = np.random.default_rng(0)
series = pd.Series(np.sin(np.arange(120) * 2 * np.pi / 12) + rng.normal(0, 0.3, 120))
ax = plot_autocorrelation(series, lags=36)
```

### `plot_seasonal_decomposition(data, *, model='additive', period=None, ax=None, output_path=None, **kwargs)`

Decompose a Series with `statsmodels` `seasonal_decompose` and plot four stacked panels: original, trend, seasonal and residual. `model` is `'additive'` or `'multiplicative'` and `period` is the number of observations per cycle. `period` can be omitted only when the Series has a `DatetimeIndex` (or `PeriodIndex`) with a frequency; otherwise statsmodels raises `ValueError`. The function lays out its own four panels (given an `ax`, it adds them to that axes' figure rather than drawing on it), reads only `figsize` from `**kwargs`, and returns the top panel.

**Example:**
```python
from metainformant.visualization.analysis.timeseries import plot_seasonal_decomposition
import numpy as np
import pandas as pd

index = pd.date_range('2020-01-01', periods=48, freq='MS')
series = pd.Series(np.sin(np.arange(48) * 2 * np.pi / 12) + np.arange(48) * 0.05, index=index)
ax = plot_seasonal_decomposition(series, period=12)
```

### `plot_forecast(data, forecast, *, confidence_intervals=None, ax=None, output_path=None, **kwargs)`

Plot a historical Series (blue) and a forecast Series (red, dashed) on one axes. `confidence_intervals` is an optional DataFrame with `lower` and `upper` columns, row-aligned with `forecast`, drawn as a shaded band labelled "95% CI" (the label is fixed and the function does not compute intervals). Extra keyword arguments are passed to `Axes.plot` for both lines. `color` (blue for the history, red for the forecast) and `linestyle` (solid for the history, dashed for the forecast) are defaults that a keyword argument overrides, and the override then applies to both lines; `label` is set by the function, and the confidence band stays red.

**Example:**
```python
from metainformant.visualization.analysis.timeseries import plot_forecast
import numpy as np
import pandas as pd

index = pd.date_range('2021-01-01', periods=48, freq='MS')
series = pd.Series(np.sin(np.arange(48) * 2 * np.pi / 12) + np.arange(48) * 0.05, index=index)
history, future = series[:36], series[36:]
bounds = pd.DataFrame({'lower': future - 0.4, 'upper': future + 0.4})
ax = plot_forecast(history, future, confidence_intervals=bounds)
```

### `plot_trend_analysis(data, *, ax=None, output_path=None, **kwargs)`

Plot a Series together with a least-squares straight-line trend (only a linear trend is fitted). The slope is per observation (one index step), and is shown in the legend and in a text box with the trend direction. Extra keyword arguments are passed to `Axes.plot` for the Series (`label` is set by the function; `alpha` defaults to 0.7 and a keyword argument overrides it).

**Example:**
```python
from metainformant.visualization.analysis.timeseries import plot_trend_analysis
import numpy as np
import pandas as pd

index = pd.date_range('2021-01-01', periods=48, freq='MS')
series = pd.Series(np.sin(np.arange(48) * 2 * np.pi / 12) + np.arange(48) * 0.05, index=index)
ax = plot_trend_analysis(series)
```
