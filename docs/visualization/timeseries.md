# Time Series Plots

Time series visualization functions including time series plots, autocorrelation plots, seasonal decomposition plots, forecast plots, and trend analysis.

## Functions

### `time_series_plot(time_points, values, *, ax=None, title='Time Series', xlabel='Time', ylabel='Value', **kwargs)`

Create a time series plot.

**Example:**
```python
from metainformant.visualization import time_series_plot
import numpy as np

time = np.arange(0, 10, 0.1)
values = np.sin(time)
ax = time_series_plot(time, values)
```

### `autocorrelation_plot(data, *, max_lag=None, ax=None, title='Autocorrelation Plot', **kwargs)`

Create an autocorrelation plot.

### `seasonal_decomposition_plot(data, period, *, ax=None, title='Seasonal Decomposition', **kwargs)`

Create a seasonal decomposition plot.

### `forecast_plot(time_points, observed, forecast, confidence_intervals=None, *, ax=None, title='Forecast Plot', **kwargs)`

Create a forecast plot with observed and predicted values.

### `trend_plot(time_points, values, *, trend_type='linear', ax=None, title='Trend Analysis', **kwargs)`

Create a trend plot with trend line.

