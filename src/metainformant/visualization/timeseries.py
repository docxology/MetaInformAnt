"""Time series visualization functions.

This module provides visualization functions for time series data including
time series plots, autocorrelation plots, seasonal decomposition plots,
forecast plots, and trend analysis.
"""

from __future__ import annotations

from typing import Sequence

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Use non-interactive backend by default for tests/headless
matplotlib.use("Agg", force=True)


def time_series_plot(
    time_points: Sequence[float],
    values: Sequence[float],
    *,
    ax: plt.Axes | None = None,
    title: str = "Time Series",
    xlabel: str = "Time",
    ylabel: str = "Value",
    **kwargs
) -> plt.Axes:
    """Create a time series plot.

    Args:
        time_points: Time points
        values: Values at each time point
        ax: Matplotlib axes (creates new if None)
        title: Plot title
        xlabel: X-axis label
        ylabel: Y-axis label
        **kwargs: Additional arguments for plot

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import time_series_plot
        >>> import numpy as np
        >>> time = np.arange(0, 10, 0.1)
        >>> values = np.sin(time)
        >>> ax = time_series_plot(time, values)
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 6))

    ax.plot(time_points, values, **kwargs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)

    return ax


def autocorrelation_plot(
    data: Sequence[float],
    *,
    max_lag: int | None = None,
    ax: plt.Axes | None = None,
    title: str = "Autocorrelation Plot",
    **kwargs
) -> plt.Axes:
    """Create an autocorrelation plot.

    Args:
        data: Time series data
        max_lag: Maximum lag to plot (if None, uses len(data) // 2)
        ax: Matplotlib axes (creates new if None)
        title: Plot title
        **kwargs: Additional arguments

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import autocorrelation_plot
        >>> import numpy as np
        >>> data = np.sin(np.arange(0, 100, 0.1))
        >>> ax = autocorrelation_plot(data)
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 6))

    data = np.array(data)
    n = len(data)
    
    if max_lag is None:
        max_lag = n // 2
    
    # Calculate autocorrelation
    autocorr = []
    for lag in range(max_lag + 1):
        if lag == 0:
            corr = 1.0
        else:
            corr = np.corrcoef(data[:-lag], data[lag:])[0, 1]
        autocorr.append(corr)
    
    lags = range(max_lag + 1)
    ax.plot(lags, autocorr, marker='o', **kwargs)
    ax.axhline(y=0, color='k', linestyle='-', linewidth=0.5)
    ax.axhline(y=0.95, color='r', linestyle='--', alpha=0.5, label='95% confidence')
    ax.axhline(y=-0.95, color='r', linestyle='--', alpha=0.5)
    
    ax.set_xlabel("Lag")
    ax.set_ylabel("Autocorrelation")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend()

    return ax


def seasonal_decomposition_plot(
    data: Sequence[float],
    period: int,
    *,
    ax: plt.Axes | None = None,
    title: str = "Seasonal Decomposition",
    **kwargs
) -> plt.Axes:
    """Create a seasonal decomposition plot.

    Args:
        data: Time series data
        period: Seasonal period
        ax: Matplotlib axes (creates new if None)
        title: Plot title
        **kwargs: Additional arguments

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import seasonal_decomposition_plot
        >>> import numpy as np
        >>> t = np.arange(0, 100)
        >>> data = np.sin(2 * np.pi * t / 12) + 0.1 * t + np.random.normal(0, 0.1, 100)
        >>> ax = seasonal_decomposition_plot(data, period=12)
    """
    if ax is None:
        fig, axes = plt.subplots(4, 1, figsize=(12, 10))
        ax = axes[0]  # Return first subplot
    else:
        # Create subplots within the provided axes
        fig = ax.figure
        axes = fig.subplots(4, 1)

    data = np.array(data)
    n = len(data)

    # Simple moving average for trend
    window = min(period, n // 4)
    if window % 2 == 0:
        window += 1
    
    trend = np.convolve(data, np.ones(window) / window, mode='same')
    
    # Detrend
    detrended = data - trend
    
    # Seasonal component (average over periods)
    n_periods = n // period
    seasonal = np.zeros(period)
    for i in range(period):
        seasonal[i] = np.mean([detrended[i + j * period] for j in range(n_periods) if i + j * period < n])
    
    # Extend seasonal to match data length
    seasonal_extended = np.tile(seasonal, n_periods + 1)[:n]
    
    # Residual
    residual = detrended - seasonal_extended
    
    # Plot components
    time = np.arange(n)
    axes[0].plot(time, data, label='Original', **kwargs)
    axes[0].set_ylabel("Original")
    axes[0].legend()
    axes[0].set_title(title)
    
    axes[1].plot(time, trend, label='Trend', **kwargs)
    axes[1].set_ylabel("Trend")
    axes[1].legend()
    
    axes[2].plot(time, seasonal_extended, label='Seasonal', **kwargs)
    axes[2].set_ylabel("Seasonal")
    axes[2].legend()
    
    axes[3].plot(time, residual, label='Residual', **kwargs)
    axes[3].set_ylabel("Residual")
    axes[3].set_xlabel("Time")
    axes[3].legend()
    
    plt.tight_layout()

    return axes[0]


def forecast_plot(
    time_points: Sequence[float],
    observed: Sequence[float],
    forecast: Sequence[float],
    confidence_intervals: tuple[Sequence[float], Sequence[float]] | None = None,
    *,
    ax: plt.Axes | None = None,
    title: str = "Forecast Plot",
    **kwargs
) -> plt.Axes:
    """Create a forecast plot with observed and predicted values.

    Args:
        time_points: Time points
        observed: Observed values
        forecast: Forecasted values
        confidence_intervals: Optional (lower, upper) confidence intervals
        ax: Matplotlib axes (creates new if None)
        title: Plot title
        **kwargs: Additional arguments

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import forecast_plot
        >>> import numpy as np
        >>> time = np.arange(0, 20)
        >>> observed = np.sin(time[:10])
        >>> forecast = np.sin(time[10:])
        >>> ax = forecast_plot(time, observed, forecast)
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 6))

    time_points = np.array(time_points)
    observed = np.array(observed)
    forecast = np.array(forecast)

    # Split time points
    n_observed = len(observed)
    time_observed = time_points[:n_observed]
    time_forecast = time_points[n_observed:]

    # Plot observed
    ax.plot(time_observed, observed, 'o-', label='Observed', color='blue', **kwargs)

    # Plot forecast
    ax.plot(time_forecast, forecast, 's-', label='Forecast', color='red', **kwargs)

    # Plot confidence intervals if provided
    if confidence_intervals:
        lower, upper = confidence_intervals
        ax.fill_between(time_forecast, lower, upper, alpha=0.3, color='red', label='Confidence Interval')

    # Add vertical line at forecast start
    if len(time_observed) > 0:
        ax.axvline(x=time_observed[-1], color='gray', linestyle='--', alpha=0.7, label='Forecast Start')

    ax.set_xlabel("Time")
    ax.set_ylabel("Value")
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3)

    return ax


def trend_plot(
    time_points: Sequence[float],
    values: Sequence[float],
    *,
    trend_type: str = "linear",
    ax: plt.Axes | None = None,
    title: str = "Trend Analysis",
    **kwargs
) -> plt.Axes:
    """Create a trend plot with trend line.

    Args:
        time_points: Time points
        values: Values at each time point
        trend_type: Type of trend ('linear', 'polynomial')
        ax: Matplotlib axes (creates new if None)
        title: Plot title
        **kwargs: Additional arguments

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import trend_plot
        >>> import numpy as np
        >>> time = np.arange(0, 10)
        >>> values = 0.5 * time + np.random.normal(0, 0.1, 10)
        >>> ax = trend_plot(time, values, trend_type='linear')
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 6))

    time_points = np.array(time_points)
    values = np.array(values)

    # Plot data
    ax.scatter(time_points, values, alpha=0.6, label='Data', **kwargs)

    # Fit trend
    if trend_type == "linear":
        coeffs = np.polyfit(time_points, values, 1)
        trend_line = np.polyval(coeffs, time_points)
        ax.plot(time_points, trend_line, 'r-', label='Linear Trend', linewidth=2)
    elif trend_type == "polynomial":
        coeffs = np.polyfit(time_points, values, 2)
        trend_line = np.polyval(coeffs, time_points)
        ax.plot(time_points, trend_line, 'r-', label='Polynomial Trend', linewidth=2)
    else:
        raise ValueError(f"Unknown trend type: {trend_type}")

    ax.set_xlabel("Time")
    ax.set_ylabel("Value")
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3)

    return ax

