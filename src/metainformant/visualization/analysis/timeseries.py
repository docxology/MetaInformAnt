"""Time series analysis visualization functions.

This module provides specialized plotting functions for time series data
including basic plots, autocorrelation analysis, seasonal decomposition,
forecasting visualizations, and trend analysis.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.axes import Axes

from metainformant.core import logging, paths, validation

logger = logging.get_logger(__name__)

# Optional imports with graceful fallbacks
try:
    from statsmodels.tsa.seasonal import seasonal_decompose
    from statsmodels.graphics.tsaplots import plot_acf
    HAS_STATSMODELS = True
except ImportError:
    seasonal_decompose = None
    plot_acf = None
    HAS_STATSMODELS = False


def plot_time_series(
    data: pd.DataFrame,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs
) -> Axes:
    """Create a basic time series plot.

    Args:
        data: DataFrame with datetime index and time series columns
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib plot().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If data is not properly formatted
    """
    validation.validate_type(data, pd.DataFrame, "data")

    if data.empty:
        raise ValueError("Time series data DataFrame cannot be empty")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop('figsize', (12, 6)))

    # Plot each column as a separate line
    for column in data.columns:
        ax.plot(data.index, data[column],
               label=column, **kwargs)

    ax.set_xlabel('Time')
    ax.set_ylabel('Value')
    ax.set_title('Time Series Plot')

    if len(data.columns) > 1:
        ax.legend()

    # Format x-axis for time series
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right')

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Time series plot saved to {output_path}")

    return ax


def plot_autocorrelation(
    data: pd.Series,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs
) -> Axes:
    """Create an autocorrelation function (ACF) plot.

    Args:
        data: Pandas Series with time series data
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to ACF plotting.

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If data is not a proper time series
        ImportError: If statsmodels is not available
    """
    validation.validate_type(data, pd.Series, "data")

    if data.empty:
        raise ValueError("Time series data Series cannot be empty")

    if not HAS_STATSMODELS:
        raise ImportError("statsmodels required for autocorrelation plotting")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop('figsize', (10, 6)))

    # Use statsmodels ACF plot
    plot_acf(data.values, ax=ax, **kwargs)

    ax.set_title('Autocorrelation Function (ACF)')

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Autocorrelation plot saved to {output_path}")

    return ax


def plot_seasonal_decomposition(
    data: pd.Series,
    *,
    model: str = 'additive',
    period: int | None = None,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs
) -> Axes:
    """Create a seasonal decomposition plot.

    Args:
        data: Pandas Series with time series data
        model: Decomposition model ('additive' or 'multiplicative')
        period: Period for seasonal decomposition (auto-detected if None)
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to decomposition.

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If data is not suitable for decomposition
        ImportError: If statsmodels is not available
    """
    validation.validate_type(data, pd.Series, "data")

    if data.empty:
        raise ValueError("Time series data Series cannot be empty")

    if not HAS_STATSMODELS:
        raise ImportError("statsmodels required for seasonal decomposition")

    # Create decomposition
    decomposition = seasonal_decompose(data, model=model, period=period)

    if ax is None:
        fig, axes = plt.subplots(4, 1, figsize=kwargs.pop('figsize', (12, 10)),
                                sharex=True)
    else:
        # If single ax provided, create subplots anyway for decomposition
        fig = ax.figure
        axes = []
        for i in range(4):
            axes.append(fig.add_subplot(4, 1, i+1, sharex=axes[0] if axes else None))

    # Plot components
    axes[0].plot(data.index, data.values, label='Original')
    axes[0].set_title('Original Time Series')
    axes[0].legend()

    axes[1].plot(data.index, decomposition.trend, label='Trend', color='orange')
    axes[1].set_title('Trend Component')
    axes[1].legend()

    axes[2].plot(data.index, decomposition.seasonal, label='Seasonal', color='green')
    axes[2].set_title('Seasonal Component')
    axes[2].legend()

    axes[3].plot(data.index, decomposition.resid, label='Residual', color='red')
    axes[3].set_title('Residual Component')
    axes[3].legend()

    plt.tight_layout()

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Seasonal decomposition plot saved to {output_path}")

    return axes[0]  # Return first axis for consistency


def plot_forecast(
    data: pd.Series,
    forecast: pd.Series,
    *,
    confidence_intervals: pd.DataFrame | None = None,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs
) -> Axes:
    """Create a forecast plot with optional confidence intervals.

    Args:
        data: Pandas Series with historical time series data
        forecast: Pandas Series with forecast values
        confidence_intervals: DataFrame with 'lower' and 'upper' columns for CI
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to plotting.

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If data and forecast don't align properly
    """
    validation.validate_type(data, pd.Series, "data")
    validation.validate_type(forecast, pd.Series, "forecast")

    if data.empty or forecast.empty:
        raise ValueError("Data and forecast Series cannot be empty")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop('figsize', (12, 6)))

    # Plot historical data
    ax.plot(data.index, data.values, label='Historical', color='blue', **kwargs)

    # Plot forecast
    ax.plot(forecast.index, forecast.values, label='Forecast',
           color='red', linestyle='--', **kwargs)

    # Plot confidence intervals if provided
    if confidence_intervals is not None:
        validation.validate_type(confidence_intervals, pd.DataFrame, "confidence_intervals")
        if 'lower' in confidence_intervals.columns and 'upper' in confidence_intervals.columns:
            ax.fill_between(forecast.index,
                          confidence_intervals['lower'],
                          confidence_intervals['upper'],
                          alpha=0.3, color='red', label='95% CI')

    ax.set_xlabel('Time')
    ax.set_ylabel('Value')
    ax.set_title('Time Series Forecast')
    ax.legend()

    # Format x-axis
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right')

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Forecast plot saved to {output_path}")

    return ax


def plot_trend_analysis(
    data: pd.Series,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs
) -> Axes:
    """Create a trend analysis plot with fitted trend line.

    Args:
        data: Pandas Series with time series data
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to plotting.

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If data is not suitable for trend analysis
    """
    validation.validate_type(data, pd.Series, "data")

    if data.empty:
        raise ValueError("Time series data Series cannot be empty")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop('figsize', (12, 6)))

    # Plot original data
    ax.plot(data.index, data.values, label='Original', alpha=0.7, **kwargs)

    # Fit linear trend
    x_numeric = np.arange(len(data))
    coefficients = np.polyfit(x_numeric, data.values, 1)
    trend_line = np.polyval(coefficients, x_numeric)

    # Plot trend line
    ax.plot(data.index, trend_line, label=f'Trend (slope: {coefficients[0]:.4f})',
           color='red', linewidth=2)

    # Add trend information as text
    slope = coefficients[0]
    trend_direction = "increasing" if slope > 0 else "decreasing"
    ax.text(0.02, 0.98, f'Trend: {trend_direction}\nSlope: {slope:.4f}',
           transform=ax.transAxes, verticalalignment='top',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax.set_xlabel('Time')
    ax.set_ylabel('Value')
    ax.set_title('Trend Analysis')
    ax.legend()

    # Format x-axis
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right')

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Trend analysis plot saved to {output_path}")

    return ax






