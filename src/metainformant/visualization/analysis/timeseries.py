"""Time series analysis visualization functions.

This module provides specialized plotting functions for time series data
including basic plots, autocorrelation analysis, seasonal decomposition,
forecasting visualizations, and trend analysis.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes

from metainformant.core.data import validation
from metainformant.core.io import paths
from metainformant.core.utils import logging
from metainformant.visualization.config.conventions import save_figure_deterministic
from metainformant.visualization.config.figures import panel_grid, root_figure

logger = logging.get_logger(__name__)

# Optional imports with graceful fallbacks
try:
    from statsmodels.graphics.tsaplots import plot_acf
    from statsmodels.tsa.seasonal import seasonal_decompose

    HAS_STATSMODELS = True
except ImportError:
    seasonal_decompose = None
    plot_acf = None
    HAS_STATSMODELS = False


def _save_plot(ax: Axes, output_path: str | Path, label: str) -> str:
    """Ensure the output directory exists and save the plotted figure deterministically.

    Consolidates the repeated ensure_directory / save_figure_deterministic / logger
    triple used by every plot function in this module. Saves ``ax.figure`` so the
    plotted figure is written even when a different pyplot figure is current.
    """
    paths.ensure_directory(Path(output_path).parent)
    save_figure_deterministic(root_figure(ax), output_path, dpi=300, bbox_inches="tight")
    logger.info(f"{label} saved to {output_path}")
    return str(output_path)


def plot_time_series(
    data: pd.DataFrame, *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs: Any
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
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (12, 6)))

    # Plot each column as a separate line
    for column in data.columns:
        ax.plot(data.index, data[column], label=column, **kwargs)

    ax.set_xlabel("Time")
    ax.set_ylabel("Value")
    ax.set_title("Time Series Plot")

    if len(data.columns) > 1:
        ax.legend()

    # Format x-axis for time series
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")

    if output_path:
        _save_plot(ax, output_path, "Time series plot")

    return ax


def plot_autocorrelation(
    data: pd.Series, *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs: Any
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
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (10, 6)))

    # Use statsmodels ACF plot
    plot_acf(data.values, ax=ax, **kwargs)

    ax.set_title("Autocorrelation Function (ACF)")

    if output_path:
        _save_plot(ax, output_path, "Autocorrelation plot")

    return ax


def plot_seasonal_decomposition(
    data: pd.Series,
    *,
    model: str = "additive",
    period: int | None = None,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs: Any,
) -> Axes:
    """Create a seasonal decomposition plot with four stacked panels: original, trend, seasonal and residual.

    Args:
        data: Pandas Series with time series data
        model: Decomposition model ('additive' or 'multiplicative')
        period: Period for seasonal decomposition. If None, statsmodels infers it from the frequency of the
            series' datetime index, and raises ValueError when the index has no frequency.
        ax: Optional matplotlib axes that should hold the four panels: they are laid out inside the grid slot of
            ``ax``, which is hidden. Creates a new figure if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to ``statsmodels.tsa.seasonal.seasonal_decompose`` (for example
            ``two_sided`` or ``extrapolate_trend``). ``figsize`` sizes a new figure.

    Returns:
        The axes of the first (original series) panel

    Raises:
        ValueError: If data is not suitable for decomposition
        ImportError: If statsmodels is not available
    """
    validation.validate_type(data, pd.Series, "data")

    if data.empty:
        raise ValueError("Time series data Series cannot be empty")

    if not HAS_STATSMODELS:
        raise ImportError("statsmodels required for seasonal decomposition")

    figsize = kwargs.pop("figsize", (12, 10))

    # Create decomposition
    decomposition = seasonal_decompose(data, model=model, period=period, **kwargs)

    axes = panel_grid(ax, 4, max_cols=1, figsize=figsize, sharex=True)

    # Plot components
    axes[0].plot(data.index, data.values, label="Original")
    axes[0].set_title("Original Time Series")
    axes[0].legend()

    axes[1].plot(data.index, decomposition.trend, label="Trend", color="orange")
    axes[1].set_title("Trend Component")
    axes[1].legend()

    axes[2].plot(data.index, decomposition.seasonal, label="Seasonal", color="green")
    axes[2].set_title("Seasonal Component")
    axes[2].legend()

    axes[3].plot(data.index, decomposition.resid, label="Residual", color="red")
    axes[3].set_title("Residual Component")
    axes[3].legend()

    if ax is None:
        root_figure(axes[0]).tight_layout()

    if output_path:
        _save_plot(axes[0], output_path, "Seasonal decomposition plot")

    return axes[0]  # Return first axis for consistency


def plot_forecast(
    data: pd.Series,
    forecast: pd.Series,
    *,
    confidence_intervals: pd.DataFrame | None = None,
    confidence_level: float | None = None,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs: Any,
) -> Axes:
    """Create a forecast plot with optional confidence intervals.

    Args:
        data: Pandas Series with historical time series data
        forecast: Pandas Series with forecast values
        confidence_intervals: DataFrame with 'lower' and 'upper' columns, one row per forecast value
        confidence_level: Nominal coverage of ``confidence_intervals`` as a fraction (0.95 for a 95% interval),
            used only in the legend. The level is not known to this function, so without it the band is labelled
            "Confidence interval" rather than being given a level nobody supplied.
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to plotting.

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If data and forecast are empty, ``confidence_intervals`` lacks 'lower' or 'upper' or does not
            have one row per forecast value, or ``confidence_level`` is outside (0, 1)
    """
    validation.validate_type(data, pd.Series, "data")
    validation.validate_type(forecast, pd.Series, "forecast")

    if data.empty or forecast.empty:
        raise ValueError("Data and forecast Series cannot be empty")

    if confidence_intervals is not None:
        validation.validate_type(confidence_intervals, pd.DataFrame, "confidence_intervals")
        missing = [column for column in ("lower", "upper") if column not in confidence_intervals.columns]
        if missing:
            raise ValueError(f"confidence_intervals must have 'lower' and 'upper' columns; missing {missing}")
        if len(confidence_intervals) != len(forecast):
            raise ValueError(
                f"confidence_intervals has {len(confidence_intervals)} rows for {len(forecast)} forecast values"
            )
    if confidence_level is not None and not 0.0 < confidence_level < 1.0:
        raise ValueError(f"confidence_level must be a fraction between 0 and 1, got {confidence_level}")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (12, 6)))

    # Caller keyword arguments override these per-line defaults instead of colliding with them
    history_kwargs: dict[str, Any] = {"color": "blue", **kwargs}
    forecast_kwargs: dict[str, Any] = {"color": "red", "linestyle": "--", **kwargs}

    # Plot historical data
    ax.plot(data.index, data.values, label="Historical", **history_kwargs)

    # Plot forecast
    ax.plot(forecast.index, forecast.values, label="Forecast", **forecast_kwargs)

    # Plot confidence intervals if provided
    if confidence_intervals is not None:
        ax.fill_between(
            forecast.index,
            confidence_intervals["lower"],
            confidence_intervals["upper"],
            alpha=0.3,
            color="red",
            label="Confidence interval" if confidence_level is None else f"{confidence_level:.0%} CI",
        )

    ax.set_xlabel("Time")
    ax.set_ylabel("Value")
    ax.set_title("Time Series Forecast")
    ax.legend()

    # Format x-axis
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")

    if output_path:
        _save_plot(ax, output_path, "Forecast plot")

    return ax


def plot_trend_analysis(
    data: pd.Series, *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs: Any
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
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (12, 6)))

    # Plot original data (caller keyword arguments override the default alpha)
    original_kwargs: dict[str, Any] = {"alpha": 0.7, **kwargs}
    ax.plot(data.index, data.values, label="Original", **original_kwargs)

    # Fit linear trend
    x_numeric = np.arange(len(data))
    coefficients = np.polyfit(x_numeric, data.values, 1)
    trend_line = np.polyval(coefficients, x_numeric)

    # Plot trend line
    ax.plot(data.index, trend_line, label=f"Trend (slope: {coefficients[0]:.4f})", color="red", linewidth=2)

    # Add trend information as text. A slope that is zero up to the rounding of the fit is flat: a constant
    # series has a slope of exactly 0.0 (or noise around it), which is neither increasing nor decreasing.
    slope = coefficients[0]
    flat_tolerance = 1e-12 * max(float(np.ptp(data.values)), 1.0)
    if abs(slope) * max(len(data) - 1, 1) <= flat_tolerance:
        trend_direction = "flat"
    else:
        trend_direction = "increasing" if slope > 0 else "decreasing"
    ax.text(
        0.02,
        0.98,
        f"Trend: {trend_direction}\nSlope: {slope:.4f}",
        transform=ax.transAxes,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
    )

    ax.set_xlabel("Time")
    ax.set_ylabel("Value")
    ax.set_title("Trend Analysis")
    ax.legend()

    # Format x-axis
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")

    if output_path:
        _save_plot(ax, output_path, "Trend analysis plot")

    return ax
