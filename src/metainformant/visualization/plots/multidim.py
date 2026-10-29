"""Multi-dimensional data visualization functions.

This module provides specialized plotting functions for visualizing
multi-dimensional datasets including pairwise relationships, parallel coordinates,
radar charts, and 3D scatter plots.
"""

from __future__ import annotations

from math import pi
from pathlib import Path
from typing import Any, cast

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.colors import is_color_like

from metainformant.core.data import validation
from metainformant.core.io import paths
from metainformant.core.utils import logging
from metainformant.visualization.config import palettes
from metainformant.visualization.config.conventions import save_figure_deterministic
from metainformant.visualization.config.figures import root_figure

logger = logging.get_logger(__name__)

# Optional imports with graceful fallbacks
try:
    import seaborn as sns

    HAS_SEABORN = True
except ImportError:
    sns = None
    HAS_SEABORN = False


def _save_fig(ax: Axes, output_path: str | Path | None, label: str) -> None:
    """Persist the plotted figure (``ax.figure``) when *output_path* is requested."""
    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        save_figure_deterministic(root_figure(ax), output_path, dpi=300, bbox_inches="tight")
        logger.info(f"{label} saved to {output_path}")


def plot_pairwise_relationships(
    data: pd.DataFrame, *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs: Any
) -> Axes:
    """Create a pairwise relationships plot (scatter plot matrix).

    Args:
        data: DataFrame with numeric columns for pairwise plotting
        ax: Not supported: seaborn's pairplot builds a grid of axes in a figure of its own and cannot draw onto
            one axes, so passing ``ax`` raises ValueError instead of being ignored.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to seaborn pairplot.

    Returns:
        The first axes of the grid (top left)

    Raises:
        ValueError: If data doesn't have enough numeric columns, or ``ax`` is given
        ImportError: If seaborn is not available
    """
    validation.validate_type(data, pd.DataFrame, "data")

    if ax is not None:
        raise ValueError("plot_pairwise_relationships builds its own figure and cannot draw onto ax; pass ax=None")

    if not HAS_SEABORN:
        raise ImportError("seaborn required for pairwise relationships plotting")

    # Select only numeric columns
    numeric_data = data.select_dtypes(include=[np.number])
    if numeric_data.shape[1] < 2:
        raise ValueError("Data must have at least 2 numeric columns for pairwise plotting")

    # Create pairplot (this creates its own figure)
    g = sns.pairplot(numeric_data, **kwargs)

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        save_figure_deterministic(g.figure, output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Pairwise relationships plot saved to {output_path}")

    # Return the first axes from the grid
    return cast("Axes", g.fig.axes[0] if g.fig.axes else None)


def plot_parallel_coordinates(
    data: pd.DataFrame, *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs: Any
) -> Axes:
    """Create a parallel coordinates plot: one line per row across the min-max scaled numeric columns.

    Every numeric column is scaled to 0-1 (a column with no variation is drawn at 0.5). The lines are coloured by
    group, one colour and one legend entry per distinct value of the grouping column.

    Args:
        data: DataFrame with numeric columns for parallel coordinates
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: ``color`` names the grouping column; if it is not a column but is a matplotlib colour, every line
            gets that colour, otherwise the lines are drawn plain. Without ``color`` the first non-numeric column
            groups the lines, if there is one. ``figsize`` sizes a new figure. Any other keyword argument is
            passed to matplotlib plot() for every line (``alpha``, ``linewidth``, ...).

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If data doesn't have enough numeric columns
    """
    validation.validate_type(data, pd.DataFrame, "data")

    # Select only numeric columns
    numeric_data = data.select_dtypes(include=[np.number])
    if numeric_data.shape[1] < 2:
        raise ValueError("Data must have at least 2 numeric columns for parallel coordinates")

    # Decide how the lines are coloured: by the values of a column, one colour for all, or the default cycle
    color = kwargs.pop("color", None)
    group_col = None
    if isinstance(color, str) and color in data.columns:
        group_col, color = color, None
    elif color is None:
        non_numeric = [column for column in data.columns if column not in numeric_data.columns]
        group_col = non_numeric[0] if non_numeric else None
    elif not is_color_like(color):
        logger.warning(f"color {color!r} is neither a column of the data nor a colour; drawing the lines plain")
        color = None

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (12, 6)))

    # Scale each numeric column to 0-1; a constant column has no range and sits mid-axis, a missing value stays missing
    lows, highs = numeric_data.min(), numeric_data.max()
    spans = highs - lows
    scaled = (numeric_data - lows) / spans.where(spans > 0)
    for column in spans.index[spans == 0]:
        scaled[column] = scaled[column].where(numeric_data[column].isna(), 0.5)

    line_kwargs: dict[str, Any] = {"alpha": 0.7, **kwargs}
    if color is not None:
        line_kwargs["color"] = color

    positions = np.arange(scaled.shape[1])
    if group_col is None:
        for _, row in scaled.iterrows():
            ax.plot(positions, row.to_numpy(), **line_kwargs)
    else:
        groups = list(pd.unique(data[group_col]))  # first-seen order; missing values form one group
        group_colors = dict(zip(range(len(groups)), palettes.categorical(len(groups))))
        for index, group in enumerate(groups):
            member = (data[group_col] == group) if pd.notna(group) else data[group_col].isna()
            lines = ax.plot(positions, scaled[member].to_numpy().T, **{**line_kwargs, "color": group_colors[index]})
            lines[0].set_label(str(group))  # one legend entry per group, not one per row
        if groups:
            ax.legend(title=str(group_col))

    # Set x-ticks to column names
    ax.set_xticks(range(len(numeric_data.columns)))
    ax.set_xticklabels(numeric_data.columns, rotation=45, ha="right")
    ax.set_xlabel("Variables")
    ax.set_ylabel("Normalized Values")
    ax.set_title("Parallel Coordinates Plot")

    # Add grid for better readability
    ax.grid(True, alpha=0.3)

    _save_fig(ax, output_path, "Parallel coordinates plot")

    return ax


def plot_radar_chart(
    data: pd.DataFrame, *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs: Any
) -> Axes:
    """Create a radar chart (spider chart).

    Args:
        data: DataFrame with numeric columns (one row per category/group)
        ax: Optional polar matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib plot() for every row (``linewidth``, ``alpha``,
            ``marker``, ``color``, ...); ``figsize`` sizes a new figure.

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If data format is incompatible, or ``ax`` is not a polar axes
    """
    validation.validate_type(data, pd.DataFrame, "data")

    # Select only numeric columns
    numeric_data = data.select_dtypes(include=[np.number])
    if numeric_data.shape[1] < 3:
        raise ValueError("Data must have at least 3 numeric columns for radar chart")

    if ax is None:
        fig = plt.figure(figsize=kwargs.pop("figsize", (8, 8)))
        ax = fig.add_subplot(111, polar=True)
    elif ax.name != "polar":
        raise ValueError("A radar chart needs a polar axes (create it with projection='polar')")

    # Get categories (column names)
    categories = list(numeric_data.columns)
    n_categories = len(categories)

    # Calculate angles for each category
    angles = [n / float(n_categories) * 2 * pi for n in range(n_categories)]
    angles += angles[:1]  # Close the circle

    # Plot each row as a separate line
    colors = plt.cm.tab10(np.linspace(0, 1, len(numeric_data)))

    for i, (idx, row) in enumerate(numeric_data.iterrows()):
        values = list(row.values)
        values += values[:1]  # Close the circle

        line_kwargs: dict[str, Any] = {
            "marker": "o",
            "linestyle": "-",
            "linewidth": 2,
            "label": str(idx),
            "color": colors[i % len(colors)],
            "alpha": 0.8,
            **kwargs,
        }
        ax.plot(angles, values, **line_kwargs)
        ax.fill(angles, values, alpha=0.25, color=line_kwargs["color"])

    # Set category labels
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories)
    ax.set_title("Radar Chart")

    if len(numeric_data) > 1:
        ax.legend(loc="upper right", bbox_to_anchor=(1.1, 1.1))

    _save_fig(ax, output_path, "Radar chart")

    return ax


def plot_3d_scatter(
    data: pd.DataFrame,
    *,
    x_col: str | None = None,
    y_col: str | None = None,
    z_col: str | None = None,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs: Any,
) -> Axes:
    """Create a 3D scatter plot.

    Args:
        data: DataFrame with at least 3 numeric columns
        x_col: Column name for x-axis (first numeric column if None)
        y_col: Column name for y-axis (second numeric column if None)
        z_col: Column name for z-axis (third numeric column if None)
        ax: Optional matplotlib axes to plot on. Creates new 3D axes if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib scatter. ``c`` may be the name of a numeric column
            of ``data`` that colours the points, or any colour or per-point values scatter accepts; ``figsize``
            sizes a new figure.

    Returns:
        matplotlib 3D Axes object

    Raises:
        ValueError: If data doesn't have enough numeric columns, or ``ax`` is not a 3D axes
    """
    validation.validate_type(data, pd.DataFrame, "data")

    # Select numeric columns
    numeric_data = data.select_dtypes(include=[np.number])
    if numeric_data.shape[1] < 3:
        raise ValueError("Data must have at least 3 numeric columns for 3D scatter plot")

    # Determine column names
    if x_col is None:
        x_col = numeric_data.columns[0]
    if y_col is None:
        y_col = numeric_data.columns[1]
    if z_col is None:
        z_col = numeric_data.columns[2]

    for col in [x_col, y_col, z_col]:
        if col not in numeric_data.columns:
            raise ValueError(f"Column '{col}' not found in numeric data")

    if ax is None:
        fig = plt.figure(figsize=kwargs.pop("figsize", (10, 8)))
        ax = fig.add_subplot(111, projection="3d")
    elif ax.name != "3d":
        raise ValueError("A 3D scatter plot needs a 3D axes (create it with projection='3d')")

    # ``c`` may name a column of the data (its values colour the points); anything else (a colour, or one value
    # per point) goes to scatter unchanged. A ``c`` that matched no column used to be dropped.
    colors = kwargs.pop("c", None)
    if isinstance(colors, str) and colors in data.columns:
        colors = data[colors]
    scatter_kwargs: dict[str, Any] = dict(kwargs)
    if colors is not None:
        scatter_kwargs["c"] = colors

    ax.scatter(data[x_col], data[y_col], data[z_col], **scatter_kwargs)

    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    cast("Any", ax).set_zlabel(z_col)
    ax.set_title("3D Scatter Plot")

    _save_fig(ax, output_path, "3D scatter plot")

    return ax
