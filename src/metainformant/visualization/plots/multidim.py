"""Multi-dimensional data visualization functions.

This module provides specialized plotting functions for visualizing
multi-dimensional datasets including pairwise relationships, parallel coordinates,
radar charts, and 3D scatter plots.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any
from math import pi

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.axes import Axes

from metainformant.core import logging, paths, validation

logger = logging.get_logger(__name__)

# Optional imports with graceful fallbacks
try:
    import seaborn as sns

    HAS_SEABORN = True
except ImportError:
    sns = None
    HAS_SEABORN = False


def plot_pairwise_relationships(
    data: pd.DataFrame, *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs
) -> Axes:
    """Create a pairwise relationships plot (scatter plot matrix).

    Args:
        data: DataFrame with numeric columns for pairwise plotting
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to seaborn pairplot.

    Returns:
        matplotlib Axes object (main axes from the grid)

    Raises:
        ValueError: If data doesn't have enough numeric columns
        ImportError: If seaborn is not available
    """
    validation.validate_type(data, pd.DataFrame, "data")

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
        g.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Pairwise relationships plot saved to {output_path}")

    # Return the first axes from the grid
    return g.fig.axes[0] if g.fig.axes else None


def plot_parallel_coordinates(
    data: pd.DataFrame, *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs
) -> Axes:
    """Create a parallel coordinates plot.

    Args:
        data: DataFrame with numeric columns for parallel coordinates
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to pandas parallel_coordinates.

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

    # Check if there's a categorical column for coloring
    categorical_cols = data.select_dtypes(include=["object", "category"]).columns
    color_col = kwargs.pop("color", categorical_cols[0] if len(categorical_cols) > 0 else None)

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (12, 6)))

    # Normalize data to 0-1 range for better visualization
    normalized_data = (numeric_data - numeric_data.min()) / (numeric_data.max() - numeric_data.min())

    # If color column provided, add it back
    if color_col and color_col in data.columns:
        plot_data = normalized_data.copy()
        plot_data[color_col] = data[color_col]
    else:
        plot_data = normalized_data

    # Create parallel coordinates plot
    if HAS_SEABORN and color_col:
        # Use seaborn for better coloring
        palette = kwargs.pop("palette", "husl")
        for category in plot_data[color_col].unique():
            subset = plot_data[plot_data[color_col] == category]
            ax.plot(subset.drop(columns=[color_col]).T.values, alpha=0.7, label=str(category))
        ax.legend()
    else:
        # Basic matplotlib version
        for i in range(len(plot_data)):
            ax.plot(plot_data.iloc[i].values, alpha=0.7)

    # Set x-ticks to column names
    ax.set_xticks(range(len(numeric_data.columns)))
    ax.set_xticklabels(numeric_data.columns, rotation=45, ha="right")
    ax.set_xlabel("Variables")
    ax.set_ylabel("Normalized Values")
    ax.set_title("Parallel Coordinates Plot")

    # Add grid for better readability
    ax.grid(True, alpha=0.3)

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Parallel coordinates plot saved to {output_path}")

    return ax


def plot_radar_chart(
    data: pd.DataFrame, *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs
) -> Axes:
    """Create a radar chart (spider chart).

    Args:
        data: DataFrame with numeric columns (one row per category/group)
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to plotting.

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If data format is incompatible
    """
    validation.validate_type(data, pd.DataFrame, "data")

    # Select only numeric columns
    numeric_data = data.select_dtypes(include=[np.number])
    if numeric_data.shape[1] < 3:
        raise ValueError("Data must have at least 3 numeric columns for radar chart")

    if ax is None:
        fig = plt.figure(figsize=kwargs.pop("figsize", (8, 8)))
        ax = fig.add_subplot(111, polar=True)

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

        ax.plot(angles, values, "o-", linewidth=2, label=str(idx), color=colors[i % len(colors)], alpha=0.8)
        ax.fill(angles, values, alpha=0.25, color=colors[i % len(colors)])

    # Set category labels
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories)
    ax.set_title("Radar Chart")

    if len(numeric_data) > 1:
        ax.legend(loc="upper right", bbox_to_anchor=(1.1, 1.1))

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Radar chart saved to {output_path}")

    return ax


def plot_3d_scatter(
    data: pd.DataFrame,
    *,
    x_col: str | None = None,
    y_col: str | None = None,
    z_col: str | None = None,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs,
) -> Axes:
    """Create a 3D scatter plot.

    Args:
        data: DataFrame with at least 3 numeric columns
        x_col: Column name for x-axis (first numeric column if None)
        y_col: Column name for y-axis (second numeric column if None)
        z_col: Column name for z-axis (third numeric column if None)
        ax: Optional matplotlib axes to plot on. Creates new 3D axes if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib scatter.

    Returns:
        matplotlib 3D Axes object

    Raises:
        ValueError: If data doesn't have enough numeric columns
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

    # Check if we have a color column
    color_col = kwargs.pop("c", None)
    if color_col and color_col in data.columns:
        c = data[color_col]
    else:
        c = None

    ax.scatter(data[x_col], data[y_col], data[z_col], c=c, **kwargs)

    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.set_zlabel(z_col)
    ax.set_title("3D Scatter Plot")

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"3D scatter plot saved to {output_path}")

    return ax
