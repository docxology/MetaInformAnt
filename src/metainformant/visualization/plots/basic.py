"""Basic plotting functions for fundamental chart types.

This module provides core plotting functionality for line plots, scatter plots,
heatmaps, bar charts, pie charts, area plots, and step plots.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes

from metainformant.core.data import validation
from metainformant.core.io import paths
from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def lineplot(
    x: np.ndarray,
    y: np.ndarray | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs,
) -> Axes:
    """Create a line plot.

    Args:
        x: X-axis data
        y: Y-axis data. If None, x will be used as y and indices as x.
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib plot().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If data dimensions are incompatible
    """
    validation.validate_type(x, np.ndarray, "x")
    if y is not None:
        validation.validate_type(y, np.ndarray, "y")
        if len(x) != len(y):
            raise ValueError(f"x and y arrays must have same length: {len(x)} != {len(y)}")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (8, 6)))

    y_data = y if y is not None else x
    x_data = x if y is not None else np.arange(len(x))

    ax.plot(x_data, y_data, **kwargs)

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Line plot saved to {output_path}")

    return ax


def scatter_plot(
    x: np.ndarray, y: np.ndarray, *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs
) -> Axes:
    """Create a scatter plot.

    Args:
        x: X-axis data
        y: Y-axis data
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib scatter().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If data dimensions are incompatible
    """
    validation.validate_type(x, np.ndarray, "x")
    validation.validate_type(y, np.ndarray, "y")

    if len(x) != len(y):
        raise ValueError(f"x and y arrays must have same length: {len(x)} != {len(y)}")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (8, 6)))

    ax.scatter(x, y, **kwargs)

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Scatter plot saved to {output_path}")

    return ax


def heatmap(
    data: np.ndarray, *, cmap: str = "viridis", ax: Axes | None = None, output_path: str | Path | None = None, **kwargs
) -> Axes:
    """Create a 2D heatmap.

    Args:
        data: 2D array of data to plot
        cmap: Colormap name
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib imshow().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If data is not 2D
    """
    validation.validate_type(data, np.ndarray, "data")

    if data.ndim != 2:
        raise ValueError(f"Data must be 2D array, got {data.ndim}D")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (8, 6)))

    im = ax.imshow(data, cmap=cmap, **kwargs)
    plt.colorbar(im, ax=ax)

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Heatmap saved to {output_path}")

    return ax


def bar_plot(
    x: np.ndarray | Sequence[str],
    height: np.ndarray,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs,
) -> Axes:
    """Create a bar plot.

    Args:
        x: X-axis positions or labels
        height: Heights of bars
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib bar().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If data dimensions are incompatible
    """
    validation.validate_type(height, np.ndarray, "height")

    if not isinstance(x, (np.ndarray, list)):
        raise ValueError("x must be array-like")

    if len(x) != len(height):
        raise ValueError(f"x and height arrays must have same length: {len(x)} != {len(height)}")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (8, 6)))

    ax.bar(x, height, **kwargs)

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Bar plot saved to {output_path}")

    return ax


def pie_chart(
    sizes: np.ndarray,
    labels: list[str] | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs,
) -> Axes:
    """Create a pie chart.

    Args:
        sizes: Sizes of pie slices
        labels: Labels for each slice
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib pie().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If data dimensions are incompatible
    """
    validation.validate_type(sizes, np.ndarray, "sizes")

    if labels is not None:
        validation.validate_type(labels, list, "labels")
        if len(labels) != len(sizes):
            raise ValueError(f"labels and sizes arrays must have same length: {len(labels)} != {len(sizes)}")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (8, 6)))

    ax.pie(sizes, labels=labels, **kwargs)
    ax.axis("equal")  # Equal aspect ratio ensures that pie is drawn as a circle

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Pie chart saved to {output_path}")

    return ax


def area_plot(
    x: np.ndarray, y: np.ndarray, *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs
) -> Axes:
    """Create a filled area plot.

    Args:
        x: X-axis data
        y: Y-axis data
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib fill_between().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If data dimensions are incompatible
    """
    validation.validate_type(x, np.ndarray, "x")
    validation.validate_type(y, np.ndarray, "y")

    if len(x) != len(y):
        raise ValueError(f"x and y arrays must have same length: {len(x)} != {len(y)}")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (8, 6)))

    ax.fill_between(x, y, **kwargs)

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Area plot saved to {output_path}")

    return ax


def step_plot(
    x: np.ndarray, y: np.ndarray, *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs
) -> Axes:
    """Create a step plot.

    Args:
        x: X-axis data
        y: Y-axis data
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib step().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If data dimensions are incompatible
    """
    validation.validate_type(x, np.ndarray, "x")
    validation.validate_type(y, np.ndarray, "y")

    if len(x) != len(y):
        raise ValueError(f"x and y arrays must have same length: {len(x)} != {len(y)}")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (8, 6)))

    ax.step(x, y, **kwargs)

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Step plot saved to {output_path}")

    return ax
