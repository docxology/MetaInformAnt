"""Basic plotting functions for simple visualizations.

This module provides fundamental plotting functions including line plots,
scatter plots, bar charts, pie charts, area plots, and basic heatmaps.
"""

from __future__ import annotations

from typing import Sequence

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

# Use non-interactive backend by default for tests/headless
matplotlib.use("Agg", force=True)

try:
    import seaborn as sns

    SEABORN_AVAILABLE = True
except ImportError:
    SEABORN_AVAILABLE = False
    sns = None


def lineplot(
    x: Sequence[float] | None,
    y: Sequence[float],
    *,
    label: str | None = None,
    ax: plt.Axes | None = None,
    style: str = "-",
    color: str | None = None,
) -> plt.Axes:
    """Simple line plot wrapper returning the Axes.

    If x is None, indices are used.

    Args:
        x: X-axis values (if None, uses indices)
        y: Y-axis values
        label: Label for legend
        ax: Matplotlib axes (creates new if None)
        style: Line style (e.g., '-', '--', '-.')
        color: Line color

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import lineplot
        >>> ax = lineplot(None, [1, 4, 2, 8, 5])
        >>> ax.set_xlabel("Index")
        >>> ax.set_ylabel("Value")
    """
    if ax is None:
        _, ax = plt.subplots()
    if x is None:
        x = list(range(len(y)))
    ax.plot(x, y, style, label=label, color=color)
    if label:
        ax.legend()
    return ax


def scatter_plot(
    x: Sequence[float],
    y: Sequence[float],
    *,
    ax: plt.Axes | None = None,
    color: str | Sequence[str] | None = None,
    size: float | Sequence[float] = 20,
    alpha: float = 0.7,
    xlabel: str | None = None,
    ylabel: str | None = None,
    title: str | None = None,
) -> plt.Axes:
    """Create a scatter plot.

    Args:
        x: X-axis values
        y: Y-axis values
        ax: Matplotlib axes (creates new if None)
        color: Point colors
        size: Point sizes
        alpha: Point transparency
        xlabel: X-axis label
        ylabel: Y-axis label
        title: Plot title

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import scatter_plot
        >>> ax = scatter_plot([1, 2, 3], [4, 5, 6], xlabel="X", ylabel="Y")
    """
    if ax is None:
        _, ax = plt.subplots()

    ax.scatter(x, y, c=color, s=size, alpha=alpha)

    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)

    return ax


def heatmap(
    data: Sequence[Sequence[float]] | pd.DataFrame,
    *,
    cmap: str = "viridis",
    cbar: bool = True,
    ax: plt.Axes | None = None,
    annot: bool = False,
) -> plt.Axes:
    """Create a heatmap using seaborn.

    Args:
        data: 2D array or DataFrame to visualize
        cmap: Colormap name (default: 'viridis')
        cbar: Whether to show colorbar (default: True)
        ax: Matplotlib axes to plot on (creates new if None)
        annot: Whether to annotate cells with values (default: False)

    Returns:
        Matplotlib axes object containing the heatmap

    Returns:
        Matplotlib axes object

    Example:
        >>> import numpy as np
        >>> from metainformant.visualization import heatmap
        >>> data = np.random.random((10, 10))
        >>> ax = heatmap(data, annot=True)
    """
    if ax is None:
        _, ax = plt.subplots()

    if SEABORN_AVAILABLE:
        sns.heatmap(data, cmap=cmap, cbar=cbar, ax=ax, annot=annot)
    else:
        # Fallback to matplotlib
        import numpy as np

        data_array = np.array(data)
        im = ax.imshow(data_array, cmap=cmap)
        if cbar:
            plt.colorbar(im, ax=ax)
        if annot:
            for i in range(data_array.shape[0]):
                for j in range(data_array.shape[1]):
                    ax.text(j, i, f"{data_array[i, j]:.2f}", ha="center", va="center")
    return ax


def bar_plot(
    x: Sequence[str] | Sequence[float],
    y: Sequence[float],
    *,
    ax: plt.Axes | None = None,
    color: str | Sequence[str] | None = None,
    alpha: float = 0.7,
    xlabel: str | None = None,
    ylabel: str | None = None,
    title: str | None = None,
    horizontal: bool = False,
) -> plt.Axes:
    """Create a bar plot.

    Args:
        x: Categories or x-axis values
        y: Bar heights
        ax: Matplotlib axes (creates new if None)
        color: Bar color(s)
        alpha: Bar transparency
        xlabel: X-axis label
        ylabel: Y-axis label
        title: Plot title
        horizontal: Whether to create horizontal bars

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import bar_plot
        >>> ax = bar_plot(["A", "B", "C"], [10, 20, 15])
    """
    if ax is None:
        _, ax = plt.subplots()

    if horizontal:
        ax.barh(x, y, color=color, alpha=alpha)
    else:
        ax.bar(x, y, color=color, alpha=alpha)

    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)

    return ax


def pie_chart(
    sizes: Sequence[float],
    labels: Sequence[str] | None = None,
    *,
    ax: plt.Axes | None = None,
    colors: Sequence[str] | None = None,
    autopct: str = "%1.1f%%",
    title: str | None = None,
) -> plt.Axes:
    """Create a pie chart.

    Args:
        sizes: Values for each slice
        labels: Labels for each slice
        ax: Matplotlib axes (creates new if None)
        colors: Colors for slices
        autopct: Format string for percentage labels
        title: Plot title

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import pie_chart
        >>> ax = pie_chart([30, 25, 45], ["A", "B", "C"])
    """
    if ax is None:
        _, ax = plt.subplots()

    ax.pie(sizes, labels=labels, colors=colors, autopct=autopct, startangle=90)
    ax.axis("equal")

    if title:
        ax.set_title(title)

    return ax


def area_plot(
    x: Sequence[float],
    y: Sequence[float],
    *,
    ax: plt.Axes | None = None,
    color: str | None = None,
    alpha: float = 0.5,
    xlabel: str | None = None,
    ylabel: str | None = None,
    title: str | None = None,
) -> plt.Axes:
    """Create an area plot (filled line plot).

    Args:
        x: X-axis values
        y: Y-axis values
        ax: Matplotlib axes (creates new if None)
        color: Fill color
        alpha: Fill transparency
        xlabel: X-axis label
        ylabel: Y-axis label
        title: Plot title

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import area_plot
        >>> ax = area_plot([1, 2, 3, 4], [1, 4, 2, 3])
    """
    if ax is None:
        _, ax = plt.subplots()

    ax.fill_between(x, y, alpha=alpha, color=color)

    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)

    return ax


def step_plot(
    x: Sequence[float],
    y: Sequence[float],
    *,
    ax: plt.Axes | None = None,
    where: str = "pre",
    label: str | None = None,
    color: str | None = None,
    xlabel: str | None = None,
    ylabel: str | None = None,
    title: str | None = None,
) -> plt.Axes:
    """Create a step plot.

    Args:
        x: X-axis values
        y: Y-axis values
        ax: Matplotlib axes (creates new if None)
        where: Step position ('pre', 'post', 'mid')
        label: Label for legend
        color: Line color
        xlabel: X-axis label
        ylabel: Y-axis label
        title: Plot title

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import step_plot
        >>> ax = step_plot([1, 2, 3, 4], [1, 4, 2, 3])
    """
    if ax is None:
        _, ax = plt.subplots()

    ax.step(x, y, where=where, label=label, color=color)

    if label:
        ax.legend()
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)

    return ax

