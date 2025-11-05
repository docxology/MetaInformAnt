"""Layout utilities for multi-panel figures.

This module provides utilities for creating multi-panel figures, managing
subplot grids, and composing complex figure layouts.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np


def create_subplot_grid(
    n_plots: int,
    ncols: int = 3,
    *,
    figsize: tuple[float, float] | None = None,
    sharex: bool = False,
    sharey: bool = False,
) -> tuple[plt.Figure, np.ndarray]:
    """Create a subplot grid for multiple plots.

    Args:
        n_plots: Number of plots
        ncols: Number of columns
        figsize: Figure size tuple (if None, auto-calculates)
        sharex: Whether to share x-axis
        sharey: Whether to share y-axis

    Returns:
        Tuple of (figure, axes array)

    Example:
        >>> from metainformant.visualization.layout import create_subplot_grid
        >>> fig, axes = create_subplot_grid(6, ncols=3)
    """
    nrows = (n_plots + ncols - 1) // ncols

    if figsize is None:
        figsize = (4 * ncols, 3 * nrows)

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize,
                            sharex=sharex, sharey=sharey)

    if nrows == 1 and ncols == 1:
        axes = np.array([[axes]])
    elif nrows == 1:
        axes = axes.reshape(1, -1)
    elif ncols == 1:
        axes = axes.reshape(-1, 1)

    return fig, axes


def create_multi_panel(
    n_panels: int,
    layout: str = 'grid',
    *,
    figsize: tuple[float, float] | None = None,
    **kwargs
) -> tuple[plt.Figure, list]:
    """Create a multi-panel figure layout.

    Args:
        n_panels: Number of panels
        layout: Layout type ('grid', 'vertical', 'horizontal')
        figsize: Figure size tuple
        **kwargs: Additional arguments for subplots

    Returns:
        Tuple of (figure, list of axes)

    Example:
        >>> from metainformant.visualization.layout import create_multi_panel
        >>> fig, axes = create_multi_panel(4, layout='grid')
    """
    if figsize is None:
        if layout == 'vertical':
            figsize = (6, 3 * n_panels)
        elif layout == 'horizontal':
            figsize = (6 * n_panels, 3)
        else:
            ncols = int(np.ceil(np.sqrt(n_panels)))
            nrows = int(np.ceil(n_panels / ncols))
            figsize = (4 * ncols, 3 * nrows)

    if layout == 'vertical':
        fig, axes = plt.subplots(n_panels, 1, figsize=figsize, **kwargs)
        if n_panels == 1:
            axes = [axes]
    elif layout == 'horizontal':
        fig, axes = plt.subplots(1, n_panels, figsize=figsize, **kwargs)
        if n_panels == 1:
            axes = [axes]
    else:  # grid
        ncols = int(np.ceil(np.sqrt(n_panels)))
        nrows = int(np.ceil(n_panels / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize, **kwargs)
        if nrows == 1 and ncols == 1:
            axes = [axes]
        elif nrows == 1:
            axes = axes.flatten().tolist()
        elif ncols == 1:
            axes = axes.flatten().tolist()
        else:
            axes = axes.flatten().tolist()

    return fig, axes


def add_shared_axis_labels(
    fig: plt.Figure,
    xlabel: str | None = None,
    ylabel: str | None = None,
):
    """Add shared axis labels to a multi-panel figure.

    Args:
        fig: Matplotlib figure
        xlabel: Shared x-axis label
        ylabel: Shared y-axis label

    Example:
        >>> from metainformant.visualization.layout import add_shared_axis_labels
        >>> fig, axes = create_multi_panel(4)
        >>> add_shared_axis_labels(fig, xlabel='Time', ylabel='Value')
    """
    if xlabel:
        fig.text(0.5, 0.02, xlabel, ha='center', va='bottom')
    if ylabel:
        fig.text(0.02, 0.5, ylabel, ha='left', va='center', rotation='vertical')


def hide_unused_subplots(
    axes: np.ndarray,
    n_used: int,
):
    """Hide unused subplots in a grid.

    Args:
        axes: Array of axes
        n_used: Number of subplots actually used

    Example:
        >>> from metainformant.visualization.layout import hide_unused_subplots
        >>> fig, axes = create_subplot_grid(9, ncols=3)
        >>> # Use first 6 subplots
        >>> hide_unused_subplots(axes, 6)
    """
    axes_flat = axes.flatten()
    for i in range(n_used, len(axes_flat)):
        axes_flat[i].set_visible(False)


def adjust_spacing(
    fig: plt.Figure,
    *,
    hspace: float | None = None,
    wspace: float | None = None,
    top: float | None = None,
    bottom: float | None = None,
    left: float | None = None,
    right: float | None = None,
):
    """Adjust spacing in a multi-panel figure.

    Args:
        fig: Matplotlib figure
        hspace: Height spacing between subplots
        wspace: Width spacing between subplots
        top: Top margin
        bottom: Bottom margin
        left: Left margin
        right: Right margin

    Example:
        >>> from metainformant.visualization.layout import adjust_spacing
        >>> fig, axes = create_multi_panel(4)
        >>> adjust_spacing(fig, hspace=0.3, wspace=0.2)
    """
    plt.subplots_adjust(
        hspace=hspace,
        wspace=wspace,
        top=top,
        bottom=bottom,
        left=left,
        right=right,
    )

