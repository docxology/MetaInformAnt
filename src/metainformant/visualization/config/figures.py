"""Figure-resolution helper shared by the plot functions that export an ``Axes``.

Plot functions accept a caller-supplied ``ax`` and, when ``output_path`` is
given, save the figure that axes is drawn on. ``Axes.figure`` is a
:class:`~matplotlib.figure.SubFigure` (which has no ``savefig``) when the axes
was created through ``fig.subfigures(...)``; :func:`root_figure` always
returns the savable root :class:`~matplotlib.figure.Figure`.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt

if TYPE_CHECKING:
    from matplotlib.artist import Artist
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure


def panel_grid(
    ax: Axes | None,
    n_panels: int,
    *,
    max_cols: int,
    figsize: tuple[float, float],
    sharex: bool = False,
) -> list[Axes]:
    """Return one axes per panel for a multi-panel plot, in row-major order, honouring a caller-supplied ``ax``.

    The panels fill a compact grid of at most ``max_cols`` columns, so a plot with fewer panels than the grid
    has cells is not padded with empty axes.

    * ``ax=None``: a new figure of ``figsize`` is created with exactly ``n_panels`` axes.
    * ``ax`` given: the panels are laid out inside the grid slot ``ax`` occupies (a sub-gridspec of its
      subplot spec), on the figure or subfigure ``ax`` belongs to, and ``ax`` itself is hidden because it is only
      the frame that holds them. Nothing is created on the pyplot state or stacked over ``ax``.

    Args:
        ax: Axes that should hold the panels, or None to create a figure.
        n_panels: Number of panels (at least 1).
        max_cols: Maximum number of columns.
        figsize: Size of the new figure when ``ax`` is None.
        sharex: Share the x axis between the panels.

    Raises:
        ValueError: If ``n_panels`` is below 1, or ``ax`` was not created with ``add_subplot``/``subplots`` and so
            has no grid slot to hold the panels.
    """
    if n_panels < 1:
        raise ValueError("A multi-panel plot needs at least one panel")
    ncols = min(max_cols, n_panels)
    nrows = math.ceil(n_panels / ncols)

    if ax is None:
        fig, grid = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False, sharex=sharex)
        panels = list(grid.flatten())
        for surplus in panels[n_panels:]:
            fig.delaxes(surplus)
        return panels[:n_panels]

    spec = ax.get_subplotspec()
    if spec is None:
        raise ValueError(
            "ax must come from add_subplot, subplots or subfigure.subplots so that the panels can be laid out "
            "inside its grid slot; pass ax=None to create a new figure"
        )
    subgrid = spec.subgridspec(nrows, ncols)
    panels = []
    for index in range(n_panels):
        panels.append(
            ax.figure.add_subplot(
                subgrid[index // ncols, index % ncols], sharex=panels[0] if sharex and panels else None
            )
        )
    ax.set_axis_off()
    return panels


def root_figure(artist: Artist) -> Figure:
    """Return the root :class:`~matplotlib.figure.Figure` that *artist* is drawn on.

    For an artist on an ordinary figure this is ``artist.figure`` itself. For an
    artist inside one or more nested subfigures it is the top-level figure that
    owns them, so the result can always be handed to
    :func:`~metainformant.visualization.config.conventions.save_figure_deterministic`
    or ``Figure.savefig``.

    Args:
        artist: Any matplotlib artist, typically the ``Axes`` a plot was drawn on.

    Returns:
        The root figure of *artist*.
    """
    return artist.figure.figure
