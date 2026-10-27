"""Figure-resolution helper shared by the plot functions that export an ``Axes``.

Plot functions accept a caller-supplied ``ax`` and, when ``output_path`` is
given, save the figure that axes is drawn on. ``Axes.figure`` is a
:class:`~matplotlib.figure.SubFigure` (which has no ``savefig``) when the axes
was created through ``fig.subfigures(...)``; :func:`root_figure` always
returns the savable root :class:`~matplotlib.figure.Figure`.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from matplotlib.artist import Artist
    from matplotlib.figure import Figure


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
