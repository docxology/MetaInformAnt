from __future__ import annotations

from typing import Sequence

import matplotlib

# Non-interactive for tests/headless
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.animation import FuncAnimation  # noqa: E402


def animate_time_series(values: Sequence[float], *, interval_ms: int = 200):
    """Create a simple time-series animation.

    Returns (fig, animation).
    """
    fig, ax = plt.subplots()
    line, = ax.plot([], [], lw=2)
    ax.set_xlim(0, max(1, len(values)))
    ymin = min(values) if values else 0.0
    ymax = max(values) if values else 1.0
    if ymin == ymax:
        ymax = ymin + 1.0
    ax.set_ylim(ymin, ymax)

    def init():
        line.set_data([], [])
        return (line,)

    def update(frame: int):
        x = list(range(frame + 1))
        y = list(values[: frame + 1])
        line.set_data(x, y)
        return (line,)

    anim = FuncAnimation(fig, update, init_func=init, frames=len(values), interval=interval_ms, blit=True)
    return fig, anim

from __future__ import annotations

from typing import Sequence, Iterable, Callable

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


def animate_time_series(
    series: Sequence[float],
    *,
    interval_ms: int = 100,
    repeat: bool = False,
    init_points: int = 1,
    ax: plt.Axes | None = None,
) -> tuple[plt.Figure, FuncAnimation]:
    """Animate a time series as a growing line.

    Returns the Figure and FuncAnimation for saving/displaying.
    """
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    x = list(range(len(series)))
    (line,) = ax.plot([], [], "-b")
    ax.set_xlim(0, max(1, len(series) - 1))
    ymin = min(series) if series else 0.0
    ymax = max(series) if series else 1.0
    if ymin == ymax:
        ymin -= 1.0
        ymax += 1.0
    ax.set_ylim(ymin, ymax)

    def init():
        line.set_data(x[:init_points], series[:init_points])
        return (line,)

    def update(frame_idx: int):
        i = min(frame_idx + 1, len(series))
        line.set_data(x[:i], series[:i])
        return (line,)

    anim = FuncAnimation(
        fig,
        update,
        init_func=init,
        frames=len(series),
        interval=interval_ms,
        blit=True,
        repeat=repeat,
    )
    return fig, anim




