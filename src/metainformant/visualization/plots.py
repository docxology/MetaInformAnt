from __future__ import annotations

from typing import Iterable, Sequence

import matplotlib

# Use non-interactive backend by default for tests/headless
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt  # noqa: E402
import seaborn as sns  # noqa: E402


def lineplot(x: Iterable[float] | None, y: Sequence[float], *, label: str | None = None):
    """Simple line plot; if x is None, uses index.

    Returns the matplotlib Axes instance.
    """
    ax = plt.gca()
    if x is None:
        ax.plot(list(range(len(y))), y, label=label)
    else:
        ax.plot(list(x), y, label=label)
    if label:
        ax.legend()
    return ax


def heatmap(matrix: Sequence[Sequence[float]], *, cmap: str = "viridis"):
    """2D heatmap using seaborn. Returns Axes.
    """
    ax = sns.heatmap(matrix, cmap=cmap)
    return ax

from __future__ import annotations

from typing import Iterable, Sequence

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


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
    """
    if ax is None:
        _, ax = plt.subplots()
    if x is None:
        x = list(range(len(y)))
    ax.plot(x, y, style, label=label, color=color)
    if label:
        ax.legend()
    return ax


def heatmap(
    data: Sequence[Sequence[float]] | pd.DataFrame,
    *,
    cmap: str = "viridis",
    cbar: bool = True,
    ax: plt.Axes | None = None,
    annot: bool = False,
) -> plt.Axes:
    """Heatmap using seaborn; accepts 2D sequences or DataFrame."""
    if ax is None:
        _, ax = plt.subplots()
    sns.heatmap(data, cmap=cmap, cbar=cbar, ax=ax, annot=annot)
    return ax


def pairplot_dataframe(df: pd.DataFrame, *, hue: str | None = None):
    """Pairplot for a tidy DataFrame, returns seaborn PairGrid."""
    grid = sns.pairplot(df, hue=hue)
    return grid




