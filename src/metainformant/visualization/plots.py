from __future__ import annotations

from typing import Iterable, Sequence

import matplotlib

# Use non-interactive backend by default for tests/headless
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

try:
    import seaborn as sns  # noqa: E402

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


def pairplot_dataframe(df: pd.DataFrame, *, hue: str | None = None):
    """Pairplot for a tidy DataFrame, returns seaborn PairGrid or matplotlib figure."""
    if SEABORN_AVAILABLE:
        grid = sns.pairplot(df, hue=hue)
        return grid
    else:
        # Fallback: simple scatter plot matrix with matplotlib
        import numpy as np

        numeric_cols = df.select_dtypes(include=[np.number]).columns
        n_cols = len(numeric_cols)

        fig, axes = plt.subplots(n_cols, n_cols, figsize=(n_cols * 3, n_cols * 3))
        if n_cols == 1:
            axes = np.array([[axes]])
        elif n_cols == 2:
            axes = axes.reshape(2, 2)

        for i, col1 in enumerate(numeric_cols):
            for j, col2 in enumerate(numeric_cols):
                ax = axes[i, j]
                if i == j:
                    # Diagonal: histogram
                    ax.hist(df[col1], alpha=0.7)
                    ax.set_xlabel(col1)
                else:
                    # Off-diagonal: scatter plot
                    ax.scatter(df[col2], df[col1], alpha=0.6)
                    ax.set_xlabel(col2)
                    ax.set_ylabel(col1)

        plt.tight_layout()
        return fig
