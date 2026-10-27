"""Regression tests: phenotype plots must return a real, savable ``Figure``.

``Axes.figure`` is a :class:`matplotlib.figure.SubFigure` when the caller passes
an axes created through ``fig.subfigures(...)``. The plot functions are declared
to return ``plt.Figure`` but returned ``ax.figure`` unchanged, so callers that
supplied a subfigure axes received a ``SubFigure`` (no ``savefig``) back. They
now return the root ``Figure``.

All tests use real matplotlib Agg figures and real PNG files in ``tmp_path``;
no test doubles.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import pytest
from matplotlib.axes import Axes
from matplotlib.colors import to_rgba
from matplotlib.figure import Figure

from metainformant.phenotype.visualization.plots import (
    plot_boxplot_with_swarm,
    plot_categorical_proportions,
    plot_correlation_heatmap,
    plot_violin,
)

matplotlib.use("Agg")

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"

PlotCall = Callable[[pd.DataFrame, Axes], Any]


@pytest.fixture
def sample_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "group": ["A", "A", "A", "B", "B", "B"],
            "category": ["X", "Y", "X", "Y", "Y", "X"],
            "value": [1.1, 2.2, 1.5, 5.5, 6.1, 5.8],
            "value2": [10, 20, 15, 50, 60, 55],
        }
    )


_PLOTS: list[tuple[str, PlotCall]] = [
    ("plot_boxplot_with_swarm", lambda df, ax: plot_boxplot_with_swarm(df, "group", "value", ax=ax)),
    ("plot_violin", lambda df, ax: plot_violin(df, "group", "value", ax=ax)),
    ("plot_categorical_proportions", lambda df, ax: plot_categorical_proportions(df, "group", "category", ax=ax)),
    # With a caller-supplied axes the dendrogram layout is skipped and the
    # simple heatmap branch (which returns ``ax.figure``) is used.
    ("plot_correlation_heatmap", lambda df, ax: plot_correlation_heatmap(df[["value", "value2"]].corr(), ax=ax)),
]


@pytest.mark.parametrize("plot", [pytest.param(p, id=i) for i, p in _PLOTS])
def test_subfigure_axes_return_the_savable_root_figure(tmp_path: Path, sample_df: pd.DataFrame, plot: PlotCall) -> None:
    figure = plt.figure(figsize=(10, 5))
    left, _right = figure.subfigures(1, 2)
    ax = left.subplots()
    assert type(ax.figure).__name__ == "SubFigure"

    result = plot(sample_df, ax)

    assert isinstance(result, Figure)
    assert result is figure
    output_path = tmp_path / "phenotype_plot.png"
    result.savefig(output_path)
    assert output_path.read_bytes()[:8] == PNG_SIGNATURE
    plt.close(figure)


@pytest.mark.parametrize("plot", [pytest.param(p, id=i) for i, p in _PLOTS])
def test_regular_axes_still_return_their_own_figure(sample_df: pd.DataFrame, plot: PlotCall) -> None:
    figure, ax = plt.subplots(figsize=(8, 5))

    result = plot(sample_df, ax)

    assert result is figure
    plt.close(figure)


def test_categorical_proportions_use_the_tab20_listed_colormap(sample_df: pd.DataFrame) -> None:
    """The stacked bars are coloured from ``tab20`` (a ListedColormap with ``.colors``)."""
    figure, ax = plt.subplots()

    plot_categorical_proportions(sample_df, "group", "category", ax=ax)

    tab20 = plt.get_cmap("tab20")
    assert hasattr(tab20, "colors")
    colors_by_category = {container.get_label(): container.patches[0].get_facecolor() for container in ax.containers}
    assert colors_by_category["X"] == pytest.approx(to_rgba(tab20.colors[0]))
    assert colors_by_category["Y"] == pytest.approx(to_rgba(tab20.colors[1]))
    plt.close(figure)
