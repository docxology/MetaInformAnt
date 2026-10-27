"""Regression tests: the static scatter fallback must handle a categorical colour column.

Without Plotly, ``interactive_scatter`` fell back to ``ax.scatter(..., c=data[color])``. matplotlib reads a
non-numeric ``c`` as colour specifications, so any column of category labels (``"ctrl"``, ``"treated"``)
raised ``ValueError: 'c' argument must be a color, a sequence of colors, or a sequence of numbers``.
Plotly colours such a column category by category with a legend; the fallback now does the same, while a
numeric column keeps its continuous colour scale.

The fallback is only reachable without Plotly, so these tests skip when Plotly is installed. All data are
real frames drawn onto real matplotlib Agg axes; no test doubles.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterator

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.axes import Axes
from matplotlib.colors import to_hex

from metainformant.visualization.config import palettes
from metainformant.visualization.dashboards import interactive
from metainformant.visualization.dashboards.interactive import interactive_scatter

matplotlib.use("Agg")

pytestmark = pytest.mark.skipif(interactive.HAS_PLOTLY, reason="the matplotlib fallback is only used without Plotly")

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


@pytest.fixture(autouse=True)
def _close_figures() -> Iterator[None]:
    yield
    plt.close("all")


@pytest.fixture
def frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "x": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            "y": [2.0, 1.0, 4.0, 3.0, 6.0, 5.0],
            "score": [0.1, 0.5, 0.9, 0.3, 0.7, 0.2],
            "group": ["ctrl", "treated", "ctrl", "dose", "treated", "ctrl"],
        }
    )


def _legend_labels(ax: Axes) -> list[str]:
    legend = ax.get_legend()
    assert legend is not None, "a categorical colour column needs a legend"
    return [text.get_text() for text in legend.get_texts()]


def _drawn_points(ax: Axes) -> int:
    return sum(len(collection.get_offsets()) for collection in ax.collections)


def test_numeric_colour_column_keeps_the_continuous_colour_scale(frame: pd.DataFrame) -> None:
    ax = interactive_scatter(frame, "x", "y", color="score")

    assert len(ax.collections) == 1
    assert np.asarray(ax.collections[0].get_array()).tolist() == frame["score"].tolist()
    assert ax.get_legend() is None


def test_string_colour_column_draws_one_coloured_series_per_category(frame: pd.DataFrame) -> None:
    ax = interactive_scatter(frame, "x", "y", color="group")

    categories = ["ctrl", "treated", "dose"]  # order of first appearance, as Plotly orders them
    assert _legend_labels(ax) == categories
    assert ax.get_legend().get_title().get_text() == "group"
    assert len(ax.collections) == len(categories)
    for collection, category in zip(ax.collections, categories):
        expected = frame.loc[frame["group"] == category, ["x", "y"]].to_numpy().tolist()
        assert np.asarray(collection.get_offsets()).tolist() == expected
    colours = [to_hex(collection.get_facecolor()[0]).lower() for collection in ax.collections]
    assert colours == [colour.lower() for colour in palettes.categorical(len(categories))]
    assert len(set(colours)) == len(categories)


def test_categorical_dtype_keeps_every_row_and_labels_missing_values(frame: pd.DataFrame) -> None:
    frame["group"] = pd.Categorical(
        ["ctrl", None, "ctrl", "dose", "treated", "ctrl"], categories=["dose", "treated", "ctrl", "never_used"]
    )

    ax = interactive_scatter(frame, "x", "y", color="group")

    assert _legend_labels(ax) == ["ctrl", "dose", "treated", "missing"]
    assert _drawn_points(ax) == len(frame)
    missing = ax.collections[-1]
    assert np.asarray(missing.get_offsets()).tolist() == [[2.0, 1.0]]
    assert to_hex(missing.get_facecolor()[0]).lower() == "#9e9e9e"


def test_more_categories_than_the_palette_holds_still_get_distinct_colours() -> None:
    n = 12
    frame = pd.DataFrame(
        {"x": np.arange(n, dtype=float), "y": np.arange(n, dtype=float), "group": [f"g{i}" for i in range(n)]}
    )

    ax = interactive_scatter(frame, "x", "y", color="group")

    assert len(_legend_labels(ax)) == n
    assert len({to_hex(collection.get_facecolor()[0]) for collection in ax.collections}) == n


def test_categorical_colour_column_can_be_exported(tmp_path: Path, frame: pd.DataFrame) -> None:
    output = tmp_path / "scatter.png"

    interactive_scatter(frame, "x", "y", color="group", output_path=output)

    assert output.read_bytes().startswith(PNG_SIGNATURE)


def test_empty_frame_with_a_categorical_colour_column_draws_an_empty_plot(frame: pd.DataFrame) -> None:
    ax = interactive_scatter(frame.iloc[0:0], "x", "y", color="group")

    assert _drawn_points(ax) == 0
    assert ax.get_legend() is None
