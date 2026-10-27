"""Regression tests: deterministic figure export gives every format only the metadata it accepts.

``save_figure_deterministic`` built one PNG-style mapping (``Software``, ``Creator``,
``Creation Time``, ``Date``) and handed it to every matplotlib backend, but each backend
validates its own metadata keys. The SVG writer raised
``ValueError: Unknown metadata key(s) passed to SVG writer: Software,Creation Time``, so no
SVG could be saved at all, and the Pillow-backed raster writers (JPEG, TIFF, WebP) refuse
metadata outright (``ValueError: metadata not supported for format 'jpeg'``).

Accepting the keys is not enough to make SVG reproducible: the writer stamps ``Date`` with the
current time and derives the ids of clip paths and tick markers from a fresh ``uuid4`` unless
``svg.hashsalt`` is fixed, so two saves of one figure differ in the date and in every id.

All tests draw real matplotlib Agg figures, write real files in ``tmp_path`` and compare the saved
bytes; no test doubles are involved.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterator

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.figure import Figure

from metainformant.visualization.config.conventions import DETERMINISTIC_RCPARAMS, save_figure_deterministic
from metainformant.visualization.plots import basic

matplotlib.use("Agg")

RASTER_SIGNATURES = {
    "jpg": (b"\xff\xd8\xff",),
    "tiff": (b"II*\x00", b"MM\x00*"),
}


@pytest.fixture
def figure() -> Iterator[Figure]:
    """A figure with markers, ticks, a legend and a title: every one of them draws salted SVG ids."""
    fig, ax = plt.subplots()
    ax.plot([0.0, 1.0, 2.0, 3.0], [1.0, 3.0, 2.0, 5.0], marker="o", label="series")
    ax.set_title("deterministic export")
    ax.legend()
    yield fig
    plt.close(fig)


def test_svg_save_no_longer_raises(tmp_path: Path, figure: Figure) -> None:
    path = tmp_path / "figure.svg"

    save_figure_deterministic(figure, path)

    text = path.read_text(encoding="utf-8")
    assert text.startswith("<?xml")
    assert "<svg" in text
    assert text.rstrip().endswith("</svg>")


@pytest.mark.parametrize("suffix", ["png", "pdf", "svg"])
def test_two_saves_are_byte_identical_whatever_the_clock_says(
    tmp_path: Path, figure: Figure, monkeypatch: pytest.MonkeyPatch, suffix: str
) -> None:
    """SOURCE_DATE_EPOCH is matplotlib's own clock override: two values must leave no trace in the file."""
    first, second = tmp_path / f"first.{suffix}", tmp_path / f"second.{suffix}"
    monkeypatch.setenv("SOURCE_DATE_EPOCH", "100")
    save_figure_deterministic(figure, first)
    monkeypatch.setenv("SOURCE_DATE_EPOCH", "200")
    save_figure_deterministic(figure, second)

    assert first.stat().st_size > 0
    assert first.read_bytes() == second.read_bytes()


def test_svg_carries_neither_a_creator_nor_a_date(tmp_path: Path, figure: Figure) -> None:
    path = tmp_path / "figure.svg"
    save_figure_deterministic(figure, path)

    text = path.read_text(encoding="utf-8")

    assert "<dc:date>" not in text
    assert "<dc:creator>" not in text
    assert "matplotlib.org" not in text


def test_pdf_carries_neither_a_producer_nor_a_date(tmp_path: Path, figure: Figure) -> None:
    path = tmp_path / "figure.pdf"
    save_figure_deterministic(figure, path)

    data = path.read_bytes()

    assert b"/CreationDate" not in data
    assert b"/ModDate" not in data
    assert b"Matplotlib" not in data


def test_png_carries_no_producer_text(tmp_path: Path, figure: Figure) -> None:
    path = tmp_path / "figure.png"
    save_figure_deterministic(figure, path)

    assert b"atplotlib" not in path.read_bytes()


def test_svg_ids_do_not_depend_on_the_deterministic_rcparams(tmp_path: Path, figure: Figure) -> None:
    """Text drawn as paths (the matplotlib default) and as text elements are both reproducible."""
    for index, rc in enumerate(({}, DETERMINISTIC_RCPARAMS)):
        with plt.rc_context(rc):
            first, second = tmp_path / f"first{index}.svg", tmp_path / f"second{index}.svg"
            save_figure_deterministic(figure, first)
            save_figure_deterministic(figure, second)
        assert first.read_bytes() == second.read_bytes(), f"rc overrides: {rc}"


def test_svg_hashsalt_is_pinned_only_for_the_duration_of_the_save(tmp_path: Path, figure: Figure) -> None:
    with plt.rc_context({"svg.hashsalt": None}):  # matplotlib's default: a fresh uuid4 per save
        save_figure_deterministic(figure, tmp_path / "figure.svg")

        assert plt.rcParams["svg.hashsalt"] is None


def test_a_caller_supplied_svg_hashsalt_is_respected(tmp_path: Path, figure: Figure) -> None:
    ours, direct = tmp_path / "ours.svg", tmp_path / "direct.svg"
    with plt.rc_context({"svg.hashsalt": "caller-salt"}):
        save_figure_deterministic(figure, ours)
        figure.savefig(direct, dpi=300, metadata={"Creator": None, "Date": None})
        assert plt.rcParams["svg.hashsalt"] == "caller-salt"

    assert ours.read_bytes() == direct.read_bytes()


@pytest.mark.parametrize("suffix", sorted(RASTER_SIGNATURES))
def test_raster_formats_without_metadata_support_save_reproducibly(tmp_path: Path, figure: Figure, suffix: str) -> None:
    first, second = tmp_path / f"first.{suffix}", tmp_path / f"second.{suffix}"

    save_figure_deterministic(figure, first, dpi=100)
    save_figure_deterministic(figure, second, dpi=100)

    assert first.read_bytes().startswith(RASTER_SIGNATURES[suffix])
    assert first.read_bytes() == second.read_bytes()


@pytest.mark.parametrize("use_subfigure", [False, True], ids=["figure-axes", "subfigure-axes"])
def test_plot_functions_export_reproducible_svg_through_root_figure(tmp_path: Path, use_subfigure: bool) -> None:
    """Plot functions pass ``root_figure(ax)`` to the helper; that route must work for SVG too."""
    saved = []
    for run in range(2):
        fig = plt.figure(figsize=(8, 4))
        ax = fig.subfigures(1, 2)[0].subplots() if use_subfigure else fig.subplots()
        path = tmp_path / f"line{run}.svg"

        basic.lineplot(np.arange(5), np.arange(5) ** 2, ax=ax, output_path=path)

        saved.append(path.read_bytes())
        plt.close(fig)

    assert saved[0].startswith(b"<?xml")
    assert saved[0] == saved[1]
