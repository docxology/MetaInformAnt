"""Regression tests: plot exports must work for an ``Axes`` that lives in a ``SubFigure``.

``Axes.figure`` is a :class:`matplotlib.figure.SubFigure` (which has no
``savefig``) when the caller passes an axes created through
``fig.subfigures(...)``. Every ``output_path`` export used to hand that
SubFigure straight to ``save_figure_deterministic`` and crashed with
``AttributeError: 'SubFigure' object has no attribute 'savefig'``. The exports
now resolve the root ``Figure`` first.

All tests use real matplotlib Agg figures and real PNG files in ``tmp_path``;
no test doubles.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable

import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import pytest
from matplotlib.axes import Axes

from metainformant.visualization.analysis import (
    dimred,
    information,
    quality_assessment,
    quality_omics,
    quality_sequencing,
    statistical,
    timeseries,
)
from metainformant.visualization.genomics import expression, genomics, networks, trees
from metainformant.visualization.plots import basic, general, multidim, specialized

matplotlib.use("Agg")

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"

ExportCall = Callable[[Axes, Path], Any]


def _subfigure_axes(polar: bool = False) -> Axes:
    """Return a real axes that lives inside the left half of a two-subfigure figure."""
    figure = plt.figure(figsize=(8, 4))
    left, _right = figure.subfigures(1, 2)
    if polar:
        return left.subplots(subplot_kw={"projection": "polar"})
    return left.subplots()


def _assert_png(path: Path) -> None:
    assert path.exists(), f"{path.name} was not written"
    assert path.read_bytes()[:8] == PNG_SIGNATURE, f"{path.name} is not a PNG"


_GENES = ["g1", "g2", "g3", "g4", "g5"]
_SAMPLES = ["s1", "s2", "s3", "s4", "s5"]
_MATRIX = np.linspace(-2.0, 2.0, 25).reshape(5, 5)
_EXPRESSION = pd.DataFrame(_MATRIX, index=_GENES, columns=_SAMPLES)
_TREE = nx.DiGraph([("root", "A"), ("root", "B"), ("A", "A1"), ("A", "A2"), ("B", "B1")])
_STAR_TREE = nx.DiGraph([("root", "leaf1"), ("root", "leaf2")])
_GRAPH = nx.complete_graph(6)
_COMMUNITIES = {0: 0, 1: 0, 2: 1, 3: 1, 4: 2, 5: 2}
_GWAS = pd.DataFrame({"CHR": [1, 2, 3], "BP": [1000, 2000, 3000], "P": [0.01, 0.001, 0.05]})
_DE = pd.DataFrame(
    {"log2FoldChange": [1.0, -2.0, 0.5, -0.8], "padj": [1e-6, 1e-7, 0.01, 0.05], "baseMean": [100, 200, 50, 75]}
)

# (id, needs a polar axes, call(ax, output_path)). The first nine entries are the
# shared per-module save helpers; the rest are public entry points whose
# ``output_path`` branch exports ``ax.figure`` directly or through a helper.
_CASES: list[tuple[str, bool, ExportCall]] = [
    ("basic._save_fig", False, lambda ax, p: basic._save_fig(ax, p, "x")),
    ("multidim._save_fig", False, lambda ax, p: multidim._save_fig(ax, p, "x")),
    ("dimred._save_plot", False, lambda ax, p: dimred._save_plot(ax, p, "x")),
    ("information._save_plot", False, lambda ax, p: information._save_plot(ax, p, "x")),
    ("statistical._save_plot", False, lambda ax, p: statistical._save_plot(ax, p, "x")),
    ("timeseries._save_plot", False, lambda ax, p: timeseries._save_plot(ax, p, "x")),
    ("quality_assessment._save_figure", False, lambda ax, p: quality_assessment._save_figure(ax, p, "x")),
    ("quality_omics._save_figure", False, lambda ax, p: quality_omics._save_figure(ax, p, "x")),
    ("quality_sequencing._save_figure", False, lambda ax, p: quality_sequencing._save_figure(ax, p, "x")),
    ("basic.lineplot", False, lambda ax, p: basic.lineplot(np.arange(5), np.arange(5) ** 2, ax=ax, output_path=p)),
    (
        "statistical.histogram",
        False,
        lambda ax, p: statistical.histogram(np.linspace(0.0, 1.0, 50), ax=ax, output_path=p),
    ),
    (
        "information.plot_entropy_profile",
        False,
        lambda ax, p: information.plot_entropy_profile({"s": [1.0, 1.8, 2.5]}, ax=ax, output_path=p),
    ),
    (
        "quality_sequencing.plot_gc_distribution",
        False,
        lambda ax, p: quality_sequencing.plot_gc_distribution(list(np.linspace(30, 60, 40)), ax=ax, output_path=p),
    ),
    (
        "timeseries.plot_time_series",
        False,
        lambda ax, p: timeseries.plot_time_series(
            pd.DataFrame({"v": np.arange(10.0)}, index=pd.date_range("2026-01-01", periods=10)),
            ax=ax,
            output_path=p,
        ),
    ),
    ("general.expression_heatmap", False, lambda ax, p: general.expression_heatmap(_EXPRESSION, ax=ax, output_path=p)),
    (
        "general.pca_plot",
        False,
        lambda ax, p: general.pca_plot(
            pd.DataFrame(
                {
                    "PC1": np.linspace(-1, 1, 20),
                    "PC2": np.cos(np.linspace(0, 3, 20)),
                    "PC1_variance": [25.5] * 20,
                    "PC2_variance": [18.3] * 20,
                }
            ),
            ax=ax,
            output_path=p,
        ),
    ),
    (
        "specialized.plot_venn_diagram",
        False,
        lambda ax, p: specialized.plot_venn_diagram({"A": {1, 2, 3}, "B": {2, 3, 4}}, ax=ax, output_path=p),
    ),
    (
        "specialized.plot_chord_diagram",
        True,
        lambda ax, p: specialized.plot_chord_diagram(
            np.array([[0.0, 2.0, 1.0], [2.0, 0.0, 0.5], [1.0, 0.5, 0.0]]), ["x", "y", "z"], ax=ax, output_path=p
        ),
    ),
    (
        "specialized.plot_circular_barplot",
        True,
        lambda ax, p: specialized.plot_circular_barplot(
            np.array([3.0, 5.0, 2.0, 4.0]), ["a", "b", "c", "d"], ax=ax, output_path=p
        ),
    ),
    (
        "specialized.plot_network_circular_layout",
        True,
        lambda ax, p: specialized.plot_network_circular_layout(nx.cycle_graph(6), ax=ax, output_path=p),
    ),
    (
        "specialized.plot_alluvial_diagram",
        False,
        lambda ax, p: specialized.plot_alluvial_diagram(
            pd.DataFrame({"s1": [1.0, 2.0, 3.0], "s2": [2.0, 1.0, 2.5]}), ["s1", "s2"], ax=ax, output_path=p
        ),
    ),
    (
        "expression.plot_expression_heatmap",
        False,
        lambda ax, p: expression.plot_expression_heatmap(_EXPRESSION, ax=ax, output_path=p),
    ),
    (
        "expression.plot_enrichment_barplot",
        False,
        lambda ax, p: expression.plot_enrichment_barplot(
            pd.DataFrame({"Term": ["Pathway 1", "Pathway 2"], "padj": [1e-5, 1e-4], "Count": [8, 6]}),
            ax=ax,
            output_path=p,
        ),
    ),
    (
        "expression.plot_differential_expression",
        False,
        lambda ax, p: expression.plot_differential_expression(_DE, ax=ax, output_path=p),
    ),
    ("genomics.manhattan_plot", False, lambda ax, p: genomics.manhattan_plot(_GWAS, ax=ax, output_path=p)),
    (
        "genomics.volcano_plot",
        False,
        lambda ax, p: genomics.volcano_plot(_DE[["log2FoldChange", "padj"]], ax=ax, output_path=p),
    ),
    (
        "genomics.regional_plot",
        False,
        lambda ax, p: genomics.regional_plot(
            pd.DataFrame({"CHR": ["1"] * 3, "BP": [1000000, 1010000, 1020000], "P": [0.01, 0.001, 0.05]}),
            chr="1",
            start=900000,
            end=1100000,
            ax=ax,
            output_path=p,
        ),
    ),
    (
        "genomics.circular_manhattan_plot",
        True,
        lambda ax, p: genomics.circular_manhattan_plot(
            pd.DataFrame({"CHR": [1, 2, 1], "BP": [1000000, 2000000, 3000000], "P": [0.01, 0.001, 0.05]}),
            ax=ax,
            output_path=p,
        ),
    ),
    ("genomics.chromosome_ideogram", False, lambda ax, p: genomics.chromosome_ideogram(ax=ax, output_path=p)),
    (
        "genomics.coverage_plot",
        False,
        lambda ax, p: genomics.coverage_plot(
            np.array([10, 20, 15]), np.array([1000000, 1001000, 1002000]), ax=ax, output_path=p
        ),
    ),
    (
        "genomics.variant_plot",
        False,
        lambda ax, p: genomics.variant_plot(pd.DataFrame({"POS": [1000000, 1010000, 1020000]}), ax=ax, output_path=p),
    ),
    ("networks.plot_network_basic", False, lambda ax, p: networks.plot_network_basic(_GRAPH, ax=ax, output_path=p)),
    (
        "networks.plot_network_circular",
        False,
        lambda ax, p: networks.plot_network_circular(_GRAPH, ax=ax, output_path=p),
    ),
    (
        "networks.plot_network_force_directed",
        False,
        lambda ax, p: networks.plot_network_force_directed(_GRAPH, ax=ax, output_path=p),
    ),
    (
        "networks.plot_network_hierarchical",
        False,
        lambda ax, p: networks.plot_network_hierarchical(_TREE, ax=ax, output_path=p),
    ),
    (
        "networks.plot_community_network",
        False,
        lambda ax, p: networks.plot_community_network(_GRAPH, _COMMUNITIES, ax=ax, output_path=p),
    ),
    ("trees.plot_phylo_tree", False, lambda ax, p: trees.plot_phylo_tree(_TREE, ax=ax, output_path=p)),
    ("trees.circular_tree_plot", True, lambda ax, p: trees.circular_tree_plot(_STAR_TREE, ax=ax, output_path=p)),
    ("trees.unrooted_tree_plot", False, lambda ax, p: trees.unrooted_tree_plot(_TREE, ax=ax, output_path=p)),
    (
        "trees.tree_annotation_plot",
        False,
        lambda ax, p: trees.tree_annotation_plot(
            nx.DiGraph([("root", "leaf")]),
            {"leaf": {"color": "green", "label": "Annotated leaf"}},
            ax=ax,
            output_path=p,
        ),
    ),
    ("trees.tree_comparison_plot", False, lambda ax, p: trees.tree_comparison_plot(_TREE, _TREE, ax=ax, output_path=p)),
]


@pytest.mark.parametrize(("polar", "call"), [pytest.param(p, c, id=i) for i, p, c in _CASES])
def test_export_from_subfigure_axes_writes_png(tmp_path: Path, polar: bool, call: ExportCall) -> None:
    """A caller-supplied SubFigure axes plus ``output_path`` must export the root figure."""
    ax = _subfigure_axes(polar=polar)
    assert type(ax.figure).__name__ == "SubFigure"
    output_path = tmp_path / "subfigure_axes.png"

    call(ax, output_path)

    _assert_png(output_path)
    plt.close("all")


def test_export_from_regular_axes_is_unchanged(tmp_path: Path) -> None:
    """Axes on an ordinary Figure keep exporting that same Figure."""
    figure, ax = plt.subplots()
    output_path = tmp_path / "regular_axes.png"

    returned = basic.lineplot(np.arange(5), np.arange(5) ** 2, ax=ax, output_path=output_path)

    assert returned is ax
    _assert_png(output_path)
    plt.close(figure)


class TestRootFigure:
    """``root_figure`` resolves the savable root Figure for any artist."""

    def test_plain_axes_resolve_to_their_figure(self) -> None:
        from metainformant.visualization.config.figures import root_figure

        figure, ax = plt.subplots()

        assert root_figure(ax) is figure
        plt.close(figure)

    def test_subfigure_axes_resolve_to_root_figure(self) -> None:
        from metainformant.visualization.config.figures import root_figure

        figure = plt.figure()
        left, right = figure.subfigures(1, 2)
        left_ax = left.subplots()
        right_ax = right.subplots()

        assert root_figure(left_ax) is figure
        assert root_figure(right_ax) is figure
        assert hasattr(root_figure(left_ax), "savefig")
        plt.close(figure)

    def test_nested_subfigure_axes_resolve_to_root_figure(self) -> None:
        from metainformant.visualization.config.figures import root_figure

        figure = plt.figure()
        outer = figure.subfigures(1, 2)[1]
        inner = outer.subfigures(1, 2)[0]
        nested_ax = inner.subplots()

        assert nested_ax.figure is inner
        assert root_figure(nested_ax) is figure
        plt.close(figure)

    def test_non_axes_artist_resolves_to_root_figure(self) -> None:
        from metainformant.visualization.config.figures import root_figure

        figure = plt.figure()
        ax = figure.subfigures(1, 2)[0].subplots()
        (line,) = ax.plot([0, 1], [0, 1])

        assert root_figure(line) is figure
        plt.close(figure)
