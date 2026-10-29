"""Regression tests for the plot, genomics, network, animation and dashboard defects found in the round-4 sweep.

Each assertion reads what actually landed on the artists (offsets, sizes, patch geometry, tick/legend/text labels,
axes counts, saved bytes), so a plot that is only returned does not pass. Real matplotlib Agg figures and real data;
no test doubles. The tests that run a child interpreter pin the child's ``PYTHONPATH`` to the package under test.

Defects covered (each test fails on the code before the fix):

* genomics: an empty ``regional_plot`` ignored ``output_path``; Manhattan plots sorted chromosomes as strings
  ("1", "10", "2") and drew invalid p-values (negative ones as a peak of 300); ``chromosome_ideogram`` ignored
  ``**kwargs`` and could only draw GRCh38;
* expression: per-sample z-scores labelled as gene z-scores, a hard-coded top 20, a ``s`` argument overridden;
* multidim: parallel coordinates did not colour by group (one legend entry per row), turned constant columns into NaN
  and ignored ``**kwargs``; radar/3-D scatter ignored ``**kwargs`` or ``c``; the pairplot ignored ``ax``;
* networks: unseeded spring layouts, hash-ordered community colours, a real community named -1 drawn as unassigned;
* animations: limits cut data, ``animate_evolution`` showed empty or wrong text for most frames,
  ``animate_clustering`` raised for list labels, ``animate_network`` for nodes added later;
* specialized: the built-in Venn drew inclusive counts in the wrong regions, the Sankey fallback labelled targets with
  source names, and the number of sets raised ``ImportError`` or ``ValueError`` depending on an optional install;
* dashboards: stamped, non-deterministic exports, blank/mis-titled panels, ignored ``size``/``colorscale``.
"""

from __future__ import annotations

import importlib.util
import math
import os
import subprocess
import sys
from pathlib import Path
from typing import Any, Iterator, List

import matplotlib
import matplotlib.patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

import metainformant
from metainformant.visualization.analysis.information import plot_information_network
from metainformant.visualization.dashboards import composite, interactive
from metainformant.visualization.genomics import expression, genomics, networks
from metainformant.visualization.plots import animations, cross_species, multidim, specialized

nx = pytest.importorskip("networkx")

matplotlib.use("Agg")

HAS_PLOTLY = importlib.util.find_spec("plotly") is not None
HAS_MPL_VENN = importlib.util.find_spec("matplotlib_venn") is not None
HAS_SEABORN = importlib.util.find_spec("seaborn") is not None


@pytest.fixture(autouse=True)
def _close_figures() -> Iterator[None]:
    yield
    plt.close("all")


def _offsets(ax: Any, index: int = 0) -> np.ndarray:
    return np.asarray(ax.collections[index].get_offsets(), dtype=float)


def _run_with_hash_seed(program: str, seed: int = 1) -> str:
    """Run *program* in a fresh interpreter under a string-hash seed; its PYTHONPATH is pinned to this package."""
    package_root = str(Path(metainformant.__file__).resolve().parent.parent)
    env = {**os.environ, "PYTHONHASHSEED": str(seed), "PYTHONPATH": package_root, "MPLBACKEND": "Agg"}
    done = subprocess.run([sys.executable, "-c", program], env=env, capture_output=True, text=True, check=True)
    return done.stdout.strip()


# ---------------------------------------------------------------------------
# genomics
# ---------------------------------------------------------------------------


def test_regional_plot_writes_the_file_for_a_region_without_data(tmp_path: Path) -> None:
    """The empty-region branch returned before the ``output_path`` block, so no file was ever written."""
    data = pd.DataFrame({"CHR": ["1"], "BP": [100], "P": [0.01]})
    output = tmp_path / "region.png"

    ax = genomics.regional_plot(data, chr="2", start=0, end=1000, output_path=output)

    assert output.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    assert [text.get_text() for text in ax.texts] == ["No data in region 2:0-1000"]


def test_manhattan_plots_order_chromosomes_naturally() -> None:
    frame = pd.DataFrame(
        {
            "CHR": ["chr1", "chr10", "chr2", "chr1", "chr10", "chr2"],
            "BP": [10, 10, 10, 20, 20, 20],
            "P": [0.5, 0.01, 0.2, 0.3, 0.04, 0.9],
        }
    )

    ax = genomics.manhattan_plot(frame)

    assert [label.get_text() for label in ax.get_xticklabels()] == ["chr1", "chr2", "chr10"]
    assert list(ax.get_xticks()) == sorted(ax.get_xticks())  # laid out left to right in that order

    mixed = pd.DataFrame({"CHR": [2, "X", 1, "Y", 10], "BP": [5, 5, 5, 5, 5], "P": [0.1] * 5})
    assert [label.get_text() for label in genomics.manhattan_plot(mixed).get_xticklabels()] == [
        "1",
        "2",
        "10",
        "X",
        "Y",
    ]  # integers and strings together raised TypeError


def test_circular_manhattan_labels_chromosomes_in_natural_order() -> None:
    frame = pd.DataFrame({"CHR": ["1", "10", "2"] * 2, "BP": [1, 1, 1, 9, 9, 9], "P": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]})

    ax = genomics.circular_manhattan_plot(frame)

    by_angle = sorted(ax.texts, key=lambda text: text.get_position()[0])
    assert [text.get_text() for text in by_angle] == ["1", "2", "10"]


@pytest.mark.parametrize("bad", [-0.1, 1.5])
def test_p_values_outside_zero_to_one_are_refused(bad: float) -> None:
    """A negative p-value was clipped to 1e-300 and drawn at a height of 300; one above 1 was drawn below zero."""
    frame = pd.DataFrame({"CHR": ["1", "1"], "BP": [1, 2], "P": [0.2, bad]})

    with pytest.raises(ValueError, match="between 0 and 1"):
        genomics.manhattan_plot(frame)
    with pytest.raises(ValueError, match="between 0 and 1"):
        genomics.regional_plot(frame, chr="1", start=0, end=10)


def test_chromosome_ideogram_draws_grch38_lengths_by_default() -> None:
    ax = genomics.chromosome_ideogram()

    widths = [patch.get_width() for patch in ax.patches]
    assert len(widths) == 24 and widths[0] == 1.0  # chromosome 1 is the longest
    assert widths[-1] == pytest.approx(57227415 / 248956422)  # chromosome Y
    assert [label.get_text() for label in ax.get_yticklabels()] == [str(i) for i in range(1, 23)] + ["X", "Y"]
    assert ax.get_title() == "Human Chromosome Lengths (GRCh38)"


def test_chromosome_ideogram_accepts_any_genome_and_caller_keywords() -> None:
    """The lengths of a non-human genome could not be given, and ``color``/``alpha`` were silently ignored."""
    ax = genomics.chromosome_ideogram(chromosome_lengths={"LG1": 40_000_000, "LG2": 20_000_000}, color="red", alpha=0.2)

    assert [patch.get_width() for patch in ax.patches] == [1.0, 0.5]
    assert [label.get_text() for label in ax.get_yticklabels()] == ["LG1", "LG2"]
    assert {tuple(patch.get_facecolor()[:3]) for patch in ax.patches} == {(1.0, 0.0, 0.0)}
    assert {patch.get_alpha() for patch in ax.patches} == {0.2}
    assert ax.get_title() == "Chromosome Lengths"
    with pytest.raises(ValueError, match="positive"):
        genomics.chromosome_ideogram(chromosome_lengths={"LG1": 0})
    with pytest.raises(ValueError, match="cannot be empty"):
        genomics.chromosome_ideogram(chromosome_lengths={})


# ---------------------------------------------------------------------------
# expression
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not HAS_SEABORN, reason="the seaborn heatmap path")
def test_expression_heatmap_shows_per_gene_z_scores() -> None:
    """The heatmap was z-scored down each column (per sample) while its title and the documented layout (genes as
    rows) are about genes, and a constant gene would have been NaN."""
    values = np.array([[1.0, 2.0, 3.0, 4.0], [10.0, 40.0, 20.0, 50.0], [5.0, 5.0, 5.0, 5.0]])
    frame = pd.DataFrame(values, index=["g1", "g2", "flat"], columns=["s1", "s2", "s3", "s4"])

    ax = expression.plot_expression_heatmap(frame)

    drawn = np.asarray(ax.collections[0].get_array(), dtype=float).reshape(3, 4)
    expected = (values[:2] - values[:2].mean(axis=1, keepdims=True)) / values[:2].std(axis=1, ddof=1, keepdims=True)
    np.testing.assert_allclose(drawn[:2], expected)
    np.testing.assert_allclose(drawn[2], 0.0)  # a gene that does not vary has no pattern to show
    assert ax.get_title() == "Gene Expression Heatmap (per-gene Z-score)"


def _enrichment(n: int) -> pd.DataFrame:
    return pd.DataFrame(
        {"Term": [f"term {i}" for i in range(n)], "p.adjust": np.linspace(1e-6, 0.04, n), "Count": range(n)}
    )


def test_enrichment_barplot_top_n_is_a_documented_argument() -> None:
    ax = expression.plot_enrichment_barplot(_enrichment(25))
    assert len(ax.patches) == 20 and "top 20 of 25" in ax.get_title()  # the cut-off is stated on the figure

    five = expression.plot_enrichment_barplot(_enrichment(25), top_n=5)
    assert len(five.patches) == 5 and "top 5 of 25" in five.get_title()
    assert expression.plot_enrichment_barplot(_enrichment(8)).get_title() == "Pathway Enrichment Results"
    with pytest.raises(ValueError, match="top_n"):
        expression.plot_enrichment_barplot(_enrichment(5), top_n=0)


def _de_frame(base_mean: Any) -> pd.DataFrame:
    return pd.DataFrame(
        {"log2FoldChange": [2.0, -1.5, 0.2, 0.1], "padj": [1e-5, 0.01, 0.5, 0.9], "baseMean": base_mean}
    )


def test_differential_expression_honours_the_point_size_argument() -> None:
    scaled = expression.plot_differential_expression(_de_frame([10.0, 100.0, 1000.0, 55.0]))
    sized = expression.plot_differential_expression(_de_frame([10.0, 100.0, 1000.0, 55.0]), s=77)

    assert sorted(scaled.collections[0].get_sizes())[0] == 20 and sorted(scaled.collections[0].get_sizes())[-1] == 100
    assert set(sized.collections[0].get_sizes()) == {77.0}  # ``s`` used to be overwritten by the baseMean scaling


def test_differential_expression_with_a_constant_mean_column_has_finite_sizes() -> None:
    ax = expression.plot_differential_expression(_de_frame([5.0, 5.0, 5.0, 5.0]))

    assert set(ax.collections[0].get_sizes()) == {20.0}  # a zero range used to divide by zero (NaN sizes)


# ---------------------------------------------------------------------------
# multidim
# ---------------------------------------------------------------------------


def _groups_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "v1": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            "v2": [6.0, 5.0, 4.0, 3.0, 2.0, 1.0],
            "v3": [7.0] * 6,
            "species": ["A", "A", "B", "B", "C", "C"],
        }
    )


def test_parallel_coordinates_colours_by_group_with_one_legend_entry_per_group() -> None:
    ax = multidim.plot_parallel_coordinates(_groups_frame())

    colours = [line.get_color() for line in ax.lines]
    assert len(ax.lines) == 6
    assert colours[0] == colours[1] and colours[2] == colours[3] and colours[4] == colours[5]
    assert len(set(colours)) == 3  # rows used to cycle through the default colours whatever their group
    assert [text.get_text() for text in ax.get_legend().get_texts()] == ["A", "B", "C"]  # one entry per row before
    assert ax.get_legend().get_title().get_text() == "species"


def test_parallel_coordinates_puts_a_constant_column_mid_axis_and_forwards_keywords() -> None:
    ax = multidim.plot_parallel_coordinates(_groups_frame(), linewidth=4, color="red")

    assert [float(line.get_ydata()[2]) for line in ax.lines] == [0.5] * 6  # a zero range was NaN
    assert {line.get_linewidth() for line in ax.lines} == {4.0}  # extra keywords were ignored
    assert {line.get_color() for line in ax.lines} == {"red"}  # a colour (not a column) was ignored
    assert ax.get_legend() is None


def test_parallel_coordinates_names_its_grouping_column_through_color() -> None:
    frame = _groups_frame().assign(other=["x", "y", "x", "y", "x", "y"])

    ax = multidim.plot_parallel_coordinates(frame, color="other")

    assert [text.get_text() for text in ax.get_legend().get_texts()] == ["x", "y"]


def test_radar_chart_forwards_keywords_and_needs_a_polar_axes() -> None:
    frame = pd.DataFrame({"a": [0.8, 0.6], "b": [0.7, 0.9], "c": [0.9, 0.5]})

    ax = multidim.plot_radar_chart(frame, linewidth=5, alpha=0.4)

    assert {line.get_linewidth() for line in ax.lines} == {5.0} and {line.get_alpha() for line in ax.lines} == {0.4}
    _, cartesian = plt.subplots()
    with pytest.raises(ValueError, match="polar"):
        multidim.plot_radar_chart(frame, ax=cartesian)


def test_3d_scatter_uses_a_colour_given_as_c() -> None:
    frame = pd.DataFrame({"x": [1.0, 2.0, 3.0], "y": [3.0, 1.0, 2.0], "z": [0.0, 1.0, 4.0]})

    literal = multidim.plot_3d_scatter(frame, c="red")
    by_column = multidim.plot_3d_scatter(frame, c="z")

    assert {tuple(face[:3]) for face in literal.collections[0].get_facecolor()} == {(1.0, 0.0, 0.0)}  # c was dropped
    np.testing.assert_allclose(by_column.collections[0].get_array(), frame["z"])
    _, flat = plt.subplots()
    with pytest.raises(ValueError, match="3D"):
        multidim.plot_3d_scatter(frame, ax=flat)


@pytest.mark.skipif(not HAS_SEABORN, reason="seaborn pairplot")
def test_pairplot_refuses_an_axes_and_exports_deterministically(tmp_path: Path) -> None:
    frame = pd.DataFrame({"a": np.arange(8.0), "b": np.arange(8.0) ** 2})
    _, ax = plt.subplots()
    with pytest.raises(ValueError, match="cannot draw onto ax"):
        multidim.plot_pairwise_relationships(frame, ax=ax)

    output = tmp_path / "pairs.pdf"
    multidim.plot_pairwise_relationships(frame, output_path=output)
    assert b"/CreationDate" not in output.read_bytes()  # g.savefig stamped the clock into every PDF


# ---------------------------------------------------------------------------
# networks
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("plot", [networks.plot_network_basic, networks.plot_network_force_directed])
def test_force_directed_layouts_are_reproducible(plot: Any) -> None:
    """``spring_layout`` was unseeded, so the same graph was drawn differently on every call."""
    graph = nx.gnp_random_graph(25, 0.15, seed=1)

    first = _offsets(plot(graph))
    second = _offsets(plot(graph))

    np.testing.assert_array_equal(first, second)
    other = _offsets(plot(graph, layout_kwargs={"seed": 7}))
    assert not np.array_equal(first, other)  # an explicit seed still wins


def test_community_network_layout_and_information_network_are_reproducible() -> None:
    graph = nx.gnp_random_graph(25, 0.15, seed=1)
    communities = {node: node % 3 for node in graph}
    np.testing.assert_array_equal(
        _offsets(networks.plot_community_network(graph, communities)),
        _offsets(networks.plot_community_network(graph, communities)),
    )

    nodes = [f"n{i}" for i in range(12)]
    edges = [(nodes[i], nodes[(i * 5 + 1) % 12], 1.0 + i) for i in range(12) if nodes[i] != nodes[(i * 5 + 1) % 12]]
    np.testing.assert_array_equal(
        _offsets(plot_information_network(nodes, edges)), _offsets(plot_information_network(nodes, edges))
    )


def test_community_colours_are_assigned_in_sorted_order_in_every_process() -> None:
    """A set of community names iterates in hash order, so a community changed colour from run to run."""
    program = (
        "import networkx as nx\n"
        "from metainformant.visualization.genomics.networks import plot_community_network\n"
        "graph = nx.path_graph(['n0', 'n1', 'n2'])\n"
        "ax = plot_community_network(graph, {'n0': 'beta', 'n1': 'alpha', 'n2': 'gamma'})\n"
        "print([t.get_text() for t in ax.get_legend().get_texts()])\n"
    )

    assert _run_with_hash_seed(program) == "['Community alpha', 'Community beta', 'Community gamma']"


def test_a_community_called_minus_one_is_a_community_and_only_missing_nodes_are_unassigned() -> None:
    graph = nx.path_graph(["a", "b", "c"])

    ax = networks.plot_community_network(graph, {"a": -1, "b": 0})

    faces = ax.collections[0].get_facecolor()
    assert tuple(faces[2][:3]) == pytest.approx(matplotlib.colors.to_rgb("gray"))  # c is not in the dictionary
    assert tuple(faces[0][:3]) != pytest.approx(matplotlib.colors.to_rgb("gray"))  # a is community -1, not gray
    assert [text.get_text() for text in ax.get_legend().get_texts()] == ["Community -1", "Community 0", "Unassigned"]


# ---------------------------------------------------------------------------
# animations (each frame is drawn by calling the animation's own frame function)
# ---------------------------------------------------------------------------


def test_time_series_and_trajectory_axes_contain_all_the_data() -> None:
    """The limits were ``min * 1.1`` and ``max * 1.1``: a positive minimum lay outside the axes."""
    fig, _ = animations.animate_time_series(np.array([[2.0, 3.0, 5.0, 4.0]]))
    low, high = fig.axes[0].get_ylim()
    assert low < 2.0 and high > 5.0

    fig, _ = animations.animate_trajectory([np.array([[1.0, 2.0], [2.0, 3.0], [4.0, 5.0]])])
    (x_low, x_high), (y_low, y_high) = fig.axes[0].get_xlim(), fig.axes[0].get_ylim()
    assert x_low < 1.0 and x_high > 4.0 and y_low < 2.0 and y_high > 5.0


def test_time_series_animation_accepts_a_dataframe_with_column_labels() -> None:
    """The docstring promised a DataFrame; it was rejected."""
    frame = pd.DataFrame({"alpha": [1.0, 2.0, 3.0], "beta": [3.0, 2.0, 1.0]})

    fig, anim = animations.animate_time_series(frame)

    anim._func(2)
    lines = fig.axes[0].lines
    assert [line.get_label() for line in lines] == ["alpha", "beta"]
    np.testing.assert_array_equal(lines[0].get_ydata(), [1.0, 2.0, 3.0])
    np.testing.assert_array_equal(lines[1].get_ydata(), [3.0, 2.0, 1.0])


def test_evolution_animation_reveals_each_generation_in_turn() -> None:
    """Frame ``k`` mapped to generation ``k``, so after the first character every middle generation was blank."""
    sequences = ["AAAA", "AAGA", "ACGA"]

    fig, anim = animations.animate_evolution(sequences)

    generation_text, mutation_text = fig.axes[0].texts
    seen = []
    for frame in range(12):
        anim._func(frame)
        seen.append((generation_text.get_text(), mutation_text.get_text()))
    assert seen[0] == ("Generation 1: A", "")
    assert seen[3] == ("Generation 1: AAAA", "")
    assert seen[5] == ("Generation 2: AA", "")  # shown so far matches generation 1
    assert seen[6] == ("Generation 2: AAG", "Mutations: 3:A->G")
    assert seen[11] == ("Generation 3: ACGA", "Mutations: 2:A->C, 3:A->G")
    with pytest.raises(ValueError, match="at least one character"):
        animations.animate_evolution(["", ""])


def test_clustering_animation_recolours_by_each_frames_labels_with_one_colour_scale() -> None:
    data = np.array([[0.0, 0.0], [1.0, 1.0], [2.0, 0.5], [3.0, 2.0]])
    labels = [[0, 0, 1, 1], [0, 1, 1, 2], [2, 2, 2, 2]]  # plain lists, as documented

    fig, anim = animations.animate_clustering(data, labels, cmap="viridis", alpha=0.5)

    scatter = fig.axes[0].collections[0]
    assert scatter.get_cmap().name == "viridis" and scatter.get_alpha() == 0.5  # no keyword collision
    for frame, expected in enumerate(labels):
        anim._func(frame)
        np.testing.assert_array_equal(scatter.get_array(), expected)
    assert scatter.get_clim() == (0, 2)  # the scale spans every frame, so a label keeps its colour
    assert fig.axes[0].get_title() == "Clustering Animation - Iteration 3"


def test_clustering_animation_checks_one_label_per_point() -> None:
    with pytest.raises(ValueError, match="one label for each of the 3 data points"):
        animations.animate_clustering(np.zeros((3, 2)), [[0, 1, 0], [0, 1]])


def test_network_animation_survives_nodes_that_appear_later() -> None:
    """The layout came from the first graph only, so ``nx.draw`` raised for a node added in a later frame."""
    first = nx.Graph([("a", "b")])
    second = nx.Graph([("a", "b"), ("b", "c")])

    fig, anim = animations.animate_network([first, second])

    anim._func(1)
    assert len(_offsets(fig.axes[0])) == 3
    assert {text.get_text() for text in fig.axes[0].texts} == {"a", "b", "c"}


# ---------------------------------------------------------------------------
# specialized
# ---------------------------------------------------------------------------


def _region_counts(ax: Any, sets: List[set]) -> List[tuple]:
    """For every numeric label: which circles contain it, and the number it shows."""
    circles = [patch for patch in ax.patches if hasattr(patch, "radius")]
    found = []
    for text in ax.texts:
        head = text.get_text().split("\n")[0]
        if not head.isdigit():
            continue
        x, y = text.get_position()
        membership = tuple(math.hypot(x - c.center[0], y - c.center[1]) < c.radius for c in circles)
        found.append((membership, int(head)))
    return found


@pytest.mark.skipif(HAS_MPL_VENN, reason="matplotlib-venn draws the diagram, not the built-in fallback")
def test_builtin_venn3_puts_each_exclusive_count_in_its_own_region() -> None:
    a, b, c = {1, 2, 3, 4, 5}, {4, 5, 6, 7}, {5, 7, 8, 9, 10}

    ax = specialized.plot_venn_diagram({"A": a, "B": b, "C": c})

    found = _region_counts(ax, [a, b, c])
    assert len(found) == 7 and len({membership for membership, _ in found}) == 7  # one label in each of the 7 regions
    for membership, shown in found:
        universe = a | b | c
        exact = {x for x in universe if tuple(x in s for s in (a, b, c)) == membership}
        assert shown == len(exact), (membership, shown)  # exclusive counts: whole-set totals were printed before
    assert {text.get_text() for text in ax.texts} >= {"A", "B", "C"}


@pytest.mark.skipif(HAS_MPL_VENN, reason="matplotlib-venn draws the diagram, not the built-in fallback")
def test_builtin_venn2_puts_the_intersection_in_the_overlap() -> None:
    a, b = {1, 2, 3, 4, 5}, {4, 5, 6, 7}

    ax = specialized.plot_venn_diagram({"A": a, "B": b})

    assert sorted(_region_counts(ax, [a, b])) == [((False, True), 2), ((True, False), 3), ((True, True), 2)]


@pytest.mark.parametrize("count", [1, 4])
def test_the_number_of_sets_is_a_value_error_whether_or_not_matplotlib_venn_is_installed(count: int) -> None:
    with pytest.raises(ValueError, match="two or three sets"):
        specialized.plot_venn_diagram({f"S{i}": {i} for i in range(count)})


@pytest.mark.skipif(HAS_PLOTLY, reason="plotly draws the Sankey diagram, not the matplotlib fallback")
def test_sankey_fallback_labels_targets_with_target_names_and_scales_by_share() -> None:
    flows = [("in", "mid", 10.0), ("mid", "out", 7.0), ("mid", "waste", 3.0)]

    fig = specialized.plot_sankey_diagram(flows)

    left, right = fig.axes
    assert [label.get_text() for label in left.get_yticklabels()] == ["in", "mid"]
    assert [label.get_text() for label in right.get_yticklabels()] == ["mid", "out", "waste"]  # sources were repeated
    widths = [line.get_linewidth() for line in left.lines]
    assert widths == sorted(widths, reverse=True) and max(widths) == 10.0

    thousand = specialized.plot_sankey_diagram([(s, t, v * 1000) for s, t, v in flows])
    assert [line.get_linewidth() for line in thousand.axes[0].lines] == widths  # the same figure in any unit


# ---------------------------------------------------------------------------
# chord and alluvial diagrams
# ---------------------------------------------------------------------------


def test_chord_layout_gives_each_node_and_each_end_an_angle_proportional_to_its_entries() -> None:
    matrix = np.array([[0.0, 2.0, 1.0], [2.0, 0.0, 0.5], [1.0, 0.5, 0.0]])
    gap = 0.1
    unit = (2 * math.pi - 3 * gap) / matrix.sum()

    spans, ends = specialized._chord_layout(matrix, gap)

    assert [end - start for start, end in spans] == pytest.approx([3.0 * unit, 2.5 * unit, 1.5 * unit])
    assert ends[(0, 1)][1] - ends[(0, 1)][0] == pytest.approx(2.0 * unit)  # the ends of one chord are as wide...
    assert ends[(1, 0)][1] - ends[(1, 0)][0] == pytest.approx(2.0 * unit)  # ...as each other in a symmetric matrix
    assert spans[1][0] == pytest.approx(spans[0][1] + gap) and spans[2][1] + gap == pytest.approx(2 * math.pi)


def test_chord_diagram_draws_one_ribbon_per_pair_inside_the_circle() -> None:
    """The ribbons were a Bezier path with the wrong number of control points drawn in (angle, radius) coordinates."""
    matrix = np.array([[0.0, 2.0, 1.0], [2.0, 0.0, 0.5], [1.0, 0.5, 0.0]])

    ax = specialized.plot_chord_diagram(matrix, ["x", "y", "z"])

    polygons = [np.asarray(patch.get_xy(), dtype=float) for patch in ax.patches]
    assert len(polygons) == 3 + 3  # three ribbons (0-1, 0-2, 1-2) and three node arcs
    ribbons, nodes = polygons[:3], polygons[3:]
    assert max(np.hypot(*ribbon.T).max() for ribbon in ribbons) <= 1.0 + 1e-9  # chords stay inside the ring
    assert all(0.99 < np.hypot(*node.T).min() and np.hypot(*node.T).max() < 1.07 for node in nodes)
    assert [text.get_text() for text in ax.texts] == ["x", "y", "z"]
    assert ax.get_aspect() == 1.0


def test_chord_diagram_uses_both_triangles_of_an_asymmetric_matrix() -> None:
    """Only ``matrix[i, j]`` with ``i < j`` was read, so a relationship stored below the diagonal vanished."""
    lower_only = np.array([[0.0, 0.0], [2.0, 0.0]])
    spans, ends = specialized._chord_layout(np.array([[0.0, 3.0], [1.0, 0.0]]), 0.0)

    ax = specialized.plot_chord_diagram(lower_only)

    assert len(ax.patches) == 1 + 2  # the one ribbon and the two node arcs
    widths = (ends[(0, 1)][1] - ends[(0, 1)][0], ends[(1, 0)][1] - ends[(1, 0)][0])
    assert widths[0] == pytest.approx(3 * widths[1])  # a different width at each end, one per direction


def test_chord_diagram_on_a_polar_axes_and_input_validation() -> None:
    matrix = np.array([[0.0, 1.0], [1.0, 0.0]])
    figure = plt.figure()
    polar = figure.add_subplot(111, projection="polar")

    ax = specialized.plot_chord_diagram(matrix, ax=polar)

    radii = np.concatenate([np.asarray(patch.get_xy())[:, 1] for patch in ax.patches])
    assert radii.max() < 1.07  # (angle, radius) coordinates: the ring is at radius 1 whatever the projection
    with pytest.raises(ValueError, match="non-negative"):
        specialized.plot_chord_diagram(np.array([[0.0, -1.0], [-1.0, 0.0]]))
    with pytest.raises(ValueError, match="no positive entry"):
        specialized.plot_chord_diagram(np.zeros((2, 2)))
    with pytest.raises(ValueError, match="colors"):
        specialized.plot_chord_diagram(matrix, colors=["red"])


def _alluvial_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "s1": ["x", "x", "y", "y", "y"],
            "s2": ["p", "q", "p", "q", "q"],
            "s3": ["u", "u", "v", "v", "u"],
        }
    )


@pytest.mark.skipif(importlib.util.find_spec("alluvial") is not None, reason="the alluvial package draws the diagram")
def test_builtin_alluvial_draws_a_bar_per_category_and_a_ribbon_per_transition() -> None:
    """The fallback drew one coloured line per row through the raw values: no bars, no ribbons, no flow."""
    ax = specialized.plot_alluvial_diagram(_alluvial_frame(), ["s1", "s2", "s3"])

    bars = [patch for patch in ax.patches if isinstance(patch, matplotlib.patches.Rectangle)]
    ribbons = [patch for patch in ax.patches if isinstance(patch, matplotlib.patches.PathPatch)]
    assert sorted(bar.get_height() for bar in bars) == [2.0, 2.0, 2.0, 3.0, 3.0, 3.0]  # x,y / p,q / u,v by row count
    assert {text.get_text() for text in ax.texts} == {"x", "y", "p", "q", "u", "v"}
    assert len(ribbons) == 8  # four transitions between s1 and s2, four between s2 and s3
    heights = sorted(float(ribbon.get_path().vertices[0][1] - ribbon.get_path().vertices[7][1]) for ribbon in ribbons)
    assert heights == pytest.approx([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0])  # rows making each transition
    for ribbon in ribbons:  # a ribbon is as tall at its target end as at its source end
        v = ribbon.get_path().vertices
        assert v[0][1] - v[7][1] == pytest.approx(v[3][1] - v[4][1])


@pytest.mark.skipif(importlib.util.find_spec("alluvial") is not None, reason="the alluvial package draws the diagram")
def test_builtin_alluvial_accounts_for_every_row_at_both_ends_of_every_bar() -> None:
    ax = specialized.plot_alluvial_diagram(_alluvial_frame(), ["s1", "s2", "s3"])
    bars = [patch for patch in ax.patches if isinstance(patch, matplotlib.patches.Rectangle)]
    ribbons = [patch.get_path().vertices for patch in ax.patches if isinstance(patch, matplotlib.patches.PathPatch)]

    for bar in bars:
        low, high = bar.get_y(), bar.get_y() + bar.get_height()
        right_edge, left_edge = bar.get_x() + bar.get_width(), bar.get_x()
        leaving = sum(
            v[0][1] - v[7][1] for v in ribbons if abs(v[0][0] - right_edge) < 1e-9 and low - 1e-9 <= v[7][1] < high
        )
        arriving = sum(
            v[3][1] - v[4][1] for v in ribbons if abs(v[3][0] - left_edge) < 1e-9 and low - 1e-9 <= v[4][1] < high
        )
        stage = round(bar.get_x() + bar.get_width() / 2)
        assert leaving == (pytest.approx(bar.get_height()) if stage < 2 else 0.0)  # the last stage has nothing leaving
        assert arriving == (pytest.approx(bar.get_height()) if stage > 0 else 0.0)  # the first has nothing arriving


@pytest.mark.skipif(importlib.util.find_spec("alluvial") is not None, reason="the alluvial package draws the diagram")
def test_builtin_alluvial_keeps_rows_with_a_missing_stage_in_their_own_category() -> None:
    frame = pd.DataFrame({"a": ["x", "x", None], "b": ["p", "q", "q"]})

    ax = specialized.plot_alluvial_diagram(frame, ["a", "b"])

    assert {text.get_text() for text in ax.texts} == {"x", "missing", "p", "q"}
    assert (
        sum(bar.get_height() for bar in ax.patches if isinstance(bar, matplotlib.patches.Rectangle)) == 6.0
    )  # 3 rows x 2


# ---------------------------------------------------------------------------
# composite dashboards and interactive fallbacks
# ---------------------------------------------------------------------------


def test_genomic_overview_refuses_unusable_values_instead_of_leaving_empty_panels() -> None:
    with pytest.raises(ValueError, match="DataFrame"):
        composite.genomic_overview({"expression": np.zeros((3, 3))})
    with pytest.raises(ValueError, match="2-D array"):
        composite.genomic_overview({"pca_data": np.arange(5.0)})
    with pytest.raises(ValueError, match="same shape"):
        composite.genomic_overview({"log2fc": np.zeros(4), "pvalues": np.ones(5)})


def test_genomic_overview_says_when_the_heatmap_shows_only_the_first_rows_and_columns() -> None:
    big = pd.DataFrame(np.arange(30 * 25.0).reshape(30, 25))
    small = pd.DataFrame(np.arange(25.0).reshape(5, 5))

    big_title = composite.genomic_overview({"expression": big}).axes[0].get_title()
    small_title = composite.genomic_overview({"expression": small}).axes[0].get_title()

    assert big_title == "Expression Heatmap (first 20 genes x 20 samples)"  # not "top 20": nothing was ranked
    assert small_title == "Expression Heatmap"


def test_dashboard_exports_are_deterministic(tmp_path: Path) -> None:
    """``fig.savefig`` stamped a creation date and producer into every PDF; the shared exporter strips them."""
    frame = pd.DataFrame(
        {"x": [1.0, 2.0, 3.0], "y": [3.0, 1.0, 2.0], "log2FoldChange": [2.0, -2.0, 0.1], "pvalue": [0.01, 0.02, 0.9]}
    )
    outputs = {
        "multi_panel": lambda p: composite.multi_panel(
            [{"func": lambda ax: ax.plot([1, 2]), "title": "t"}], ncols=1, output_path=p
        ),
        "genomic_overview": lambda p: composite.genomic_overview({"gc_content": np.arange(10.0)}, output_path=p),
        "qc_summary": lambda p: composite.qc_summary({"gc_content": np.arange(10.0)}, output_path=p),
        "scatter": lambda p: interactive.interactive_scatter(frame, "x", "y", output_path=p),
        "heatmap": lambda p: interactive.interactive_heatmap(np.arange(9.0).reshape(3, 3), output_path=p),
        "volcano": lambda p: interactive.interactive_volcano(frame, output_path=p),
    }
    if HAS_PLOTLY:
        outputs = {
            name: call for name, call in outputs.items() if name in ("multi_panel", "genomic_overview", "qc_summary")
        }

    for name, call in outputs.items():
        path = tmp_path / f"{name}.pdf"
        call(path)
        content = path.read_bytes()
        assert content.startswith(b"%PDF") and b"/CreationDate" not in content, name


@pytest.mark.skipif(HAS_PLOTLY, reason="plotly draws the interactive figure, not the static fallback")
def test_static_scatter_draws_the_size_column_and_colormap_names_are_honoured() -> None:
    frame = pd.DataFrame({"x": [1.0, 2.0, 3.0], "y": [3.0, 1.0, 2.0], "w": [10.0, 20.0, 30.0], "g": ["a", "b", "a"]})

    ax = interactive.interactive_scatter(frame, "x", "y", size="w")
    assert list(ax.collections[0].get_sizes()) == [20.0, 110.0, 200.0]  # the size column was ignored

    grouped = interactive.interactive_scatter(frame, "x", "y", color="g", size="w")
    sizes = {coll.get_label(): list(coll.get_sizes()) for coll in grouped.collections}
    assert sizes == {"a": [20.0, 200.0], "b": [110.0]}

    heat = interactive.interactive_heatmap(np.arange(9.0).reshape(3, 3), colorscale="Viridis")
    assert heat.images[0].get_cmap().name == "viridis"  # colorscale was ignored
    unknown = interactive.interactive_heatmap(np.arange(9.0).reshape(3, 3), colorscale="NoSuchScale")
    assert unknown.images[0].get_cmap().name == "RdBu_r"


@pytest.mark.skipif(HAS_PLOTLY, reason="plotly draws the interactive figure, not the static fallback")
def test_static_volcano_legend_lists_only_categories_that_have_genes() -> None:
    frame = pd.DataFrame({"log2FoldChange": [2.0, -2.0, 0.1], "pvalue": [0.001, 0.001, 0.9], "gene": ["a", "b", "c"]})

    ax = interactive.interactive_volcano(frame)

    labels = [text.get_text() for text in ax.get_legend().get_texts()]
    assert labels == [
        "Not Significant",
        "Up-regulated",
        "Down-regulated",
    ]  # "Significant" is unreachable and had an entry


# ---------------------------------------------------------------------------
# cross-species
# ---------------------------------------------------------------------------


def _nine_species() -> pd.DataFrame:
    rng = np.random.default_rng(0)
    names = [f"sp_{i}" for i in range(9)]
    matrix = rng.uniform(0.1, 1.6, (9, 9))
    matrix = (matrix + matrix.T) / 2
    np.fill_diagonal(matrix, 0.0)
    return pd.DataFrame(matrix, index=names, columns=names)


def test_family_violin_is_byte_identical_between_calls(tmp_path: Path) -> None:
    """seaborn's stripplot draws its jitter from the global NumPy random state, so every call differed."""
    divergence = _nine_species()
    families = {name: f"F{i // 3}" for i, name in enumerate(divergence.index)}
    state = np.random.get_state()

    cross_species.plot_family_violin(divergence, tmp_path / "a.png", families)
    cross_species.plot_family_violin(divergence, tmp_path / "b.png", families)

    assert (tmp_path / "a.png").read_bytes() == (tmp_path / "b.png").read_bytes()
    after = np.random.get_state()
    assert (
        state[1].tolist() == after[1].tolist() and state[2] == after[2]
    )  # and the caller's random stream is untouched


def test_combined_summary_dendrogram_is_monochrome_like_the_stand_alone_one(tmp_path: Path) -> None:
    """Branches below an arbitrary distance of 0.4 were coloured, implying clusters nothing in the data defines."""
    divergence = _nine_species()
    combined, alone = tmp_path / "combined.svg", tmp_path / "alone.svg"

    cross_species.plot_combined_summary(divergence, combined)
    cross_species.plot_dendrogram(divergence, alone)

    assert "#7f8c8d" not in combined.read_text()  # the colour used for the links above the threshold
    assert "#7f8c8d" not in alone.read_text()
