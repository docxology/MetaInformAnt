"""Regression tests for the visualization-package defects fixed in round 3.

Every assertion reads the data that actually landed on the artists (line data,
collection offsets and sizes, colour maps, patch geometry, axes counts), so a
plot that is merely returned without the requested properties fails.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

from matplotlib.collections import PathCollection  # noqa: E402
from matplotlib.colors import to_rgb  # noqa: E402

from metainformant.visualization.analysis.information import (  # noqa: E402
    plot_entropy_profile,
    plot_information_landscape,
    plot_information_network,
    plot_mutual_information_matrix,
    plot_renyi_spectra,
)
from metainformant.visualization.analysis.timeseries import plot_forecast, plot_trend_analysis  # noqa: E402
from metainformant.visualization.dashboards.composite import multi_panel  # noqa: E402
from metainformant.visualization.genomics.genomics import circular_manhattan_plot, regional_plot  # noqa: E402
from metainformant.visualization.plots.specialized import plot_circular_barplot  # noqa: E402


def _offsets(collection: PathCollection) -> np.ndarray:
    return np.asarray(collection.get_offsets(), dtype=float)


# ---------------------------------------------------------------------------
# Caller keyword arguments must override defaults instead of colliding with them
# ---------------------------------------------------------------------------


class TestInformationKeywordOverrides:
    """``cmap``/``marker``/... used to raise ``TypeError: multiple values for keyword argument``."""

    def test_mutual_information_matrix_honours_cmap(self):
        ax = plot_mutual_information_matrix(np.eye(3), cmap="magma")
        assert ax.images[0].get_cmap().name == "magma"

    def test_mutual_information_matrix_default_cmap_is_viridis(self):
        ax = plot_mutual_information_matrix(np.eye(3))
        assert ax.images[0].get_cmap().name == "viridis"

    def test_entropy_profile_honours_marker(self):
        ax = plot_entropy_profile({"a": [1.0, 2.0, 3.0], "b": [2.0, 1.0, 0.5]}, marker="s")
        assert [line.get_marker() for line in ax.lines] == ["s", "s"]

    def test_entropy_profile_default_marker_is_circle(self):
        ax = plot_entropy_profile({"a": [1.0, 2.0, 3.0]})
        assert ax.lines[0].get_marker() == "o"

    def test_renyi_spectra_honours_marker(self):
        ax = plot_renyi_spectra({"a": [3.0, 2.0, 1.0]}, [0.5, 1.0, 2.0], marker="^")
        assert ax.lines[0].get_marker() == "^"
        assert ax.lines[0].get_xdata().tolist() == [0.5, 1.0, 2.0]

    def test_renyi_spectra_default_marker_is_circle(self):
        ax = plot_renyi_spectra({"a": [3.0, 2.0, 1.0]}, [0.5, 1.0, 2.0])
        assert ax.lines[0].get_marker() == "o"

    def test_information_landscape_honours_cmap_and_levels(self):
        ax = plot_information_landscape(np.arange(16.0).reshape(4, 4), cmap="magma", levels=[0, 5, 10, 15])
        contour = ax.collections[0]
        assert contour.get_cmap().name == "magma"
        assert contour.levels.tolist() == [0.0, 5.0, 10.0, 15.0]

    def test_information_landscape_defaults(self):
        ax = plot_information_landscape(np.arange(16.0).reshape(4, 4))
        contour = ax.collections[0]
        assert contour.get_cmap().name == "viridis"
        assert len(contour.levels) > 5  # 20 requested levels, not the 3-level override above

    def test_information_network_honours_drawing_and_layout_kwargs(self):
        nx = pytest.importorskip("networkx")
        nodes = ["a", "b", "c"]
        edges = [("a", "b", 1.0), ("b", "c", 2.0)]

        ax = plot_information_network(
            nodes, edges, node_color="red", node_size=500, alpha=0.3, layout_kwargs={"seed": 3}
        )

        nodes_artist = ax.collections[0]
        assert np.allclose(nodes_artist.get_facecolor()[0][:3], to_rgb("red"))
        assert set(nodes_artist.get_sizes().tolist()) == {500.0}  # a scalar size applies to every node
        assert nodes_artist.get_alpha() == pytest.approx(0.3)
        # layout_kwargs reach the layout (not networkx.draw): positions equal the seeded spring layout
        graph = nx.DiGraph()
        graph.add_nodes_from(nodes)
        graph.add_weighted_edges_from(edges)
        expected = nx.spring_layout(graph, seed=3)
        np.testing.assert_allclose(_offsets(nodes_artist), [expected[n] for n in nodes])

    def test_information_network_default_node_sizes_follow_degree(self):
        pytest.importorskip("networkx")
        ax = plot_information_network(["a", "b", "c"], [("a", "b", 1.0), ("b", "c", 2.0)])
        # degree * 100 + 300 for the path a-b-c (degrees 1, 2, 1)
        assert ax.collections[0].get_sizes().tolist() == [400.0, 500.0, 400.0]


class TestTimeSeriesKeywordOverrides:
    """``color``/``alpha``/``linestyle`` used to collide with the hard-coded defaults."""

    @staticmethod
    def _series():
        history = pd.Series(np.arange(10.0), index=pd.date_range("2020-01-01", periods=10))
        forecast = pd.Series(np.arange(10.0, 14.0), index=pd.date_range("2020-01-11", periods=4))
        return history, forecast

    def test_forecast_honours_color(self):
        history, forecast = self._series()
        ax = plot_forecast(history, forecast, color="green")
        assert [line.get_color() for line in ax.lines[:2]] == ["green", "green"]

    def test_forecast_default_colors_and_linestyle(self):
        history, forecast = self._series()
        ax = plot_forecast(history, forecast)
        assert [line.get_color() for line in ax.lines[:2]] == ["blue", "red"]
        assert [line.get_linestyle() for line in ax.lines[:2]] == ["-", "--"]

    def test_forecast_honours_linestyle(self):
        history, forecast = self._series()
        ax = plot_forecast(history, forecast, linestyle=":")
        assert [line.get_linestyle() for line in ax.lines[:2]] == [":", ":"]

    def test_forecast_honours_alpha(self):
        history, forecast = self._series()
        ax = plot_forecast(history, forecast, alpha=0.5)
        assert [line.get_alpha() for line in ax.lines[:2]] == [0.5, 0.5]

    def test_trend_analysis_honours_alpha(self):
        history, _ = self._series()
        ax = plot_trend_analysis(history, alpha=0.3)
        assert ax.lines[0].get_alpha() == pytest.approx(0.3)

    def test_trend_analysis_default_alpha(self):
        history, _ = self._series()
        ax = plot_trend_analysis(history)
        assert ax.lines[0].get_alpha() == pytest.approx(0.7)


class TestCircularBarplotKeywordOverrides:
    def test_alpha_is_honoured(self):
        ax = plot_circular_barplot(np.array([1.0, 2.0, 3.0]), alpha=0.2)
        assert [patch.get_alpha() for patch in ax.patches] == [pytest.approx(0.2)] * 3

    def test_default_alpha(self):
        ax = plot_circular_barplot(np.array([1.0, 2.0, 3.0]))
        assert [patch.get_alpha() for patch in ax.patches] == [pytest.approx(0.7)] * 3

    def test_width_and_bottom_are_honoured(self):
        ax = plot_circular_barplot(np.array([1.0, 2.0, 3.0]), width=0.1, bottom=0.5)
        assert [patch.get_width() for patch in ax.patches] == [pytest.approx(0.1)] * 3
        assert [patch.get_y() for patch in ax.patches] == [pytest.approx(0.5)] * 3


# ---------------------------------------------------------------------------
# Sankey node indices
# ---------------------------------------------------------------------------


class TestSankeyLinks:
    """Plotly is not installed in the test environment, so the index building is a pure function."""

    FLOWS = [("a", "x", 5.0), ("b", "x", 3.0), ("c", "y", 2.0), ("a", "y", 1.0), ("d", "z", 4.0)]

    def test_indices_decode_to_the_original_flows(self):
        from metainformant.visualization.plots.specialized import _build_sankey_links

        labels, sources, targets, values = _build_sankey_links(self.FLOWS)

        assert [(labels[s], labels[t]) for s, t in zip(sources, targets)] == [(s, t) for s, t, _ in self.FLOWS]
        assert values == [5.0, 3.0, 2.0, 1.0, 4.0]

    def test_labels_are_unique_and_in_first_seen_order(self):
        from metainformant.visualization.plots.specialized import _build_sankey_links

        labels, _, _, _ = _build_sankey_links(self.FLOWS)

        assert labels == ["a", "x", "b", "c", "y", "d", "z"]

    def test_random_flows_round_trip(self):
        from metainformant.visualization.plots.specialized import _build_sankey_links

        rng = np.random.default_rng(7)
        names = [f"node{i}" for i in range(12)]
        flows = [(names[rng.integers(12)], names[rng.integers(12)], float(rng.uniform(1, 9))) for _ in range(60)]

        labels, sources, targets, values = _build_sankey_links(flows)

        assert len(labels) == len(set(labels))
        assert [(labels[s], labels[t]) for s, t in zip(sources, targets)] == [(s, t) for s, t, _ in flows]
        assert values == [v for _, _, v in flows]

    def test_empty_flows(self):
        from metainformant.visualization.plots.specialized import _build_sankey_links

        assert _build_sankey_links([]) == ([], [], [], [])


# ---------------------------------------------------------------------------
# Genomics plots
# ---------------------------------------------------------------------------


class TestCircularManhattanStyling:
    """``s=`` and ``alpha=`` used to be popped inside the chromosome loop, so only chromosome 1 got them."""

    @staticmethod
    def _frame():
        return pd.DataFrame(
            {
                "CHR": [1, 1, 2, 2, 3, 3],
                "BP": [10, 20, 10, 20, 10, 20],
                "P": [1e-3, 1e-2, 1e-4, 1e-1, 1e-5, 0.5],
            }
        )

    def test_size_and_alpha_apply_to_every_chromosome(self):
        ax = circular_manhattan_plot(self._frame(), s=50, alpha=0.2)

        assert len(ax.collections) == 3
        assert [c.get_sizes().tolist() for c in ax.collections] == [[50.0]] * 3
        assert [c.get_alpha() for c in ax.collections] == [pytest.approx(0.2)] * 3

    def test_defaults_apply_to_every_chromosome(self):
        ax = circular_manhattan_plot(self._frame())

        assert [c.get_sizes().tolist() for c in ax.collections] == [[2.0]] * 3
        assert [c.get_alpha() for c in ax.collections] == [pytest.approx(0.7)] * 3


class TestHierarchicalNetworkLayout:
    """``subset_key`` was given node->level, which networkx rejects, so the spring layout always won."""

    def test_nodes_are_placed_in_one_column_per_level(self):
        nx = pytest.importorskip("networkx")
        from metainformant.visualization.genomics.networks import plot_network_hierarchical

        graph = nx.DiGraph([("r", "a"), ("r", "b"), ("a", "c"), ("a", "d"), ("b", "e")])
        ax = plot_network_hierarchical(graph)

        x = dict(zip(graph.nodes(), _offsets(ax.collections[0])[:, 0]))
        assert x["a"] == pytest.approx(x["b"])
        assert x["c"] == pytest.approx(x["d"])
        assert x["d"] == pytest.approx(x["e"])
        assert x["r"] < x["a"] < x["c"]

    def test_unreachable_nodes_go_to_the_last_column(self):
        nx = pytest.importorskip("networkx")
        from metainformant.visualization.genomics.networks import plot_network_hierarchical

        graph = nx.DiGraph([("r", "a"), ("r", "b")])
        graph.add_edge("p", "q")  # a second root: not reachable from the first root
        graph.add_edge("q", "p")  # ... and a cycle, so it has no in-degree-zero node of its own
        ax = plot_network_hierarchical(graph)

        x = dict(zip(graph.nodes(), _offsets(ax.collections[0])[:, 0]))
        assert x["r"] < x["a"] == pytest.approx(x["b"])
        assert x["p"] == pytest.approx(x["q"])
        assert x["p"] > x["a"]


class TestRegionalRecombination:
    """The recombination-rate placeholder is replaced by real data when given and omitted otherwise."""

    @staticmethod
    def _data():
        return pd.DataFrame(
            {
                "CHR": ["1"] * 5,
                "BP": [100, 200, 300, 400, 500],
                "P": [1e-3, 1e-4, 1e-6, 1e-2, 0.5],
            }
        )

    @staticmethod
    def _twin(ax):
        others = [other for other in ax.figure.axes if other is not ax]
        assert len(others) == 1
        return others[0]

    def test_supplied_rates_are_drawn_on_a_secondary_axis(self):
        recombination = pd.DataFrame({"BP": [0, 300, 600], "RATE": [1.0, 5.0, 2.0]})

        ax = regional_plot(self._data(), chr="1", start=0, end=600, recombination_rate=recombination)

        twin = self._twin(ax)
        assert len(twin.lines) == 1
        assert twin.lines[0].get_xdata().tolist() == [0, 300, 600]
        assert twin.lines[0].get_ydata().tolist() == [1.0, 5.0, 2.0]
        assert twin.get_ylabel() == "Recombination rate (cM/Mb)"

    def test_rate_positions_are_relative_to_the_region_start_like_the_p_values(self):
        recombination = pd.DataFrame({"BP": [150, 300, 450], "RATE": [1.0, 5.0, 2.0]})

        ax = regional_plot(self._data(), chr="1", start=100, end=500, recombination_rate=recombination)

        twin = self._twin(ax)
        assert twin.lines[0].get_xdata().tolist() == [50, 200, 350]
        scatter = [c for c in ax.collections if isinstance(c, PathCollection)][0]
        assert _offsets(scatter)[:, 0].tolist() == [0, 100, 200, 300, 400]

    def test_only_rows_inside_the_region_and_chromosome_are_drawn(self):
        recombination = pd.DataFrame(
            {
                "CHR": ["1", "1", "1", "2"],
                "BP": [50, 300, 900, 300],
                "RATE": [9.0, 5.0, 8.0, 7.0],
            }
        )

        ax = regional_plot(self._data(), chr="1", start=100, end=500, recombination_rate=recombination)

        twin = self._twin(ax)
        assert twin.lines[0].get_xdata().tolist() == [200]
        assert twin.lines[0].get_ydata().tolist() == [5.0]

    def test_no_secondary_axis_without_data(self):
        ax = regional_plot(self._data(), chr="1", start=0, end=600)

        assert len(ax.figure.axes) == 1  # nothing is invented when no recombination map is supplied

    def test_no_secondary_axis_when_no_rate_falls_in_the_region(self):
        recombination = pd.DataFrame({"BP": [10_000, 20_000], "RATE": [1.0, 2.0]})

        ax = regional_plot(self._data(), chr="1", start=0, end=600, recombination_rate=recombination)

        assert len(ax.figure.axes) == 1

    def test_rates_without_the_required_columns_are_rejected(self):
        with pytest.raises(ValueError, match="RATE"):
            regional_plot(
                self._data(), chr="1", start=0, end=600, recombination_rate=pd.DataFrame({"BP": [1, 2], "cM": [1, 2]})
            )


# ---------------------------------------------------------------------------
# Dashboards
# ---------------------------------------------------------------------------


class TestMultiPanel:
    @staticmethod
    def _draw_line(ax, ys=(0, 1, 4)):
        ax.plot(range(len(ys)), list(ys))
        return ax

    def test_single_panel_with_default_columns(self):
        fig = multi_panel([{"func": self._draw_line, "title": "only"}])

        visible = [ax for ax in fig.axes if ax.get_visible()]
        assert len(visible) == 1
        assert visible[0].get_title() == "only"
        assert visible[0].lines[0].get_ydata().tolist() == [0, 1, 4]

    def test_single_panel_in_one_column(self):
        fig = multi_panel([{"func": self._draw_line}], ncols=1)

        assert [ax.get_visible() for ax in fig.axes] == [True]

    def test_three_panels_leave_one_hidden_axis(self):
        panels = [{"func": self._draw_line, "kwargs": {"ys": (0, i + 1)}} for i in range(3)]

        fig = multi_panel(panels, ncols=2)

        assert [ax.get_visible() for ax in fig.axes] == [True, True, True, False]
        assert [ax.lines[0].get_ydata().tolist() for ax in fig.axes[:3]] == [[0, 1], [0, 2], [0, 3]]
