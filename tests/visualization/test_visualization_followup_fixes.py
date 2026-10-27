"""Regression tests for the follow-up visualization fixes.

* The matplotlib Sankey fallback ordered its rows from ``set`` iteration, so the figure changed with the
  interpreter's hash seed.
* The circular Manhattan plot silently dropped every chromosome whose variants all share one position.
* The hierarchical network layout measured levels from the first root only, so every other root (and its
  descendants) was pushed into the last column.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

from matplotlib.collections import PathCollection  # noqa: E402

from metainformant.visualization.genomics.genomics import circular_manhattan_plot  # noqa: E402
from metainformant.visualization.plots import specialized  # noqa: E402
from metainformant.visualization.plots.specialized import plot_sankey_diagram  # noqa: E402


def _offsets(collection) -> np.ndarray:
    return np.asarray(collection.get_offsets(), dtype=float)


@pytest.mark.skipif(specialized.HAS_PLOTLY, reason="with Plotly installed the Sankey is a Plotly figure")
class TestSankeyFallbackOrder:
    """Rows follow the order in which nodes first appear in ``flows``, whatever the hash seed."""

    FLOWS = [
        ("zeta", "t_c", 4.0),
        ("alpha", "t_a", 3.0),
        ("mid", "t_c", 2.0),
        ("alpha", "t_b", 1.0),
        ("omega", "t_a", 5.0),
        ("beta", "t_d", 6.0),
        ("gamma", "t_b", 7.0),
        ("delta", "t_e", 8.0),
    ]
    SOURCES = ["zeta", "alpha", "mid", "omega", "beta", "gamma", "delta"]
    TARGETS = ["t_c", "t_a", "t_b", "t_d", "t_e"]

    def test_row_labels_are_in_first_seen_order(self):
        ax = plot_sankey_diagram(self.FLOWS).axes[0]

        assert [text.get_text() for text in ax.get_yticklabels()] == self.SOURCES
        assert list(ax.get_yticks()) == list(range(len(self.SOURCES)))

    def test_every_flow_runs_from_the_row_of_its_source_to_the_row_of_its_target(self):
        ax = plot_sankey_diagram(self.FLOWS).axes[0]

        markers = [c for c in ax.collections if isinstance(c, PathCollection)]
        assert len(markers) == 2 * len(self.FLOWS)  # a source marker then a target marker per flow
        for k, (source, target, _) in enumerate(self.FLOWS):
            assert _offsets(markers[2 * k])[0].tolist() == [0, self.SOURCES.index(source)]
            assert _offsets(markers[2 * k + 1])[0].tolist() == [1, self.TARGETS.index(target)]

    def test_two_calls_draw_the_same_figure(self):
        first = plot_sankey_diagram(self.FLOWS).axes[0]
        second = plot_sankey_diagram(self.FLOWS).axes[0]

        assert [t.get_text() for t in first.get_yticklabels()] == [t.get_text() for t in second.get_yticklabels()]


class TestCircularManhattanSinglePositionChromosomes:
    """A chromosome whose variants share one position used to be skipped (its position range was zero)."""

    @staticmethod
    def _frame() -> pd.DataFrame:
        return pd.DataFrame(
            {
                "CHR": [1, 1, 2, 3, 3, 3],
                "BP": [10, 20, 500, 700, 700, 700],
                "P": [1e-3, 1e-2, 1e-4, 1e-5, 1e-6, 0.5],
            }
        )

    def test_every_chromosome_is_drawn(self):
        ax = circular_manhattan_plot(self._frame())

        assert [len(_offsets(c)) for c in ax.collections] == [2, 1, 3]

    def test_single_position_chromosomes_sit_in_the_middle_of_their_sector(self):
        ax = circular_manhattan_plot(self._frame())

        sector = 2 * np.pi / 3
        theta = [_offsets(c)[:, 0] for c in ax.collections]
        assert theta[0] == pytest.approx([0.0, sector])  # chromosome 1 still spans its sector
        assert theta[1] == pytest.approx([1.5 * sector])  # chromosome 2: one variant
        assert theta[2] == pytest.approx([2.5 * sector] * 3)  # chromosome 3: three variants at one position

    def test_radii_are_minus_log10_p(self):
        ax = circular_manhattan_plot(self._frame())

        radii = [_offsets(c)[:, 1] for c in ax.collections]
        assert radii[0] == pytest.approx([3.0, 2.0])
        assert radii[1] == pytest.approx([4.0])
        assert radii[2] == pytest.approx([5.0, 6.0, -np.log10(0.5)])

    def test_size_and_alpha_still_apply_to_the_newly_drawn_chromosomes(self):
        ax = circular_manhattan_plot(self._frame(), s=40, alpha=0.3)

        assert [c.get_sizes().tolist() for c in ax.collections] == [[40.0]] * 3
        assert [c.get_alpha() for c in ax.collections] == [pytest.approx(0.3)] * 3


class TestHierarchicalNetworkMultipleRoots:
    """Levels are the distance from the nearest root, so every root starts the first column."""

    @staticmethod
    def _columns(graph) -> dict:
        from metainformant.visualization.genomics.networks import plot_network_hierarchical

        ax = plot_network_hierarchical(graph)
        return dict(zip(graph.nodes(), _offsets(ax.collections[0])[:, 0]))

    def test_every_root_is_in_the_first_column(self):
        nx = pytest.importorskip("networkx")
        graph = nx.DiGraph([("r1", "a"), ("r1", "b"), ("r2", "c"), ("c", "d")])

        x = self._columns(graph)

        assert x["r1"] == pytest.approx(x["r2"])
        assert x["a"] == pytest.approx(x["b"])
        assert x["b"] == pytest.approx(x["c"])
        assert x["r1"] < x["a"] < x["d"]

    def test_a_node_below_two_roots_takes_the_nearest_level(self):
        nx = pytest.importorskip("networkx")
        # n is two steps from r1 (r1 -> m -> n) but one step from r2, so it sits with the first-level nodes
        graph = nx.DiGraph([("r1", "m"), ("m", "n"), ("r2", "n")])

        x = self._columns(graph)

        assert x["r1"] == pytest.approx(x["r2"])
        assert x["m"] == pytest.approx(x["n"])
        assert x["r1"] < x["m"]

    def test_nodes_no_root_reaches_still_go_to_the_last_column(self):
        nx = pytest.importorskip("networkx")
        graph = nx.DiGraph([("r1", "a"), ("r2", "b")])
        graph.add_edge("p", "q")
        graph.add_edge("q", "p")  # a cycle without a root of its own

        x = self._columns(graph)

        assert x["r1"] == pytest.approx(x["r2"])
        assert x["a"] == pytest.approx(x["b"])
        assert x["p"] == pytest.approx(x["q"])
        assert x["r1"] < x["a"] < x["p"]

    def test_single_root_layout_is_unchanged(self):
        nx = pytest.importorskip("networkx")
        graph = nx.DiGraph([("r", "a"), ("r", "b"), ("a", "c"), ("a", "d"), ("b", "e")])

        x = self._columns(graph)

        assert x["a"] == pytest.approx(x["b"])
        assert x["c"] == pytest.approx(x["d"]) == pytest.approx(x["e"])
        assert x["r"] < x["a"] < x["c"]
