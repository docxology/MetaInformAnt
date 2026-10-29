"""Multi-omics plots draw the data they are given, label it with its own names and never invent a value.

Regression tests for ``metainformant.multiomics.visualization.visualization``:

* ``plot_integrated_omics_heatmap`` labelled the rows of every matrix after the first with the *first* feature names
  again (``feature_names[: data.shape[0]]`` per matrix), so with two matrices the second one's rows carried the first
  one's names;
* ``plot_pathway_enrichment_integration`` read a missing ``pvalue`` as ``1.0``, so a pathway that an omic reported
  without a p-value was drawn as "not enriched" (-log10 p = 0) instead of as missing;
* the network panel of ``create_interactive_multiomics_dashboard`` was a fixed three-point zig-zag
  (``x=[0, 1, 2], y=[0, 1, 0]``) for every ``network_data``.

The dashboard needs plotly, which is optional (the test for the figure is skipped without it); the network layout it
draws is a plain function that is tested directly. Real computation on the Agg backend; no test doubles
(real-implementation policy).
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pytest

import metainformant
from metainformant.multiomics.visualization.visualization import (
    _network_layout,
    create_interactive_multiomics_dashboard,
    plot_integrated_omics_heatmap,
    plot_multiomics_network,
    plot_pathway_enrichment_integration,
)

matplotlib.use("Agg")


@pytest.fixture(autouse=True)
def _close_figures() -> Any:
    plt.close("all")
    yield
    plt.close("all")


def _plotted_array(ax: Any) -> np.ma.MaskedArray:
    """The values a heatmap drew, from a seaborn mesh or a plain image."""
    artist = ax.collections[0] if ax.collections else ax.images[0]
    return np.ma.masked_invalid(np.ma.asarray(artist.get_array()))


class TestIntegratedHeatmapLabels:
    @staticmethod
    def _omics() -> dict[str, np.ndarray]:
        rng = np.random.default_rng(1)
        return {"rna": rng.normal(size=(3, 4)), "protein": rng.normal(size=(2, 4))}

    def test_every_matrix_gets_its_own_row_names(self) -> None:
        ax = plot_integrated_omics_heatmap(self._omics(), feature_names=["a", "b", "c", "d", "e"])

        assert [label.get_text() for label in ax.get_yticklabels()] == ["a", "b", "c", "d", "e"]  # was a b c a b

    def test_the_stacked_matrix_has_the_rows_of_all_matrices(self) -> None:
        omics = self._omics()

        ax = plot_integrated_omics_heatmap(omics, feature_names=["a", "b", "c", "d", "e"])

        assert _plotted_array(ax).shape == (5, 4)
        np.testing.assert_allclose(_plotted_array(ax)[3:], omics["protein"])

    @pytest.mark.parametrize("n_names", [2, 3, 4, 6])
    def test_a_name_list_of_the_wrong_length_raises(self, n_names: int) -> None:
        with pytest.raises(ValueError, match=rf"feature_names has {n_names} entries but the stacked matrix has 5 rows"):
            plot_integrated_omics_heatmap(self._omics(), feature_names=[f"f{i}" for i in range(n_names)])

    def test_the_error_lists_the_rows_of_each_matrix(self) -> None:
        with pytest.raises(ValueError, match=r"rna: 3, protein: 2"):
            plot_integrated_omics_heatmap(self._omics(), feature_names=["only one"])

    def test_without_names_nothing_is_checked_or_invented(self) -> None:
        ax = plot_integrated_omics_heatmap(self._omics())

        assert _plotted_array(ax).shape == (5, 4)


class TestPathwayEnrichmentHeatmap:
    RESULTS = {
        "rna": [{"pathway": "glycolysis", "pvalue": 1e-3}, {"pathway": "tca", "pvalue": 0.5}],
        "protein": [{"pathway": "glycolysis"}, {"term": "urea", "pvalue": 0.0}],  # no p-value for glycolysis
    }

    def _drawn(self) -> tuple[dict[tuple[str, str], float], np.ma.MaskedArray]:
        ax = plot_pathway_enrichment_integration(self.RESULTS)
        values = _plotted_array(ax)
        rows = [label.get_text() for label in ax.get_yticklabels()]
        columns = [label.get_text() for label in ax.get_xticklabels()]
        cells = {
            (rows[i], columns[j]): float(values[i, j])
            for i in range(len(rows))
            for j in range(len(columns))
            if not values.mask[i, j]
        }
        return cells, values

    def test_a_pathway_reported_without_a_p_value_is_blank_not_p_equal_one(self) -> None:
        cells, _ = self._drawn()

        assert ("glycolysis", "protein") not in cells  # was drawn as -log10(1.0) = 0.0
        assert cells[("glycolysis", "rna")] == pytest.approx(3.0)

    def test_values_are_minus_log10_p_and_a_zero_p_value_is_capped(self) -> None:
        cells, _ = self._drawn()

        assert cells[("tca", "rna")] == pytest.approx(-np.log10(0.5))
        assert cells[("urea", "protein")] == pytest.approx(300.0)

    def test_pathways_an_omic_did_not_report_are_blank(self) -> None:
        cells, values = self._drawn()

        assert ("urea", "rna") not in cells and ("tca", "protein") not in cells
        # 3 pathways x 2 omics = 6 cells; only rna/glycolysis, rna/tca and protein/urea have a p-value
        assert len(cells) == 3 and int(values.mask.sum()) == 3

    def test_every_drawn_value_is_finite(self) -> None:
        cells, _ = self._drawn()

        assert all(np.isfinite(v) for v in cells.values())


class TestNetworkLayout:
    ADJACENCY = {"TP53": ["MDM2", "BAX"], "MDM2": ["TP53"], "BAX": [], "LONELY": []}

    def test_an_adjacency_mapping_gives_its_nodes_and_edges(self) -> None:
        nodes, positions, edges = _network_layout(self.ADJACENCY)

        assert sorted(nodes) == ["BAX", "LONELY", "MDM2", "TP53"]  # isolated nodes are drawn too
        assert sorted(tuple(sorted(edge)) for edge in edges) == [("BAX", "TP53"), ("MDM2", "TP53")]
        assert set(positions) == set(nodes)

    def test_a_networkx_graph_gives_the_same_network(self) -> None:
        graph = nx.Graph()
        graph.add_nodes_from(self.ADJACENCY)
        graph.add_edges_from([("TP53", "MDM2"), ("TP53", "BAX")])

        from_graph = _network_layout(graph)
        from_mapping = _network_layout(self.ADJACENCY)

        assert sorted(from_graph[0]) == sorted(from_mapping[0])
        assert sorted(tuple(sorted(e)) for e in from_graph[2]) == sorted(tuple(sorted(e)) for e in from_mapping[2])

    def test_the_layout_is_deterministic_and_finite(self) -> None:
        first = _network_layout(self.ADJACENCY)[1]
        second = _network_layout(self.ADJACENCY)[1]

        assert first == second
        assert all(np.isfinite(x) and np.isfinite(y) for x, y in first.values())

    def test_different_networks_are_laid_out_differently(self) -> None:
        # the old panel drew the same three points for every network
        small = _network_layout({"a": ["b"], "b": []})[1]
        larger = _network_layout({"a": ["b", "c", "d"], "b": [], "c": [], "d": []})[1]

        assert len(small) == 2 and len(larger) == 4

    @pytest.mark.parametrize("bad", [None, 7, "a-b", [("a", "b")]])
    def test_other_input_is_refused(self, bad: Any) -> None:
        with pytest.raises(ValueError, match=r"network_data must be a networkx graph or a"):
            _network_layout(bad)

    def test_an_empty_network_is_refused(self) -> None:
        with pytest.raises(ValueError, match=r"no nodes to draw"):
            _network_layout({})


class TestDashboardNetworkPanel:
    def test_the_network_traces_hold_the_network_itself(self) -> None:
        pytest.importorskip("plotly")
        adjacency = {"TP53": ["MDM2", "BAX"], "MDM2": ["TP53"], "BAX": []}

        figure = create_interactive_multiomics_dashboard({"network_data": adjacency})

        edges, nodes = figure.data[0], figure.data[1]
        assert len(edges.x) == 3 * 2 and edges.x[2] is None  # two edges, each as (x0, x1, gap)
        assert sorted(nodes.text) == ["BAX", "MDM2", "TP53"]
        assert list(nodes.x) != [0, 1, 2]  # the old fixed zig-zag

    def test_a_network_that_cannot_be_drawn_is_refused_not_replaced(self) -> None:
        pytest.importorskip("plotly")

        with pytest.raises(ValueError, match=r"network_data must be"):
            create_interactive_multiomics_dashboard({"network_data": 42})


class TestPlotsAreReproducible:
    """Set iteration order and the unseeded spring layout used to change the figures from run to run."""

    @staticmethod
    def _network() -> tuple[nx.Graph, dict[str, str]]:
        graph = nx.cycle_graph(8)
        graph.add_edges_from([(0, 4), (2, 6)])
        return graph, {n: ("rna" if n % 2 else "protein") for n in graph.nodes()}

    def test_the_layout_is_the_same_every_time(self) -> None:
        graph, omics = self._network()

        first = plot_multiomics_network(graph, omics).collections[0].get_offsets().data.copy()
        plt.close("all")
        second = plot_multiomics_network(graph, omics).collections[0].get_offsets().data.copy()

        np.testing.assert_array_equal(first, second)

    def test_the_legend_does_not_trigger_matplotlibs_single_rgb_colour_warning(self) -> None:
        import warnings

        graph, omics = self._network()

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            plot_multiomics_network(graph, omics)

        assert [str(w.message) for w in caught if "single numeric RGB" in str(w.message)] == []

    def test_another_layout_seed_gives_another_layout(self) -> None:
        graph, omics = self._network()

        default = plot_multiomics_network(graph, omics).collections[0].get_offsets().data.copy()
        plt.close("all")
        other = plot_multiomics_network(graph, omics, layout_seed=7).collections[0].get_offsets().data.copy()

        assert not np.allclose(default, other)

    def test_drawing_does_not_consume_the_global_numpy_random_state(self) -> None:
        graph, omics = self._network()
        np.random.seed(123)
        before = np.random.get_state()

        plot_multiomics_network(graph, omics)

        after = np.random.get_state()
        assert before[0] == after[0] and np.array_equal(before[1], after[1]) and before[2:] == after[2:]

    def test_the_pathways_of_the_enrichment_heatmap_come_in_sorted_order(self) -> None:
        results = {
            "rna": [{"pathway": name, "pvalue": 0.01} for name in ("zeta", "alpha", "mu")],
            "protein": [{"pathway": name, "pvalue": 0.02} for name in ("beta", "alpha")],
        }

        ax = plot_pathway_enrichment_integration(results)

        assert [label.get_text() for label in ax.get_yticklabels()] == ["alpha", "beta", "mu", "zeta"]

    def test_the_enrichment_heatmap_does_not_depend_on_string_hashing(self) -> None:
        program = r"""
import json
import matplotlib
matplotlib.use("Agg")
from metainformant.multiomics.visualization.visualization import plot_pathway_enrichment_integration
names = [f"PW{i:02d}" for i in range(15)]
results = {"rna": [{"pathway": n, "pvalue": 0.01 * (i + 1)} for i, n in enumerate(names)],
           "protein": [{"pathway": n, "pvalue": 0.5} for n in names[::2]]}
ax = plot_pathway_enrichment_integration(results)
print(json.dumps([t.get_text() for t in ax.get_yticklabels()]))
"""
        package_root = str(Path(metainformant.__file__).resolve().parents[1])  # the very package this process uses
        outputs = []
        for seed in ("0", "2", "3", "5"):
            env = {**os.environ, "PYTHONHASHSEED": seed, "MPLBACKEND": "Agg"}
            env["PYTHONPATH"] = os.pathsep.join(filter(None, [package_root, os.environ.get("PYTHONPATH")]))
            completed = subprocess.run(
                [sys.executable, "-c", program], env=env, capture_output=True, text=True, timeout=300
            )
            assert completed.returncode == 0, completed.stderr[-2000:]
            outputs.append(json.loads(completed.stdout))

        assert all(labels == outputs[0] for labels in outputs)
        assert outputs[0] == sorted(outputs[0])
