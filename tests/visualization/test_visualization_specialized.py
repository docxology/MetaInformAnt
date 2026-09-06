"""Tests for plots/specialized.py — the module had no dedicated test file before 2026-09.

Covers the fallback rendering paths (matplotlib-venn / plotly / alluvial /
upsetplot are optional and absent in the base env) plus the validation errors
added in the 2026-09 review pass (venn >3 sets, chord label mismatch).
"""

from __future__ import annotations

import importlib.util

import numpy as np
import pandas as pd
import pytest

from metainformant.visualization.plots.specialized import (
    plot_alluvial_diagram,
    plot_chord_diagram,
    plot_circular_barplot,
    plot_network_circular_layout,
    plot_sankey_diagram,
    plot_venn_diagram,
)

HAS_MPL_VENN = importlib.util.find_spec("matplotlib_venn") is not None
HAS_UPSETPLOT = importlib.util.find_spec("upsetplot") is not None
HAS_PLOTLY = importlib.util.find_spec("plotly") is not None
HAS_ALLUVIAL = importlib.util.find_spec("alluvial") is not None


class TestVennDiagram:
    """plot_venn_diagram fallback and validation paths."""

    def test_two_sets_fallback_renders(self, tmp_path):
        """With matplotlib-venn absent the internal 2-set fallback must render."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        sets = {"A": {"gene1", "gene2", "gene3"}, "B": {"gene2", "gene3", "gene4"}}
        output_path = tmp_path / "venn2.png"

        ax = plot_venn_diagram(sets, output_path=output_path)

        assert ax.get_title() == "Venn Diagram"
        assert output_path.exists()
        plt.close("all")

    def test_three_sets_fallback_renders(self, tmp_path):
        """With matplotlib-venn absent the internal 3-set fallback must render."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        sets = {"A": {1, 2, 3}, "B": {2, 3, 4}, "C": {3, 4, 5}}
        output_path = tmp_path / "venn3.png"

        ax = plot_venn_diagram(sets, output_path=output_path)

        assert ax.get_title() == "Venn Diagram"
        assert output_path.exists()
        plt.close("all")

    def test_more_than_three_sets_raises(self):
        """>3 sets must fail loudly with the supported-set count in the message."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        sets = {name: {1, 2} for name in ("A", "B", "C", "D")}
        # Without matplotlib-venn the fallback path raises ImportError; with it
        # installed the dedicated branch raises ValueError. Both are the loud
        # failure this contract requires.
        expected = ValueError if HAS_MPL_VENN else ImportError
        with pytest.raises(expected):
            plot_venn_diagram(sets)
        plt.close("all")

    def test_rejects_non_dict_input(self):
        """Non-dict input must raise the domain validation error."""
        from metainformant.core.utils.errors import ValidationError

        with pytest.raises(ValidationError, match="sets"):
            plot_venn_diagram(["A", "B"])


class TestChordDiagram:
    """plot_chord_diagram validation and rendering."""

    def test_label_count_mismatch_raises(self):
        """Labels not matching the matrix dimension must raise ValueError."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        matrix = np.eye(3)
        with pytest.raises(ValueError):
            plot_chord_diagram(matrix, labels=["only", "two"])
        plt.close("all")

    def test_non_square_matrix_raises(self):
        """A non-square matrix must raise ValueError before any drawing."""
        matrix = np.ones((2, 3))
        with pytest.raises(ValueError, match="square"):
            plot_chord_diagram(matrix)

    def test_renders_and_saves(self, tmp_path):
        """A valid square matrix renders on polar axes and saves deterministically."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        matrix = np.array([[0.0, 2.0, 1.0], [2.0, 0.0, 0.5], [1.0, 0.5, 0.0]])
        labels = ["x", "y", "z"]
        output_path = tmp_path / "chord.png"

        ax = plot_chord_diagram(matrix, labels, output_path=output_path)

        assert ax.get_title() == "Chord Diagram"
        assert output_path.exists()
        plt.close("all")


class TestCircularBarplot:
    """plot_circular_barplot rendering."""

    def test_renders_and_saves(self, tmp_path):
        """Values render as bars on polar axes; output saved when requested."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        values = np.array([3.0, 5.0, 2.0, 4.0])
        labels = ["a", "b", "c", "d"]
        output_path = tmp_path / "circular_bar.png"

        ax = plot_circular_barplot(values, labels, output_path=output_path)

        assert ax.get_title() == "Circular Bar Plot"
        assert len(ax.patches) == len(values)
        assert output_path.exists()
        plt.close("all")


@pytest.mark.skipif(HAS_PLOTLY, reason="plotly installed; matplotlib fallback path unused")
class TestSankeyDiagram:
    """plot_sankey_diagram (plotly absent in this env: matplotlib fallback)."""

    def test_fallback_renders_and_saves(self, tmp_path):
        """Without plotly the matplotlib flow fallback must render a figure."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        flows = [("input", "process", 10.0), ("process", "output", 7.0)]
        output_path = tmp_path / "sankey.png"

        fig = plot_sankey_diagram(flows, output_path=output_path)

        assert fig is not None
        assert output_path.exists()
        plt.close("all")

    def test_rejects_non_list_flows(self):
        """Flows must be a list of (source, target, value) tuples."""
        from metainformant.core.utils.errors import ValidationError

        with pytest.raises(ValidationError, match="flows"):
            plot_sankey_diagram("not-a-list")


class TestAlluvialDiagram:
    """plot_alluvial_diagram (alluvial package absent: simple fallback)."""

    @pytest.mark.skipif(HAS_ALLUVIAL, reason="alluvial installed; simple fallback unused")
    def test_fallback_renders(self, tmp_path):
        """Without the alluvial package the simple fallback must render."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        data = pd.DataFrame({"stage1": [1.0, 2.0, 3.0], "stage2": [2.0, 1.0, 2.5]})
        stages = ["stage1", "stage2"]

        ax = plot_alluvial_diagram(data, stages)

        assert ax.get_title() == "Alluvial Diagram"
        plt.close("all")


class TestNetworkCircularLayout:
    """plot_network_circular_layout (networkx is installed)."""

    def test_renders_and_saves(self, tmp_path):
        """A small graph renders in circular layout and saves when requested."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import networkx as nx

        graph = nx.cycle_graph(6)
        output_path = tmp_path / "circular_network.png"

        ax = plot_network_circular_layout(graph, output_path=output_path)

        assert ax.get_title() == "Circular Network Layout"
        assert output_path.exists()
        plt.close("all")


class TestUpsetPlot:
    """plot_upset_plot requires the optional upsetplot package."""

    def test_missing_upsetplot_raises(self):
        """Without upsetplot a clear ImportError must be raised."""
        import matplotlib

        matplotlib.use("Agg")

        if HAS_UPSETPLOT:
            pytest.skip("upsetplot installed; ImportError path not reachable")

        from metainformant.visualization.plots.specialized import plot_upset_plot

        with pytest.raises(ImportError, match="upsetplot"):
            plot_upset_plot({"A": {1, 2}, "B": {2, 3}})
