"""Regression tests for ``create_population_summary_plot``.

It drew a population that did not report a metric as a bar of height 0 (a measured zero), and it ordered the
panels by iterating a ``set`` of metric names, so the figure depended on the interpreter's hash seed.
"""

# ruff: noqa: E402 - Agg backend must be forced before pyplot-backed modules are imported.

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import pytest  # noqa: E402

from metainformant.dna.population.visualization_core import create_population_summary_plot  # noqa: E402

# Metric names deliberately not in alphabetical order (and enough of them that a lucky set order is
# vanishingly unlikely)
METRICS = ["watterson_theta", "tajima_d", "segregating_sites", "pi", "heterozygosity", "fst", "diversity", "ne"]


def _panel(fig, metric: str):
    title = f'{metric.replace("_", " ").title()} by Population'
    (axis,) = [ax for ax in fig.axes if ax.get_title() == title]
    return axis


class TestPanelOrder:
    def test_panels_are_ordered_alphabetically_by_metric(self):
        data = {
            "pop1": {m: float(i) for i, m in enumerate(METRICS)},
            "pop2": {m: float(i) for i, m in enumerate(METRICS)},
        }

        fig = create_population_summary_plot(data)

        assert [ax.get_title() for ax in fig.axes] == [
            f'{m.replace("_", " ").title()} by Population' for m in sorted(METRICS)
        ]

    def test_order_does_not_depend_on_which_population_reports_a_metric_first(self):
        first = {
            "pop1": {"zeta": 1.0, "alpha": 2.0, "delta": 4.0},
            "pop2": {"mid": 3.0, "beta": 5.0, "omega": 6.0},
        }
        second = {
            "pop2": {"omega": 6.0, "beta": 5.0, "mid": 3.0},
            "pop1": {"delta": 4.0, "alpha": 2.0, "zeta": 1.0},
        }

        titles_first = [ax.get_title() for ax in create_population_summary_plot(first).axes]
        titles_second = [ax.get_title() for ax in create_population_summary_plot(second).axes]

        expected = [f"{m} by Population" for m in ("Alpha", "Beta", "Delta", "Mid", "Omega", "Zeta")]
        assert titles_first == titles_second == expected

    def test_populations_keep_the_order_they_were_given_in(self):
        data = {"pop_b": {"diversity": 0.02}, "pop_a": {"diversity": 0.01}, "pop_c": {"diversity": 0.03}}

        fig = create_population_summary_plot(data)

        axis = fig.axes[0]
        assert [t.get_text() for t in axis.get_xticklabels()] == ["pop_b", "pop_a", "pop_c"]
        assert [p.get_height() for p in axis.patches] == pytest.approx([0.02, 0.01, 0.03])


class TestMissingMetrics:
    DATA = {
        "popA": {"diversity": 0.01, "fst": 0.05},
        "popB": {"diversity": 0.02},
        "popC": {"fst": 0.04, "diversity": 0.03},
    }

    def test_a_population_without_the_metric_has_no_bar(self):
        fig = create_population_summary_plot(self.DATA)

        fst = _panel(fig, "fst")
        assert [p.get_height() for p in fst.patches] == pytest.approx([0.05, 0.04])  # popA and popC only
        assert [p.get_x() + p.get_width() / 2 for p in fst.patches] == pytest.approx([0.0, 2.0])

    def test_the_missing_value_is_annotated_not_drawn_as_zero(self):
        fig = create_population_summary_plot(self.DATA)

        fst = _panel(fig, "fst")
        labels = {round(t.get_position()[0]): t.get_text() for t in fst.texts}
        assert labels[1] == "n/a"  # popB sits at x = 1
        assert labels[0] == "0.0500" and labels[2] == "0.0400"

    def test_every_population_stays_on_the_axis(self):
        fig = create_population_summary_plot(self.DATA)

        fst = _panel(fig, "fst")
        assert [t.get_text() for t in fst.get_xticklabels()] == ["popA", "popB", "popC"]

    def test_none_and_nan_values_are_missing_too(self):
        data = {"p1": {"fst": 0.1}, "p2": {"fst": None}, "p3": {"fst": float("nan")}, "p4": {"fst": 0.2}}

        fig = create_population_summary_plot(data)

        axis = fig.axes[0]
        assert [p.get_height() for p in axis.patches] == pytest.approx([0.1, 0.2])
        assert sorted(t.get_text() for t in axis.texts) == ["0.1000", "0.2000", "n/a", "n/a"]

    def test_complete_data_is_unchanged(self):
        data = {"pop1": {"diversity": 0.01, "fst": 0.05}, "pop2": {"diversity": 0.015, "fst": 0.03}}

        fig = create_population_summary_plot(data)

        diversity = _panel(fig, "diversity")
        assert [p.get_height() for p in diversity.patches] == pytest.approx([0.01, 0.015])
        assert sorted(t.get_text() for t in diversity.texts) == ["0.0100", "0.0150"]
        assert not any(t.get_text() == "n/a" for ax in fig.axes for t in ax.texts)

    def test_negative_values_are_real_values(self):
        data = {"pop1": {"tajima_d": -1.5}, "pop2": {"tajima_d": 0.5}}

        fig = create_population_summary_plot(data)

        axis = fig.axes[0]
        assert [p.get_height() for p in axis.patches] == pytest.approx([-1.5, 0.5])

    def test_figure_is_returned_closed_and_saved_when_a_path_is_given(self, tmp_path):
        out = tmp_path / "summary.png"

        fig = create_population_summary_plot(self.DATA, output_file=str(out))

        assert out.exists() and out.stat().st_size > 0
        assert len(plt.get_fignums()) == 0
        assert len(fig.axes) == 2  # diversity and fst

    def test_empty_inputs_still_return_none(self):
        assert create_population_summary_plot({}) is None
        assert create_population_summary_plot({"pop1": {}}) is None
