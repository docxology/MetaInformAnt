"""Regression tests for dna.population plotting defects fixed in round 3.

* ``plot_population_diversity`` ignored its documented ``output_path`` alias.
* ``plot_pca_results`` returned ``None`` for the bare coordinate array its docstring promises to accept.
* ``plot_summary_statistics_grid`` / ``plot_neutrality_test_summary`` silently dropped every entry beyond
  the sixth / fourth.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402 - Agg backend must be forced before pyplot is imported
import numpy as np  # noqa: E402
import pytest  # noqa: E402

from metainformant.dna.population import visualization_core as vizcore, visualization_stats as vizstats  # noqa: E402


class TestPopulationDiversityOutputPath:
    def test_output_path_alias_writes_the_figure(self, tmp_path):
        out = tmp_path / "diversity_alias.png"

        fig = vizcore.plot_population_diversity({"pop1": 0.01, "pop2": 0.02}, output_path=str(out))

        assert fig is not None
        assert out.exists() and out.stat().st_size > 0
        assert len(plt.get_fignums()) == 0  # the figure is closed, not leaked into pyplot

    def test_output_file_still_writes_the_figure(self, tmp_path):
        out = tmp_path / "diversity_file.png"

        vizcore.plot_population_diversity({"pop1": 0.01, "pop2": 0.02}, output_file=str(out))

        assert out.exists() and out.stat().st_size > 0

    def test_no_file_is_written_without_a_path(self, tmp_path):
        fig = vizcore.plot_population_diversity({"pop1": 0.01})

        assert fig is not None
        assert list(tmp_path.iterdir()) == []

    def test_bars_carry_the_supplied_values(self, tmp_path):
        fig = vizcore.plot_population_diversity({"pop1": 0.01, "pop2": 0.02}, output_path=str(tmp_path / "d.png"))

        heights = [patch.get_height() for patch in fig.axes[0].patches]
        assert heights == pytest.approx([0.01, 0.02])


class TestPcaResultsBareArray:
    COORDS = np.array([[0.0, 0.0], [1.0, 1.0], [2.0, 0.5], [0.5, 2.0]])

    def test_bare_array_is_plotted_as_a_pc1_pc2_scatter(self):
        fig = vizstats.plot_pca_results(self.COORDS)

        assert fig is not None
        ax = fig.axes[0]
        np.testing.assert_allclose(ax.collections[0].get_offsets(), self.COORDS)
        assert (ax.get_xlabel(), ax.get_ylabel()) == ("PC1", "PC2")
        assert len(plt.get_fignums()) == 0

    def test_only_the_first_two_components_are_drawn(self):
        coords = np.arange(18.0).reshape(6, 3)

        fig = vizstats.plot_pca_results(coords)

        assert fig is not None
        np.testing.assert_allclose(fig.axes[0].collections[0].get_offsets(), coords[:, :2])

    def test_list_of_rows_is_accepted_like_an_array(self):
        fig = vizstats.plot_pca_results(self.COORDS.tolist())

        assert fig is not None
        np.testing.assert_allclose(fig.axes[0].collections[0].get_offsets(), self.COORDS)

    def test_bare_array_is_saved_when_a_path_is_given(self, tmp_path):
        out = tmp_path / "pca.png"

        fig = vizstats.plot_pca_results(self.COORDS, output_path=str(out))

        assert fig is not None
        assert out.exists() and out.stat().st_size > 0

    def test_sample_names_annotate_a_bare_array(self):
        fig = vizstats.plot_pca_results(self.COORDS, sample_names=["a", "b", "c", "d"])

        assert [text.get_text() for text in fig.axes[0].texts] == ["a", "b", "c", "d"]

    @pytest.mark.parametrize("unusable", [np.arange(5.0), np.zeros((4, 1)), np.zeros((0, 2)), [[1.0, 2.0], [3.0]]])
    def test_unusable_arrays_still_return_none(self, unusable):
        assert vizstats.plot_pca_results(unusable) is None
        assert len(plt.get_fignums()) == 0


def _stats(n_metrics: int) -> dict[str, dict[str, float]]:
    return {
        "popA": {f"stat_{j}": float(j + 1) for j in range(n_metrics)},
        "popB": {f"stat_{j}": float(10 + j) for j in range(n_metrics)},
    }


class TestSummaryStatisticsGrid:
    @staticmethod
    def _drawn(fig) -> dict[str, list[float]]:
        return {ax.get_title(): [patch.get_height() for patch in ax.patches] for ax in fig.axes if ax.axison}

    def test_every_metric_beyond_the_sixth_is_plotted(self):
        fig = vizstats.plot_summary_statistics_grid(_stats(8))

        drawn = self._drawn(fig)
        assert sorted(drawn) == sorted(f"Stat {j}" for j in range(8))
        for j in range(8):
            assert drawn[f"Stat {j}"] == pytest.approx([j + 1, 10 + j])

    def test_seven_metrics_use_a_fourth_row_with_one_unused_axis(self):
        fig = vizstats.plot_summary_statistics_grid(_stats(7))

        assert len(fig.axes) == 8
        assert sum(1 for ax in fig.axes if ax.axison) == 7

    def test_six_or_fewer_metrics_keep_the_three_by_two_layout(self):
        for n_metrics in (2, 6):
            fig = vizstats.plot_summary_statistics_grid(_stats(n_metrics))
            assert len(fig.axes) == 6
            assert sum(1 for ax in fig.axes if ax.axison) == n_metrics

    def test_taller_grid_keeps_the_per_row_figure_height(self):
        height_3_rows = vizstats.plot_summary_statistics_grid(_stats(6)).get_size_inches()[1]
        height_5_rows = vizstats.plot_summary_statistics_grid(_stats(10)).get_size_inches()[1]

        assert height_5_rows == pytest.approx(height_3_rows * 5 / 3)

    def test_grid_is_saved(self, tmp_path):
        out = tmp_path / "grid.png"

        vizstats.plot_summary_statistics_grid(_stats(8), output_path=str(out))

        assert out.exists() and out.stat().st_size > 0


class TestNeutralityTestSummary:
    @staticmethod
    def _summary(n_tests: int) -> dict[str, dict[str, float]]:
        return {f"test_{i}": {"statistic": float(i) - 3.0} for i in range(n_tests)}

    @staticmethod
    def _drawn(fig) -> dict[str, float]:
        return {ax.get_title(): ax.patches[0].get_height() for ax in fig.axes if ax.axison}

    def test_every_test_beyond_the_fourth_is_plotted(self):
        fig = vizstats.plot_neutrality_test_summary(self._summary(7))

        drawn = self._drawn(fig)
        assert sorted(drawn) == sorted(f"test_{i}" for i in range(7))
        for i in range(7):
            assert drawn[f"test_{i}"] == pytest.approx(float(i) - 3.0)

    def test_five_tests_use_a_third_row(self):
        fig = vizstats.plot_neutrality_test_summary(self._summary(5))

        assert len(fig.axes) == 6
        assert sum(1 for ax in fig.axes if ax.axison) == 5

    def test_four_or_fewer_tests_keep_the_two_by_two_layout(self):
        for n_tests in (1, 4):
            fig = vizstats.plot_neutrality_test_summary(self._summary(n_tests))
            assert len(fig.axes) == 4
            assert sum(1 for ax in fig.axes if ax.axison) == n_tests

    def test_p_values_still_appear_in_the_titles(self):
        summary = {f"t{i}": {"statistic": 1.0, "p_value": 0.01234} for i in range(5)}

        fig = vizstats.plot_neutrality_test_summary(summary)

        titles = sorted(ax.get_title() for ax in fig.axes if ax.axison)
        assert titles == sorted(f"t{i} (p=0.0123)" for i in range(5))
