"""Regression tests for the analysis-plot defects found in the round-4 sweep of ``visualization``.

Every test reads the data that landed on the artists (scatter offsets, line data, bar patches, tick and legend
labels, axes counts and positions), so a figure that is merely returned does not pass. Real matplotlib Agg figures
and real data; no test doubles.

Defects covered (each test fails on the code before the fix):

* ``qq_plot`` used ``linspace(0.01, 0.99, n)`` as plotting positions, squeezing the tails of every sample;
* ``leverage_plot`` drew ``sum(x_i^2) / sum(x^2)`` (not leverage; it summed to 1 instead of the parameter count);
* ``biplot`` could not take its ``scaling_factor`` (it was also passed to ``scatter``, which rejects it);
* ``plot_forecast`` dropped confidence intervals without ``lower``/``upper`` silently and labelled every band
  "95% CI"; ``plot_trend_analysis`` called a flat series "decreasing";
* ``plot_seasonal_decomposition`` ignored ``ax`` and its documented ``**kwargs``;
* the multi-panel QC plots left empty panels for absent sections, ignored ``ax`` and (batch colours) assigned colours
  in hash order; ``plot_quality_metrics`` drew absent statistics as bars of 0 and gave GC bars a fixed width;
  ``plot_multiomics_quality_overview`` read a score of 0.0 as missing and drew a line chart on Cartesian axes;
* caller keyword arguments collided with hard-coded defaults in six sequencing/integrity plots; an empty
  ``plot_overrepresented_sequences`` ignored ``output_path``; ``plot_error_profiles`` hid zeros on its log axis.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path
from typing import Any, Iterator, List, Sequence

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

import metainformant
from metainformant.visualization.analysis import quality_assessment, quality_omics, quality_sequencing
from metainformant.visualization.analysis.dimred import biplot
from metainformant.visualization.analysis.statistical import leverage_plot, qq_plot
from metainformant.visualization.analysis.timeseries import (
    plot_forecast,
    plot_seasonal_decomposition,
    plot_trend_analysis,
)

matplotlib.use("Agg")


@pytest.fixture(autouse=True)
def _close_figures() -> Iterator[None]:
    yield
    plt.close("all")


def _offsets(ax: Any, index: int = 0) -> np.ndarray:
    return np.asarray(ax.collections[index].get_offsets(), dtype=float)


def run_under_hash_seeds(program: str, seeds: Sequence[int] = (1,)) -> List[str]:
    """Run *program* in a fresh interpreter per string-hash seed and return each one's stdout.

    Hash seed 1 puts ``{"zeta", "alpha", "mu"}`` and ``{"beta", "alpha", "gamma"}`` in an order that is neither
    first-seen nor sorted, so code that iterates a set of strings gives a different answer under it. (Importing the
    visualization package takes several seconds, so one seed is used.)

    The child imports the same ``metainformant`` as this process: its ``PYTHONPATH`` is pinned to the package root
    of the imported module, so a test run against a copy of the source never reaches the installed one.
    """
    package_root = str(Path(metainformant.__file__).resolve().parent.parent)
    outputs = []
    for seed in seeds:
        env = {**os.environ, "PYTHONHASHSEED": str(seed), "PYTHONPATH": package_root, "MPLBACKEND": "Agg"}
        done = subprocess.run([sys.executable, "-c", program], env=env, capture_output=True, text=True, check=True)
        outputs.append(done.stdout.strip())
    return outputs


# ---------------------------------------------------------------------------
# statistical: qq_plot, leverage_plot
# ---------------------------------------------------------------------------


def test_qq_plot_of_a_normal_sample_is_a_straight_line_into_the_tails() -> None:
    """``linspace(0.01, 0.99, n)`` paired the largest of 1000 observations with the 99th percentile (2.33 sd) instead
    of the 99.95th (3.29 sd), so a perfectly normal sample bent away from the line in both tails."""
    stats = pytest.importorskip("scipy.stats")
    n = 1000
    positions = (np.arange(1, n + 1) - 0.5) / n
    data = stats.norm.ppf(positions)  # exactly the quantiles of a normal sample

    ax = qq_plot(data)

    theoretical, observed = _offsets(ax).T
    np.testing.assert_allclose(theoretical, stats.norm.ppf(positions, loc=data.mean(), scale=data.std()))
    away_from_centre = np.abs(theoretical) > 1e-6
    ratio = observed[away_from_centre] / theoretical[away_from_centre]
    assert np.ptp(ratio) < 1e-9  # one slope (1 / sd) everywhere, tails included
    assert theoretical.max() == pytest.approx(3.29 * data.std(), rel=0.01)


def test_qq_plot_against_a_uniform_sample_is_straight() -> None:
    stats = pytest.importorskip("scipy.stats")
    n = 200
    data = (np.arange(1, n + 1) - 0.5) / n

    ax = qq_plot(data, distribution="uniform")

    theoretical, observed = _offsets(ax).T
    np.testing.assert_allclose(theoretical, stats.uniform.ppf(data, loc=data.min(), scale=np.ptp(data)))
    assert np.ptp(observed - theoretical) < 0.02  # a uniform sample lies on y = x


def test_qq_plot_refuses_constant_data_and_unknown_distributions_before_drawing() -> None:
    pytest.importorskip("scipy")
    with pytest.raises(ValueError, match="not constant"):
        qq_plot(np.full(20, 3.0))
    before = len(plt.get_fignums())
    with pytest.raises(ValueError, match="Unsupported distribution"):
        qq_plot(np.arange(10.0), distribution="cauchy")
    assert len(plt.get_fignums()) == before  # no figure was left behind


def test_leverage_plot_draws_the_hat_values_of_a_simple_regression() -> None:
    x = np.array([1.0, 2.0, 3.0, 4.0, 10.0])
    y = np.array([1.1, 1.9, 3.2, 3.9, 9.0])
    n = len(x)

    ax = leverage_plot(x, y)  # a vector is one feature

    leverage, studentized = _offsets(ax).T
    expected_leverage = 1 / n + (x - x.mean()) ** 2 / np.sum((x - x.mean()) ** 2)
    np.testing.assert_allclose(leverage, expected_leverage)
    assert leverage.sum() == pytest.approx(2.0)  # the number of fitted parameters, not 1
    slope, intercept = np.polyfit(x, y, 1)
    residuals = y - (slope * x + intercept)
    variance = np.sum(residuals**2) / (n - 2)
    np.testing.assert_allclose(studentized, residuals / np.sqrt(variance * (1 - expected_leverage)))
    assert ax.get_xlabel() == "Leverage"


def test_leverage_plot_of_several_features_matches_the_hat_matrix() -> None:
    rng = np.random.default_rng(5)
    features = rng.normal(size=(30, 3))
    y = features @ np.array([1.0, -2.0, 0.5]) + rng.normal(scale=0.3, size=30)
    design = np.column_stack([np.ones(30), features])

    ax = leverage_plot(features, y)

    leverage = _offsets(ax)[:, 0]
    np.testing.assert_allclose(leverage, np.diag(design @ np.linalg.pinv(design)), atol=1e-10)
    assert leverage.sum() == pytest.approx(4.0)  # intercept plus three features


def test_leverage_plot_needs_more_samples_than_parameters() -> None:
    with pytest.raises(ValueError, match="more samples"):
        leverage_plot(np.array([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]), np.arange(3.0))  # 3 samples, 3 parameters


# ---------------------------------------------------------------------------
# dimred: biplot
# ---------------------------------------------------------------------------


def test_biplot_scaling_factor_scales_the_feature_arrows_and_alpha_does_not_collide() -> None:
    decomposition = pytest.importorskip("sklearn.decomposition")
    rng = np.random.default_rng(3)
    frame = pd.DataFrame(rng.normal(size=(40, 3)), columns=["gc", "length", "depth"])
    model = decomposition.PCA(n_components=2).fit(frame)
    loadings = model.components_[:2].T

    scaled = biplot(frame, model, scaling_factor=2.0, alpha=0.3)

    labels = {text.get_text(): text.get_position() for text in scaled.texts}
    for index, name in enumerate(frame.columns):
        np.testing.assert_allclose(labels[name], loadings[index] * 2.0 * 1.1)
    assert scaled.collections[0].get_alpha() == pytest.approx(0.3)

    default = biplot(frame, model)
    scores = model.transform(frame.to_numpy())
    expected = scores.std() / loadings.std() * 0.8
    default_labels = {text.get_text(): text.get_position() for text in default.texts}
    np.testing.assert_allclose(default_labels["gc"], loadings[0] * expected * 1.1)


# ---------------------------------------------------------------------------
# timeseries
# ---------------------------------------------------------------------------


def _weekly_series(length: int = 70) -> pd.Series:
    index = pd.date_range("2024-01-01", periods=length, freq="D")
    return pd.Series(np.sin(2 * np.pi * np.arange(length) / 7) + np.linspace(0, 3, length), index=index)


def test_forecast_confidence_interval_is_not_labelled_with_a_level_nobody_gave() -> None:
    history = pd.Series(np.arange(10.0))
    forecast = pd.Series(np.arange(10.0, 15.0), index=range(10, 15))
    interval = pd.DataFrame({"lower": forecast - 1, "upper": forecast + 1})

    unlabelled = plot_forecast(history, forecast, confidence_intervals=interval)
    labelled = plot_forecast(history, forecast, confidence_intervals=interval, confidence_level=0.9)

    assert unlabelled.get_legend_handles_labels()[1] == ["Historical", "Forecast", "Confidence interval"]
    assert labelled.get_legend_handles_labels()[1] == ["Historical", "Forecast", "90% CI"]
    band = unlabelled.collections[0].get_paths()[0].vertices
    assert band[:, 1].min() == 9.0 and band[:, 1].max() == 15.0  # the supplied bounds are what is drawn


def test_forecast_interval_must_have_one_row_per_forecast_value_and_a_valid_level() -> None:
    history = pd.Series(np.arange(10.0))
    forecast = pd.Series(np.arange(10.0, 15.0))

    with pytest.raises(ValueError, match="rows"):
        plot_forecast(history, forecast, confidence_intervals=pd.DataFrame({"lower": [1.0], "upper": [2.0]}))
    with pytest.raises(ValueError, match="confidence_level"):
        plot_forecast(history, forecast, confidence_level=95)


@pytest.mark.parametrize(
    ("values", "word"),
    [
        (np.full(30, 4.0), "flat"),  # a slope of exactly 0.0 used to read "decreasing"
        (np.linspace(0.0, 5.0, 30), "increasing"),
        (np.linspace(5.0, 0.0, 30), "decreasing"),
    ],
)
def test_trend_direction_names_a_flat_series_flat(values: np.ndarray, word: str) -> None:
    series = pd.Series(values, index=pd.date_range("2024-01-01", periods=len(values)))

    ax = plot_trend_analysis(series)

    assert [text.get_text() for text in ax.texts][0].startswith(f"Trend: {word}")


def test_seasonal_decomposition_draws_its_four_panels_inside_the_axes_it_is_given() -> None:
    """``ax`` was ignored: four subplots were stacked over the whole figure and ``ax`` stayed an empty frame."""
    series = _weekly_series()
    figure, host = plt.subplots()

    first = plot_seasonal_decomposition(series, period=7, ax=host)

    panels = [axes for axes in figure.axes if axes is not host]
    assert len(panels) == 4 and first is panels[0] and not host.axison
    host_box = host.get_position()
    for panel in panels:
        box = panel.get_position()
        assert host_box.x0 - 1e-9 <= box.x0 and box.x1 <= host_box.x1 + 1e-9
        assert host_box.y0 - 1e-9 <= box.y0 and box.y1 <= host_box.y1 + 1e-9
    assert [panel.get_title() for panel in panels] == [
        "Original Time Series",
        "Trend Component",
        "Seasonal Component",
        "Residual Component",
    ]


def test_seasonal_decomposition_passes_its_keyword_arguments_to_statsmodels() -> None:
    seasonal = pytest.importorskip("statsmodels.tsa.seasonal")
    series = _weekly_series()

    ax = plot_seasonal_decomposition(series, period=7, two_sided=False)

    expected = seasonal.seasonal_decompose(series, period=7, two_sided=False).trend
    trend = ax.figure.axes[1].lines[0].get_ydata()
    np.testing.assert_allclose(trend, expected.to_numpy(), equal_nan=True)
    assert np.isnan(trend[:6]).all()  # a one-sided filter has no trend for the first period


# ---------------------------------------------------------------------------
# panel_grid
# ---------------------------------------------------------------------------


def test_panel_grid_is_compact_and_refuses_an_axes_without_a_grid_slot() -> None:
    from metainformant.visualization.config.figures import panel_grid  # new helper: imported here so the other tests

    # in this module still collect (and fail on their own) against code that predates it

    panels = panel_grid(None, 5, max_cols=3, figsize=(6, 4))

    assert len(panels) == 5 and len(panels[0].figure.axes) == 5  # no padding to the 2 x 3 grid
    free_axes = plt.figure().add_axes((0.1, 0.1, 0.8, 0.8))  # made with add_axes: it has no subplot spec
    with pytest.raises(ValueError, match="grid slot"):
        panel_grid(free_axes, 2, max_cols=2, figsize=(1, 1))
    with pytest.raises(ValueError, match="at least one"):
        panel_grid(None, 0, max_cols=2, figsize=(1, 1))


def test_panel_grid_inside_a_subfigure_axes() -> None:
    from metainformant.visualization.config.figures import panel_grid

    figure = plt.figure(figsize=(8, 4))
    left, _right = figure.subfigures(1, 2)
    host = left.subplots()

    panels = panel_grid(host, 3, max_cols=2, figsize=(1, 1))

    assert len(panels) == 3 and all(panel.figure is left for panel in panels) and not host.axison


# ---------------------------------------------------------------------------
# quality_assessment
# ---------------------------------------------------------------------------


def test_batch_effects_qc_draws_one_panel_per_section_and_no_empty_ones() -> None:
    data = {"batch_sizes": {"b1": 5, "b2": 7}, "correction_metrics": {"kBET": 0.3, "ASW": 0.5}}

    first = quality_assessment.plot_batch_effects_qc(data)

    panels = first.figure.axes
    assert [panel.get_title() for panel in panels] == ["Batch Size Distribution", "Batch Correction Effectiveness"]
    assert [len(panel.patches) for panel in panels] == [2, 2]  # every panel holds data, none is blank


def test_batch_effects_qc_refuses_unusable_input_instead_of_leaving_blank_panels() -> None:
    with pytest.raises(ValueError, match="none of the recognised sections"):
        quality_assessment.plot_batch_effects_qc({"unrelated": 1})
    with pytest.raises(ValueError, match="needs the keys"):
        quality_assessment.plot_batch_effects_qc({"batch_pca": {"pc1": [1.0], "pc2": [2.0]}})  # no "batches"
    with pytest.raises(ValueError, match="must be a dict"):
        quality_assessment.plot_batch_effects_qc({"batch_sizes": [1, 2]})


def test_batch_effects_qc_draws_onto_a_given_axes() -> None:
    figure, host = plt.subplots()

    first = quality_assessment.plot_batch_effects_qc(
        {"batch_sizes": {"b1": 5}, "silhouette_scores": {"b1": 0.4}, "batch_variance": {"PC1": 0.6, "PC2": 0.4}},
        ax=host,
    )

    assert not host.axison and len(figure.axes) == 4 and first is figure.axes[1]


def test_batch_pca_colours_follow_first_seen_batch_order_in_every_process() -> None:
    """A ``set`` of batch names iterates in hash order, so the colour of a batch changed from run to run."""
    program = (
        "from metainformant.visualization.analysis.quality_assessment import plot_batch_effects_qc\n"
        "batches = ['zeta', 'alpha', 'mu', 'zeta', 'alpha', 'mu']\n"
        "ax = plot_batch_effects_qc({'batch_pca': {'pc1': [0, 1, 2, 3, 4, 5], 'pc2': [5, 4, 3, 2, 1, 0], "
        "'batches': batches}})\n"
        "print([t.get_text() for t in ax.get_legend().get_texts()])\n"
    )

    (output,) = run_under_hash_seeds(program)

    assert output == "['zeta', 'alpha', 'mu']"


def test_error_profile_marks_the_zero_a_log_axis_cannot_show() -> None:
    ax = quality_assessment.plot_error_profiles(
        {"substitution": np.array([0.0, 0.001, 0.01]), "indel": np.array([1e-4] * 3)}
    )

    assert ax.get_yscale() == "log"
    assert "zero at axis floor" in ax.get_ylabel()
    clean = quality_assessment.plot_error_profiles({"substitution": np.array([0.001, 0.01])})
    assert "zero" not in clean.get_ylabel()  # no mark when nothing is zero


def test_data_integrity_metrics_accepts_caller_colour_and_alpha() -> None:
    ax = quality_assessment.plot_data_integrity_metrics(
        {"completeness": 0.95, "missing_rate": 0.2}, color="tab:purple", alpha=0.3
    )

    purple = matplotlib.colors.to_rgb("tab:purple")
    assert [tuple(patch.get_facecolor()[:3]) for patch in ax.patches] == [
        purple,
        purple,
    ]  # not the traffic-light colours
    assert [patch.get_alpha() for patch in ax.patches] == [0.3, 0.3]


# ---------------------------------------------------------------------------
# quality_omics
# ---------------------------------------------------------------------------


def test_vcf_quality_metrics_has_one_panel_per_section_and_a_labelled_titv_reference() -> None:
    data = {
        "allele_frequencies": [0.1, 0.2, 0.2, 0.4],
        "titv_by_qual": {"qualities": [20, 30], "titv_ratios": [1.9, 2.0]},
    }

    first = quality_omics.plot_vcf_quality_metrics(data)

    allele_panel, titv_panel = first.figure.axes
    assert (allele_panel.get_title(), titv_panel.get_title()) == (
        "Allele Frequency Spectrum",
        "Transition/Transversion Ratio",
    )
    reference = [line for line in titv_panel.lines if line.get_label().startswith("Expected")]
    assert [line.get_ydata()[0] for line in reference] == [2.1] and reference[0].get_label() == "Expected (2.1)"

    organism = quality_omics.plot_vcf_quality_metrics(data, expected_titv=3.0)
    assert [
        line.get_ydata()[0] for line in organism.figure.axes[1].lines if line.get_label().startswith("Expected")
    ] == [3.0]
    none = quality_omics.plot_vcf_quality_metrics(data, expected_titv=None)
    assert not [line for line in none.figure.axes[1].lines if line.get_label().startswith("Expected")]


def test_vcf_quality_metrics_refuses_a_section_it_cannot_draw() -> None:
    with pytest.raises(ValueError, match="needs the keys"):
        quality_omics.plot_vcf_quality_metrics({"qual_distribution": {"qualities": [1, 2]}})
    with pytest.raises(ValueError, match="none of the recognised"):
        quality_omics.plot_vcf_quality_metrics({})


def test_singlecell_qc_leaves_out_cells_without_counts_from_the_complexity_panel() -> None:
    """Genes per UMI divided by zero counts (a RuntimeWarning, which this suite turns into an error)."""
    metrics = {"n_counts": np.array([0.0, 100.0, 200.0, 400.0]), "n_genes": np.array([0.0, 50.0, 80.0, 120.0])}

    first = quality_omics.plot_singlecell_qc_metrics(metrics)

    titles = [panel.get_title() for panel in first.figure.axes if panel.get_title()]  # the heatmap adds a colour bar
    assert titles == [
        "Library Size Distribution",
        "Genes Detected per Cell",
        "Library Complexity (1 cells with no counts left out)",
        "QC Metric Correlations",
    ]
    complexity = first.figure.axes[2]
    assert sum(patch.get_height() for patch in complexity.patches) == 3  # three cells have a ratio


def test_singlecell_qc_with_one_metric_has_one_panel() -> None:
    first = quality_omics.plot_singlecell_qc_metrics({"percent_mito": np.array([1.0, 2.0, 3.0])})

    assert len(first.figure.axes) == 1
    with pytest.raises(ValueError, match="none of the supported"):
        quality_omics.plot_singlecell_qc_metrics({"other": np.arange(3.0)})


def test_protein_structure_quality_draws_only_the_present_sections() -> None:
    first = quality_omics.plot_protein_structure_quality(
        {"ramachandran_stats": {"favoured": 0.9, "outliers": 0.1}, "overall_quality": 0.88}
    )

    ramachandran, overall = first.figure.axes
    assert ramachandran.get_title() == "Ramachandran Statistics"  # a bar chart of regions, not a phi/psi plot
    assert [text.get_text() for text in overall.texts] == ["Overall Quality: 0.88"]
    assert len(first.figure.axes) == 2


def test_multiomics_overview_is_a_radar_chart_and_a_zero_score_is_a_score() -> None:
    reports = {
        "Transcriptomics": {"overall_quality": 0.0, "completeness": 0.9},  # the worst quality, not a missing value
        "Proteomics": {"quality_score": 0.8},
        "Metabolomics": {"a": 0.5, "b": 0.7, "label": "ignored"},  # no score key: mean of the numeric values
        "Epigenomics": {},  # nothing to score
    }

    ax = quality_omics.plot_multiomics_quality_overview(reports)

    assert ax.name == "polar"
    closed = ax.lines[0].get_ydata()
    assert list(closed) == [0.0, 0.8, pytest.approx(0.6), 0.0, 0.0]  # four layers, the first repeated to close
    assert [text.get_text() for text in ax.texts] == ["0.00", "0.80", "0.60", "n/a"]
    _, cartesian = plt.subplots()
    with pytest.raises(ValueError, match="polar"):
        quality_omics.plot_multiomics_quality_overview(reports, ax=cartesian)
    with pytest.raises(ValueError, match="cannot be empty"):
        quality_omics.plot_multiomics_quality_overview({})


# ---------------------------------------------------------------------------
# quality_sequencing
# ---------------------------------------------------------------------------


def test_quality_metrics_does_not_draw_missing_statistics_as_zero() -> None:
    first = quality_sequencing.plot_quality_metrics({"basic_statistics": {"num_reads": 5000, "mean_length": 150.5}})

    (panel,) = first.figure.axes
    first.figure.canvas.draw()  # categorical tick labels are formatted at draw time
    assert [label.get_text() for label in panel.get_xticklabels()] == ["num_reads", "mean_length"]
    assert [patch.get_height() for patch in panel.patches] == [5000, 150.5]  # three more bars of 0 used to be drawn
    with pytest.raises(ValueError, match="none of the statistics"):
        quality_sequencing.plot_quality_metrics({"basic_statistics": {"unrelated": 1}})


def test_quality_metrics_gc_bars_are_as_wide_as_their_bins() -> None:
    """The bar width was fixed at 2 whatever the bins were, leaving gaps or overlapping bars."""
    gc = {"bins": [0, 10, 30, 100], "counts": [4, 9, 2]}

    first = quality_sequencing.plot_quality_metrics({"gc_content_distribution": gc})

    (panel,) = first.figure.axes
    assert [patch.get_width() for patch in panel.patches] == [10, 20, 70]
    assert [patch.get_x() + patch.get_width() / 2 for patch in panel.patches] == [5, 20, 65]
    with pytest.raises(ValueError, match="one count per bin"):
        quality_sequencing.plot_quality_metrics({"gc_content_distribution": {"bins": [0, 10, 30], "counts": [1, 2, 3]}})


def test_quality_metrics_draws_present_sections_in_a_compact_grid() -> None:
    qc = {
        "per_base_quality": {"positions": [1, 2, 3], "mean_qualities": [30.0, 31.0, 29.0]},
        "sequence_length_distribution": {"lengths": [100, 101], "counts": [5, 6]},
        "basic_statistics": {"num_reads": 11},
    }

    first = quality_sequencing.plot_quality_metrics(qc)

    assert [panel.get_title() for panel in first.figure.axes] == [
        "Per-Base Quality",
        "Read Length Distribution",
        "Basic Statistics",
    ]


@pytest.mark.parametrize(
    "call",
    [
        lambda: quality_sequencing.plot_adapter_content({"Illumina": [1.0, 4.0]}, alpha=0.25),
        lambda: quality_sequencing.plot_sequence_duplication_levels({"1": 80.0, "2": 12.0}, alpha=0.25),
        lambda: quality_sequencing.plot_kmer_profiles({"ACG": 5, "CGT": 3}, alpha=0.25),
        lambda: quality_sequencing.plot_overrepresented_sequences(
            [{"sequence": "ACGT", "percentage": 1.0}], alpha=0.25
        ),
        lambda: quality_sequencing.plot_gc_distribution([40.0, 50.0, 60.0], alpha=0.25, edgecolor="red"),
        lambda: quality_sequencing.plot_length_distribution([100, 120, 150], alpha=0.25, edgecolor="red"),
        lambda: quality_assessment.plot_data_integrity_metrics({"completeness": 0.9}, alpha=0.25),
    ],
    ids=["adapter", "duplication", "kmer", "overrepresented", "gc", "length", "integrity"],
)
def test_caller_keywords_override_the_hard_coded_bar_defaults(call: Any) -> None:
    """``alpha``/``edgecolor`` were also passed explicitly, so the call raised ``TypeError: multiple values``."""
    ax = call()

    assert ax.patches and {patch.get_alpha() for patch in ax.patches} == {0.25}


def test_per_base_quality_boxplot_accepts_its_own_widths() -> None:
    ax = quality_sequencing.plot_per_base_quality_boxplot(
        {"1-5": [30.0, 31.0, 29.0], "6-10": [28.0, 27.0, 30.0]}, widths=1
    )

    assert [np.ptp(box.get_path().vertices[:, 0]) for box in ax.patches] == [1.0, 1.0]  # widths=1, not the default 2


def test_overrepresented_sequences_writes_the_file_even_when_there_are_none(tmp_path: Path) -> None:
    """The empty case returned before saving, so the requested ``output_path`` was silently never written."""
    output = tmp_path / "overrepresented.png"

    ax = quality_sequencing.plot_overrepresented_sequences([], output_path=output)

    assert output.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    assert [text.get_text() for text in ax.texts] == ["No overrepresented sequences found"]


def test_overrepresented_sequences_truncates_only_what_it_cuts() -> None:
    long = "ACGT" * 8
    ax = quality_sequencing.plot_overrepresented_sequences(
        [{"sequence": "ACGTACGT", "percentage": 2.0, "count": 4}, {"sequence": long, "percentage": 1.0, "count": 2}]
    )

    assert [label.get_text() for label in ax.get_yticklabels()] == ["ACGTACGT", long[:20] + "..."]
