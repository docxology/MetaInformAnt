"""Regression tests: one failing plot of ``analyze_life_course`` must not drop the others.

The three summary plots (sequence lengths, domains, outcomes) sat in a single ``try/except Exception``.
The first failure silently dropped every plot after it, and it was only logged: ``results`` showed an
empty or short ``visualizations`` list and nothing said why. Every plot is now handled on its own and
failures are recorded in ``results["visualization_errors"]``.

A failure is provoked without test doubles: a directory squatting on one plot's file name makes writing
that single file fail with ``IsADirectoryError``.
"""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, List

import pytest

from metainformant.core import io
from metainformant.life_events.core.events import Event, EventSequence
from metainformant.life_events.workflow.workflow import _generate_analysis_visualizations, analyze_life_course

# Plot type -> file name, in the order the plots are made.
PLOT_FILES: Dict[str, str] = {
    "sequence_length_distribution": "sequence_lengths.png",
    "domain_distribution": "domain_distribution.png",
    "outcome_distribution": "outcome_distribution.png",
}
OUTCOMES = [0, 1, 0, 1]


def make_sequences() -> List[EventSequence]:
    plans = [["x", "y"], ["x", "z", "y"], ["z"], ["x", "y", "z"]]
    return [
        EventSequence(f"p{i}", [Event(name, datetime(2000 + j, 1, 1), "life") for j, name in enumerate(plan)])
        for i, plan in enumerate(plans)
    ]


def squat_on(directory: Path, *plot_types: str) -> None:
    for plot_type in plot_types:
        (directory / PLOT_FILES[plot_type]).mkdir()


@pytest.mark.parametrize("broken", list(PLOT_FILES))
def test_a_failing_plot_does_not_drop_the_other_plots(tmp_path: Path, broken: str):
    squat_on(tmp_path, broken)

    visualizations = _generate_analysis_visualizations(make_sequences(), OUTCOMES, tmp_path, {})

    assert [entry["type"] for entry in visualizations] == [plot for plot in PLOT_FILES if plot != broken]
    for entry in visualizations:
        assert Path(entry["file"]).is_file() and Path(entry["file"]).stat().st_size > 0


@pytest.mark.parametrize("broken", list(PLOT_FILES))
def test_a_failing_plot_is_recorded_in_the_results(tmp_path: Path, broken: str):
    squat_on(tmp_path, broken)
    results: Dict = {}

    _generate_analysis_visualizations(make_sequences(), OUTCOMES, tmp_path, results)

    assert len(results["visualization_errors"]) == 1
    error = results["visualization_errors"][0]
    assert error["type"] == broken
    assert error["error_type"] == "IsADirectoryError"
    assert PLOT_FILES[broken] in error["error"]


def test_the_failure_is_still_logged(tmp_path: Path, caplog: pytest.LogCaptureFixture):
    squat_on(tmp_path, "domain_distribution")

    with caplog.at_level(logging.WARNING):
        _generate_analysis_visualizations(make_sequences(), OUTCOMES, tmp_path, {})

    assert "domain_distribution" in caplog.text and "IsADirectoryError" in caplog.text


def test_every_failing_plot_is_recorded_and_nothing_raises(tmp_path: Path):
    squat_on(tmp_path, *PLOT_FILES)
    results: Dict = {}

    visualizations = _generate_analysis_visualizations(make_sequences(), OUTCOMES, tmp_path, results)

    assert visualizations == []
    assert [error["type"] for error in results["visualization_errors"]] == list(PLOT_FILES)


def test_no_errors_are_recorded_when_every_plot_succeeds(tmp_path: Path):
    results: Dict = {}

    visualizations = _generate_analysis_visualizations(make_sequences(), OUTCOMES, tmp_path, results)

    assert [entry["type"] for entry in visualizations] == list(PLOT_FILES)
    assert results["visualization_errors"] == []


def test_a_plot_that_is_not_attempted_is_not_an_error(tmp_path: Path):
    results: Dict = {}

    visualizations = _generate_analysis_visualizations(make_sequences(), None, tmp_path, results)

    assert [entry["type"] for entry in visualizations] == ["sequence_length_distribution", "domain_distribution"]
    assert results["visualization_errors"] == []


def test_analyze_life_course_keeps_going_and_reports_a_failed_plot(tmp_path: Path):
    squat_on(tmp_path, "sequence_length_distribution")

    results = analyze_life_course(make_sequences(), outcomes=OUTCOMES, output_dir=tmp_path, embedding_dim=8, epochs=1)

    assert [entry["type"] for entry in results["visualizations"]] == ["domain_distribution", "outcome_distribution"]
    assert [error["type"] for error in results["visualization_errors"]] == ["sequence_length_distribution"]
    # The rest of the analysis is unaffected by the failed plot.
    assert len(results["predictions"]) == len(OUTCOMES)
    assert Path(results["embeddings"]).is_file()
    # The record is part of the results file, not only of the returned dictionary.
    saved = io.load_json(tmp_path / "analysis_results.json")
    assert saved["visualization_errors"] == results["visualization_errors"]
    assert saved["visualization_errors"][0]["error_type"] == "IsADirectoryError"


def test_analyze_life_course_reports_no_plot_errors_on_success(tmp_path: Path):
    results = analyze_life_course(make_sequences(), outcomes=OUTCOMES, output_dir=tmp_path, embedding_dim=8, epochs=1)

    assert results["visualization_errors"] == []
    assert [entry["type"] for entry in results["visualizations"]] == list(PLOT_FILES)


def test_an_exception_raised_by_the_plot_function_itself_is_recorded(tmp_path: Path):
    """Not only file errors: the outcome plot cannot draw text labels and raises from inside the plot function."""
    results: Dict = {}

    visualizations = _generate_analysis_visualizations(make_sequences(), ["a", "b", "a", "b"], tmp_path, results)

    produced = {entry["type"] for entry in visualizations}
    failed = {error["type"]: error for error in results["visualization_errors"]}
    assert failed["outcome_distribution"]["error_type"] in {"TypeError", "ValueError"}
    assert not (tmp_path / PLOT_FILES["outcome_distribution"]).exists()
    assert "domain_distribution" in produced  # the plots around the failing one are still made
    assert produced | set(failed) == set(PLOT_FILES)  # every plot is either written or recorded


def test_a_plot_that_returns_no_figure_is_recorded(tmp_path: Path):
    """A plot function that has nothing to draw returns None; that used to vanish without a trace.

    The length plot of an empty list of sequences is such a case (without matplotlib, every plot function
    returns None as well).
    """
    results: Dict = {}

    visualizations = _generate_analysis_visualizations([], None, tmp_path, results)

    assert visualizations == []
    assert [(error["type"], error["error_type"]) for error in results["visualization_errors"]] == [
        ("sequence_length_distribution", "NoFigure")
    ]
    assert "no figure" in results["visualization_errors"][0]["error"]
    assert not list(tmp_path.glob("*.png"))
