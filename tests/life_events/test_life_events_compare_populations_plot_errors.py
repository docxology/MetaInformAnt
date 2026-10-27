"""Regression tests: a failed comparison plot of ``compare_populations`` is recorded in the results.

Defect: the one plot of ``compare_populations`` sat in a ``try/except Exception`` that only logged the failure
("Error generating comparison visualizations: ..."). ``results["visualizations"]`` then came back empty with no
trace of why, unlike ``analyze_life_course``, which records every failed plot in ``visualization_errors``.
``compare_populations`` now does the same and still carries on to return the statistical comparison.

A failure is provoked without test doubles: a directory squatting on the plot's file name makes writing it fail
with ``IsADirectoryError``.
"""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path
from typing import List

import pytest

from metainformant.core import io
from metainformant.life_events.core.events import Event, EventSequence
from metainformant.life_events.workflow.workflow import compare_populations

PLOT_FILE = "population_comparison.png"


def _sequences(plans: List[List[str]], prefix: str) -> List[EventSequence]:
    return [
        EventSequence(f"{prefix}{i}", [Event(name, datetime(2000 + j, 1, 1), "life") for j, name in enumerate(plan)])
        for i, plan in enumerate(plans)
    ]


# Group one: lengths 2 and 3, event types {x, y, z}. Group two: lengths 1 and 1, event types {y, w}.
GROUP_ONE = _sequences([["x", "y"], ["x", "z", "y"]], "a")
GROUP_TWO = _sequences([["y"], ["w"]], "b")


def test_a_failing_plot_is_recorded_in_the_results(tmp_path: Path) -> None:
    (tmp_path / PLOT_FILE).mkdir()

    results = compare_populations(GROUP_ONE, GROUP_TWO, output_dir=tmp_path)

    assert results["visualizations"] == []
    assert len(results["visualization_errors"]) == 1
    error = results["visualization_errors"][0]
    assert error["type"] == "population_comparison"
    assert error["error_type"] == "IsADirectoryError"
    assert PLOT_FILE in error["error"]


def test_the_comparison_is_still_returned_when_the_plot_fails(tmp_path: Path) -> None:
    (tmp_path / PLOT_FILE).mkdir()

    results = compare_populations(GROUP_ONE, GROUP_TWO, group_names=("cases", "controls"), output_dir=tmp_path)

    comparison = results["comparison"]
    # Mean lengths (2 + 3) / 2 = 2.5 and (1 + 1) / 2 = 1; event types {x, y, z} and {y, w} share only y,
    # so the Jaccard similarity is 1 / 4.
    assert comparison["length_comparison"]["group1_mean"] == 2.5
    assert comparison["length_comparison"]["group2_mean"] == 1.0
    assert comparison["event_types"]["overlap"] == 1
    assert comparison["event_types"]["jaccard_similarity"] == pytest.approx(0.25)
    assert comparison["common_event_types"] == ["y"]
    assert (results["group1_name"], results["group2_name"]) == ("cases", "controls")


def test_the_failure_is_part_of_the_saved_results_file(tmp_path: Path) -> None:
    (tmp_path / PLOT_FILE).mkdir()

    results = compare_populations(GROUP_ONE, GROUP_TWO, output_dir=tmp_path)

    saved = io.load_json(tmp_path / "comparison_results.json")
    assert saved["visualization_errors"] == results["visualization_errors"]
    assert saved["visualization_errors"][0]["error_type"] == "IsADirectoryError"


def test_the_failure_is_still_logged(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    (tmp_path / PLOT_FILE).mkdir()

    with caplog.at_level(logging.WARNING):
        compare_populations(GROUP_ONE, GROUP_TWO, output_dir=tmp_path)

    assert "population_comparison" in caplog.text and "IsADirectoryError" in caplog.text


def test_no_error_is_recorded_when_the_plot_is_written(tmp_path: Path) -> None:
    results = compare_populations(GROUP_ONE, GROUP_TWO, output_dir=tmp_path)

    assert results["visualization_errors"] == []
    assert [entry["type"] for entry in results["visualizations"]] == ["population_comparison"]
    plot = Path(results["visualizations"][0]["file"])
    assert plot == tmp_path / PLOT_FILE
    assert plot.is_file() and plot.stat().st_size > 0
    assert io.load_json(tmp_path / "comparison_results.json")["visualization_errors"] == []


def test_invalid_input_is_still_rejected(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="non-empty"):
        compare_populations([], GROUP_TWO, output_dir=tmp_path)
    with pytest.raises(ValueError, match="non-empty"):
        compare_populations(GROUP_ONE, [], output_dir=tmp_path)
