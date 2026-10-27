"""Regression tests: ``analyze_life_course`` states its mean duration honestly and survives sequences without events.

Defect: the result held ``"analysis_timestamp": statistics.mean([...durations...])``. The key promised a timestamp but
held a mean duration (no documentation, test or caller treats it as a time), and the mean was taken over an empty list
whenever no sequence had any events, so ``analyze_life_course`` raised ``statistics.StatisticsError`` for such input
-- for instance a cohort of people with no recorded events yet. The value now lives under ``mean_duration`` and, with
no events anywhere, the summaries are explicitly unavailable (None / empty) instead of an exception.

Oracles: durations, means and counts are worked out by hand from the timestamps below.
"""

from __future__ import annotations

import json
import statistics
from typing import Optional

import pytest

from metainformant.phenotype.analysis.life_course import Event, EventSequence, analyze_life_course


def sequence(person: str, *events: tuple[float, str], outcomes: Optional[object] = None) -> EventSequence:
    metadata = {} if outcomes is None else {"outcomes": outcomes}
    return EventSequence(person, [Event(timestamp=t, event_type=e) for t, e in events], metadata=metadata)


# Durations (last timestamp minus first): 20 and 40; the third person has no events.
LONG = sequence("p1", (0, "a"), (5, "b"), (20, "c"))
LONGER = sequence("p2", (10, "a"), (30, "a"), (50, "d"))
NO_EVENTS = sequence("p3")


def test_the_mean_duration_has_an_honest_name_and_the_old_key_is_gone() -> None:
    result = analyze_life_course([LONG, LONGER, NO_EVENTS])

    assert "analysis_timestamp" not in result
    assert result["mean_duration"] == 30.0  # (20 + 40) / 2
    assert result["mean_duration"] == result["duration_stats"]["mean"]


def test_sequences_without_events_do_not_dilute_the_mean_duration() -> None:
    with_empty = analyze_life_course([LONG, LONGER, NO_EVENTS, sequence("p4")])
    without_empty = analyze_life_course([LONG, LONGER])

    assert with_empty["mean_duration"] == without_empty["mean_duration"] == 30.0  # not 60 / 4 = 15
    assert with_empty["n_sequences"] == 4
    assert with_empty["n_sequences_with_events"] == 2
    assert without_empty["n_sequences_with_events"] == 2


def test_the_other_summaries_of_a_normal_cohort_are_unchanged() -> None:
    result = analyze_life_course([LONG, LONGER, NO_EVENTS])

    assert result["duration_stats"] == {
        "mean": 30.0,
        "median": 30.0,
        "min": 20.0,
        "max": 40.0,
        "std": statistics.stdev([20.0, 40.0]),
    }
    assert result["event_frequency"] == {"a": 3, "b": 1, "c": 1, "d": 1}
    # Complexity: 3 distinct types in 3 events = 1.0 and 2 distinct in 3 events = 2/3 (the empty one is left out).
    assert result["complexity_stats"]["mean"] == pytest.approx((1.0 + 2 / 3) / 2)
    assert result["trajectory_patterns"]["n_sequences"] == 3


def test_a_cohort_without_any_events_gives_an_explicit_unavailable_result_not_an_exception() -> None:
    """Was ``statistics.StatisticsError: mean requires at least one data point``."""
    result = analyze_life_course([sequence("q1"), sequence("q2"), sequence("q3")])

    assert result["n_sequences"] == 3
    assert result["n_sequences_with_events"] == 0
    assert result["mean_duration"] is None
    assert result["duration_stats"] is None
    assert result["complexity_stats"] is None
    assert result["event_frequency"] == {}
    assert result["trajectory_patterns"]["patterns"] == {"empty": 3}
    assert "analysis_timestamp" not in result


def test_a_single_sequence_without_events_is_handled_too() -> None:
    result = analyze_life_course([sequence("only")])

    assert (result["n_sequences"], result["n_sequences_with_events"], result["mean_duration"]) == (1, 0, None)


def test_a_single_sequence_with_one_event_has_zero_duration_not_an_unavailable_one() -> None:
    result = analyze_life_course([sequence("one", (7, "a"))])

    assert result["mean_duration"] == 0.0
    assert result["duration_stats"]["std"] == 0
    assert result["n_sequences_with_events"] == 1


def test_outcome_analysis_still_works_when_there_are_no_events() -> None:
    cohort = [sequence("r1", outcomes=["relapse"]), sequence("r2", outcomes=[]), sequence("r3", outcomes=["relapse"])]

    result = analyze_life_course(cohort, outcomes=["relapse"])

    analysis = result["outcome_analysis"]["relapse"]
    assert (analysis["n_with_outcome"], analysis["n_without_outcome"]) == (2, 1)
    assert analysis["prevalence"] == pytest.approx(2 / 3)
    # Nobody has a duration or a complexity, so those comparisons are unavailable, not invented
    assert analysis["trajectory_comparison"]["duration"] is None
    assert analysis["trajectory_comparison"]["complexity"] is None
    assert analysis["event_type_association"] == {}  # no event types to test


def test_the_unavailable_result_is_plain_json() -> None:
    result = analyze_life_course([sequence("q1"), sequence("q2")], outcomes=["relapse"])

    json.dumps(result, allow_nan=False)


def test_an_empty_list_of_sequences_is_still_reported_as_an_error() -> None:
    assert analyze_life_course([]) == {"error": "No sequences provided"}
