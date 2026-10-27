"""Regression tests: ``workflow.event_importance(method="frequency")`` ranks by frequency.

Defect: the ``"frequency"`` branch first stored each event type's number of occurrences and then overwrote every
score with ``1 - mean relative position``, so the "frequency" importance was really an earliness score: the most
frequent event type could rank last, and a rare event that happened to come first ranked first. The docstring
listed a ``"temporal"`` method as well, but it was not implemented (it returned ``{}``). The earliness score now
lives under ``"temporal"``, ``"frequency"`` counts occurrences, and an unknown method raises instead of returning
an empty result.

Oracles: occurrence counts, relative positions and transition counts are worked out by hand for the sequences
below (and the counts are re-derived with ``collections.Counter``).
"""

from __future__ import annotations

from collections import Counter
from datetime import datetime, timedelta

import pytest

from metainformant.life_events.core.events import Event, EventSequence
from metainformant.life_events.workflow.workflow import event_importance


def _sequence(person: str, event_types: list[str]) -> EventSequence:
    start = datetime(2020, 1, 1)
    return EventSequence(person, [Event(t, start + timedelta(days=i), "life") for i, t in enumerate(event_types)])


# 'a' is by far the most frequent event type but comes late; 'x' and 'b' come first but are rare.
SEQUENCES = [
    _sequence("p1", ["x", "y", "a", "a", "a"]),
    _sequence("p2", ["b", "a", "a"]),
    _sequence("p3", ["b", "c"]),
]
COUNTS = {"a": 5, "b": 2, "x": 1, "y": 1, "c": 1}


def test_the_hand_counts_match_the_sequences() -> None:
    counted = Counter(event.event_type for sequence in SEQUENCES for event in sequence.events)

    assert dict(counted) == COUNTS


def test_frequency_importance_is_the_normalised_occurrence_count() -> None:
    scores = event_importance(SEQUENCES, method="frequency")

    # min-max scaling of the counts 5, 2, 1, 1, 1: (count - 1) / (5 - 1)
    assert scores == pytest.approx({"a": 1.0, "b": 0.25, "x": 0.0, "y": 0.0, "c": 0.0})


def test_frequency_is_the_default_method() -> None:
    assert event_importance(SEQUENCES) == event_importance(SEQUENCES, method="frequency")


def test_the_most_frequent_event_type_ranks_first_by_frequency() -> None:
    scores = event_importance(SEQUENCES, method="frequency")

    ranking = sorted(scores, key=lambda event_type: (-scores[event_type], event_type))
    assert ranking[:2] == ["a", "b"]  # by earliness 'a' would have ranked fourth of five


def test_unnormalised_frequency_scores_are_the_counts() -> None:
    scores = event_importance(SEQUENCES, method="frequency", normalize=False)

    assert scores == {event_type: float(count) for event_type, count in COUNTS.items()}


def test_equally_frequent_event_types_get_the_neutral_score() -> None:
    sequences = [_sequence("q1", ["u", "v"]), _sequence("q2", ["v", "u"])]

    assert event_importance(sequences, method="frequency") == {"u": 0.5, "v": 0.5}


def test_temporal_importance_is_one_minus_the_mean_relative_position() -> None:
    scores = event_importance(SEQUENCES, method="temporal", normalize=False)

    # Relative position of an event is index / (length - 1).
    #   x: position 0 -> 1.0                         y: 1/4 -> 0.75
    #   a: positions 2/4, 3/4, 4/4 (p1) and 1/2, 2/2 (p2) -> mean 0.75 -> 0.25
    #   b: 0 in p2 and 0 in p3 -> 1.0                c: 1 in p3 -> 0.0
    assert scores == pytest.approx({"x": 1.0, "y": 0.75, "a": 0.25, "b": 1.0, "c": 0.0})


def test_temporal_importance_is_what_frequency_used_to_return() -> None:
    """The earliness ranking is not lost: it differs from the frequency ranking on this data."""
    temporal = event_importance(SEQUENCES, method="temporal")
    frequency = event_importance(SEQUENCES, method="frequency")

    assert temporal != frequency
    assert temporal["a"] == pytest.approx(0.25)  # 'a' is common but late
    assert frequency["a"] == pytest.approx(1.0)


def test_transition_importance_is_unchanged() -> None:
    scores = event_importance(SEQUENCES, method="transition", normalize=False)

    # Transitions: x->y, y->a, a->a (x3), b->a, b->c. A type scores its incoming plus outgoing transitions
    # (a self-transition counts on both sides): a 5 + 3, y 1 + 1, b 0 + 2, x 0 + 1, c 1 + 0.
    assert scores == {"x": 1.0, "y": 2.0, "a": 8.0, "b": 2.0, "c": 1.0}


@pytest.mark.parametrize("method", ["popularity", "", "Frequency"])
def test_an_unknown_method_raises_instead_of_returning_nothing(method: str) -> None:
    with pytest.raises(ValueError, match="Unknown method"):
        event_importance(SEQUENCES, method=method)


def test_an_unknown_method_raises_even_without_sequences() -> None:
    with pytest.raises(ValueError, match="Unknown method"):
        event_importance([], method="popularity")


@pytest.mark.parametrize("method", ["frequency", "temporal", "transition"])
def test_no_sequences_give_no_scores(method: str) -> None:
    assert event_importance([], method=method) == {}
