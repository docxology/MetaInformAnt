"""Regression tests: ``biological_embedding(..., "temporal")`` works on real event sequences.

Defect: the temporal embedding took ``np.diff`` of the events' ``datetime`` timestamps, which gives an array of
``timedelta`` objects; ``np.std`` then raised ``TypeError`` ("loop of ufunc does not support argument 0 of type
datetime.timedelta"). The call failed for every ``EventSequence`` with at least two events, so the ``'temporal'``
embedding type could never be used. The intervals are now converted to days before the statistics are taken.

Oracles: the pooled intervals and the six statistics are worked out by hand below and re-derived with the standard
library's ``statistics`` module, independently of NumPy.
"""

from __future__ import annotations

import statistics
from datetime import datetime, timedelta

import numpy as np
import pytest

from metainformant.life_events.core.events import Event, EventSequence
from metainformant.life_events.models.embeddings import biological_embedding
from metainformant.phenotype.analysis.life_course import Event as NumericEvent, EventSequence as NumericEventSequence

T0 = datetime(2020, 1, 1)


def _sequence(person: str, start: datetime, steps: list[tuple[str, float]]) -> EventSequence:
    """Events of the given types, each ``days`` after the previous one (the first is at ``start``)."""
    events = []
    moment = start
    for event_type, days in steps:
        moment = moment + timedelta(days=days)
        events.append(Event(event_type, moment, "life"))
    return EventSequence(person, events)


def _expected_vector(intervals: list[float], embedding_dim: int) -> np.ndarray:
    """Mean, population SD, min, max, count and median of the intervals, scaled by the largest entry."""
    raw = [
        statistics.fmean(intervals),
        statistics.pstdev(intervals),
        min(intervals),
        max(intervals),
        float(len(intervals)),
        statistics.median(intervals),
    ]
    vector = np.array(raw + [0.0] * (embedding_dim - 6))
    return vector / np.max(np.abs(vector))


# Two sequences share the event pattern (a, b, c): intervals 10, 30 days and 5, 30 days.
SEQ_ONE = _sequence("p1", T0, [("a", 0), ("b", 10), ("c", 30)])
SEQ_TWO = _sequence("p2", datetime(2021, 3, 1), [("a", 0), ("b", 5), ("c", 30)])
# One sequence with its own pattern (x, y), 36 hours apart.
SEQ_THREE = _sequence("p3", datetime(2022, 1, 1), [("x", 0), ("y", 1.5)])


def test_temporal_embedding_no_longer_raises_for_event_sequences() -> None:
    result = biological_embedding([SEQ_ONE], "temporal")

    assert result["type"] == "temporal"
    assert list(result["embeddings"]) == [str(("a", "b", "c"))]


def test_embedding_matches_the_hand_computed_interval_statistics() -> None:
    result = biological_embedding([SEQ_ONE, SEQ_TWO, SEQ_THREE], "temporal", embedding_dim=10)

    # Pooled intervals of pattern (a, b, c): 10 and 30 (p1), 5 and 30 (p2).
    #   mean 18.75, population SD sqrt(129.6875) = 11.3880, min 5, max 30, count 4, median 20; largest entry 30.
    by_hand = np.array([18.75, 11.387, 5.0, 30.0, 4.0, 20.0, 0, 0, 0, 0]) / 30.0
    np.testing.assert_allclose(result["embeddings"][str(("a", "b", "c"))], by_hand, atol=1e-4)
    np.testing.assert_allclose(
        result["embeddings"][str(("a", "b", "c"))], _expected_vector([10.0, 30.0, 5.0, 30.0], 10), rtol=1e-12
    )
    # Pattern (x, y): one interval of 1.5 days -> [1.5, 0, 1.5, 1.5, 1, 1.5] scaled by 1.5.
    np.testing.assert_allclose(
        result["embeddings"][str(("x", "y"))], [1.0, 0.0, 1.0, 1.0, 1 / 1.5, 1.0, 0, 0, 0, 0], rtol=1e-12
    )


def test_intervals_are_measured_in_days_including_fractions() -> None:
    result = biological_embedding([SEQ_THREE], "temporal")

    assert result["temporal_patterns"][("x", "y")] == [1.5]  # 36 hours


def test_pooled_intervals_keep_the_order_of_the_sequences() -> None:
    result = biological_embedding([SEQ_ONE, SEQ_TWO], "temporal")

    assert result["temporal_patterns"] == {("a", "b", "c"): [10.0, 30.0, 5.0, 30.0]}


def test_only_the_first_three_event_types_form_the_pattern_but_every_interval_counts() -> None:
    long_sequence = _sequence("p4", T0, [("e1", 0), ("e2", 1), ("e3", 2), ("e4", 3), ("e5", 4)])  # days 0, 1, 3, 6, 10

    result = biological_embedding([long_sequence], "temporal", embedding_dim=8)

    assert result["temporal_patterns"] == {("e1", "e2", "e3"): [1.0, 2.0, 3.0, 4.0]}
    np.testing.assert_allclose(
        result["embeddings"][str(("e1", "e2", "e3"))], _expected_vector([1.0, 2.0, 3.0, 4.0], 8), rtol=1e-12
    )


def test_vector_has_embedding_dim_entries_zero_padding_and_unit_maximum() -> None:
    for dim in (6, 7, 30):
        vector = biological_embedding([SEQ_ONE, SEQ_TWO], "temporal", embedding_dim=dim)["embeddings"][
            str(("a", "b", "c"))
        ]
        assert vector.shape == (dim,)
        assert np.max(np.abs(vector)) == pytest.approx(1.0)
        assert np.all(vector[6:] == 0.0)


def test_default_dimension_is_thirty() -> None:
    result = biological_embedding([SEQ_ONE], "temporal")

    assert result["embedding_dim"] == 30
    assert result["embeddings"][str(("a", "b", "c"))].shape == (30,)


def test_the_result_is_reproducible() -> None:
    first = biological_embedding([SEQ_ONE, SEQ_TWO, SEQ_THREE], "temporal")
    second = biological_embedding([SEQ_ONE, SEQ_TWO, SEQ_THREE], "temporal")

    assert first["embeddings"].keys() == second["embeddings"].keys()
    for key in first["embeddings"]:
        np.testing.assert_array_equal(first["embeddings"][key], second["embeddings"][key])


def test_sequences_with_fewer_than_two_events_contribute_nothing() -> None:
    single = EventSequence("solo", [Event("a", T0, "life")])

    assert biological_embedding([single], "temporal")["embeddings"] == {}
    assert biological_embedding([], "temporal")["embeddings"] == {}
    only_solo_and_pair = biological_embedding([single, SEQ_THREE], "temporal")
    assert list(only_solo_and_pair["embeddings"]) == [str(("x", "y"))]


def test_simultaneous_events_give_a_zero_interval_not_an_error() -> None:
    """Two events at the same instant: the interval is 0, so the count (1) is the only non-zero statistic."""
    same_time = EventSequence("tie", [Event("a", T0, "life"), Event("b", T0, "life")])

    vector = biological_embedding([same_time], "temporal", embedding_dim=8)["embeddings"][str(("a", "b"))]

    expected = np.zeros(8)
    expected[4] = 1.0  # the count of intervals (1) is the only non-zero statistic; it becomes the maximum
    np.testing.assert_array_equal(vector, expected)


def test_numeric_timestamps_are_subtracted_as_given() -> None:
    """Sequences whose events carry plain numbers (the phenotype life-course classes) work as well."""
    sequence = NumericEventSequence(
        "n1", [NumericEvent(timestamp=0.0, event_type="a"), NumericEvent(timestamp=4.0, event_type="b")]
    )

    result = biological_embedding([sequence], "temporal", embedding_dim=6)

    assert result["temporal_patterns"] == {("a", "b"): [4.0]}
    np.testing.assert_allclose(result["embeddings"][str(("a", "b"))], [1.0, 0.0, 1.0, 1.0, 0.25, 1.0], rtol=1e-12)


@pytest.mark.parametrize("dim", [0, 3, 5])
def test_an_embedding_dimension_below_the_six_statistics_is_rejected(dim: int) -> None:
    """The vector holds six statistics; a smaller ``embedding_dim`` used to yield longer vectors than it claimed."""
    with pytest.raises(ValueError, match="at least 6"):
        biological_embedding([SEQ_ONE], "temporal", embedding_dim=dim)


def test_the_other_embedding_types_are_unaffected() -> None:
    assert biological_embedding([SEQ_ONE], "event_type")["type"] == "event_type"
    assert biological_embedding([SEQ_ONE], "domain")["type"] == "domain"
    with pytest.raises(ValueError, match="Unknown embedding type"):
        biological_embedding([SEQ_ONE], "calendar")
