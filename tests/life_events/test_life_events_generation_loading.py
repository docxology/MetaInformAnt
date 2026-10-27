"""Regression tests for ``generate_synthetic_life_events`` and ``load_sequences_from_json``.

Defects covered:

* An unknown ``outcome_relationship`` (``"complex"``, ``"health_focused"``, ...) silently produced
  random outcomes instead of failing.
* ``temporal_patterns`` was documented ("realistic", "random", "clustered") but never read.
* ``load_sequences_from_json`` dropped each sequence's ``metadata`` and each event's ``person_id``,
  ``description`` and ``confidence``.
"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import pytest
from scipy import stats

from metainformant.life_events.core.events import Event, EventDatabase, EventSequence
from metainformant.life_events.core.utils import (
    generate_cohort_sequences,
    generate_synthetic_life_events,
    load_sequences_from_json,
)

IMPLEMENTED_RELATIONSHIPS = ("random", "event_based", "temporal")
NEGATIVE_EVENTS = {"diagnosis", "divorce", "financial_difficulty"}


# --- outcome_relationship --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name", ["complex", "health_focused", "education_focused", "Random", "event-based", "", "temporal "]
)
def test_unknown_outcome_relationship_is_rejected(name):
    with pytest.raises(ValueError, match="outcome_relationship"):
        generate_synthetic_life_events(n_sequences=3, generate_outcomes=True, outcome_relationship=name, random_state=1)


def test_unknown_outcome_relationship_is_rejected_even_when_no_outcomes_are_requested():
    with pytest.raises(ValueError, match="outcome_relationship"):
        generate_synthetic_life_events(n_sequences=3, generate_outcomes=False, outcome_relationship="complex")


def test_the_error_names_the_implemented_relationships():
    with pytest.raises(ValueError) as error:
        generate_synthetic_life_events(n_sequences=1, generate_outcomes=True, outcome_relationship="complex")

    message = str(error.value)
    assert "'complex'" in message
    for name in IMPLEMENTED_RELATIONSHIPS:
        assert name in message


def test_event_based_outcomes_follow_the_documented_rule():
    sequences, outcomes = generate_synthetic_life_events(
        n_sequences=60, generate_outcomes=True, outcome_relationship="event_based", random_state=5
    )

    expected = [int(any(event.event_type in NEGATIVE_EVENTS for event in seq.events)) for seq in sequences]
    assert outcomes == expected
    assert set(outcomes) == {0, 1}


def test_temporal_outcomes_follow_the_documented_rule():
    sequences, outcomes = generate_synthetic_life_events(
        n_sequences=60, generate_outcomes=True, outcome_relationship="temporal", random_state=5
    )

    assert outcomes == [int(len(seq.events) > 12) for seq in sequences]
    assert set(outcomes) == {0, 1}


def test_random_outcomes_are_binary_and_independent_of_both_rules():
    sequences, outcomes = generate_synthetic_life_events(
        n_sequences=80, generate_outcomes=True, outcome_relationship="random", random_state=5
    )

    assert set(outcomes) == {0, 1}
    event_rule = [int(any(event.event_type in NEGATIVE_EVENTS for event in seq.events)) for seq in sequences]
    length_rule = [int(len(seq.events) > 12) for seq in sequences]
    assert outcomes != event_rule
    assert outcomes != length_rule


@pytest.mark.parametrize("relationship", IMPLEMENTED_RELATIONSHIPS)
def test_outcomes_are_reproducible_for_a_seed(relationship):
    first = generate_synthetic_life_events(
        n_sequences=20, generate_outcomes=True, outcome_relationship=relationship, random_state=9
    )
    second = generate_synthetic_life_events(
        n_sequences=20, generate_outcomes=True, outcome_relationship=relationship, random_state=9
    )

    assert first[1] == second[1]
    assert [seq.to_dict() for seq in first[0]] == [seq.to_dict() for seq in second[0]]


# --- temporal_patterns ---------------------------------------------------------------------------------


def cohort(pattern: Optional[str] = None, seed: int = 3) -> List[EventSequence]:
    options: Dict[str, Any] = {} if pattern is None else {"temporal_patterns": pattern}
    sequences, _ = generate_synthetic_life_events(
        n_sequences=40, min_events_per_sequence=6, max_events_per_sequence=14, random_state=seed, **options
    )
    return sequences


def timestamps(sequence: EventSequence) -> List[datetime]:
    return [event.timestamp for event in sequence.events]


def earliest(sequence: EventSequence, event_type: str) -> Optional[datetime]:
    times = [event.timestamp for event in sequence.events if event.event_type == event_type]
    return min(times) if times else None


def school_before_college(sequences: List[EventSequence]) -> List[bool]:
    outcomes = []
    for sequence in sequences:
        start, finish = earliest(sequence, "started_school"), earliest(sequence, "graduated_college")
        if start is not None and finish is not None:
            outcomes.append(start < finish)
    return outcomes


def long_gaps(sequence: EventSequence, span_seconds: float, fraction: float) -> int:
    times = timestamps(sequence)
    return sum(1 for a, b in zip(times, times[1:]) if (b - a).total_seconds() > fraction * span_seconds)


@pytest.mark.parametrize("pattern", ["seasonal", "Realistic", "", None])
def test_unknown_temporal_pattern_is_rejected(pattern):
    options: Dict[str, Any] = {"temporal_patterns": pattern}
    with pytest.raises(ValueError, match="temporal_patterns"):
        generate_synthetic_life_events(n_sequences=2, random_state=1, **options)


def test_realistic_is_the_default_layout():
    default = cohort()
    explicit = cohort("realistic")

    assert [seq.to_dict() for seq in default] == [seq.to_dict() for seq in explicit]


@pytest.mark.parametrize("pattern", ["random", "clustered"])
def test_other_patterns_keep_the_events_and_only_move_their_timestamps(pattern):
    realistic = cohort("realistic")
    other = cohort(pattern)

    assert len(realistic) == len(other)
    changed = 0
    for original, retimed in zip(realistic, other):
        assert original.person_id == retimed.person_id
        assert sorted((e.event_type, e.domain) for e in original.events) == sorted(
            (e.event_type, e.domain) for e in retimed.events
        )
        start, end = timestamps(original)[0], timestamps(original)[-1]
        assert all(start <= stamp <= end for stamp in timestamps(retimed))  # inside the realistic span
        assert timestamps(retimed) == sorted(timestamps(retimed))
        changed += timestamps(original) != timestamps(retimed)
    assert changed >= 35  # nearly every timeline is laid out differently


@pytest.mark.parametrize("pattern", ["random", "clustered"])
def test_other_patterns_are_reproducible_for_a_seed(pattern):
    first = [seq.to_dict() for seq in cohort(pattern, seed=11)]
    again = [seq.to_dict() for seq in cohort(pattern, seed=11)]
    other_seed = [seq.to_dict() for seq in cohort(pattern, seed=12)]

    assert first == again
    assert first != other_seed


def test_random_pattern_removes_the_life_stage_ordering():
    realistic = school_before_college(cohort("realistic"))
    shuffled = school_before_college(cohort("random"))

    assert len(realistic) >= 20 and len(shuffled) >= 20
    # In life-stage chains school precedes graduating from college (only unrelated filler events break this) ...
    assert sum(realistic) / len(realistic) >= 0.85
    # ... while independent timestamps make the order a coin flip.
    assert 0.2 <= sum(shuffled) / len(shuffled) <= 0.8


def test_random_pattern_spreads_events_uniformly_over_the_timeline():
    realistic = cohort("realistic")
    shuffled = cohort("random")

    positions = []
    for original, retimed in zip(realistic, shuffled):
        start, end = timestamps(original)[0], timestamps(original)[-1]
        span = (end - start).total_seconds()
        positions.extend((stamp - start).total_seconds() / span for stamp in timestamps(retimed))

    assert len(positions) > 250
    assert stats.kstest(positions, "uniform").pvalue > 0.001


def test_clustered_pattern_places_events_in_a_few_bursts():
    realistic = cohort("realistic")
    clustered = cohort("clustered")

    realistic_gap_counts = []
    for original, retimed in zip(realistic, clustered):
        span = (timestamps(original)[-1] - timestamps(original)[0]).total_seconds()
        assert long_gaps(retimed, span, 0.15) <= 2  # at most three bursts, so at most two quiet gaps
        realistic_gap_counts.append(long_gaps(original, span, 0.15))
    assert max(realistic_gap_counts) >= 3  # the realistic layout does not have this structure


def test_non_realistic_patterns_do_not_change_which_events_are_generated_or_their_outcomes():
    baseline, baseline_outcomes = generate_synthetic_life_events(
        n_sequences=30, generate_outcomes=True, outcome_relationship="event_based", random_state=4
    )
    for pattern in ("random", "clustered"):
        sequences, outcomes = generate_synthetic_life_events(
            n_sequences=30,
            generate_outcomes=True,
            outcome_relationship="event_based",
            temporal_patterns=pattern,
            random_state=4,
        )
        assert outcomes == baseline_outcomes
        assert [sorted(e.event_type for e in seq.events) for seq in sequences] == [
            sorted(e.event_type for e in seq.events) for seq in baseline
        ]


# --- load_sequences_from_json --------------------------------------------------------------------------


def detailed_sequence() -> EventSequence:
    return EventSequence(
        "p1",
        [
            Event(
                "degree",
                datetime(2010, 6, 1, 8, 30),
                "education",
                attributes={"school": "Redwoods", "years": 4},
                person_id="p1",
                description="Bachelor's degree",
                confidence=0.8,
            ),
            Event("diagnosis", datetime(2015, 1, 2), "health", confidence=0.25, description="self reported"),
        ],
        metadata={"cohort": "A", "site": 3, "tags": ["pilot"]},
    )


def test_load_sequences_from_json_keeps_sequence_metadata_and_event_details(tmp_path: Path):
    original = detailed_sequence()
    path = tmp_path / "sequences.json"
    path.write_text(json.dumps([original.to_dict()]))

    (loaded,) = load_sequences_from_json(path)

    assert loaded.metadata == {"cohort": "A", "site": 3, "tags": ["pilot"]}
    first, second = loaded.events
    assert (first.person_id, first.description, first.confidence) == ("p1", "Bachelor's degree", 0.8)
    assert first.attributes == {"school": "Redwoods", "years": 4}
    assert (second.person_id, second.description, second.confidence) == (None, "self reported", 0.25)
    assert loaded.to_dict() == original.to_dict()


def test_load_sequences_from_json_reads_a_database_file_back_unchanged(tmp_path: Path):
    cohorts = generate_cohort_sequences(
        n_cohorts=2, n_sequences_per_cohort=3, min_events_per_sequence=3, max_events_per_sequence=6, random_state=2
    )
    assert isinstance(cohorts, dict)
    originals = [seq for group in cohorts.values() for seq in group] + [detailed_sequence()]
    path = tmp_path / "database.json"
    EventDatabase(originals).save_json(path)

    loaded = load_sequences_from_json(path)

    assert [seq.to_dict() for seq in loaded] == [seq.to_dict() for seq in originals]
    assert {seq.metadata["cohort"] for seq in loaded} == {"cohort_0", "cohort_1", "A"}


def test_load_sequences_from_json_uses_defaults_for_missing_optional_fields(tmp_path: Path):
    path = tmp_path / "minimal.json"
    path.write_text(
        json.dumps(
            [
                {
                    "id": "p9",
                    "events": [{"event_type": "degree", "timestamp": "2010-06-01T00:00:00", "domain": "education"}],
                }
            ]
        )
    )

    (loaded,) = load_sequences_from_json(path)

    event = loaded.events[0]
    assert loaded.person_id == "p9"
    assert loaded.metadata == {}
    assert (event.person_id, event.description, event.confidence, event.attributes) == (None, "", 1.0, {})


def test_load_sequences_from_json_treats_null_optional_fields_as_missing(tmp_path: Path):
    path = tmp_path / "nulls.json"
    path.write_text(
        json.dumps(
            {
                "person_id": "p1",
                "metadata": None,
                "events": [
                    {
                        "event_type": "degree",
                        "timestamp": 1262304000,
                        "domain": "education",
                        "attributes": None,
                        "person_id": None,
                        "description": None,
                        "confidence": None,
                    }
                ],
            }
        )
    )

    (loaded,) = load_sequences_from_json(path)

    event = loaded.events[0]
    assert loaded.metadata == {}
    assert (event.person_id, event.description, event.confidence, event.attributes) == (None, "", 1.0, {})
