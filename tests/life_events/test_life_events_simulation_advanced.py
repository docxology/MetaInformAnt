"""Tests for simulation functions in life_events module.

``generate_synthetic_life_events`` used to accept ``transition_probabilities``, ``seasonal_patterns``,
``rare_event_probability``, ``temporal_noise``, ``missing_data_probability`` and
``cooccurrence_patterns`` through ``**kwargs`` and ignore all of them, and the tests that passed them
asserted nothing about them. ``rare_event_probability``, ``missing_data_probability`` and
``temporal_noise`` are implemented now; any other keyword argument raises ``TypeError``.
"""

from __future__ import annotations

from datetime import datetime
from typing import Any, List, Tuple

import numpy as np
import pytest

from metainformant.life_events.core.events import Event, EventSequence
from metainformant.life_events.core.utils import (
    add_temporal_noise,
    generate_cohort_sequences,
    generate_event_chain,
    generate_synthetic_life_events,
)

NEGATIVE_EVENTS = {"diagnosis", "divorce", "financial_difficulty"}
SECONDS_PER_YEAR = 365.25 * 24 * 3600
# The rare events documented for ``rare_event_probability``: event type -> life domain.
RARE_EVENT_DOMAINS = {
    "serious_accident": "health",
    "bereavement": "family",
    "natural_disaster": "housing",
    "windfall": "finance",
    "business_failure": "career",
}


def test_generate_event_chain(tmp_path):
    """Test event chain generation with transition probabilities."""
    chain_rules = {
        "education:degree": {"occupation:job_change": 0.8, "income:raise": 0.2},
        "occupation:job_change": {"address:move": 0.5, "income:raise": 0.5},
    }

    event_types_by_domain = {
        "education": ["degree"],
        "occupation": ["job_change"],
        "income": ["raise"],
        "address": ["move"],
    }

    events = generate_event_chain(
        chain_rules=chain_rules,
        start_event="degree",
        start_domain="education",
        n_events=5,
        start_timestamp=datetime(2010, 1, 1).timestamp(),
        time_span=365 * 5 * 86400,
        event_types_by_domain=event_types_by_domain,
        random_state=42,
    )

    assert len(events) == 5
    assert all(isinstance(e, Event) for e in events)
    assert events[0].event_type == "degree"
    assert events[0].domain == "education"


def test_add_temporal_noise(tmp_path):
    """Test temporal noise addition to sequences."""
    from datetime import datetime

    seq = EventSequence(
        "p1",
        [
            Event("degree", datetime(2010, 1, 1), "education"),
            Event("job_change", datetime(2015, 1, 1), "occupation"),
        ],
    )

    noisy_seq = add_temporal_noise(seq, noise_level=0.5, max_days_shift=7, missing_probability=0.0, random_state=42)

    assert len(noisy_seq.events) <= len(seq.events)
    assert noisy_seq.person_id == seq.person_id


def test_generate_synthetic_life_events(tmp_path):
    """Test realistic life events generation with a rare-event rate and event-based outcomes."""
    sequences, outcomes = generate_synthetic_life_events(
        n_sequences=20,
        min_events_per_sequence=5,
        max_events_per_sequence=15,
        rare_event_probability=0.1,
        generate_outcomes=True,
        outcome_relationship="event_based",
        random_state=42,
    )

    assert len(sequences) == 20
    assert outcomes is not None
    assert len(outcomes) == 20
    assert all(len(seq.events) >= 5 for seq in sequences)
    assert all(len(seq.events) <= 15 for seq in sequences)

    # The rare-event rate is an actual effect: about one event in ten is rare ...
    events = [event for seq in sequences for event in seq.events]
    rare_share = sum(event.event_type in RARE_EVENT_DOMAINS for event in events) / len(events)
    assert 0.03 <= rare_share <= 0.2
    # ... and the event-based outcome follows the documented rule on the sequences that were returned.
    assert outcomes == [int(any(event.event_type in NEGATIVE_EVENTS for event in seq.events)) for seq in sequences]
    assert set(outcomes) <= {0, 1}


@pytest.mark.parametrize(
    "option",
    [
        {"transition_probabilities": {"education": {"occupation": 0.7}, "occupation": {"income": 0.6}}},
        {"seasonal_patterns": True},
        {"cooccurrence_patterns": {"job_change": ["move"], "degree": ["job_change"]}},
        {"not_an_option": 1},
    ],
)
def test_unsupported_options_raise_type_error_instead_of_being_ignored(option):
    name = next(iter(option))

    with pytest.raises(TypeError, match=name):
        generate_synthetic_life_events(n_sequences=3, random_state=1, **option)


def test_the_options_that_were_silently_ignored_before_are_rejected_together():
    with pytest.raises(TypeError, match="transition_probabilities|seasonal_patterns"):
        generate_synthetic_life_events(
            n_sequences=20,
            transition_probabilities={"education": {"occupation": 0.7}},
            seasonal_patterns=True,
            rare_event_probability=0.1,
            random_state=42,
        )


def test_generate_cohort_sequences(tmp_path):
    """Test cohort sequence generation."""
    cohorts = generate_cohort_sequences(
        n_cohorts=2, n_sequences_per_cohort=10, min_events_per_sequence=5, max_events_per_sequence=10, random_state=42
    )

    assert len(cohorts) == 2
    assert all(len(seqs) == 10 for seqs in cohorts.values())
    assert all("cohort" in seq.metadata for cohort_seqs in cohorts.values() for seq in cohort_seqs)


# --- the implemented options ----------------------------------------------------------------------------


def cohort(n: int = 300, seed: int = 11, **options: Any) -> Tuple[List[EventSequence], Any]:
    return generate_synthetic_life_events(
        n_sequences=n, min_events_per_sequence=3, max_events_per_sequence=15, random_state=seed, **options
    )


def describe(sequences: List[EventSequence]) -> List[List[Tuple[str, str, datetime]]]:
    return [[(e.event_type, e.domain, e.timestamp) for e in seq.events] for seq in sequences]


def n_events(sequences: List[EventSequence]) -> int:
    return sum(len(seq.events) for seq in sequences)


def is_subsequence(part: List[Any], whole: List[Any]) -> bool:
    remaining = iter(whole)
    return all(item in remaining for item in part)


def mean_displacement_years(original: List[EventSequence], noisy: List[EventSequence]) -> float:
    """Mean distance a timestamp moved: sorted timestamps are the optimal matching of two point sets."""
    total, count = 0.0, 0
    for before, after in zip(original, noisy):
        a = np.sort([e.timestamp.timestamp() for e in before.events])
        b = np.sort([e.timestamp.timestamp() for e in after.events])
        total += float(np.abs(a - b).sum())
        count += a.size
    return total / count / SECONDS_PER_YEAR


def mean_span_years(sequences: List[EventSequence]) -> float:
    return float(np.mean([(seq.events[-1].timestamp - seq.events[0].timestamp).days / 365.25 for seq in sequences]))


def test_generate_realistic_with_noise(tmp_path):
    """Temporal noise moves timestamps and missing data removes about the requested share of the events."""
    baseline, _ = generate_synthetic_life_events(n_sequences=10, random_state=42)
    sequences, outcomes = generate_synthetic_life_events(
        n_sequences=10, temporal_noise=0.2, missing_data_probability=0.1, random_state=42
    )

    assert len(sequences) == 10
    assert outcomes is None
    assert all(len(seq.events) > 0 for seq in sequences)
    # Missing data is an effect: events disappeared (about one in ten), and none appeared out of nowhere.
    assert n_events(sequences) < n_events(baseline)
    assert n_events(sequences) >= 0.6 * n_events(baseline)
    # Temporal noise is an effect: the surviving events carry different timestamps, in chronological order.
    assert all([e.timestamp for e in seq.events] == sorted(e.timestamp for e in seq.events) for seq in sequences)
    original_times = {timestamp for events in describe(baseline) for _, _, timestamp in events}
    assert any(timestamp not in original_times for events in describe(sequences) for _, _, timestamp in events)


def test_cooccurrence_patterns_are_not_supported(tmp_path):
    """Co-occurrence patterns are not implemented: asking for them is an error, not a no-op."""
    with pytest.raises(TypeError, match="cooccurrence_patterns"):
        generate_synthetic_life_events(
            n_sequences=20, cooccurrence_patterns={"job_change": ["move"], "degree": ["job_change"]}, random_state=42
        )


def test_zero_valued_options_reproduce_the_default_output():
    plain, plain_outcomes = cohort(generate_outcomes=True, outcome_relationship="event_based")
    explicit, explicit_outcomes = cohort(
        generate_outcomes=True,
        outcome_relationship="event_based",
        rare_event_probability=0.0,
        missing_data_probability=0.0,
        temporal_noise=0.0,
    )

    assert describe(plain) == describe(explicit)
    assert plain_outcomes == explicit_outcomes
    assert not any(event.event_type in RARE_EVENT_DOMAINS for seq in plain for event in seq.events)


# rare_event_probability


def test_the_documented_rare_events_are_exported():
    from metainformant.life_events.core.utils import RARE_EVENT_DOMAINS as exported

    assert exported == RARE_EVENT_DOMAINS


@pytest.mark.parametrize("probability", [0.1, 0.5])
def test_rare_event_probability_sets_the_share_of_rare_events(probability):
    baseline, _ = cohort()
    sequences, _ = cohort(rare_event_probability=probability)

    rare = sum(event.event_type in RARE_EVENT_DOMAINS for seq in sequences for event in seq.events)

    assert abs(rare / n_events(sequences) - probability) < 0.04
    assert n_events(sequences) == n_events(baseline)  # events are replaced, none are added or removed


def test_rare_events_replace_events_in_place():
    baseline, _ = cohort()
    sequences, _ = cohort(rare_event_probability=0.3)

    replaced = kept = 0
    for before, after in zip(describe(baseline), describe(sequences)):
        assert [timestamp for _, _, timestamp in before] == [timestamp for _, _, timestamp in after]
        for old, new in zip(before, after):
            if new[0] in RARE_EVENT_DOMAINS:
                replaced += 1
            else:
                assert new == old
                kept += 1
    assert replaced > 0 and kept > 0


def test_rare_events_have_their_documented_domains_and_all_types_occur():
    sequences, _ = cohort(rare_event_probability=1.0)

    events = [event for seq in sequences for event in seq.events]
    assert all(event.event_type in RARE_EVENT_DOMAINS for event in events)  # probability 1: every event is rare
    assert all(event.domain == RARE_EVENT_DOMAINS[event.event_type] for event in events)
    assert {event.event_type for event in events} == set(RARE_EVENT_DOMAINS)


def test_event_based_outcomes_describe_the_sequence_with_its_rare_events():
    sequences, outcomes = cohort(rare_event_probability=0.4, generate_outcomes=True, outcome_relationship="event_based")

    assert outcomes == [int(any(event.event_type in NEGATIVE_EVENTS for event in seq.events)) for seq in sequences]


# missing_data_probability


@pytest.mark.parametrize("probability", [0.2, 0.6])
def test_missing_data_probability_sets_the_share_of_missing_events(probability):
    baseline, _ = cohort()
    sequences, _ = cohort(missing_data_probability=probability)

    assert abs(n_events(sequences) / n_events(baseline) - (1 - probability)) < 0.04


def test_missing_data_only_removes_events_of_the_complete_sequences():
    baseline, _ = cohort()
    sequences, _ = cohort(missing_data_probability=0.4)

    for complete, observed in zip(describe(baseline), describe(sequences)):
        assert is_subsequence(observed, complete)
    assert any(len(observed) < len(complete) for complete, observed in zip(describe(baseline), describe(sequences)))


@pytest.mark.parametrize("relationship", ["random", "event_based", "temporal"])
def test_missing_data_does_not_change_the_outcomes(relationship):
    _, complete_outcomes = cohort(generate_outcomes=True, outcome_relationship=relationship)
    sequences, outcomes = cohort(
        missing_data_probability=0.5, generate_outcomes=True, outcome_relationship=relationship
    )

    # Outcomes describe the complete record, so hiding events hides information from the sequence.
    assert outcomes == complete_outcomes
    assert len(sequences) == len(outcomes)


def test_missing_data_probability_one_removes_every_event():
    sequences, outcomes = cohort(n=20, missing_data_probability=1.0, generate_outcomes=True)

    assert all(len(seq.events) == 0 for seq in sequences)
    assert outcomes is not None and len(outcomes) == 20


def test_missing_events_hide_the_diagnosis_that_defines_the_outcome():
    """A positive event-based outcome whose diagnosis went missing is the point of the option."""
    sequences, outcomes = cohort(
        missing_data_probability=0.5, generate_outcomes=True, outcome_relationship="event_based"
    )

    hidden = [
        outcome == 1 and not any(event.event_type in NEGATIVE_EVENTS for event in seq.events)
        for seq, outcome in zip(sequences, outcomes)
    ]
    assert any(hidden)


# temporal_noise


def test_zero_temporal_noise_leaves_timestamps_unchanged():
    baseline, _ = cohort()
    sequences, _ = cohort(temporal_noise=0.0)

    assert describe(baseline) == describe(sequences)


def test_temporal_noise_moves_timestamps_but_keeps_the_events():
    baseline, _ = cohort()
    sequences, _ = cohort(temporal_noise=0.2)

    changed = 0
    for before, after in zip(baseline, sequences):
        assert sorted((e.event_type, e.domain) for e in before.events) == sorted(
            (e.event_type, e.domain) for e in after.events
        )
        stamps = [e.timestamp for e in after.events]
        assert stamps == sorted(stamps)
        changed += stamps != [e.timestamp for e in before.events]
    assert changed >= 0.95 * len(baseline)


def test_temporal_noise_scales_with_the_requested_level():
    baseline, _ = cohort()
    span = mean_span_years(baseline)
    displacement = {level: mean_displacement_years(baseline, cohort(temporal_noise=level)[0]) for level in (0.05, 0.2)}

    assert displacement[0.05] < displacement[0.2]
    for level, moved in displacement.items():
        # Gaussian jitter with a standard deviation of level * span moves a timestamp by about 0.8 * level * span.
        assert 0.4 < moved / (level * span) < 0.9


def test_temporal_noise_does_not_change_outcomes_or_which_events_occur():
    baseline, baseline_outcomes = cohort(generate_outcomes=True, outcome_relationship="event_based")
    sequences, outcomes = cohort(temporal_noise=0.3, generate_outcomes=True, outcome_relationship="event_based")

    assert outcomes == baseline_outcomes
    assert [sorted(e.event_type for e in seq.events) for seq in sequences] == [
        sorted(e.event_type for e in seq.events) for seq in baseline
    ]


# combinations and reproducibility


@pytest.mark.parametrize(
    "options",
    [
        {"rare_event_probability": 0.3},
        {"missing_data_probability": 0.3},
        {"temporal_noise": 0.2},
        {"rare_event_probability": 0.3, "missing_data_probability": 0.3, "temporal_noise": 0.2},
    ],
)
def test_options_are_reproducible_for_a_seed(options):
    first, _ = cohort(seed=5, **options)
    again, _ = cohort(seed=5, **options)
    other, _ = cohort(seed=6, **options)

    assert describe(first) == describe(again)
    assert describe(first) != describe(other)


def test_options_use_independent_random_streams():
    """Switching one option on does not change what another option does for the same seed."""
    rare_only, _ = cohort(rare_event_probability=0.3)
    both, _ = cohort(rare_event_probability=0.3, missing_data_probability=0.4)

    for with_rare, observed in zip(describe(rare_only), describe(both)):
        assert is_subsequence(observed, with_rare)  # the same rare events, minus the missing ones


@pytest.mark.parametrize(
    "options",
    [
        {"rare_event_probability": -0.1},
        {"rare_event_probability": 1.5},
        {"missing_data_probability": -0.1},
        {"missing_data_probability": 2},
        {"missing_data_probability": float("nan")},
        {"temporal_noise": -0.5},
        {"temporal_noise": float("nan")},
    ],
)
def test_invalid_option_values_are_rejected(options):
    with pytest.raises(ValueError, match=next(iter(options))):
        generate_synthetic_life_events(n_sequences=3, random_state=1, **options)
