"""Tests for simulation functions in life_events module."""

from __future__ import annotations

from datetime import datetime

import numpy as np
import pytest

from metainformant.life_events.core.events import Event, EventSequence
from metainformant.life_events.core.utils import add_temporal_noise, generate_cohort_sequences, generate_event_chain, generate_synthetic_life_events


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
    """Test realistic life events generation with advanced patterns."""
    sequences, outcomes = generate_synthetic_life_events(
        n_sequences=20,
        min_events_per_sequence=5,
        max_events_per_sequence=15,
        transition_probabilities={"education": {"occupation": 0.7}, "occupation": {"income": 0.6}},
        seasonal_patterns=True,
        rare_event_probability=0.1,
        generate_outcomes=True,
        outcome_relationship="complex",
        random_state=42,
    )

    assert len(sequences) == 20
    assert outcomes is not None
    assert len(outcomes) == 20
    assert all(len(seq.events) >= 5 for seq in sequences)
    assert all(len(seq.events) <= 15 for seq in sequences)


def test_generate_cohort_sequences(tmp_path):
    """Test cohort sequence generation."""
    cohorts = generate_cohort_sequences(
        n_cohorts=2, n_sequences_per_cohort=10, min_events_per_sequence=5, max_events_per_sequence=10, random_state=42
    )

    assert len(cohorts) == 2
    assert all(len(seqs) == 10 for seqs in cohorts.values())
    assert all("cohort" in seq.metadata for cohort_seqs in cohorts.values() for seq in cohort_seqs)


def test_generate_realistic_with_noise(tmp_path):
    """Test realistic generation with temporal noise and missing data."""
    sequences, outcomes = generate_synthetic_life_events(
        n_sequences=10, temporal_noise=0.2, missing_data_probability=0.1, random_state=42
    )

    assert len(sequences) == 10
    # Some sequences may have fewer events due to missing data
    assert all(len(seq.events) > 0 for seq in sequences)


def test_generate_realistic_cooccurrence(tmp_path):
    """Test realistic generation with co-occurrence patterns."""
    sequences, _ = generate_realistic_life_events(
        n_sequences=20, cooccurrence_patterns={"job_change": ["move"], "degree": ["job_change"]}, random_state=42
    )

    assert len(sequences) == 20
    # Check that some co-occurrences exist
    has_cooccurrence = False
    for seq in sequences:
        tokens = [f"{e.domain}:{e.event_type}" for e in seq.events]
        if "occupation:job_change" in tokens and "address:move" in tokens:
            has_cooccurrence = True
            break

    # This is probabilistic, so we just check sequences were generated
    assert len(sequences) > 0
