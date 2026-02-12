"""Tests for life_events utility functions."""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import pytest

from metainformant.life_events.core.events import Event, EventSequence
from metainformant.life_events.core.utils import get_event_statistics, validate_sequence


def test_load_sequences_from_json_database_format(tmp_path):
    """Test loading sequences from database format JSON."""
    from metainformant.life_events.core.utils import load_sequences_from_json

    data = {
        "sequences": [
            {
                "person_id": "person_001",
                "events": [
                    {
                        "event_type": "degree",
                        "timestamp": "2010-06-01T00:00:00",
                        "domain": "education",
                        "attributes": {},
                    }
                ],
                "metadata": {},
            }
        ],
        "metadata": {},
    }

    json_file = tmp_path / "events.json"
    with open(json_file, "w") as f:
        json.dump(data, f)

    sequences = load_sequences_from_json(json_file)
    assert len(sequences) == 1
    assert sequences[0].person_id == "person_001"


def test_load_sequences_from_json_list_format(tmp_path):
    """Test loading sequences from list format JSON."""
    from metainformant.life_events.core.utils import load_sequences_from_json

    data = [
        {
            "person_id": "person_001",
            "events": [
                {"event_type": "degree", "timestamp": "2010-06-01T00:00:00", "domain": "education", "attributes": {}}
            ],
            "metadata": {},
        }
    ]

    json_file = tmp_path / "events.json"
    with open(json_file, "w") as f:
        json.dump(data, f)

    sequences = load_sequences_from_json(json_file)
    assert len(sequences) == 1


def test_load_sequences_from_json_not_found():
    """Test loading from non-existent file."""
    from metainformant.life_events.core.utils import load_sequences_from_json

    with pytest.raises(FileNotFoundError):
        load_sequences_from_json("nonexistent.json")


def test_validate_sequence_valid():
    """Test validation of valid sequence."""
    sequence = EventSequence(person_id="person_001", events=[Event("degree", datetime(2010, 6, 1), "education")])

    is_valid, errors = validate_sequence(sequence)
    assert is_valid
    assert len([e for e in errors if not e.startswith("warning")]) == 0


def test_validate_sequence_empty_person_id():
    """Test validation of sequence with empty person_id."""
    sequence = EventSequence(person_id="", events=[Event("degree", datetime(2010, 6, 1), "education")])

    is_valid, errors = validate_sequence(sequence)
    assert not is_valid
    assert any("person_id cannot be empty" in e for e in errors)


def test_validate_sequence_no_events():
    """Test validation of sequence with no events."""
    sequence = EventSequence(person_id="person_001", events=[])

    is_valid, errors = validate_sequence(sequence)
    assert not is_valid
    assert any("no events" in e.lower() for e in errors)


def test_convert_sequences_to_tokens():
    """Test converting sequences to token format."""
    from metainformant.life_events.core.utils import convert_sequences_to_tokens

    sequences = [
        EventSequence(
            person_id="person_001",
            events=[
                Event("degree", datetime(2010, 6, 1), "education"),
                Event("job_change", datetime(2015, 3, 1), "occupation"),
            ],
        )
    ]

    tokens = convert_sequences_to_tokens(sequences)
    assert len(tokens) == 1
    assert tokens[0] == ["education:degree", "occupation:job_change"]


def test_get_event_statistics():
    """Test event statistics computation."""
    sequences = [
        EventSequence(
            person_id="person_001",
            events=[
                Event("degree", datetime(2010, 6, 1), "education"),
                Event("job_change", datetime(2015, 3, 1), "occupation"),
            ],
        ),
        EventSequence(person_id="person_002", events=[Event("diagnosis", datetime(2020, 1, 15), "health")]),
    ]

    stats = get_event_statistics(sequences)

    assert stats["total_sequences"] == 2
    assert stats["total_events"] == 3
    assert stats["domains"]["education"] == 1
    assert stats["domains"]["occupation"] == 1
    assert stats["domains"]["health"] == 1
    assert "temporal" in stats


def test_get_event_statistics_empty():
    """Test statistics with empty sequences."""
    stats = get_event_statistics([])

    assert stats["total_sequences"] == 0
    assert stats["total_events"] == 0
