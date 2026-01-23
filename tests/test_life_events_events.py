"""Tests for life_events event data structures."""

from __future__ import annotations

from datetime import datetime

import pytest

from metainformant.life_events import Event, EventDatabase, EventSequence


def test_event_creation():
    """Test creating Event objects."""
    event = Event(
        event_type="diagnosis", timestamp=datetime(2020, 1, 15), domain="health", attributes={"condition": "diabetes"}
    )

    assert event.event_type == "diagnosis"
    assert event.domain == "health"
    assert event.attributes["condition"] == "diabetes"
    assert isinstance(event.timestamp, datetime)


def test_event_validation():
    """Test Event validation."""
    # Empty event_type should raise error
    with pytest.raises(ValueError, match="event_type cannot be empty"):
        Event("", datetime(2020, 1, 1), "health")

    # Empty domain should raise error
    with pytest.raises(ValueError, match="domain cannot be empty"):
        Event("diagnosis", datetime(2020, 1, 1), "")


def test_event_to_from_dict():
    """Test Event serialization."""
    event = Event(
        event_type="job_change",
        timestamp=datetime(2015, 3, 1),
        domain="occupation",
        attributes={"company": "TechCorp"},
        person_id="person_001",
    )

    event_dict = event.to_dict()
    assert event_dict["event_type"] == "job_change"
    assert event_dict["domain"] == "occupation"

    # Reconstruct from dict
    event2 = Event.from_dict(event_dict)
    assert event2.event_type == event.event_type
    assert event2.domain == event.domain
    assert event2.attributes == event.attributes


def test_event_sequence_creation():
    """Test creating EventSequence objects."""
    events = [
        Event("degree", datetime(2010, 6, 1), "education", {"degree": "BS"}),
        Event("job_change", datetime(2015, 3, 1), "occupation", {"company": "TechCorp"}),
        Event("diagnosis", datetime(2020, 1, 15), "health", {"condition": "diabetes"}),
    ]

    sequence = EventSequence(person_id="person_001", events=events)

    assert len(sequence.events) == 3
    assert sequence.person_id == "person_001"
    # Events should be sorted by timestamp
    assert sequence.events[0].event_type == "degree"


def test_event_sequence_filtering():
    """Test EventSequence filtering methods."""
    events = [
        Event("degree", datetime(2010, 6, 1), "education"),
        Event("job_change", datetime(2015, 3, 1), "occupation"),
        Event("diagnosis", datetime(2020, 1, 15), "health"),
    ]

    sequence = EventSequence(person_id="person_001", events=events)

    # Filter by domain
    health_seq = sequence.filter_by_domain("health")
    assert len(health_seq.events) == 1
    assert health_seq.events[0].domain == "health"

    # Filter by time
    recent_seq = sequence.filter_by_time(start_time=datetime(2015, 1, 1), end_time=datetime(2025, 1, 1))
    assert len(recent_seq.events) == 2


def test_event_sequence_get_methods():
    """Test EventSequence getter methods."""
    events = [
        Event("degree", datetime(2010, 6, 1), "education"),
        Event("job_change", datetime(2015, 3, 1), "occupation"),
        Event("diagnosis", datetime(2020, 1, 15), "health"),
    ]

    sequence = EventSequence(person_id="person_001", events=events)

    event_types = sequence.get_event_types()
    assert "degree" in event_types
    assert "job_change" in event_types
    assert "diagnosis" in event_types

    domains = sequence.get_domains()
    assert "education" in domains
    assert "occupation" in domains
    assert "health" in domains


def test_event_sequence_to_dataframe():
    """Test EventSequence DataFrame conversion."""
    events = [
        Event("degree", datetime(2010, 6, 1), "education", {"degree": "BS"}),
        Event("job_change", datetime(2015, 3, 1), "occupation"),
    ]

    sequence = EventSequence(person_id="person_001", events=events)
    df = sequence.to_dataframe()

    assert len(df) == 2
    assert "person_id" in df.columns
    assert "event_type" in df.columns
    assert "domain" in df.columns
    assert df.iloc[0]["event_type"] == "degree"


def test_event_sequence_to_from_dict():
    """Test EventSequence serialization."""
    events = [
        Event("degree", datetime(2010, 6, 1), "education"),
        Event("job_change", datetime(2015, 3, 1), "occupation"),
    ]

    sequence = EventSequence(person_id="person_001", events=events, metadata={"source": "test"})

    seq_dict = sequence.to_dict()
    assert seq_dict["person_id"] == "person_001"
    assert len(seq_dict["events"]) == 2
    assert seq_dict["metadata"]["source"] == "test"

    # Reconstruct
    sequence2 = EventSequence.from_dict(seq_dict)
    assert len(sequence2.events) == 2
    assert sequence2.metadata["source"] == "test"


def test_event_database_creation():
    """Test creating EventDatabase."""
    seq1 = EventSequence(person_id="person_001", events=[Event("degree", datetime(2010, 6, 1), "education")])
    seq2 = EventSequence(person_id="person_002", events=[Event("job_change", datetime(2015, 3, 1), "occupation")])

    database = EventDatabase(sequences=[seq1, seq2])
    assert len(database.sequences) == 2


def test_event_database_add_sequence():
    """Test adding sequences to database."""
    database = EventDatabase()
    seq = EventSequence(person_id="person_001", events=[Event("degree", datetime(2010, 6, 1), "education")])

    database.add_sequence(seq)
    assert len(database.sequences) == 1


def test_event_database_filter_by_domain():
    """Test EventDatabase domain filtering."""
    seq1 = EventSequence(
        person_id="person_001",
        events=[
            Event("degree", datetime(2010, 6, 1), "education"),
            Event("diagnosis", datetime(2020, 1, 15), "health"),
        ],
    )

    database = EventDatabase(sequences=[seq1])
    health_db = database.filter_by_domain("health")

    assert len(health_db.sequences) == 1
    assert len(health_db.sequences[0].events) == 1
    assert health_db.sequences[0].events[0].domain == "health"


def test_event_database_statistics():
    """Test EventDatabase statistics."""
    seq1 = EventSequence(
        person_id="person_001",
        events=[
            Event("degree", datetime(2010, 6, 1), "education"),
            Event("job_change", datetime(2015, 3, 1), "occupation"),
        ],
    )
    seq2 = EventSequence(person_id="person_002", events=[Event("diagnosis", datetime(2020, 1, 15), "health")])

    database = EventDatabase(sequences=[seq1, seq2])
    stats = database.get_statistics()

    assert stats["n_sequences"] == 2
    assert stats["n_events"] == 3
    assert stats["domains"]["education"] == 1
    assert stats["domains"]["occupation"] == 1
    assert stats["domains"]["health"] == 1


def test_event_database_to_dataframe():
    """Test EventDatabase DataFrame conversion."""
    seq1 = EventSequence(person_id="person_001", events=[Event("degree", datetime(2010, 6, 1), "education")])
    seq2 = EventSequence(person_id="person_002", events=[Event("job_change", datetime(2015, 3, 1), "occupation")])

    database = EventDatabase(sequences=[seq1, seq2])
    df = database.to_dataframe()

    assert len(df) == 2
    assert set(df["person_id"]) == {"person_001", "person_002"}


def test_event_database_to_from_dict():
    """Test EventDatabase serialization."""
    seq1 = EventSequence(person_id="person_001", events=[Event("degree", datetime(2010, 6, 1), "education")])

    database = EventDatabase(sequences=[seq1], metadata={"source": "test"})
    db_dict = database.to_dict()

    assert len(db_dict["sequences"]) == 1
    assert db_dict["metadata"]["source"] == "test"

    # Reconstruct
    database2 = EventDatabase.from_dict(db_dict)
    assert len(database2.sequences) == 1
    assert database2.metadata["source"] == "test"
