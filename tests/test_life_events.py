"""Tests for life events module - event sequence data structures.

All tests follow NO_MOCKING_POLICY and use real implementations.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import pytest

from metainformant.life_events.core.events import Event, EventDatabase, EventSequence


class TestEvent:
    """Test Event dataclass functionality."""

    def test_event_creation_basic(self):
        """Test basic event creation."""
        event = Event(
            event_type="job_change",
            timestamp=datetime(2020, 1, 1, 12, 0, 0),
            domain="occupation"
        )
        assert event.event_type == "job_change"
        assert event.domain == "occupation"
        assert isinstance(event.timestamp, datetime)
        assert event.attributes == {}

    def test_event_creation_with_attributes(self):
        """Test event creation with additional attributes."""
        event = Event(
            event_type="diagnosis",
            timestamp=datetime(2020, 2, 1),
            domain="health",
            attributes={"condition": "diabetes", "severity": "mild"}
        )
        assert event.attributes["condition"] == "diabetes"
        assert event.attributes["severity"] == "mild"

    def test_event_creation_with_numeric_timestamp(self):
        """Test event creation with numeric timestamp."""
        event = Event(
            event_type="move",
            timestamp=1577836800.0,  # Unix timestamp
            domain="address"
        )
        assert isinstance(event.timestamp, datetime)
        assert event.timestamp.timestamp() == 1577836800.0

    def test_event_validation_empty_type(self):
        """Test that empty event_type raises ValueError."""
        with pytest.raises(ValueError, match="event_type cannot be empty"):
            Event(
                event_type="",
                timestamp=datetime.now(),
                domain="occupation"
            )

    def test_event_validation_empty_domain(self):
        """Test that empty domain raises ValueError."""
        with pytest.raises(ValueError, match="domain cannot be empty"):
            Event(
                event_type="job_change",
                timestamp=datetime.now(),
                domain=""
            )

    def test_event_to_dict(self):
        """Test event serialization to dictionary."""
        event = Event(
            event_type="job_change",
            timestamp=datetime(2020, 1, 1, 12, 0, 0),
            domain="occupation",
            attributes={"company": "ABC Corp"}
        )
        data = event.to_dict()
        assert data["event_type"] == "job_change"
        assert data["domain"] == "occupation"
        assert data["attributes"]["company"] == "ABC Corp"
        assert isinstance(data["timestamp"], str)  # ISO format

    def test_event_from_dict(self):
        """Test event deserialization from dictionary."""
        data = {
            "event_type": "job_change",
            "timestamp": "2020-01-01T12:00:00",
            "domain": "occupation",
            "attributes": {"company": "ABC Corp"}
        }
        event = Event.from_dict(data)
        assert event.event_type == "job_change"
        assert isinstance(event.timestamp, datetime)
        assert event.attributes["company"] == "ABC Corp"

    def test_event_from_dict_numeric_timestamp(self):
        """Test event from dict with numeric timestamp."""
        data = {
            "event_type": "move",
            "timestamp": 1577836800.0,
            "domain": "address"
        }
        event = Event.from_dict(data)
        assert isinstance(event.timestamp, datetime)


class TestEventSequence:
    """Test EventSequence container functionality."""

    def test_sequence_creation_basic(self):
        """Test basic sequence creation."""
        seq = EventSequence(person_id="person1")
        assert seq.person_id == "person1"
        assert seq.events == []
        assert seq.metadata == {}

    def test_sequence_with_events(self):
        """Test sequence with events."""
        events = [
            Event("job_change", datetime(2020, 1, 1), "occupation"),
            Event("move", datetime(2020, 2, 1), "address"),
        ]
        seq = EventSequence(person_id="person1", events=events)
        assert len(seq.events) == 2
        assert seq.events[0].event_type == "job_change"

    def test_sequence_auto_sorts_by_time(self):
        """Test that events are automatically sorted by timestamp."""
        events = [
            Event("event2", datetime(2020, 2, 1), "occupation"),
            Event("event1", datetime(2020, 1, 1), "occupation"),
        ]
        seq = EventSequence(person_id="person1", events=events)
        assert seq.events[0].event_type == "event1"
        assert seq.events[1].event_type == "event2"

    def test_sequence_filter_by_domain(self):
        """Test filtering sequence by domain."""
        events = [
            Event("job_change", datetime(2020, 1, 1), "occupation"),
            Event("diagnosis", datetime(2020, 2, 1), "health"),
            Event("move", datetime(2020, 3, 1), "address"),
        ]
        seq = EventSequence(person_id="person1", events=events)
        health_seq = seq.filter_by_domain("health")
        assert len(health_seq.events) == 1
        assert health_seq.events[0].event_type == "diagnosis"
        assert health_seq.person_id == "person1"

    def test_sequence_filter_by_time(self):
        """Test filtering sequence by time range."""
        events = [
            Event("event1", datetime(2020, 1, 1), "occupation"),
            Event("event2", datetime(2020, 2, 1), "occupation"),
            Event("event3", datetime(2020, 3, 1), "occupation"),
        ]
        seq = EventSequence(person_id="person1", events=events)
        filtered = seq.filter_by_time(
            start_time=datetime(2020, 1, 15),
            end_time=datetime(2020, 2, 15)
        )
        assert len(filtered.events) == 1
        assert filtered.events[0].event_type == "event2"

    def test_sequence_filter_by_time_numeric(self):
        """Test filtering with numeric timestamps."""
        events = [
            Event("event1", 1577836800.0, "occupation"),  # 2020-01-01
            Event("event2", 1580515200.0, "occupation"),  # 2020-02-01
        ]
        seq = EventSequence(person_id="person1", events=events)
        filtered = seq.filter_by_time(start_time=1578000000.0, end_time=1580600000.0)
        assert len(filtered.events) == 1
        assert filtered.events[0].event_type == "event2"

    def test_sequence_get_event_types(self):
        """Test getting unique event types."""
        events = [
            Event("job_change", datetime(2020, 1, 1), "occupation"),
            Event("job_change", datetime(2020, 2, 1), "occupation"),
            Event("move", datetime(2020, 3, 1), "address"),
        ]
        seq = EventSequence(person_id="person1", events=events)
        types = seq.get_event_types()
        assert len(types) == 2
        assert "job_change" in types
        assert "move" in types

    def test_sequence_get_domains(self):
        """Test getting unique domains."""
        events = [
            Event("event1", datetime(2020, 1, 1), "occupation"),
            Event("event2", datetime(2020, 2, 1), "health"),
            Event("event3", datetime(2020, 3, 1), "occupation"),
        ]
        seq = EventSequence(person_id="person1", events=events)
        domains = seq.get_domains()
        assert len(domains) == 2
        assert "occupation" in domains
        assert "health" in domains

    def test_sequence_to_dataframe(self):
        """Test converting sequence to DataFrame."""
        events = [
            Event("job_change", datetime(2020, 1, 1), "occupation", {"company": "ABC"}),
            Event("move", datetime(2020, 2, 1), "address"),
        ]
        seq = EventSequence(person_id="person1", events=events)
        df = seq.to_dataframe()
        assert len(df) == 2
        assert "person_id" in df.columns
        assert "event_type" in df.columns
        assert "timestamp" in df.columns
        assert "domain" in df.columns
        assert "company" in df.columns  # From attributes
        assert df.iloc[0]["company"] == "ABC"

    def test_sequence_to_dict(self):
        """Test sequence serialization."""
        events = [
            Event("job_change", datetime(2020, 1, 1), "occupation"),
        ]
        seq = EventSequence(person_id="person1", events=events, metadata={"source": "test"})
        data = seq.to_dict()
        assert data["person_id"] == "person1"
        assert len(data["events"]) == 1
        assert data["metadata"]["source"] == "test"

    def test_sequence_from_dict(self):
        """Test sequence deserialization."""
        data = {
            "person_id": "person1",
            "events": [
                {
                    "event_type": "job_change",
                    "timestamp": "2020-01-01T12:00:00",
                    "domain": "occupation",
                    "attributes": {}
                }
            ],
            "metadata": {"source": "test"}
        }
        seq = EventSequence.from_dict(data)
        assert seq.person_id == "person1"
        assert len(seq.events) == 1
        assert seq.metadata["source"] == "test"


class TestEventDatabase:
    """Test EventDatabase collection functionality."""

    def test_database_creation_empty(self):
        """Test creating empty database."""
        db = EventDatabase()
        assert db.sequences == []
        assert db.metadata == {}

    def test_database_creation_with_sequences(self):
        """Test creating database with sequences."""
        seq1 = EventSequence(person_id="person1", events=[
            Event("event1", datetime(2020, 1, 1), "occupation")
        ])
        seq2 = EventSequence(person_id="person2", events=[
            Event("event2", datetime(2020, 2, 1), "health")
        ])
        db = EventDatabase(sequences=[seq1, seq2])
        assert len(db.sequences) == 2

    def test_database_add_sequence(self):
        """Test adding sequence to database."""
        db = EventDatabase()
        seq = EventSequence(person_id="person1")
        db.add_sequence(seq)
        assert len(db.sequences) == 1
        assert db.sequences[0].person_id == "person1"

    def test_database_filter_by_domain(self):
        """Test filtering database by domain."""
        seq1 = EventSequence(person_id="person1", events=[
            Event("job_change", datetime(2020, 1, 1), "occupation"),
            Event("diagnosis", datetime(2020, 2, 1), "health"),
        ])
        seq2 = EventSequence(person_id="person2", events=[
            Event("move", datetime(2020, 3, 1), "address"),
        ])
        db = EventDatabase(sequences=[seq1, seq2])
        health_db = db.filter_by_domain("health")
        assert len(health_db.sequences) == 2  # Both sequences, but filtered
        assert len(health_db.sequences[0].events) == 1
        assert health_db.sequences[0].events[0].domain == "health"

    def test_database_get_statistics(self):
        """Test getting database statistics."""
        seq1 = EventSequence(person_id="person1", events=[
            Event("job_change", datetime(2020, 1, 1), "occupation"),
            Event("diagnosis", datetime(2020, 2, 1), "health"),
        ])
        seq2 = EventSequence(person_id="person2", events=[
            Event("move", datetime(2020, 3, 1), "address"),
        ])
        db = EventDatabase(sequences=[seq1, seq2])
        stats = db.get_statistics()
        assert stats["n_sequences"] == 2
        assert stats["n_events"] == 3
        assert stats["avg_events_per_sequence"] == 1.5
        assert stats["domains"]["occupation"] == 1
        assert stats["domains"]["health"] == 1
        assert stats["domains"]["address"] == 1

    def test_database_get_statistics_empty(self):
        """Test statistics for empty database."""
        db = EventDatabase()
        stats = db.get_statistics()
        assert stats["n_sequences"] == 0
        assert stats["n_events"] == 0
        assert stats["domains"] == {}
        assert stats["event_types"] == {}

    def test_database_to_dataframe(self):
        """Test converting database to DataFrame."""
        seq1 = EventSequence(person_id="person1", events=[
            Event("event1", datetime(2020, 1, 1), "occupation"),
        ])
        seq2 = EventSequence(person_id="person2", events=[
            Event("event2", datetime(2020, 2, 1), "health"),
        ])
        db = EventDatabase(sequences=[seq1, seq2])
        df = db.to_dataframe()
        assert len(df) == 2
        assert "person_id" in df.columns
        assert set(df["person_id"].values) == {"person1", "person2"}

    def test_database_to_dataframe_empty(self):
        """Test converting empty database to DataFrame."""
        db = EventDatabase()
        df = db.to_dataframe()
        assert len(df) == 0

    def test_database_to_dict(self):
        """Test database serialization."""
        seq = EventSequence(person_id="person1", events=[
            Event("event1", datetime(2020, 1, 1), "occupation"),
        ])
        db = EventDatabase(sequences=[seq], metadata={"source": "test"})
        data = db.to_dict()
        assert len(data["sequences"]) == 1
        assert data["metadata"]["source"] == "test"

    def test_database_from_dict(self):
        """Test database deserialization."""
        data = {
            "sequences": [
                {
                    "person_id": "person1",
                    "events": [
                        {
                            "event_type": "job_change",
                            "timestamp": "2020-01-01T12:00:00",
                            "domain": "occupation",
                            "attributes": {}
                        }
                    ],
                    "metadata": {}
                }
            ],
            "metadata": {"source": "test"}
        }
        db = EventDatabase.from_dict(data)
        assert len(db.sequences) == 1
        assert db.metadata["source"] == "test"
        assert db.sequences[0].person_id == "person1"


