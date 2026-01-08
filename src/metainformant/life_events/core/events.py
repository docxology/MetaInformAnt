"""Life events data structures and management.

This module provides classes for representing and managing life course events,
event sequences, and collections of event data.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

from metainformant.core import logging

logger = logging.get_logger(__name__)


@dataclass
class Event:
    """Represents a single life event with temporal and domain information.

    Attributes:
        event_type: Type of event (e.g., 'job_change', 'marriage', 'diagnosis')
        timestamp: When the event occurred
        domain: Life domain this event belongs to (e.g., 'occupation', 'health', 'family')
        attributes: Additional event-specific attributes
    """
    event_type: str
    timestamp: datetime
    domain: str
    attributes: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """Convert event to dictionary for serialization."""
        return {
            'event_type': self.event_type,
            'timestamp': self.timestamp.isoformat(),
            'domain': self.domain,
            'attributes': self.attributes
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> Event:
        """Create event from dictionary."""
        return cls(
            event_type=data['event_type'],
            timestamp=datetime.fromisoformat(data['timestamp']),
            domain=data['domain'],
            attributes=data.get('attributes', {})
        )

    def __str__(self) -> str:
        return f"{self.event_type} ({self.domain}) at {self.timestamp}"


class EventSequence:
    """Represents a sequence of life events for a single individual.

    Provides methods for filtering, querying, and analyzing temporal event sequences.
    """

    def __init__(self, person_id: str, events: List[Event]):
        """Initialize event sequence.

        Args:
            person_id: Unique identifier for the individual
            events: List of events in chronological order
        """
        self.person_id = person_id
        self.events = sorted(events, key=lambda e: e.timestamp)

    def filter_by_domain(self, domain: str) -> EventSequence:
        """Filter events by domain.

        Args:
            domain: Domain to filter by

        Returns:
            New EventSequence with filtered events
        """
        filtered_events = [e for e in self.events if e.domain == domain]
        return EventSequence(self.person_id, filtered_events)

    def filter_by_time_range(self, start: datetime, end: datetime) -> EventSequence:
        """Filter events by time range.

        Args:
            start: Start timestamp (inclusive)
            end: End timestamp (inclusive)

        Returns:
            New EventSequence with filtered events
        """
        filtered_events = [
            e for e in self.events
            if start <= e.timestamp <= end
        ]
        return EventSequence(self.person_id, filtered_events)

    def get_event_types(self) -> List[str]:
        """Get list of unique event types."""
        return list(set(e.event_type for e in self.events))

    def get_domains(self) -> List[str]:
        """Get list of unique domains."""
        return list(set(e.domain for e in self.events))

    def to_dataframe(self) -> pd.DataFrame:
        """Convert sequence to pandas DataFrame."""
        data = []
        for event in self.events:
            row = {
                'person_id': self.person_id,
                'event_type': event.event_type,
                'timestamp': event.timestamp,
                'domain': event.domain,
                **event.attributes
            }
            data.append(row)
        return pd.DataFrame(data)

    def to_dict(self) -> Dict[str, Any]:
        """Convert sequence to dictionary for serialization."""
        return {
            'person_id': self.person_id,
            'events': [e.to_dict() for e in self.events]
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> EventSequence:
        """Create sequence from dictionary."""
        events = [Event.from_dict(e_data) for e_data in data['events']]
        return cls(data['person_id'], events)

    def __len__(self) -> int:
        return len(self.events)

    def __getitem__(self, index: int) -> Event:
        return self.events[index]

    def __str__(self) -> str:
        return f"EventSequence(person_id={self.person_id}, events={len(self.events)})"


class EventDatabase:
    """Database of event sequences with batch processing capabilities.

    Manages collections of EventSequence objects and provides methods for
    querying, filtering, and batch operations.
    """

    def __init__(self, sequences: List[EventSequence]):
        """Initialize database.

        Args:
            sequences: List of EventSequence objects
        """
        self.sequences = sequences
        self._build_index()

    def _build_index(self):
        """Build internal indices for efficient querying."""
        self.person_index = {seq.person_id: seq for seq in self.sequences}

    def get_sequence(self, person_id: str) -> Optional[EventSequence]:
        """Get sequence for a specific person.

        Args:
            person_id: Person identifier

        Returns:
            EventSequence if found, None otherwise
        """
        return self.person_index.get(person_id)

    def filter_by_domain(self, domain: str) -> EventDatabase:
        """Filter all sequences to only include events from specified domain.

        Args:
            domain: Domain to filter by

        Returns:
            New EventDatabase with filtered sequences
        """
        filtered_sequences = [
            seq.filter_by_domain(domain)
            for seq in self.sequences
        ]
        # Remove empty sequences
        filtered_sequences = [seq for seq in filtered_sequences if len(seq) > 0]
        return EventDatabase(filtered_sequences)

    def filter_by_time_range(self, start: datetime, end: datetime) -> EventDatabase:
        """Filter all sequences to events within time range.

        Args:
            start: Start timestamp
            end: End timestamp

        Returns:
            New EventDatabase with filtered sequences
        """
        filtered_sequences = [
            seq.filter_by_time_range(start, end)
            for seq in self.sequences
        ]
        # Remove empty sequences
        filtered_sequences = [seq for seq in filtered_sequences if len(seq) > 0]
        return EventDatabase(filtered_sequences)

    def get_statistics(self) -> Dict[str, Any]:
        """Get database statistics.

        Returns:
            Dictionary with database statistics
        """
        total_events = sum(len(seq) for seq in self.sequences)
        event_types = set()
        domains = set()

        for seq in self.sequences:
            event_types.update(seq.get_event_types())
            domains.update(seq.get_domains())

        return {
            'num_sequences': len(self.sequences),
            'total_events': total_events,
            'unique_event_types': len(event_types),
            'unique_domains': len(domains),
            'avg_events_per_sequence': total_events / len(self.sequences) if self.sequences else 0
        }

    def to_dataframe(self) -> pd.DataFrame:
        """Convert entire database to pandas DataFrame."""
        dfs = [seq.to_dataframe() for seq in self.sequences]
        return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()

    def save_json(self, path: str | Path) -> None:
        """Save database to JSON file.

        Args:
            path: Output file path
        """
        data = {
            'sequences': [seq.to_dict() for seq in self.sequences],
            'statistics': self.get_statistics()
        }

        with open(path, 'w') as f:
            json.dump(data, f, indent=2, default=str)

        logger.info(f"Saved {len(self.sequences)} sequences to {path}")

    @classmethod
    def load_json(cls, path: str | Path) -> EventDatabase:
        """Load database from JSON file.

        Args:
            path: Input file path

        Returns:
            EventDatabase instance
        """
        with open(path, 'r') as f:
            data = json.load(f)

        sequences = [EventSequence.from_dict(seq_data) for seq_data in data['sequences']]
        return cls(sequences)

    def __len__(self) -> int:
        return len(self.sequences)

    def __getitem__(self, index: int) -> EventSequence:
        return self.sequences[index]

    def __iter__(self):
        return iter(self.sequences)

    def __str__(self) -> str:
        stats = self.get_statistics()
        return f"EventDatabase(sequences={stats['num_sequences']}, events={stats['total_events']})"
