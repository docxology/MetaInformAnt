"""Event sequence data structures for life course analysis.

This module provides core data structures for representing temporal event sequences,
including individual events, event sequences, and collections of sequences.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Dict, List, Optional, Sequence

import numpy as np
import pandas as pd


@dataclass
class Event:
    """Individual life event with temporal and domain information.
    
    Attributes:
        event_type: Type/category of event (e.g., "job_change", "diagnosis")
        timestamp: When the event occurred (datetime or numeric)
        domain: Domain category (health, education, occupation, income, address)
        attributes: Additional event-specific attributes as dictionary
        person_id: Optional identifier for the person (if not in sequence)
    """
    event_type: str
    timestamp: datetime | float | int
    domain: str
    attributes: Dict[str, Any] = field(default_factory=dict)
    person_id: Optional[str] = None
    
    def __post_init__(self) -> None:
        """Validate event data."""
        if not self.event_type:
            raise ValueError("event_type cannot be empty")
        if not self.domain:
            raise ValueError("domain cannot be empty")
        if self.domain not in ["health", "education", "occupation", "income", "address", "other"]:
            # Allow other domains but warn
            pass
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert event to dictionary representation.
        
        Returns:
            Dictionary representation of event
        """
        return {
            "event_type": self.event_type,
            "timestamp": self.timestamp.isoformat() if isinstance(self.timestamp, datetime) else self.timestamp,
            "domain": self.domain,
            "attributes": self.attributes,
            "person_id": self.person_id,
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "Event":
        """Create Event from dictionary.
        
        Args:
            data: Dictionary with event data
            
        Returns:
            Event object
        """
        timestamp = data["timestamp"]
        if isinstance(timestamp, str):
            try:
                timestamp = datetime.fromisoformat(timestamp)
            except (ValueError, AttributeError):
                # Try numeric conversion
                timestamp = float(timestamp)
        
        return cls(
            event_type=data["event_type"],
            timestamp=timestamp,
            domain=data["domain"],
            attributes=data.get("attributes", {}),
            person_id=data.get("person_id"),
        )


@dataclass
class EventSequence:
    """Container for temporal event sequence of a single person.
    
    Attributes:
        person_id: Unique identifier for the person
        events: List of Event objects in temporal order
        metadata: Additional metadata about the sequence
    """
    person_id: str
    events: List[Event] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    def __post_init__(self) -> None:
        """Sort events by timestamp."""
        self.sort_by_time()
    
    def sort_by_time(self) -> None:
        """Sort events by timestamp."""
        def get_timestamp(event: Event) -> float:
            if isinstance(event.timestamp, datetime):
                return event.timestamp.timestamp()
            return float(event.timestamp)
        
        self.events.sort(key=get_timestamp)
    
    def filter_by_domain(self, domain: str) -> "EventSequence":
        """Create new sequence with only events from specified domain.
        
        Args:
            domain: Domain to filter by
            
        Returns:
            New EventSequence with filtered events
        """
        filtered_events = [e for e in self.events if e.domain == domain]
        return EventSequence(
            person_id=self.person_id,
            events=filtered_events,
            metadata=self.metadata.copy()
        )
    
    def filter_by_time(
        self,
        start_time: Optional[datetime | float] = None,
        end_time: Optional[datetime | float] = None
    ) -> "EventSequence":
        """Create new sequence with events in time range.
        
        Args:
            start_time: Start of time range (inclusive)
            end_time: End of time range (inclusive)
            
        Returns:
            New EventSequence with filtered events
        """
        def get_timestamp(event: Event) -> float:
            if isinstance(event.timestamp, datetime):
                return event.timestamp.timestamp()
            return float(event.timestamp)
        
        def in_range(event: Event) -> bool:
            ts = get_timestamp(event)
            if start_time is not None:
                start_ts = start_time.timestamp() if isinstance(start_time, datetime) else float(start_time)
                if ts < start_ts:
                    return False
            if end_time is not None:
                end_ts = end_time.timestamp() if isinstance(end_time, datetime) else float(end_time)
                if ts > end_ts:
                    return False
            return True
        
        filtered_events = [e for e in self.events if in_range(e)]
        return EventSequence(
            person_id=self.person_id,
            events=filtered_events,
            metadata=self.metadata.copy()
        )
    
    def get_event_types(self) -> List[str]:
        """Get unique event types in sequence.
        
        Returns:
            List of unique event types
        """
        return sorted(list(set(e.event_type for e in self.events)))
    
    def get_domains(self) -> List[str]:
        """Get unique domains in sequence.
        
        Returns:
            List of unique domains
        """
        return sorted(list(set(e.domain for e in self.events)))
    
    def to_dataframe(self) -> pd.DataFrame:
        """Convert sequence to pandas DataFrame.
        
        Returns:
            DataFrame with columns: person_id, event_type, timestamp, domain, attributes
        """
        rows = []
        for event in self.events:
            row = {
                "person_id": self.person_id,
                "event_type": event.event_type,
                "timestamp": event.timestamp,
                "domain": event.domain,
                **event.attributes
            }
            rows.append(row)
        
        return pd.DataFrame(rows)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert sequence to dictionary representation.
        
        Returns:
            Dictionary representation of sequence
        """
        return {
            "person_id": self.person_id,
            "events": [e.to_dict() for e in self.events],
            "metadata": self.metadata,
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "EventSequence":
        """Create EventSequence from dictionary.
        
        Args:
            data: Dictionary with sequence data
            
        Returns:
            EventSequence object
        """
        events = [Event.from_dict(e) for e in data.get("events", [])]
        return cls(
            person_id=data["person_id"],
            events=events,
            metadata=data.get("metadata", {})
        )


class EventDatabase:
    """Collection of multiple event sequences.
    
    Provides batch processing, statistics, and domain filtering across
    multiple event sequences.
    
    Attributes:
        sequences: List of EventSequence objects
        metadata: Database-level metadata
    """
    
    def __init__(
        self,
        sequences: Optional[Sequence[EventSequence]] = None,
        metadata: Optional[Dict[str, Any]] = None
    ):
        """Initialize EventDatabase.
        
        Args:
            sequences: List of event sequences
            metadata: Database metadata
        """
        self.sequences: List[EventSequence] = list(sequences) if sequences else []
        self.metadata: Dict[str, Any] = metadata or {}
    
    def add_sequence(self, sequence: EventSequence) -> None:
        """Add event sequence to database.
        
        Args:
            sequence: EventSequence to add
        """
        self.sequences.append(sequence)
    
    def filter_by_domain(self, domain: str) -> "EventDatabase":
        """Create new database with sequences filtered by domain.
        
        Args:
            domain: Domain to filter by
            
        Returns:
            New EventDatabase with filtered sequences
        """
        filtered_sequences = [seq.filter_by_domain(domain) for seq in self.sequences]
        return EventDatabase(sequences=filtered_sequences, metadata=self.metadata.copy())
    
    def get_statistics(self) -> Dict[str, Any]:
        """Compute database-level statistics.
        
        Returns:
            Dictionary with statistics about the database
        """
        if not self.sequences:
            return {
                "n_sequences": 0,
                "n_events": 0,
                "domains": {},
                "event_types": {},
            }
        
        total_events = sum(len(seq.events) for seq in self.sequences)
        domain_counts: Dict[str, int] = {}
        event_type_counts: Dict[str, int] = {}
        
        for seq in self.sequences:
            for event in seq.events:
                domain_counts[event.domain] = domain_counts.get(event.domain, 0) + 1
                event_type_counts[event.event_type] = event_type_counts.get(event.event_type, 0) + 1
        
        return {
            "n_sequences": len(self.sequences),
            "n_events": total_events,
            "avg_events_per_sequence": total_events / len(self.sequences) if self.sequences else 0.0,
            "domains": domain_counts,
            "event_types": event_type_counts,
        }
    
    def to_dataframe(self) -> pd.DataFrame:
        """Convert all sequences to single DataFrame.
        
        Returns:
            Combined DataFrame of all sequences
        """
        dfs = [seq.to_dataframe() for seq in self.sequences]
        if not dfs:
            return pd.DataFrame()
        return pd.concat(dfs, ignore_index=True)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert database to dictionary representation.
        
        Returns:
            Dictionary representation of database
        """
        return {
            "sequences": [seq.to_dict() for seq in self.sequences],
            "metadata": self.metadata,
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "EventDatabase":
        """Create EventDatabase from dictionary.
        
        Args:
            data: Dictionary with database data
            
        Returns:
            EventDatabase object
        """
        sequences = [EventSequence.from_dict(s) for s in data.get("sequences", [])]
        return cls(sequences=sequences, metadata=data.get("metadata", {}))

