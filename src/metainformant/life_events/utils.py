"""Utility functions for life course analysis.

This module provides helper functions for common patterns in life course analysis,
including data loading, validation, and convenience functions.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, List

from ..core import io
from .events import Event, EventDatabase, EventSequence


def load_sequences_from_json(path: str | Path) -> List[EventSequence]:
    """Load event sequences from JSON file.
    
    Supports both single sequence and database formats.
    
    Args:
        path: Path to JSON file containing event sequences
        
    Returns:
        List of EventSequence objects
        
    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file format is invalid
        
    Examples:
        >>> sequences = load_sequences_from_json("data/events.json")
        >>> len(sequences)
        10
    """
    path = Path(path)
    
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")
    
    data = io.load_json(path)
    
    # Handle database format
    if isinstance(data, dict) and "sequences" in data:
        database = EventDatabase.from_dict(data)
        return database.sequences
    
    # Handle list of sequences
    if isinstance(data, list):
        sequences = []
        for item in data:
            if isinstance(item, dict):
                if "person_id" in item or "events" in item:
                    sequences.append(EventSequence.from_dict(item))
                else:
                    # Try to create sequence from flat structure
                    raise ValueError("Invalid sequence format in list")
        return sequences
    
    # Handle single sequence
    if isinstance(data, dict) and ("person_id" in data or "events" in data):
        return [EventSequence.from_dict(data)]
    
    raise ValueError(f"Invalid file format in {path}. Expected list of sequences or database format.")


def validate_sequence(sequence: EventSequence) -> tuple[bool, List[str]]:
    """Validate event sequence for common issues.
    
    Args:
        sequence: EventSequence to validate
        
    Returns:
        Tuple of (is_valid, list_of_errors)
        
    Examples:
        >>> seq = EventSequence("p1", [Event("degree", datetime(2010, 1, 1), "education")])
        >>> is_valid, errors = validate_sequence(seq)
        >>> is_valid
        True
    """
    errors = []
    
    if not sequence.person_id:
        errors.append("person_id cannot be empty")
    
    if not sequence.events:
        errors.append("sequence has no events")
    
    # Check for duplicate timestamps (warning, not error)
    timestamps = []
    for event in sequence.events:
        if isinstance(event.timestamp, float):
            timestamps.append(event.timestamp)
        else:
            timestamps.append(event.timestamp.timestamp())
    
    if len(timestamps) != len(set(timestamps)):
        errors.append("warning: duplicate timestamps detected")
    
    # Check for valid domains
    valid_domains = ["health", "education", "occupation", "income", "address", "other"]
    for event in sequence.events:
        if event.domain not in valid_domains:
            errors.append(f"warning: unknown domain '{event.domain}' in event {event.event_type}")
    
    return len([e for e in errors if not e.startswith("warning")]) == 0, errors


def convert_sequences_to_tokens(sequences: List[EventSequence]) -> List[List[str]]:
    """Convert event sequences to token format for embedding/ML.
    
    Args:
        sequences: List of EventSequence objects
        
    Returns:
        List of token sequences, where each sequence is a list of "domain:event_type" strings
        
    Examples:
        >>> sequences = [EventSequence("p1", [Event("degree", datetime(2010, 1, 1), "education")])]
        >>> tokens = convert_sequences_to_tokens(sequences)
        >>> tokens[0]
        ['education:degree']
    """
    return [
        [f"{e.domain}:{e.event_type}" for e in seq.events]
        for seq in sequences
    ]


def get_event_statistics(sequences: List[EventSequence]) -> dict[str, Any]:
    """Compute comprehensive statistics about event sequences.
    
    Args:
        sequences: List of EventSequence objects
        
    Returns:
        Dictionary with statistics including counts, domains, event types, temporal patterns
        
    Examples:
        >>> sequences = [EventSequence("p1", [Event("degree", datetime(2010, 1, 1), "education")])]
        >>> stats = get_event_statistics(sequences)
        >>> "total_events" in stats
        True
    """
    if not sequences:
        return {
            "total_sequences": 0,
            "total_events": 0,
            "domains": {},
            "event_types": {},
        }
    
    domain_counts = {}
    event_type_counts = {}
    total_events = 0
    
    for seq in sequences:
        total_events += len(seq.events)
        for event in seq.events:
            domain_counts[event.domain] = domain_counts.get(event.domain, 0) + 1
            event_type_counts[event.event_type] = event_type_counts.get(event.event_type, 0) + 1
    
    # Temporal statistics
    all_timestamps = []
    for seq in sequences:
        for event in seq.events:
            if isinstance(event.timestamp, float):
                all_timestamps.append(event.timestamp)
            else:
                all_timestamps.append(event.timestamp.timestamp())
    
    temporal_stats = {}
    if all_timestamps:
        temporal_stats = {
            "earliest_event": min(all_timestamps),
            "latest_event": max(all_timestamps),
            "time_span": max(all_timestamps) - min(all_timestamps),
        }
    
    return {
        "total_sequences": len(sequences),
        "total_events": total_events,
        "avg_events_per_sequence": total_events / len(sequences) if sequences else 0.0,
        "domains": domain_counts,
        "event_types": event_type_counts,
        "temporal": temporal_stats,
    }

