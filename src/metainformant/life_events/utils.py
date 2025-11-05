"""Utility functions for life course analysis.

This module provides helper functions for common patterns in life course analysis,
including data loading, validation, and convenience functions.
"""

from __future__ import annotations

from datetime import datetime, timedelta
from pathlib import Path
from typing import Any, List, Optional

import numpy as np

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


def generate_synthetic_life_events(
    n_sequences: int = 100,
    min_events_per_sequence: int = 5,
    max_events_per_sequence: int = 30,
    start_date: Optional[datetime] = None,
    end_date: Optional[datetime] = None,
    domains: Optional[List[str]] = None,
    event_types_by_domain: Optional[dict[str, List[str]]] = None,
    generate_outcomes: bool = False,
    outcome_relationship: str = "random",
    random_state: Optional[int] = None
) -> tuple[List[EventSequence], Optional[np.ndarray]]:
    """Generate synthetic life event sequences for testing and demos.
    
    Creates realistic event sequences with temporal patterns, multiple domains,
    and optionally outcome labels that relate to event patterns.
    
    Args:
        n_sequences: Number of event sequences to generate
        min_events_per_sequence: Minimum events per sequence
        max_events_per_sequence: Maximum events per sequence
        start_date: Start date for events (defaults to 10 years ago)
        end_date: End date for events (defaults to today)
        domains: List of domains to use (defaults to all standard domains)
        event_types_by_domain: Dictionary mapping domains to event types
            (defaults to realistic event types)
        generate_outcomes: Whether to generate outcome labels
        outcome_relationship: How outcomes relate to events:
            - "random": Random binary outcomes
            - "health_focused": Higher outcome for more health events
            - "education_focused": Higher outcome for education events
            - "complex": Complex pattern based on multiple domains
        random_state: Random seed for reproducibility
        
    Returns:
        Tuple of (list of EventSequence objects, optional outcome array)
        
    Examples:
        >>> from datetime import datetime
        >>> sequences, outcomes = generate_synthetic_life_events(
        ...     n_sequences=10,
        ...     random_state=42
        ... )
        >>> len(sequences)
        10
        >>> len(outcomes)
        10
    """
    if random_state is not None:
        np.random.seed(random_state)
    
    # Default domains
    if domains is None:
        domains = ["health", "education", "occupation", "income", "address", "other"]
    
    # Default event types by domain
    if event_types_by_domain is None:
        event_types_by_domain = {
            "health": [
                "diagnosis", "hospitalization", "medication_start", "surgery",
                "therapy", "checkup", "vaccination", "recovery"
            ],
            "education": [
                "degree", "certification", "course_completion", "enrollment",
                "graduation", "training", "workshop"
            ],
            "occupation": [
                "job_change", "promotion", "layoff", "retirement",
                "career_shift", "contract_signed", "project_completion"
            ],
            "income": [
                "raise", "bonus", "salary_change", "investment",
                "windfall", "loss", "budget_adjustment"
            ],
            "address": [
                "move", "relocation", "address_change", "home_purchase",
                "rental_start", "rental_end"
            ],
            "other": [
                "marriage", "divorce", "birth", "death",
                "legal_event", "travel", "milestone"
            ]
        }
    
    # Default date range
    if end_date is None:
        end_date = datetime.now()
    if start_date is None:
        start_date = end_date - timedelta(days=365 * 10)  # 10 years ago
    
    # Convert to timestamps for easier manipulation
    start_timestamp = start_date.timestamp()
    end_timestamp = end_date.timestamp()
    time_span = end_timestamp - start_timestamp
    
    sequences = []
    outcomes = None
    
    if generate_outcomes:
        outcomes = np.zeros(n_sequences, dtype=np.float64)
    
    for i in range(n_sequences):
        person_id = f"person_{i+1:04d}"
        
        # Generate number of events for this sequence
        n_events = np.random.randint(min_events_per_sequence, max_events_per_sequence + 1)
        
        # Generate events with temporal distribution
        events = []
        event_timestamps = sorted(np.random.uniform(start_timestamp, end_timestamp, n_events))
        
        # Track domain/event type usage for outcome generation
        domain_counts = {domain: 0 for domain in domains}
        event_type_list = []
        
        for timestamp in event_timestamps:
            # Select domain (weighted towards earlier domains for realism)
            domain_weights = np.array([1.0 / (j + 1) for j in range(len(domains))])
            domain_weights = domain_weights / domain_weights.sum()
            domain = np.random.choice(domains, p=domain_weights)
            
            # Select event type from domain
            if domain in event_types_by_domain and event_types_by_domain[domain]:
                event_type = np.random.choice(event_types_by_domain[domain])
            else:
                event_type = f"event_{np.random.randint(1000)}"
            
            # Create event
            event = Event(
                event_type=event_type,
                timestamp=datetime.fromtimestamp(timestamp),
                domain=domain,
                attributes={}
            )
            events.append(event)
            
            # Track for outcome generation
            domain_counts[domain] += 1
            event_type_list.append(f"{domain}:{event_type}")
        
        # Create sequence
        sequence = EventSequence(
            person_id=person_id,
            events=events,
            metadata={"synthetic": True, "generated_at": datetime.now().isoformat()}
        )
        sequences.append(sequence)
        
        # Generate outcome if requested
        if generate_outcomes:
            if outcome_relationship == "random":
                outcomes[i] = np.random.choice([0, 1])
            elif outcome_relationship == "health_focused":
                # Higher outcome for more health events
                health_ratio = domain_counts.get("health", 0) / n_events
                outcomes[i] = 1 if health_ratio > 0.3 or np.random.random() < 0.3 else 0
            elif outcome_relationship == "education_focused":
                # Higher outcome for education events
                education_ratio = domain_counts.get("education", 0) / n_events
                outcomes[i] = 1 if education_ratio > 0.2 or np.random.random() < 0.2 else 0
            elif outcome_relationship == "complex":
                # Complex pattern: positive if has education AND health events
                has_education = domain_counts.get("education", 0) > 0
                has_health = domain_counts.get("health", 0) > 0
                has_occupation = domain_counts.get("occupation", 0) > 0
                outcomes[i] = 1 if (has_education and has_health) or (has_occupation and np.random.random() < 0.4) else 0
            else:
                # Default to random
                outcomes[i] = np.random.choice([0, 1])
    
    return sequences, outcomes

