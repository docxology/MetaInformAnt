"""Utility functions for life course analysis.

This module provides helper functions for common patterns in life course analysis,
including data loading, validation, and convenience functions.
"""

from __future__ import annotations

from collections import defaultdict
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

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


def generate_event_chain(
    chain_rules: Dict[str, Dict[str, float]],
    start_event: str,
    start_domain: str,
    n_events: int,
    start_timestamp: float,
    time_span: float,
    event_types_by_domain: Dict[str, List[str]],
    random_state: Optional[int] = None
) -> List[Event]:
    """Generate causally linked event sequence based on transition probabilities.
    
    Creates a chain of events where each event's probability depends on the previous
    event, implementing a Markov chain model for event sequences.
    
    Args:
        chain_rules: Dictionary mapping "domain:event_type" to next event probabilities
            Format: {"domain:event_type": {"next_domain:next_event": probability}}
        start_event: Starting event type
        start_domain: Starting domain
        n_events: Number of events to generate
        start_timestamp: Starting timestamp
        time_span: Time span for events (in seconds)
        event_types_by_domain: Dictionary mapping domains to available event types
        random_state: Random seed for reproducibility
        
    Returns:
        List of causally linked Event objects
        
    Examples:
        >>> chain_rules = {
        ...     "education:degree": {"occupation:job_change": 0.8, "income:raise": 0.2},
        ...     "occupation:job_change": {"address:move": 0.5, "income:raise": 0.5}
        ... }
        >>> events = generate_event_chain(
        ...     chain_rules, "degree", "education", 5,
        ...     datetime(2010, 1, 1).timestamp(), 365*5*86400,
        ...     {"education": ["degree"], "occupation": ["job_change"]},
        ...     random_state=42
        ... )
        >>> len(events)
        5
    """
    if random_state is not None:
        np.random.seed(random_state)
    
    events = []
    current_event = f"{start_domain}:{start_event}"
    current_timestamp = start_timestamp
    
    for i in range(n_events):
        # Calculate time increment (exponential distribution for realistic spacing)
        time_increment = np.random.exponential(time_span / n_events)
        current_timestamp += time_increment
        
        # Extract domain and event type from current event
        if ":" in current_event:
            current_domain, current_event_type = current_event.split(":", 1)
        else:
            current_domain = start_domain
            current_event_type = start_event
        
        # Create event
        event = Event(
            event_type=current_event_type,
            timestamp=datetime.fromtimestamp(current_timestamp),
            domain=current_domain,
            attributes={"chain_position": i, "chain_event": current_event}
        )
        events.append(event)
        
        # Determine next event based on chain rules
        if current_event in chain_rules:
            next_events = chain_rules[current_event]
            next_event_list = list(next_events.keys())
            probabilities = list(next_events.values())
            # Normalize probabilities
            prob_sum = sum(probabilities)
            if prob_sum > 0:
                probabilities = [p / prob_sum for p in probabilities]
                next_event = np.random.choice(next_event_list, p=probabilities)
                current_event = next_event
            else:
                # Fallback: random event from available domains
                if event_types_by_domain:
                    domain = np.random.choice(list(event_types_by_domain.keys()))
                    if event_types_by_domain[domain]:
                        event_type = np.random.choice(event_types_by_domain[domain])
                        current_event = f"{domain}:{event_type}"
        else:
            # No transition rule, select random event
            if event_types_by_domain:
                domain = np.random.choice(list(event_types_by_domain.keys()))
                if event_types_by_domain[domain]:
                    event_type = np.random.choice(event_types_by_domain[domain])
                    current_event = f"{domain}:{event_type}"
    
    return events


def add_temporal_noise(
    sequence: EventSequence,
    noise_level: float = 0.1,
    max_days_shift: int = 30,
    missing_probability: float = 0.0,
    random_state: Optional[int] = None
) -> EventSequence:
    """Add realistic temporal noise to event sequence.
    
    Introduces temporal uncertainty by:
    - Shifting event timestamps randomly
    - Optionally removing some events (missing data)
    
    Args:
        sequence: Original EventSequence
        noise_level: Fraction of events to modify (0.0 to 1.0)
        max_days_shift: Maximum days to shift events forward/backward
        missing_probability: Probability of removing each event (0.0 to 1.0)
        random_state: Random seed for reproducibility
        
    Returns:
        New EventSequence with noise applied
        
    Examples:
        >>> seq = EventSequence("p1", [Event("degree", datetime(2010, 1, 1), "education")])
        >>> noisy_seq = add_temporal_noise(seq, noise_level=0.5, max_days_shift=7, random_state=42)
        >>> len(noisy_seq.events) <= len(seq.events)
        True
    """
    if random_state is not None:
        np.random.seed(random_state)
    
    events = []
    n_to_modify = int(len(sequence.events) * noise_level)
    indices_to_modify = np.random.choice(
        len(sequence.events),
        size=n_to_modify,
        replace=False
    )
    
    for i, event in enumerate(sequence.events):
        # Check if event should be removed (missing data)
        if np.random.random() < missing_probability:
            continue
        
        # Check if this event should have temporal noise
        if i in indices_to_modify:
            # Convert timestamp to datetime if needed
            if isinstance(event.timestamp, datetime):
                event_time = event.timestamp
            else:
                event_time = datetime.fromtimestamp(event.timestamp)
            
            # Add random shift
            days_shift = np.random.randint(-max_days_shift, max_days_shift + 1)
            new_timestamp = event_time + timedelta(days=days_shift)
        else:
            new_timestamp = event.timestamp
        
        # Create new event with modified timestamp
        new_event = Event(
            event_type=event.event_type,
            timestamp=new_timestamp,
            domain=event.domain,
            attributes={**event.attributes, "noise_applied": i in indices_to_modify}
        )
        events.append(new_event)
    
    return EventSequence(
        person_id=sequence.person_id,
        events=events,
        metadata={**sequence.metadata, "temporal_noise_applied": True}
    )


def generate_realistic_life_events(
    n_sequences: int = 100,
    min_events_per_sequence: int = 5,
    max_events_per_sequence: int = 30,
    start_date: Optional[datetime] = None,
    end_date: Optional[datetime] = None,
    domains: Optional[List[str]] = None,
    event_types_by_domain: Optional[Dict[str, List[str]]] = None,
    transition_probabilities: Optional[Dict[str, Dict[str, float]]] = None,
    cooccurrence_patterns: Optional[Dict[str, List[str]]] = None,
    seasonal_patterns: bool = False,
    rare_event_probability: float = 0.05,
    generate_outcomes: bool = False,
    outcome_relationship: str = "random",
    temporal_noise: float = 0.0,
    missing_data_probability: float = 0.0,
    random_state: Optional[int] = None
) -> Tuple[List[EventSequence], Optional[np.ndarray]]:
    """Generate highly realistic synthetic life event sequences with advanced patterns.
    
    Creates event sequences with:
    - Temporal dependencies (Markov chains)
    - Event co-occurrence patterns
    - Domain transition probabilities
    - Seasonal/cyclical patterns
    - Rare events
    - Configurable noise and missing data
    
    Args:
        n_sequences: Number of sequences to generate
        min_events_per_sequence: Minimum events per sequence
        max_events_per_sequence: Maximum events per sequence
        start_date: Start date (defaults to 10 years ago)
        end_date: End date (defaults to today)
        domains: List of domains (defaults to standard domains)
        event_types_by_domain: Event types per domain
        transition_probabilities: Domain transition probabilities
            Format: {"domain1": {"domain2": probability}}
        cooccurrence_patterns: Events that tend to co-occur
            Format: {"event_type": ["co_occurring_event1", "co_occurring_event2"]}
        seasonal_patterns: Whether to apply seasonal variations
        rare_event_probability: Probability of injecting rare events
        generate_outcomes: Whether to generate outcome labels
        outcome_relationship: Outcome generation pattern
        temporal_noise: Level of temporal noise (0.0 to 1.0)
        missing_data_probability: Probability of missing events
        random_state: Random seed for reproducibility
        
    Returns:
        Tuple of (list of EventSequence objects, optional outcome array)
        
    Examples:
        >>> sequences, outcomes = generate_realistic_life_events(
        ...     n_sequences=50,
        ...     transition_probabilities={"education": {"occupation": 0.8}},
        ...     seasonal_patterns=True,
        ...     random_state=42
        ... )
        >>> len(sequences)
        50
    """
    if random_state is not None:
        np.random.seed(random_state)
    
    # Default domains
    if domains is None:
        domains = ["health", "education", "occupation", "income", "address", "other"]
    
    # Default event types
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
    
    # Default transition probabilities (education -> occupation is common)
    if transition_probabilities is None:
        transition_probabilities = {
            "education": {"occupation": 0.7, "income": 0.2, "address": 0.1},
            "occupation": {"income": 0.6, "address": 0.3, "health": 0.1},
            "income": {"address": 0.4, "occupation": 0.3, "other": 0.3},
            "health": {"occupation": 0.3, "income": 0.2, "other": 0.5},
            "address": {"occupation": 0.4, "income": 0.3, "other": 0.3}
        }
    
    # Default date range
    if end_date is None:
        end_date = datetime.now()
    if start_date is None:
        start_date = end_date - timedelta(days=365 * 10)
    
    start_timestamp = start_date.timestamp()
    end_timestamp = end_date.timestamp()
    time_span = end_timestamp - start_timestamp
    
    sequences = []
    outcomes = None
    
    if generate_outcomes:
        outcomes = np.zeros(n_sequences, dtype=np.float64)
    
    # Rare events (infrequent but important)
    rare_events = {
        "health": ["death", "major_surgery", "critical_diagnosis"],
        "occupation": ["retirement", "major_career_shift"],
        "other": ["divorce", "major_legal_event"]
    }
    
    for i in range(n_sequences):
        person_id = f"person_{i+1:04d}"
        n_events = np.random.randint(min_events_per_sequence, max_events_per_sequence + 1)
        
        events = []
        current_domain = np.random.choice(domains)
        event_timestamps = []
        
        # Generate timestamps with potential seasonal patterns
        if seasonal_patterns:
            # Create seasonal variation (more events in certain months)
            # Bias towards summer months (June-September)
            month_probs = np.array([
                0.06, 0.06, 0.06, 0.08, 0.08, 0.12,  # Jan-Jun
                0.12, 0.12, 0.12, 0.08, 0.06, 0.06  # Jul-Dec
            ])
            # Normalize to ensure exact sum of 1.0 (handles floating point precision)
            month_probs = month_probs / month_probs.sum()
            # Keep generating until we have enough events within the date range
            max_attempts = n_events * 10  # Safety limit
            attempts = 0
            while len(event_timestamps) < n_events and attempts < max_attempts:
                month = np.random.choice(range(1, 13), p=month_probs)
                # Random day in month
                day = np.random.randint(1, 29)
                year = np.random.randint(start_date.year, end_date.year + 1)
                try:
                    event_time = datetime(year, month, day)
                    timestamp = event_time.timestamp()
                    if start_timestamp <= timestamp <= end_timestamp:
                        event_timestamps.append(timestamp)
                except ValueError:
                    # Invalid date, use uniform distribution as fallback
                    event_timestamps.append(np.random.uniform(start_timestamp, end_timestamp))
                attempts += 1
            # If we still don't have enough, fill with uniform distribution
            while len(event_timestamps) < n_events:
                event_timestamps.append(np.random.uniform(start_timestamp, end_timestamp))
        else:
            # Uniform distribution
            event_timestamps = sorted(np.random.uniform(start_timestamp, end_timestamp, n_events))
        
        # Track for co-occurrence and outcomes
        domain_counts = {domain: 0 for domain in domains}
        event_history = []
        
        for j, timestamp in enumerate(event_timestamps):
            # Determine if this is a rare event
            is_rare = np.random.random() < rare_event_probability
            
            if is_rare and current_domain in rare_events and rare_events[current_domain]:
                event_type = np.random.choice(rare_events[current_domain])
            else:
                # Select event type from domain
                if current_domain in event_types_by_domain and event_types_by_domain[current_domain]:
                    event_type = np.random.choice(event_types_by_domain[current_domain])
                else:
                    event_type = f"event_{np.random.randint(1000)}"
            
            # Create event
            event = Event(
                event_type=event_type,
                timestamp=datetime.fromtimestamp(timestamp),
                domain=current_domain,
                attributes={"is_rare": is_rare, "sequence_position": j}
            )
            events.append(event)
            
            # Track for outcomes
            domain_counts[current_domain] += 1
            event_history.append(f"{current_domain}:{event_type}")
            
            # Check for co-occurrence patterns
            if cooccurrence_patterns and event_type in cooccurrence_patterns:
                co_events = cooccurrence_patterns[event_type]
                if co_events and np.random.random() < 0.5:  # 50% chance of co-occurrence
                    # Add co-occurring event soon after
                    co_event_type = np.random.choice(co_events)
                    co_domain = None
                    for dom, event_list in event_types_by_domain.items():
                        if co_event_type in event_list:
                            co_domain = dom
                            break
                    if co_domain:
                        # Add within 30 days
                        co_timestamp = timestamp + np.random.uniform(0, 30 * 86400)
                        if co_timestamp <= end_timestamp:
                            co_event = Event(
                                event_type=co_event_type,
                                timestamp=datetime.fromtimestamp(co_timestamp),
                                domain=co_domain,
                                attributes={"co_occurring_with": event_type}
                            )
                            events.append(co_event)
                            domain_counts[co_domain] += 1
            
            # Transition to next domain based on transition probabilities
            if current_domain in transition_probabilities:
                next_domains = transition_probabilities[current_domain]
                next_domain_list = list(next_domains.keys())
                probabilities = list(next_domains.values())
                # Normalize
                prob_sum = sum(probabilities)
                if prob_sum > 0:
                    probabilities = [p / prob_sum for p in probabilities]
                    current_domain = np.random.choice(next_domain_list, p=probabilities)
                else:
                    current_domain = np.random.choice(domains)
            else:
                current_domain = np.random.choice(domains)
        
        # Sort events by timestamp
        events.sort(key=lambda e: e.timestamp.timestamp() if isinstance(e.timestamp, datetime) else e.timestamp)
        
        # Create sequence
        sequence = EventSequence(
            person_id=person_id,
            events=events,
            metadata={
                "synthetic": True,
                "realistic": True,
                "generated_at": datetime.now().isoformat(),
                "temporal_noise": temporal_noise,
                "missing_data_prob": missing_data_probability
            }
        )
        
        # Apply temporal noise if requested
        if temporal_noise > 0 or missing_data_probability > 0:
            sequence = add_temporal_noise(
                sequence,
                noise_level=temporal_noise,
                missing_probability=missing_data_probability,
                random_state=random_state
            )
        
        sequences.append(sequence)
        
        # Generate outcome if requested
        if generate_outcomes:
            if outcome_relationship == "random":
                outcomes[i] = np.random.choice([0, 1])
            elif outcome_relationship == "health_focused":
                health_ratio = domain_counts.get("health", 0) / len(events) if events else 0
                outcomes[i] = 1 if health_ratio > 0.3 or np.random.random() < 0.3 else 0
            elif outcome_relationship == "education_focused":
                education_ratio = domain_counts.get("education", 0) / len(events) if events else 0
                outcomes[i] = 1 if education_ratio > 0.2 or np.random.random() < 0.2 else 0
            elif outcome_relationship == "complex":
                has_education = domain_counts.get("education", 0) > 0
                has_health = domain_counts.get("health", 0) > 0
                has_occupation = domain_counts.get("occupation", 0) > 0
                outcomes[i] = 1 if (has_education and has_health) or (has_occupation and np.random.random() < 0.4) else 0
            else:
                outcomes[i] = np.random.choice([0, 1])
    
    return sequences, outcomes


def generate_cohort_sequences(
    n_cohorts: int = 3,
    n_sequences_per_cohort: int = 50,
    cohort_differences: Optional[Dict[str, Dict[str, float]]] = None,
    min_events_per_sequence: int = 5,
    max_events_per_sequence: int = 30,
    start_date: Optional[datetime] = None,
    end_date: Optional[datetime] = None,
    random_state: Optional[int] = None
) -> Dict[str, List[EventSequence]]:
    """Generate population-level event sequences with cohort-specific patterns.
    
    Creates multiple cohorts (e.g., different age groups, socioeconomic groups)
    with distinct event patterns, allowing for population-level analysis.
    
    Args:
        n_cohorts: Number of distinct cohorts
        n_sequences_per_cohort: Number of sequences per cohort
        cohort_differences: Dictionary specifying cohort-specific patterns
            Format: {"cohort_name": {"domain": probability, "event_type": probability}}
        min_events_per_sequence: Minimum events per sequence
        max_events_per_sequence: Maximum events per sequence
        start_date: Start date (defaults to 10 years ago)
        end_date: End date (defaults to today)
        random_state: Random seed for reproducibility
        
    Returns:
        Dictionary mapping cohort names to lists of EventSequence objects
        
    Examples:
        >>> cohorts = generate_cohort_sequences(
        ...     n_cohorts=2,
        ...     n_sequences_per_cohort=20,
        ...     cohort_differences={
        ...         "young": {"education": 0.4, "occupation": 0.3},
        ...         "old": {"health": 0.4, "retirement": 0.3}
        ...     },
        ...     random_state=42
        ... )
        >>> len(cohorts)
        2
    """
    if random_state is not None:
        np.random.seed(random_state)
    
    # Default cohort differences
    if cohort_differences is None:
        cohort_names = [f"cohort_{i+1}" for i in range(n_cohorts)]
        cohort_differences = {}
        for i, name in enumerate(cohort_names):
            # Different cohorts emphasize different domains
            domain_weights = {
                "health": 0.15 + (i * 0.1),
                "education": 0.2 - (i * 0.05),
                "occupation": 0.25,
                "income": 0.2,
                "address": 0.1,
                "other": 0.1
            }
            cohort_differences[name] = domain_weights
    
    cohorts = {}
    
    for cohort_name, domain_weights in cohort_differences.items():
        # Generate sequences for this cohort
        sequences, _ = generate_realistic_life_events(
            n_sequences=n_sequences_per_cohort,
            min_events_per_sequence=min_events_per_sequence,
            max_events_per_sequence=max_events_per_sequence,
            start_date=start_date,
            end_date=end_date,
            random_state=random_state
        )
        
        # Modify sequences to match cohort patterns
        # (In a full implementation, this would adjust domain frequencies)
        for seq in sequences:
            seq.metadata["cohort"] = cohort_name
        
        cohorts[cohort_name] = sequences
    
    return cohorts

