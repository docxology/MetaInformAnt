"""Utility functions for life events analysis.

This module provides utility functions for generating, manipulating,
and analyzing life event sequences.
"""

from __future__ import annotations

import random
from typing import List, Dict, Any, Optional
from datetime import datetime, timedelta

from metainformant.core import logging
from metainformant.life_events.events import Event, EventSequence

logger = logging.get_logger(__name__)


def add_temporal_noise(sequence: EventSequence, noise_level: float = 0.1,
                      rng: random.Random | None = None) -> EventSequence:
    """Add temporal noise to event timestamps.

    Args:
        sequence: Input event sequence
        noise_level: Fraction of total duration to use as noise scale
        rng: Random number generator (optional)

    Returns:
        New EventSequence with added noise
    """
    if rng is None:
        rng = random.Random()

    if not sequence.events:
        return sequence

    # Calculate sequence duration
    timestamps = [event.timestamp for event in sequence.events]
    duration = max(timestamps) - min(timestamps) if len(timestamps) > 1 else 1.0
    noise_scale = duration * noise_level

    # Create new events with noise
    new_events = []
    for event in sequence.events:
        # Add Gaussian noise to timestamp
        noise = rng.gauss(0, noise_scale)
        new_timestamp = event.timestamp + noise

        # Create new event with modified timestamp
        new_event = Event(
            timestamp=new_timestamp,
            event_type=event.event_type,
            description=event.description,
            metadata=event.metadata.copy(),
            confidence=event.confidence
        )
        new_events.append(new_event)

    # Sort events by timestamp
    new_events.sort(key=lambda e: e.timestamp)

    # Create new sequence
    new_sequence = EventSequence(
        person_id=sequence.person_id,
        events=new_events,
        metadata=sequence.metadata.copy()
    )

    return new_sequence


def generate_realistic_life_events(n_events: int = 50, start_age: float = 0.0,
                                 end_age: float = 80.0, seed: int | None = None) -> EventSequence:
    """Generate a realistic sequence of life events.

    Args:
        n_events: Number of events to generate
        start_age: Starting age in years
        end_age: Ending age in years
        seed: Random seed for reproducibility

    Returns:
        EventSequence with generated events
    """
    rng = random.Random(seed)

    # Define event types and their typical age ranges
    event_types = {
        'birth': (0, 0),
        'education_start': (5, 7),
        'school_graduation': (17, 19),
        'university_start': (18, 22),
        'university_graduation': (21, 26),
        'first_job': (20, 30),
        'marriage': (25, 35),
        'first_child': (25, 40),
        'career_change': (25, 60),
        'retirement': (60, 70),
        'health_diagnosis': (30, 80),
        'relocation': (20, 70)
    }

    events = []

    # Generate events
    for i in range(n_events):
        event_type = rng.choice(list(event_types.keys()))
        age_range = event_types[event_type]
        age = rng.uniform(age_range[0], age_range[1])

        # Skip if outside sequence range
        if not (start_age <= age <= end_age):
            continue

        event = Event(
            timestamp=age,
            event_type=event_type,
            description=f"Generated {event_type} event",
            metadata={'generated': True, 'event_index': i}
        )
        events.append(event)

    # Sort by timestamp
    events.sort(key=lambda e: e.timestamp)

    return EventSequence(
        person_id="generated_sequence",
        events=events,
        metadata={'generated': True, 'n_events': len(events)}
    )


def get_event_statistics(sequences: List[EventSequence]) -> Dict[str, Any]:
    """Calculate statistics across multiple event sequences.

    Args:
        sequences: List of event sequences

    Returns:
        Dictionary with sequence statistics
    """
    if not sequences:
        return {}

    stats = {
        'n_sequences': len(sequences),
        'total_events': sum(len(seq.events) for seq in sequences),
        'events_per_sequence': [],
        'event_types': set(),
        'duration_stats': []
    }

    for seq in sequences:
        stats['events_per_sequence'].append(len(seq.events))
        stats['event_types'].update(event.event_type for event in seq.events)

        if seq.events:
            timestamps = [event.timestamp for event in seq.events]
            duration = max(timestamps) - min(timestamps)
            stats['duration_stats'].append(duration)

    # Calculate aggregates
    if stats['events_per_sequence']:
        stats['mean_events_per_sequence'] = sum(stats['events_per_sequence']) / len(stats['events_per_sequence'])

    if stats['duration_stats']:
        stats['mean_duration'] = sum(stats['duration_stats']) / len(stats['duration_stats'])

    stats['unique_event_types'] = len(stats['event_types'])
    stats['event_types'] = list(stats['event_types'])

    return stats


def load_sequences_from_json(json_path: str | Path) -> List[EventSequence]:
    """Load EventSequence objects from a JSON file.

    Args:
        json_path: Path to JSON file containing sequence data

    Returns:
        List of EventSequence objects
    """
    import json
    from pathlib import Path

    json_path = Path(json_path)
    if not json_path.exists():
        raise FileNotFoundError(f"JSON file not found: {json_path}")

    with open(json_path, 'r') as f:
        data = json.load(f)

    sequences = []

    # Handle different JSON formats
    if isinstance(data, dict):
        # Database format: {"sequences": [...], "metadata": {...}}
        if "sequences" in data:
            sequence_data = data["sequences"]
        else:
            # Assume single sequence format
            sequence_data = [data]
    elif isinstance(data, list):
        # List format: [sequence1, sequence2, ...]
        sequence_data = data
    else:
        raise ValueError("Unsupported JSON format")

    for seq_data in sequence_data:
        person_id = seq_data.get("person_id", seq_data.get("id", "unknown"))
        events_data = seq_data.get("events", [])

        events = []
        for event_data in events_data:
            event = Event(
                event_type=event_data["event_type"],
                timestamp=event_data["timestamp"],
                domain=event_data.get("domain"),
                metadata=event_data.get("metadata", {})
            )
            events.append(event)

        sequence = EventSequence(person_id=person_id, events=events)
        sequences.append(sequence)

    return sequences


def convert_sequences_to_tokens(sequences: List[EventSequence]) -> List[List[str]]:
    """Convert EventSequence objects to token sequences for ML processing.

    Args:
        sequences: List of EventSequence objects

    Returns:
        List of token lists (event types)
    """
    token_sequences = []
    for seq in sequences:
        tokens = [event.event_type for event in seq.events]
        token_sequences.append(tokens)
    return token_sequences


def validate_sequence(sequence: EventSequence) -> Dict[str, Any]:
    """Validate an EventSequence for consistency and completeness.

    Args:
        sequence: EventSequence to validate

    Returns:
        Dictionary with validation results
    """
    validation = {
        "is_valid": True,
        "issues": [],
        "warnings": [],
        "stats": {
            "n_events": len(sequence.events),
            "event_types": set(),
            "domains": set(),
            "date_range": None
        }
    }

    if not sequence.events:
        validation["issues"].append("Sequence has no events")
        validation["is_valid"] = False
        return validation

    # Check events
    timestamps = []
    for i, event in enumerate(sequence.events):
        validation["stats"]["event_types"].add(event.event_type)
        if hasattr(event, 'domain') and event.domain:
            validation["stats"]["domains"].add(event.domain)

        if not event.event_type:
            validation["issues"].append(f"Event {i} has no event_type")

        if not hasattr(event, 'timestamp') or event.timestamp is None:
            validation["issues"].append(f"Event {i} has no timestamp")
        else:
            timestamps.append(event.timestamp)

    # Check temporal ordering
    if len(timestamps) > 1:
        sorted_timestamps = sorted(timestamps)
        if timestamps != sorted_timestamps:
            validation["warnings"].append("Events are not in chronological order")

        # Calculate date range
        from datetime import datetime
        if timestamps:
            min_date = min(timestamps)
            max_date = max(timestamps)
            validation["stats"]["date_range"] = {
                "start": min_date.isoformat() if hasattr(min_date, 'isoformat') else str(min_date),
                "end": max_date.isoformat() if hasattr(max_date, 'isoformat') else str(max_date),
                "span_days": (max_date - min_date).days if hasattr(max_date, '__sub__') else None
            }

    validation["stats"]["event_types"] = list(validation["stats"]["event_types"])
    validation["stats"]["domains"] = list(validation["stats"]["domains"])
    validation["is_valid"] = len(validation["issues"]) == 0

    return validation


def generate_cohort_sequences(
    n_sequences: int,
    avg_events_per_sequence: int = 10,
    event_types: Optional[List[str]] = None,
    domains: Optional[List[str]] = None,
    start_year: int = 1980,
    end_year: int = 2020,
    **kwargs: Any
) -> List[EventSequence]:
    """Generate synthetic cohort sequences for testing and simulation.

    Args:
        n_sequences: Number of sequences to generate
        avg_events_per_sequence: Average events per sequence
        event_types: List of possible event types
        domains: List of possible domains
        start_year: Start year for event timestamps
        end_year: End year for event timestamps
        **kwargs: Additional generation parameters

    Returns:
        List of generated EventSequence objects
    """
    import random
    from datetime import datetime, timedelta

    # Default event types and domains
    if event_types is None:
        event_types = [
            "started_school", "graduated_high_school", "started_college",
            "graduated_college", "got_first_job", "job_change", "promotion",
            "started_relationship", "marriage", "divorce", "had_child",
            "moved_house", "diagnosis", "recovered", "retirement"
        ]

    if domains is None:
        domains = ["education", "career", "family", "health", "housing"]

    sequences = []

    for i in range(n_sequences):
        person_id = f"person_{i:04d}"

        # Generate random number of events
        n_events = max(1, int(random.gauss(avg_events_per_sequence, 2)))

        events = []

        # Generate events with temporal progression
        current_year = start_year
        for j in range(n_events):
            # Select random event type
            event_type = random.choice(event_types)

            # Select domain based on event type
            if "school" in event_type or "college" in event_type:
                domain = "education"
            elif "job" in event_type or "promotion" in event_type or "retirement" in event_type:
                domain = "career"
            elif "relationship" in event_type or "marriage" in event_type or "divorce" in event_type or "child" in event_type:
                domain = "family"
            elif "diagnosis" in event_type or "recovered" in event_type:
                domain = "health"
            elif "moved" in event_type:
                domain = "housing"
            else:
                domain = random.choice(domains)

            # Generate timestamp with some progression
            year_offset = int((end_year - start_year) * (j / max(1, n_events - 1)))
            event_year = start_year + year_offset
            event_month = random.randint(1, 12)
            event_day = random.randint(1, 28)  # Avoid Feb 29 issues

            timestamp = datetime(event_year, event_month, event_day)

            # Create event
            event = Event(
                event_type=event_type,
                timestamp=timestamp,
                domain=domain,
                metadata={"generated": True, "sequence_index": j}
            )
            events.append(event)

        # Sort events by timestamp
        events.sort(key=lambda x: x.timestamp)

        # Create sequence
        sequence = EventSequence(person_id=person_id, events=events)
        sequences.append(sequence)

    return sequences


def sequence_embeddings(
    sequences: List[EventSequence],
    method: str = "average",
    embedding_dict: Optional[Dict[str, np.ndarray]] = None,
    **kwargs: Any
) -> Dict[str, np.ndarray]:
    """Create sequence-level embeddings from event embeddings.

    Args:
        sequences: List of EventSequence objects
        method: Aggregation method ('average', 'sum', 'max', 'attention')
        embedding_dict: Dictionary mapping event types to embeddings
        **kwargs: Additional parameters

    Returns:
        Dictionary mapping sequence IDs to embeddings
    """
    if embedding_dict is None:
        # Use simple one-hot style embeddings if none provided
        all_event_types = set()
        for seq in sequences:
            for event in seq.events:
                all_event_types.add(event.event_type)

        embedding_dim = kwargs.get('embedding_dim', 50)
        embedding_dict = {}
        for i, event_type in enumerate(all_event_types):
            embedding = np.zeros(embedding_dim)
            embedding[i % embedding_dim] = 1.0
            embedding_dict[event_type] = embedding

    sequence_embeddings_result = {}

    for seq in sequences:
        if not seq.events:
            # Empty sequence - use zero embedding
            embedding_dim = len(list(embedding_dict.values())[0]) if embedding_dict else 50
            sequence_embeddings_result[seq.person_id] = np.zeros(embedding_dim)
            continue

        # Get embeddings for all events in sequence
        event_embeddings = []
        for event in seq.events:
            if event.event_type in embedding_dict:
                event_embeddings.append(embedding_dict[event.event_type])
            else:
                # Use zero embedding for unknown event types
                embedding_dim = len(list(embedding_dict.values())[0]) if embedding_dict else 50
                event_embeddings.append(np.zeros(embedding_dim))

        if not event_embeddings:
            embedding_dim = len(list(embedding_dict.values())[0]) if embedding_dict else 50
            sequence_embeddings_result[seq.person_id] = np.zeros(embedding_dim)
            continue

        event_embeddings = np.array(event_embeddings)

        # Aggregate embeddings
        if method == "average":
            seq_embedding = np.mean(event_embeddings, axis=0)
        elif method == "sum":
            seq_embedding = np.sum(event_embeddings, axis=0)
        elif method == "max":
            seq_embedding = np.max(event_embeddings, axis=0)
        elif method == "attention":
            # Simple attention: use first event as query
            query = event_embeddings[0]
            scores = np.dot(event_embeddings, query)
            weights = np.exp(scores) / np.sum(np.exp(scores))
            seq_embedding = np.sum(event_embeddings * weights[:, np.newaxis], axis=0)
        else:
            raise ValueError(f"Unknown aggregation method: {method}")

        sequence_embeddings_result[seq.person_id] = seq_embedding

    return sequence_embeddings_result


def generate_event_chain(
    start_event: str,
    n_events: int,
    transition_matrix: Optional[Dict[str, Dict[str, float]]] = None,
    event_types: Optional[List[str]] = None,
    **kwargs: Any
) -> List[str]:
    """Generate a causally linked sequence of events using a Markov chain.

    Args:
        start_event: Initial event type
        n_events: Number of events to generate
        transition_matrix: Dictionary of transition probabilities
        event_types: List of possible event types
        **kwargs: Additional parameters

    Returns:
        List of event types in sequence
    """
    import random

    if event_types is None:
        event_types = [
            "started_school", "graduated_high_school", "started_college",
            "graduated_college", "got_first_job", "job_change", "promotion",
            "started_relationship", "marriage", "had_child", "diagnosis",
            "recovered", "moved_house", "retirement"
        ]

    # Create default transition matrix if not provided
    if transition_matrix is None:
        transition_matrix = {}
        for event in event_types:
            transitions = {}
            # Create plausible transitions
            if "school" in event:
                transitions["graduated_high_school"] = 0.3
                transitions["started_college"] = 0.4
            elif "college" in event:
                transitions["graduated_college"] = 0.3
                transitions["got_first_job"] = 0.3
            elif "job" in event or "promotion" in event:
                transitions["job_change"] = 0.2
                transitions["promotion"] = 0.2
                transitions["retirement"] = 0.1
            elif "relationship" in event:
                transitions["marriage"] = 0.3
                transitions["had_child"] = 0.2
            elif "marriage" in event:
                transitions["had_child"] = 0.4
                transitions["divorce"] = 0.1
            elif "diagnosis" in event:
                transitions["recovered"] = 0.4
            elif "moved" in event:
                transitions["job_change"] = 0.2

            # Add some random transitions
            for other_event in random.sample(event_types, min(3, len(event_types))):
                if other_event not in transitions:
                    transitions[other_event] = 0.05

            # Normalize probabilities
            total_prob = sum(transitions.values())
            if total_prob > 0:
                for evt in transitions:
                    transitions[evt] /= total_prob

            transition_matrix[event] = transitions

    # Generate sequence
    sequence = [start_event]
    current_event = start_event

    for _ in range(n_events - 1):
        if current_event not in transition_matrix:
            # Random choice if no transitions defined
            current_event = random.choice(event_types)
        else:
            transitions = transition_matrix[current_event]
            if not transitions:
                current_event = random.choice(event_types)
            else:
                # Sample from transition probabilities
                events, probs = zip(*transitions.items())
                current_event = random.choices(events, weights=probs, k=1)[0]

        sequence.append(current_event)

    return sequence


def generate_synthetic_life_events(
    n_sequences: int = 100,
    min_events_per_sequence: int = 3,
    max_events_per_sequence: int = 20,
    generate_outcomes: bool = False,
    outcome_relationship: str = "random",
    event_types: Optional[List[str]] = None,
    domains: Optional[List[str]] = None,
    temporal_patterns: str = "realistic",
    random_state: Optional[int] = None,
    **kwargs: Any
) -> Tuple[List[EventSequence], Optional[List[Any]]]:
    """Generate synthetic life event sequences with realistic patterns.

    Args:
        n_sequences: Number of sequences to generate
        min_events_per_sequence: Minimum events per sequence
        max_events_per_sequence: Maximum events per sequence
        generate_outcomes: Whether to generate outcome labels
        outcome_relationship: How outcomes relate to sequences ("random", "event_based", "temporal")
        event_types: List of possible event types
        domains: List of possible domains
        temporal_patterns: Type of temporal patterns ("realistic", "random", "clustered")
        random_state: Random seed for reproducibility
        **kwargs: Additional parameters

    Returns:
        Tuple of (sequences, outcomes) where outcomes is None if not generated
    """
    import random
    from datetime import datetime, timedelta

    if random_state is not None:
        random.seed(random_state)

    # Default event types and domains
    if event_types is None:
        event_types = [
            "started_school", "graduated_high_school", "started_college",
            "graduated_college", "got_first_job", "job_change", "promotion",
            "started_relationship", "marriage", "divorce", "had_child",
            "diagnosis", "recovered", "moved_house", "parent_illness",
            "financial_difficulty", "inheritance", "retirement"
        ]

    if domains is None:
        domains = ["education", "career", "family", "health", "finance", "housing"]

    # Define realistic event chains and patterns
    education_chain = ["started_school", "graduated_high_school", "started_college", "graduated_college"]
    career_chain = ["got_first_job", "job_change", "promotion", "retirement"]
    family_chain = ["started_relationship", "marriage", "had_child"]
    health_chain = ["diagnosis", "recovered"]

    sequences = []
    outcomes = [] if generate_outcomes else None

    for seq_id in range(n_sequences):
        person_id = f"person_{seq_id:04d}"

        # Determine number of events
        n_events = random.randint(min_events_per_sequence, max_events_per_sequence)

        events = []
        used_events = set()

        # Generate base timeline
        start_year = random.randint(1950, 2000)
        current_year = start_year

        # Add education events (most people have these)
        if random.random() < 0.8:  # 80% have education
            for event_type in education_chain:
                if len(events) < n_events and event_type not in used_events:
                    year_offset = {"started_school": 0, "graduated_high_school": 4,
                                 "started_college": 4, "graduated_college": 8}.get(event_type, 0)
                    event_year = start_year + year_offset
                    event_month = random.randint(1, 12)
                    event_day = random.randint(1, 28)

                    event = Event(
                        event_type=event_type,
                        timestamp=datetime(event_year, event_month, event_day),
                        domain="education"
                    )
                    events.append(event)
                    used_events.add(event_type)
                    current_year = max(current_year, event_year)

        # Add career events
        if random.random() < 0.7:  # 70% have career
            career_start = current_year + random.randint(0, 2)
            for i, event_type in enumerate(career_chain):
                if len(events) < n_events and event_type not in used_events:
                    event_year = career_start + i * random.randint(2, 5)
                    event_month = random.randint(1, 12)
                    event_day = random.randint(1, 28)

                    event = Event(
                        event_type=event_type,
                        timestamp=datetime(event_year, event_month, event_day),
                        domain="career"
                    )
                    events.append(event)
                    used_events.add(event_type)
                    current_year = max(current_year, event_year)

        # Add family events
        if random.random() < 0.6:  # 60% have family events
            family_start = current_year - random.randint(0, 5)
            for event_type in family_chain:
                if len(events) < n_events and event_type not in used_events:
                    event_year = family_start + random.randint(0, 3)
                    event_month = random.randint(1, 12)
                    event_day = random.randint(1, 28)

                    event = Event(
                        event_type=event_type,
                        timestamp=datetime(event_year, event_month, event_day),
                        domain="family"
                    )
                    events.append(event)
                    used_events.add(event_type)

        # Add health events (less common)
        if random.random() < 0.3 and len(events) < n_events:  # 30% have health issues
            for event_type in health_chain:
                if len(events) < n_events and event_type not in used_events:
                    event_year = current_year + random.randint(0, 10)
                    event_month = random.randint(1, 12)
                    event_day = random.randint(1, 28)

                    event = Event(
                        event_type=event_type,
                        timestamp=datetime(event_year, event_month, event_day),
                        domain="health"
                    )
                    events.append(event)
                    used_events.add(event_type)

        # Fill remaining events with random events
        while len(events) < n_events:
            event_type = random.choice(event_types)
            if event_type not in used_events:
                event_year = current_year + random.randint(0, 5)
                event_month = random.randint(1, 12)
                event_day = random.randint(1, 28)

                # Assign domain based on event type
                domain = "education" if "school" in event_type or "college" in event_type else \
                        "career" if "job" in event_type or "promotion" in event_type or "retirement" in event_type else \
                        "family" if "relationship" in event_type or "marriage" in event_type or "child" in event_type else \
                        "health" if "diagnosis" in event_type or "recovered" in event_type else \
                        "finance" if "financial" in event_type or "inheritance" in event_type else \
                        random.choice(domains)

                event = Event(
                    event_type=event_type,
                    timestamp=datetime(event_year, event_month, event_day),
                    domain=domain
                )
                events.append(event)
                used_events.add(event_type)

        # Sort events by timestamp
        events.sort(key=lambda x: x.timestamp)

        # Create sequence
        sequence = EventSequence(person_id=person_id, events=events)
        sequences.append(sequence)

        # Generate outcome if requested
        if generate_outcomes:
            if outcome_relationship == "random":
                outcome = random.choice([0, 1])
            elif outcome_relationship == "event_based":
                # Outcome based on presence of certain events
                has_negative_events = any(e.event_type in ["diagnosis", "divorce", "financial_difficulty"]
                                        for e in events)
                outcome = 1 if has_negative_events else 0
            elif outcome_relationship == "temporal":
                # Outcome based on event density/timing
                if len(events) > 12:  # Many events might indicate complex life
                    outcome = 1
                else:
                    outcome = 0
            else:
                outcome = random.choice([0, 1])

            outcomes.append(outcome)

    return sequences, outcomes
