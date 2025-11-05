"""Life course phenotype extraction from event sequences.

This module provides functions to extract phenotypes from event sequences,
enabling temporal phenotype aggregation and integration with life_events module.

Requires the life_events module to be available. Functions will raise ImportError
if life_events module is not installed.
"""

from __future__ import annotations

from datetime import datetime
from typing import Any, Dict, List, Optional

from ..core.errors import ValidationError
from ..core.logging import get_logger

try:
    from ..life_events import Event, EventSequence
    LIFE_EVENTS_AVAILABLE = True
except ImportError:
    LIFE_EVENTS_AVAILABLE = False
    Event = None  # type: ignore
    EventSequence = None  # type: ignore

logger = get_logger(__name__)


def extract_phenotypes_from_events(
    sequence: EventSequence,
    phenotype_categories: Optional[Dict[str, List[str]]] = None
) -> Dict[str, Any]:
    """Extract phenotype traits from event sequence.
    
    Maps events to phenotypic categories and aggregates temporal information.
    
    Args:
        sequence: EventSequence to extract phenotypes from (single sequence, not list)
        phenotype_categories: Optional mapping of phenotype categories to event types
        
    Returns:
        Dictionary with extracted phenotypes including:
        - person_id: Person identifier
        - total_events: Total number of events
        - domains: List of domains present
        - event_types: List of event types
        - domain_counts: Count of events per domain
        - health_events, education_events, occupation_events: Domain-specific counts
        - health_conditions, education_achievements, occupation_changes: Domain-specific event lists
        - first_event_time, last_event_time: Temporal boundaries
        - event_span_years: Time span in years
        
    Raises:
        ImportError: If life_events module is not available
        ValidationError: If sequence is invalid or empty
        
    Examples:
        >>> from metainformant.life_events import EventSequence, Event
        >>> from datetime import datetime
        >>> from metainformant.phenotype import extract_phenotypes_from_events
        >>> events = [
        ...     Event("diabetes", datetime(2020, 1, 1), "health"),
        ...     Event("bachelors", datetime(2010, 6, 1), "education"),
        ... ]
        >>> seq = EventSequence(person_id="p1", events=events)
        >>> phenotypes = extract_phenotypes_from_events(seq)
        >>> "health_events" in phenotypes
        True
        >>> phenotypes["total_events"]
        2
    """
    if not LIFE_EVENTS_AVAILABLE:
        raise ImportError(
            "life_events module not available. Install metainformant with life_events support."
        )
    
    if sequence is None:
        raise ValidationError("sequence cannot be None")
    
    logger.debug(f"Extracting phenotypes from sequence for person {sequence.person_id}")
    
    if not sequence.events:
        logger.warning(f"Empty event sequence for person {sequence.person_id}")
        return {
            "person_id": sequence.person_id,
            "total_events": 0,
            "domains": [],
            "event_types": [],
            "domain_counts": {},
        }
    
    phenotypes = {
        "person_id": sequence.person_id,
        "total_events": len(sequence.events),
        "domains": sequence.get_domains(),
        "event_types": sequence.get_event_types(),
    }
    
    # Count events by domain
    domain_counts = {}
    for event in sequence.events:
        domain_counts[event.domain] = domain_counts.get(event.domain, 0) + 1
    phenotypes["domain_counts"] = domain_counts
    
    logger.debug(f"Found {len(phenotypes['domains'])} domains in sequence")
    
    # Extract health-related phenotypes
    health_events = [e for e in sequence.events if e.domain == "health"]
    if health_events:
        phenotypes["health_events"] = len(health_events)
        phenotypes["health_conditions"] = [e.event_type for e in health_events]
    
    # Extract education-related phenotypes
    education_events = [e for e in sequence.events if e.domain == "education"]
    if education_events:
        phenotypes["education_events"] = len(education_events)
        phenotypes["education_achievements"] = [e.event_type for e in education_events]
    
    # Extract occupation-related phenotypes
    occupation_events = [e for e in sequence.events if e.domain == "occupation"]
    if occupation_events:
        phenotypes["occupation_events"] = len(occupation_events)
        phenotypes["occupation_changes"] = [e.event_type for e in occupation_events]
    
    # Temporal aggregation
    if sequence.events:
        def get_timestamp(event):
            if isinstance(event.timestamp, datetime):
                return event.timestamp.timestamp()
            return float(event.timestamp)
        
        timestamps = [get_timestamp(e) for e in sequence.events]
        phenotypes["first_event_time"] = min(timestamps)
        phenotypes["last_event_time"] = max(timestamps)
        phenotypes["event_span_years"] = (max(timestamps) - min(timestamps)) / (365.25 * 24 * 3600)
    
    # Custom phenotype categories
    if phenotype_categories:
        for category, event_types in phenotype_categories.items():
            if not isinstance(event_types, list):
                raise ValidationError(f"phenotype_categories[{category}] must be a list, got {type(event_types).__name__}")
            matching_events = [
                e for e in sequence.events
                if e.event_type in event_types
            ]
            phenotypes[f"{category}_count"] = len(matching_events)
    
    logger.info(f"Extracted {len(phenotypes)} phenotype fields from {len(sequence.events)} events")
    return phenotypes


def aggregate_temporal_phenotypes(
    sequences: List[EventSequence],
    time_window_years: float = 5.0
) -> Dict[str, Any]:
    """Aggregate phenotypes across time windows.
    
    Groups events from multiple sequences into temporal windows and computes
    aggregate statistics per window.
    
    Args:
        sequences: List of EventSequence objects (not phenotypes dicts)
        time_window_years: Size of time windows in years (must be positive)
        
    Returns:
        Dictionary with temporal phenotype aggregations:
        - time_windows: List of window statistics
        - aggregates: Overall statistics (total_events, total_people, time_span_years)
        
    Raises:
        ImportError: If life_events module is not available
        ValidationError: If sequences is invalid or time_window_years is invalid
        
    Examples:
        >>> from metainformant.life_events import EventSequence, Event
        >>> from datetime import datetime
        >>> from metainformant.phenotype import aggregate_temporal_phenotypes
        >>> seq1 = EventSequence("p1", [Event("degree", datetime(2010, 1, 1), "education")])
        >>> seq2 = EventSequence("p2", [Event("job", datetime(2015, 1, 1), "occupation")])
        >>> result = aggregate_temporal_phenotypes([seq1, seq2], time_window_years=5.0)
        >>> "time_windows" in result
        True
        >>> result["aggregates"]["total_people"]
        2
    """
    if not LIFE_EVENTS_AVAILABLE:
        raise ImportError(
            "life_events module not available. Install metainformant with life_events support."
        )
    
    if not isinstance(sequences, list):
        raise ValidationError(f"sequences must be a list, got {type(sequences).__name__}")
    
    if time_window_years <= 0:
        raise ValidationError(f"time_window_years must be positive, got {time_window_years}")
    
    logger.debug(f"Aggregating phenotypes from {len(sequences)} sequences with {time_window_years} year windows")
    
    # Collect all events with timestamps
    all_events = []
    for seq in sequences:
        for event in seq.events:
            def get_timestamp(e):
                if isinstance(e.timestamp, datetime):
                    return e.timestamp.timestamp()
                return float(e.timestamp)
            
            all_events.append({
                "event": event,
                "timestamp": get_timestamp(event),
                "person_id": seq.person_id
            })
    
    if not all_events:
        logger.warning("No events found in any sequences")
        return {"time_windows": [], "aggregates": {"total_events": 0, "total_people": 0, "time_span_years": 0.0}}
    
    # Define time windows
    min_time = min(e["timestamp"] for e in all_events)
    max_time = max(e["timestamp"] for e in all_events)
    
    window_size_seconds = time_window_years * 365.25 * 24 * 3600
    windows = []
    
    current_start = min_time
    while current_start < max_time:
        current_end = current_start + window_size_seconds
        windows.append({
            "start": current_start,
            "end": current_end,
            "events": []
        })
        current_start = current_end
    
    # Assign events to windows
    for event_data in all_events:
        for window in windows:
            if window["start"] <= event_data["timestamp"] < window["end"]:
                window["events"].append(event_data)
                break
    
    # Aggregate statistics per window
    window_stats = []
    for window in windows:
        stats = {
            "start_time": window["start"],
            "end_time": window["end"],
            "n_events": len(window["events"]),
            "n_people": len(set(e["person_id"] for e in window["events"])),
        }
        
        # Domain counts
        domain_counts = {}
        for event_data in window["events"]:
            domain = event_data["event"].domain
            domain_counts[domain] = domain_counts.get(domain, 0) + 1
        stats["domain_counts"] = domain_counts
        
        window_stats.append(stats)
    
    total_people = len(set(e["person_id"] for e in all_events))
    time_span_years = (max_time - min_time) / (365.25 * 24 * 3600)
    
    logger.info(
        f"Aggregated {len(all_events)} events from {total_people} people "
        f"into {len(window_stats)} time windows over {time_span_years:.2f} years"
    )
    
    return {
        "time_windows": window_stats,
        "aggregates": {
            "total_events": len(all_events),
            "total_people": total_people,
            "time_span_years": time_span_years,
        }
    }


def map_events_to_traits(
    sequence: EventSequence,
    trait_mapping: Optional[Dict[str, List[str]]] = None
) -> Dict[str, Any]:
    """Map events to phenotypic trait categories.
    
    Maps events in a sequence to predefined trait categories based on
    event types. Uses default mapping if none provided.
    
    Args:
        sequence: EventSequence to map (single sequence, not list)
        trait_mapping: Optional mapping of trait names to event types.
                      If None, uses default mapping for health, education, occupation, and mobility.
        
    Returns:
        Dictionary mapping trait names to dictionaries containing:
        - count: Number of events matching this trait
        - events: List of event types matching this trait
        - timestamps: List of event timestamps (ISO format for datetime)
        
    Raises:
        ImportError: If life_events module is not available
        ValidationError: If sequence is invalid
        
    Examples:
        >>> from metainformant.life_events import EventSequence, Event
        >>> from datetime import datetime
        >>> from metainformant.phenotype import map_events_to_traits
        >>> seq = EventSequence("p1", [
        ...     Event("diagnosis", datetime(2020, 1, 1), "health"),
        ...     Event("degree", datetime(2010, 1, 1), "education"),
        ... ])
        >>> traits = map_events_to_traits(seq)
        >>> "health_issues" in traits
        True
        >>> traits["health_issues"]["count"]
        1
    """
    if not LIFE_EVENTS_AVAILABLE:
        raise ImportError(
            "life_events module not available. Install metainformant with life_events support."
        )
    
    if sequence is None:
        raise ValidationError("sequence cannot be None")
    
    logger.debug(f"Mapping events to traits for person {sequence.person_id}")
    
    if trait_mapping is None:
        # Default mapping
        trait_mapping = {
            "health_issues": ["diagnosis", "treatment", "hospitalization"],
            "education_level": ["degree", "certification", "training"],
            "career_progression": ["job_change", "promotion", "salary_change"],
            "residential_mobility": ["move", "address_change"],
        }
    
    trait_events = {}
    for trait, event_types in trait_mapping.items():
        if not isinstance(event_types, list):
            raise ValidationError(f"trait_mapping[{trait}] must be a list, got {type(event_types).__name__}")
        matching = [
            e for e in sequence.events
            if e.event_type in event_types
        ]
        trait_events[trait] = {
            "count": len(matching),
            "events": [e.event_type for e in matching],
            "timestamps": [
                e.timestamp.isoformat() if isinstance(e.timestamp, datetime) else str(e.timestamp)
                for e in matching
            ]
        }
    
    logger.info(f"Mapped events to {len(trait_events)} trait categories")
    return trait_events

