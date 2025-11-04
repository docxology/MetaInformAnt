"""Life course phenotype extraction from event sequences.

This module provides functions to extract phenotypes from event sequences,
enabling temporal phenotype aggregation and integration with life_events module.
"""

from __future__ import annotations

from datetime import datetime
from typing import Any, Dict, List, Optional

import numpy as np

try:
    from ..life_events import Event, EventSequence
    LIFE_EVENTS_AVAILABLE = True
except ImportError:
    LIFE_EVENTS_AVAILABLE = False
    Event = None  # type: ignore
    EventSequence = None  # type: ignore


def extract_phenotypes_from_events(
    sequence: EventSequence,
    phenotype_categories: Optional[Dict[str, List[str]]] = None
) -> Dict[str, Any]:
    """Extract phenotype traits from event sequence.
    
    Maps events to phenotypic categories and aggregates temporal information.
    
    Args:
        sequence: EventSequence to extract phenotypes from
        phenotype_categories: Optional mapping of phenotype categories to event types
        
    Returns:
        Dictionary with extracted phenotypes
        
    Examples:
        >>> from metainformant.life_events import EventSequence, Event
        >>> from datetime import datetime
        >>> events = [
        ...     Event("diabetes", datetime(2020, 1, 1), "health"),
        ...     Event("bachelors", datetime(2010, 6, 1), "education"),
        ... ]
        >>> seq = EventSequence(person_id="p1", events=events)
        >>> phenotypes = extract_phenotypes_from_events(seq)
        >>> "health_events" in phenotypes
        True
    """
    if not LIFE_EVENTS_AVAILABLE:
        raise ImportError("life_events module not available")
    
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
            matching_events = [
                e for e in sequence.events
                if e.event_type in event_types
            ]
            phenotypes[f"{category}_count"] = len(matching_events)
    
    return phenotypes


def aggregate_temporal_phenotypes(
    sequences: List[EventSequence],
    time_window_years: float = 5.0
) -> Dict[str, Any]:
    """Aggregate phenotypes across time windows.
    
    Groups events into temporal windows and computes aggregate statistics.
    
    Args:
        sequences: List of event sequences
        time_window_years: Size of time windows in years
        
    Returns:
        Dictionary with temporal phenotype aggregations
    """
    if not LIFE_EVENTS_AVAILABLE:
        raise ImportError("life_events module not available")
    
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
        return {"time_windows": [], "aggregates": {}}
    
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
    
    return {
        "time_windows": window_stats,
        "aggregates": {
            "total_events": len(all_events),
            "total_people": len(set(e["person_id"] for e in all_events)),
            "time_span_years": (max_time - min_time) / (365.25 * 24 * 3600),
        }
    }


def map_events_to_traits(
    sequence: EventSequence,
    trait_mapping: Optional[Dict[str, List[str]]] = None
) -> Dict[str, Any]:
    """Map events to phenotypic trait categories.
    
    Args:
        sequence: EventSequence to map
        trait_mapping: Optional mapping of trait names to event types
        
    Returns:
        Dictionary mapping traits to event information
    """
    if not LIFE_EVENTS_AVAILABLE:
        raise ImportError("life_events module not available")
    
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
        matching = [
            e for e in sequence.events
            if e.event_type in event_types
        ]
        trait_events[trait] = {
            "count": len(matching),
            "events": [e.event_type for e in matching],
            "timestamps": [
                e.timestamp.isoformat() if isinstance(e.timestamp, datetime) else e.timestamp
                for e in matching
            ]
        }
    
    return trait_events

