"""Life course and temporal phenotype analysis."""

from .life_course import (
    Event,
    EventSequence,
    aggregate_temporal_phenotypes,
    analyze_life_course,
    analyze_life_course_trajectories,
    create_life_course_report,
    extract_phenotypes_from_events,
    identify_critical_periods,
    identify_trajectory_patterns,
    map_events_to_traits,
    predict_life_course_outcomes,
)

__all__ = [
    "Event",
    "EventSequence",
    "extract_phenotypes_from_events",
    "aggregate_temporal_phenotypes",
    "map_events_to_traits",
    "analyze_life_course_trajectories",
    "identify_critical_periods",
    "predict_life_course_outcomes",
    "create_life_course_report",
    "identify_trajectory_patterns",
    "analyze_life_course",
]
