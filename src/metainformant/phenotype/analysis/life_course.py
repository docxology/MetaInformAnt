"""Life course phenotype analysis.

This module provides functionality to analyze life course trajectories,
extract temporal phenotypes from event sequences, and model life course
outcomes based on phenotypic data.
"""

from __future__ import annotations

import statistics
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Set, Tuple, Union

from metainformant.core.data import validation
from metainformant.core.utils import errors
from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


@dataclass
class Event:
    """Represents a life course event."""

    timestamp: float
    event_type: str
    description: str = ""
    metadata: Dict[str, Any] = field(default_factory=dict)
    confidence: float = 1.0

    def __post_init__(self):
        validation.validate_type(self.timestamp, (int, float), "timestamp")
        validation.validate_not_empty(self.event_type, "event_type")
        validation.validate_range(self.confidence, min_val=0.0, max_val=1.0, name="confidence")


@dataclass
class EventSequence:
    """Represents a sequence of life course events for an individual."""

    person_id: str
    events: List[Event] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        validation.validate_not_empty(self.person_id, "person_id")
        validation.validate_type(self.events, list, "events")

    @property
    def duration(self) -> float:
        """Total duration of the event sequence."""
        if not self.events:
            return 0.0
        return max(e.timestamp for e in self.events) - min(e.timestamp for e in self.events)

    @property
    def event_count(self) -> int:
        """Number of events in the sequence."""
        return len(self.events)

    def get_events_in_range(self, start_time: float, end_time: float) -> List[Event]:
        """Get events within a time range."""
        return [e for e in self.events if start_time <= e.timestamp <= end_time]

    def get_events_by_type(self, event_type: str) -> List[Event]:
        """Get events of a specific type."""
        return [e for e in self.events if e.event_type == event_type]

    def add_event(self, event: Event) -> None:
        """Add an event to the sequence."""
        validation.validate_type(event, Event, "event")
        self.events.append(event)
        # Keep events sorted by timestamp
        self.events.sort(key=lambda e: e.timestamp)


def extract_phenotypes_from_events(
    event_sequence, phenotype_categories: Optional[Dict[str, List[str]]] = None,
    trait_mapping: Optional[Dict[str, List[str]]] = None,
) -> Dict[str, Any]:
    """Extract phenotypic traits from an event sequence.

    Works with both life_events.core.events.EventSequence and the local
    EventSequence class. Returns a comprehensive phenotype dict with
    person_id, domain counts, event lists, and temporal info.

    Args:
        event_sequence: EventSequence to analyze (from life_events or local)
        phenotype_categories: Optional mapping from category names to event types for custom counting
        trait_mapping: Alias for phenotype_categories (backward compat)

    Returns:
        Dictionary of extracted phenotypes

    Raises:
        errors.ValidationError: If event_sequence is None or invalid
    """
    if event_sequence is None:
        raise errors.ValidationError("event_sequence cannot be None")

    # Accept either kwarg name
    categories = phenotype_categories or trait_mapping

    # Validate custom categories if provided
    if categories is not None:
        if not isinstance(categories, dict):
            raise errors.ValidationError("phenotype_categories must be a dict")
        for key, val in categories.items():
            if not isinstance(val, list):
                raise errors.ValidationError(f"phenotype_categories['{key}'] must be a list, got {type(val).__name__}")

    person_id = getattr(event_sequence, "person_id", "unknown")
    events = getattr(event_sequence, "events", [])

    # Gather domains and event types
    domains = sorted(set(getattr(e, "domain", "unknown") for e in events))
    event_types = sorted(set(getattr(e, "event_type", "unknown") for e in events))

    # Domain-specific counts
    domain_counts: Dict[str, int] = defaultdict(int)
    for e in events:
        domain_counts[getattr(e, "domain", "unknown")] += 1

    phenotypes: Dict[str, Any] = {
        "person_id": person_id,
        "total_events": len(events),
        "domains": domains,
        "event_types": event_types,
        "domain_counts": dict(domain_counts),
    }

    # Domain-specific event counts and lists
    for domain in domains:
        domain_events = [e for e in events if getattr(e, "domain", None) == domain]
        phenotypes[f"{domain}_events"] = len(domain_events)
        # Create descriptive list based on domain
        if domain == "health":
            phenotypes["health_conditions"] = [getattr(e, "event_type", "") for e in domain_events]
        elif domain == "education":
            phenotypes["education_achievements"] = [getattr(e, "event_type", "") for e in domain_events]
        elif domain == "occupation":
            phenotypes["occupation_changes"] = [getattr(e, "event_type", "") for e in domain_events]
        else:
            phenotypes[f"{domain}_items"] = [getattr(e, "event_type", "") for e in domain_events]

    # Temporal information
    if events:
        timestamps = [getattr(e, "timestamp", None) for e in events]
        timestamps = [t for t in timestamps if t is not None]
        if timestamps:
            phenotypes["first_event_time"] = min(timestamps)
            phenotypes["last_event_time"] = max(timestamps)
            first = min(timestamps)
            last = max(timestamps)
            # Handle both datetime and numeric timestamps
            try:
                span = (last - first).days / 365.25
            except (TypeError, AttributeError):
                span = float(last - first)
            phenotypes["event_span_years"] = span

    # Custom phenotype categories
    if categories:
        for cat_name, cat_event_types in categories.items():
            matching = [e for e in events if getattr(e, "event_type", "") in cat_event_types]
            phenotypes[f"{cat_name}_count"] = len(matching)
            phenotypes[f"{cat_name}_events"] = [getattr(e, "event_type", "") for e in matching]

    # Event-type-based backward-compatible analysis (for local EventSequence without domain)
    job_end_events = [e for e in events if getattr(e, "event_type", "") == "job_end"]
    if job_end_events:
        phenotypes["career_changes"] = len(job_end_events)

    health_events = [e for e in events if getattr(e, "event_type", "") == "health_event"]
    if health_events:
        phenotypes["health_score"] = max(0, 10 - len(health_events))

    return phenotypes


def aggregate_temporal_phenotypes(
    sequences, time_window_years: float = 5.0,
    time_windows: Optional[List[Tuple[float, float]]] = None,
    trait_categories: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """Aggregate phenotypes across time windows and individuals.

    Can be called with either explicit time_windows or with time_window_years
    for auto-generated windows based on the event data.

    Args:
        sequences: List of EventSequence objects (life_events or local)
        time_window_years: Window size in years (used when time_windows not provided)
        time_windows: Optional explicit list of (start, end) time tuples
        trait_categories: Optional list of trait categories to analyze

    Returns:
        Dictionary with time_windows and aggregates

    Raises:
        errors.ValidationError: If sequences is invalid or time_window_years <= 0
    """
    if not isinstance(sequences, list):
        raise errors.ValidationError("sequences must be a list")
    if time_window_years <= 0:
        raise errors.ValidationError("time_window_years must be positive")

    # Collect all events across sequences
    all_events = []
    people_with_events = set()
    for seq in sequences:
        events = getattr(seq, "events", [])
        if events:
            people_with_events.add(getattr(seq, "person_id", "unknown"))
        all_events.extend(events)

    total_events = len(all_events)
    total_people = len(people_with_events)

    if not all_events:
        return {
            "time_windows": [],
            "aggregates": {
                "total_events": 0,
                "total_people": 0,
                "time_span_years": 0.0,
            },
        }

    # Get timestamps
    timestamps = [getattr(e, "timestamp", None) for e in all_events]
    timestamps = [t for t in timestamps if t is not None]

    first = min(timestamps)
    last = max(timestamps)
    try:
        time_span = (last - first).days / 365.25
    except (TypeError, AttributeError):
        time_span = float(last - first)

    # Auto-generate time windows if not provided
    if time_windows is None:
        from datetime import timedelta
        generated_windows = []
        current = first
        while current < last:
            try:
                window_end = current + timedelta(days=time_window_years * 365.25)
            except TypeError:
                window_end = current + time_window_years
            if window_end > last:
                window_end = last
            generated_windows.append({"start_time": current, "end_time": window_end})

            # Count events and people in this window
            window_events = [e for e in all_events
                             if current <= getattr(e, "timestamp", current) <= window_end]
            window_people = set()
            for seq in sequences:
                seq_events = getattr(seq, "events", [])
                for e in seq_events:
                    ts = getattr(e, "timestamp", None)
                    if ts is not None and current <= ts <= window_end:
                        window_people.add(getattr(seq, "person_id", "unknown"))
                        break

            domain_counts: Dict[str, int] = defaultdict(int)
            for e in window_events:
                domain_counts[getattr(e, "domain", "unknown")] += 1

            generated_windows[-1]["n_events"] = len(window_events)
            generated_windows[-1]["n_people"] = len(window_people)
            generated_windows[-1]["domain_counts"] = dict(domain_counts)

            try:
                current = current + timedelta(days=time_window_years * 365.25)
            except TypeError:
                current = current + time_window_years

        result_windows = generated_windows
    else:
        result_windows = time_windows

    return {
        "time_windows": result_windows,
        "aggregates": {
            "total_events": total_events,
            "total_people": total_people,
            "time_span_years": time_span,
        },
    }


def map_events_to_traits(
    event_sequence, trait_mapping: Optional[Dict[str, List[str]]] = None,
    trait_definitions: Optional[Dict[str, List[str]]] = None,
) -> Dict[str, Any]:
    """Map events to phenotypic traits based on definitions.

    Args:
        event_sequence: EventSequence to analyze (from life_events or local)
        trait_mapping: Dictionary mapping trait names to event types
        trait_definitions: Alias for trait_mapping (backward compat)

    Returns:
        Dictionary mapping trait names to {count, events, timestamps} dicts

    Raises:
        errors.ValidationError: If event_sequence is None or mapping invalid
    """
    if event_sequence is None:
        raise errors.ValidationError("event_sequence cannot be None")

    mapping = trait_mapping or trait_definitions
    if mapping is None:
        # Default domain-based trait mapping
        mapping = {
            "health_issues": ["diagnosis", "medical_diagnosis", "hospitalization", "health_event", "treatment"],
            "education_level": ["degree", "education_complete", "degree_earned", "education_start"],
            "career_progression": ["job_change", "job_start", "job_end", "career_change", "promotion"],
        }

    # Validate mapping
    if not isinstance(mapping, dict):
        raise errors.ValidationError("trait_mapping must be a dict")
    for key, val in mapping.items():
        if not isinstance(val, list):
            raise errors.ValidationError(f"trait_mapping['{key}'] must be a list, got {type(val).__name__}")

    events = getattr(event_sequence, "events", [])

    trait_results: Dict[str, Any] = {}
    for trait_name, event_types in mapping.items():
        relevant_events = [e for e in events if getattr(e, "event_type", "") in event_types]
        relevant_events.sort(key=lambda e: getattr(e, "timestamp", 0))

        trait_results[trait_name] = {
            "count": len(relevant_events),
            "events": [getattr(e, "event_type", "") for e in relevant_events],
            "timestamps": [getattr(e, "timestamp", None) for e in relevant_events],
        }

    return trait_results


def analyze_life_course_trajectories(
    sequences: List[EventSequence], outcome_measures: Optional[List[str]] = None
) -> Dict[str, Any]:
    """Analyze life course trajectories and their outcomes.

    Args:
        sequences: List of EventSequence objects
        outcome_measures: Optional list of outcome measures to analyze

    Returns:
        Dictionary with trajectory analysis results
    """
    if outcome_measures is None:
        outcome_measures = ["career_success", "health_outcomes", "social_stability"]

    results = {
        "total_sequences": len(sequences),
        "trajectory_patterns": {},
        "outcome_analysis": {},
        "risk_factors": {},
    }

    # Analyze trajectory patterns
    trajectory_patterns = defaultdict(int)

    for sequence in sequences:
        # Classify trajectory based on event patterns
        if sequence.event_count > 20:
            pattern = "high_activity"
        elif sequence.event_count > 10:
            pattern = "moderate_activity"
        else:
            pattern = "low_activity"

        # Check for stability indicators
        job_stability = len(sequence.get_events_by_type("job_end")) <= 2
        relationship_stability = len(sequence.get_events_by_type("relationship_end")) <= 1

        if job_stability and relationship_stability:
            pattern += "_stable"
        else:
            pattern += "_unstable"

        trajectory_patterns[pattern] += 1

    results["trajectory_patterns"] = dict(trajectory_patterns)

    # Analyze outcomes (simplified analysis)
    for outcome in outcome_measures:
        if outcome == "career_success":
            # Simple career success based on job events
            career_scores = []
            for sequence in sequences:
                job_events = sequence.get_events_by_type("job_start")
                career_score = min(10, len(job_events) * 2)  # More jobs = higher score (simplified)
                career_scores.append(career_score)

            results["outcome_analysis"][outcome] = {
                "mean_score": statistics.mean(career_scores) if career_scores else 0,
                "distribution": career_scores,
            }

        elif outcome == "health_outcomes":
            # Health outcomes based on health events
            health_scores = []
            for sequence in sequences:
                health_events = sequence.get_events_by_type("health_event")
                health_score = max(0, 10 - len(health_events))  # Fewer health events = better health
                health_scores.append(health_score)

            results["outcome_analysis"][outcome] = {
                "mean_score": statistics.mean(health_scores) if health_scores else 0,
                "distribution": health_scores,
            }

        elif outcome == "social_stability":
            # Social stability based on relationship patterns
            stability_scores = []
            for sequence in sequences:
                relationships = sequence.get_events_by_type("relationship_start")
                breakups = sequence.get_events_by_type("relationship_end")
                stability_score = max(0, len(relationships) - len(breakups))
                stability_score = min(10, stability_score)
                stability_scores.append(stability_score)

            results["outcome_analysis"][outcome] = {
                "mean_score": statistics.mean(stability_scores) if stability_scores else 0,
                "distribution": stability_scores,
            }

    # Identify risk factors
    risk_factors = {}

    # Individuals with high health events
    high_health_risk = []
    for sequence in sequences:
        health_events = len(sequence.get_events_by_type("health_event"))
        if health_events > 5:
            high_health_risk.append(sequence.person_id)

    if high_health_risk:
        risk_factors["high_health_events"] = high_health_risk

    # Individuals with frequent job changes
    job_instability_risk = []
    for sequence in sequences:
        job_changes = len(sequence.get_events_by_type("job_end"))
        if job_changes > 3:
            job_instability_risk.append(sequence.person_id)

    if job_instability_risk:
        risk_factors["job_instability"] = job_instability_risk

    results["risk_factors"] = risk_factors

    return results


def identify_critical_periods(sequences: List[EventSequence], age_ranges: List[Tuple[float, float]]) -> Dict[str, Any]:
    """Identify critical developmental periods based on event clustering.

    Args:
        sequences: List of EventSequence objects
        age_ranges: List of (start_age, end_age) tuples defining periods

    Returns:
        Dictionary with critical period analysis
    """
    validation.validate_not_empty(sequences, "sequences")
    validation.validate_not_empty(age_ranges, "age_ranges")

    results = {
        "age_ranges": age_ranges,
        "period_importance": {},
        "event_clusters": {},
    }

    # Analyze event density in each age range
    for start_age, end_age in age_ranges:
        period_key = f"{start_age}_{end_age}"
        total_events = 0
        event_types = defaultdict(int)

        for sequence in sequences:
            period_events = sequence.get_events_in_range(start_age, end_age)
            total_events += len(period_events)

            for event in period_events:
                event_types[event.event_type] += 1

        period_duration = end_age - start_age
        events_per_year = total_events / (period_duration * len(sequences)) if period_duration > 0 else 0

        results["period_importance"][period_key] = {
            "total_events": total_events,
            "events_per_year": events_per_year,
            "event_types": dict(event_types),
            "sequences_with_events": sum(1 for s in sequences if s.get_events_in_range(start_age, end_age)),
        }

    # Identify event clusters (periods with high event density)
    event_clusters = []
    importance_scores = [data["events_per_year"] for data in results["period_importance"].values()]

    if importance_scores:
        mean_importance = statistics.mean(importance_scores)
        std_importance = statistics.stdev(importance_scores) if len(importance_scores) > 1 else 0

        for period_key, data in results["period_importance"].items():
            if data["events_per_year"] > mean_importance + std_importance:
                event_clusters.append(
                    {
                        "period": period_key,
                        "importance_score": data["events_per_year"],
                        "significance": "high",
                    }
                )

    results["event_clusters"] = event_clusters

    return results


def predict_life_course_outcomes(sequences: List[EventSequence], prediction_horizon: float = 5.0) -> Dict[str, Any]:
    """Predict future life course outcomes based on current trajectories.

    Args:
        sequences: List of EventSequence objects with historical data
        prediction_horizon: Years to predict into the future

    Returns:
        Dictionary with outcome predictions
    """
    validation.validate_not_empty(sequences, "sequences")
    validation.validate_range(prediction_horizon, min_val=0.1, name="prediction_horizon")

    results = {
        "prediction_horizon": prediction_horizon,
        "outcome_predictions": {},
        "confidence_intervals": {},
    }

    # Simple prediction model based on current trajectory patterns
    predictions = {}

    for sequence in sequences:
        person_id = sequence.person_id

        # Predict career trajectory
        job_events = sequence.get_events_by_type("job_start")
        career_growth_rate = len(job_events) / max(1, sequence.duration / 365.25)  # jobs per year

        # Simple linear extrapolation
        predicted_jobs = int(career_growth_rate * prediction_horizon)
        predictions[f"{person_id}_career"] = {
            "current_jobs": len(job_events),
            "predicted_additional_jobs": predicted_jobs,
            "confidence": 0.7,  # Simplified confidence
        }

        # Predict health trajectory
        health_events = sequence.get_events_by_type("health_event")
        health_rate = len(health_events) / max(1, sequence.duration / 365.25)

        predicted_health_events = int(health_rate * prediction_horizon)
        predictions[f"{person_id}_health"] = {
            "current_health_events": len(health_events),
            "predicted_health_events": predicted_health_events,
            "risk_level": (
                "high" if predicted_health_events > 2 else "moderate" if predicted_health_events > 0 else "low"
            ),
        }

        # Predict social trajectory
        social_events = sequence.get_events_by_type("social_event")
        social_rate = len(social_events) / max(1, sequence.duration / 365.25)

        predicted_social_events = int(social_rate * prediction_horizon)
        predictions[f"{person_id}_social"] = {
            "current_social_events": len(social_events),
            "predicted_social_events": predicted_social_events,
        }

    results["outcome_predictions"] = predictions

    # Calculate confidence intervals (simplified)
    confidence_intervals = {}
    for key, prediction in predictions.items():
        if "predicted" in key:
            base_value = prediction.get(
                "current_jobs", prediction.get("current_health_events", prediction.get("current_social_events", 0))
            )
            predicted_value = prediction.get(
                "predicted_additional_jobs",
                prediction.get("predicted_health_events", prediction.get("predicted_social_events", 0)),
            )
            total = base_value + predicted_value

            # Simple confidence interval
            margin = int(total * 0.2)  # 20% margin
            confidence_intervals[key] = {
                "predicted": total,
                "lower_bound": max(0, total - margin),
                "upper_bound": total + margin,
            }

    results["confidence_intervals"] = confidence_intervals

    return results


def create_life_course_report(sequences: List[EventSequence], output_path: Optional[str | Path] = None) -> str:
    """Generate a comprehensive life course analysis report.

    Args:
        output_path: Optional path to save the report
        sequences: List of EventSequence objects to analyze

    Returns:
        Formatted report string
    """
    report_lines = []
    report_lines.append("=" * 70)
    report_lines.append("LIFE COURSE PHENOTYPE ANALYSIS REPORT")
    report_lines.append("=" * 70)
    report_lines.append("")

    # Basic statistics
    total_sequences = len(sequences)
    total_events = sum(len(s.events) for s in sequences)
    avg_events_per_person = total_events / total_sequences if total_sequences > 0 else 0

    report_lines.append(f"Total Individuals: {total_sequences}")
    report_lines.append(f"Total Events: {total_events}")
    report_lines.append(f"Average Events per Individual: {avg_events_per_person:.1f}")
    report_lines.append("")

    # Event type distribution
    event_types = defaultdict(int)
    for sequence in sequences:
        for event in sequence.events:
            event_types[event.event_type] += 1

    report_lines.append("Event Type Distribution:")
    for event_type, count in sorted(event_types.items(), key=lambda x: x[1], reverse=True):
        percentage = (count / total_events * 100) if total_events > 0 else 0
        report_lines.append(f"  {event_type}: {count} ({percentage:.1f}%)")
    report_lines.append("")

    # Trajectory analysis
    trajectory_analysis = analyze_life_course_trajectories(sequences)

    report_lines.append("Trajectory Patterns:")
    for pattern, count in trajectory_analysis["trajectory_patterns"].items():
        percentage = (count / total_sequences * 100) if total_sequences > 0 else 0
        report_lines.append(f"  {pattern}: {count} individuals ({percentage:.1f}%)")
    report_lines.append("")

    # Outcome analysis
    if trajectory_analysis["outcome_analysis"]:
        report_lines.append("Outcome Measures:")
        for outcome, data in trajectory_analysis["outcome_analysis"].items():
            mean_score = data.get("mean_score", 0)
            report_lines.append(f"  {outcome}: {mean_score:.1f} (mean score)")
        report_lines.append("")

    # Risk factors
    if trajectory_analysis["risk_factors"]:
        report_lines.append("Identified Risk Factors:")
        for risk_type, individuals in trajectory_analysis["risk_factors"].items():
            report_lines.append(f"  {risk_type}: {len(individuals)} individuals")
        report_lines.append("")

    report = "\n".join(report_lines)

    if output_path:
        output_path = Path(output_path)
        with open(output_path, "w") as f:
            f.write(report)
        logger.info(f"Life course report saved to {output_path}")

    return report


def identify_trajectory_patterns(sequences: List[EventSequence]) -> Dict[str, Any]:
    """Identify trajectory patterns in event sequences using transition analysis.

    Classifies sequences by activity level and event type diversity, then
    identifies common transition patterns between event types.

    Args:
        sequences: List of EventSequence objects

    Returns:
        Dictionary with pattern classifications and transition frequencies
    """
    if not sequences:
        return {"patterns": {}, "transitions": {}}

    pattern_counts: Dict[str, int] = defaultdict(int)
    transition_counts: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))

    for seq in sequences:
        if not seq.events:
            pattern_counts["empty"] += 1
            continue

        unique_types = set(e.event_type for e in seq.events)
        n_events = len(seq.events)

        # Classify by activity and diversity
        if n_events > 20 and len(unique_types) > 5:
            pattern = "complex_high_activity"
        elif n_events > 20:
            pattern = "repetitive_high_activity"
        elif n_events > 10 and len(unique_types) > 3:
            pattern = "diverse_moderate"
        elif n_events > 10:
            pattern = "focused_moderate"
        elif len(unique_types) > 3:
            pattern = "diverse_low"
        else:
            pattern = "simple_low"

        pattern_counts[pattern] += 1

        # Track transitions
        sorted_events = sorted(seq.events, key=lambda e: e.timestamp)
        for i in range(len(sorted_events) - 1):
            from_type = sorted_events[i].event_type
            to_type = sorted_events[i + 1].event_type
            transition_counts[from_type][to_type] += 1

    # Normalize transitions to probabilities
    transition_probs: Dict[str, Dict[str, float]] = {}
    for from_type, targets in transition_counts.items():
        total = sum(targets.values())
        transition_probs[from_type] = {to: count / total for to, count in targets.items()}

    return {
        "patterns": dict(pattern_counts),
        "transitions": transition_probs,
        "n_sequences": len(sequences),
    }


def analyze_life_course(sequences: List[EventSequence], outcomes: List[str] | None = None) -> Dict[str, Any]:
    """Comprehensive analysis of life course trajectories.

    Args:
        sequences: List of event sequences to analyze
        outcomes: Optional list of outcome measures to analyze

    Returns:
        Dictionary containing comprehensive life course analysis
    """
    if not sequences:
        return {"error": "No sequences provided"}

    results = {
        "n_sequences": len(sequences),
        "analysis_timestamp": statistics.mean([seq.duration for seq in sequences if seq.events]) if sequences else 0,
        "components": {},
    }

    # Basic statistics
    durations = [seq.duration for seq in sequences if seq.events]
    if durations:
        results["duration_stats"] = {
            "mean": statistics.mean(durations),
            "median": statistics.median(durations),
            "min": min(durations),
            "max": max(durations),
            "std": statistics.stdev(durations) if len(durations) > 1 else 0,
        }

    # Event frequency analysis
    all_events = []
    for seq in sequences:
        all_events.extend([event.event_type for event in seq.events])

    if all_events:
        event_counts = defaultdict(int)
        for event in all_events:
            event_counts[event] += 1

        results["event_frequency"] = dict(event_counts)

    # Sequence complexity analysis
    complexities = []
    for seq in sequences:
        if seq.events:
            # Simple complexity measure: number of unique event types / total events
            unique_events = len(set(event.event_type for event in seq.events))
            total_events = len(seq.events)
            complexity = unique_events / total_events if total_events > 0 else 0
            complexities.append(complexity)

    if complexities:
        results["complexity_stats"] = {
            "mean": statistics.mean(complexities),
            "median": statistics.median(complexities),
            "min": min(complexities),
            "max": max(complexities),
        }

    # Outcome analysis (placeholder)
    if outcomes:
        results["outcome_analysis"] = {}
        for outcome in outcomes:
            # This would integrate with predictive models
            results["outcome_analysis"][outcome] = {
                "n_with_outcome": sum(1 for seq in sequences if outcome in seq.metadata.get("outcomes", [])),
                "prevalence": sum(1 for seq in sequences if outcome in seq.metadata.get("outcomes", []))
                / len(sequences),
            }

    # Trajectory patterns
    results["trajectory_patterns"] = identify_trajectory_patterns(sequences)

    return results
