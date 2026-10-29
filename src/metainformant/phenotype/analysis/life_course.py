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
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy import stats as scipy_stats

from metainformant.core.data import validation
from metainformant.core.utils import errors, logging

logger = logging.get_logger(__name__)


@dataclass
class Event:
    """Represents a life course event."""

    timestamp: float
    event_type: str
    description: str = ""
    metadata: Dict[str, Any] = field(default_factory=dict)
    confidence: float = 1.0

    def __post_init__(self) -> None:
        validation.validate_type(self.timestamp, (int, float), "timestamp")
        validation.validate_not_empty(self.event_type, "event_type")
        validation.validate_range(self.confidence, min_val=0.0, max_val=1.0, name="confidence")


@dataclass
class EventSequence:
    """Represents a sequence of life course events for an individual."""

    person_id: str
    events: List[Event] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
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
    event_sequence: Any,
    phenotype_categories: Optional[Dict[str, List[str]]] = None,
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
        raw_ts = [getattr(e, "timestamp", None) for e in events]
        timestamps: list[Any] = [t for t in raw_ts if t is not None]
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
    sequences: List[Any],
    time_window_years: float = 5.0,
    time_windows: Optional[List[Tuple[float, float]]] = None,
    trait_categories: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """Aggregate phenotypes across time windows and individuals.

    Can be called with either explicit time_windows or with time_window_years
    for auto-generated windows based on the event data.

    Windows are half-open, ``[start, end)``, except the last one, which includes its end: an event on the
    boundary between two windows belongs to exactly one of them, so the window counts add up to the number of
    events (a closed interval on both sides counted it in both).

    Args:
        sequences: List of EventSequence objects (life_events or local)
        time_window_years: Window size in years (used when time_windows not provided)
        time_windows: Optional explicit list of (start, end) time tuples
        trait_categories: Optional list of trait categories (event domains) to analyze: only events whose
            ``domain`` is in the list are counted, in the totals and in every window. ``None`` keeps all events.

    Returns:
        Dictionary with ``time_windows`` (the generated window dictionaries, each with ``start_time``,
        ``end_time``, ``n_events``, ``n_people`` and ``domain_counts``; or the explicit windows as given),
        ``window_aggregates`` (those same per-window counts for either kind of window) and ``aggregates``
        (``total_events``, ``total_people``, ``time_span_years``). Events without a timestamp are counted in
        the totals but belong to no window.

    Raises:
        errors.ValidationError: If sequences is invalid or time_window_years <= 0
    """
    if not isinstance(sequences, list):
        raise errors.ValidationError("sequences must be a list")
    if time_window_years <= 0:
        raise errors.ValidationError("time_window_years must be positive")

    category_filter = set(trait_categories) if trait_categories is not None else None

    def keep(event: Any) -> bool:
        return category_filter is None or getattr(event, "domain", None) in category_filter

    # Collect the events (with their person) across sequences
    all_events: list[Any] = []
    event_people: list[Any] = []
    people_with_events = set()
    for seq in sequences:
        person = getattr(seq, "person_id", "unknown")
        for event in getattr(seq, "events", []):
            if not keep(event):
                continue
            people_with_events.add(person)
            all_events.append(event)
            event_people.append(person)

    total_events = len(all_events)
    total_people = len(people_with_events)

    if not all_events:
        return {
            "time_windows": [],
            "window_aggregates": [],
            "aggregates": {
                "total_events": 0,
                "total_people": 0,
                "time_span_years": 0.0,
            },
        }

    # Get timestamps
    stamped: list[Tuple[Any, Any, Any]] = [
        (getattr(e, "timestamp", None), e, person) for e, person in zip(all_events, event_people)
    ]
    stamped = [(ts, e, person) for ts, e, person in stamped if ts is not None]
    timestamps: list[Any] = [ts for ts, _, _ in stamped]

    if timestamps:
        first = min(timestamps)
        last = max(timestamps)
        try:
            time_span = (last - first).days / 365.25
        except (TypeError, AttributeError):
            time_span = float(last - first)
    else:
        time_span = 0.0

    def window_counts(start: Any, end: Any, closed_end: bool) -> Dict[str, Any]:
        inside = [(e, person) for ts, e, person in stamped if start <= ts and (ts <= end if closed_end else ts < end)]
        domain_counts: Dict[str, int] = defaultdict(int)
        for e, _ in inside:
            domain_counts[getattr(e, "domain", "unknown")] += 1
        return {
            "n_events": len(inside),
            "n_people": len({person for _, person in inside}),
            "domain_counts": dict(domain_counts),
        }

    # Auto-generate time windows if not provided
    if time_windows is None:
        from datetime import timedelta

        def advance(point: Any) -> Any:
            try:
                return point + timedelta(days=time_window_years * 365.25)
            except TypeError:
                return point + time_window_years

        edges = []
        if timestamps:
            current = first
            while current < last:
                edges.append((current, min(advance(current), last)))
                current = advance(current)
            if not edges:  # every event shares one timestamp: a single, degenerate window holds them all
                edges.append((first, last))

        generated_windows = [
            {"start_time": start, "end_time": end, **window_counts(start, end, closed_end=index == len(edges) - 1)}
            for index, (start, end) in enumerate(edges)
        ]
        result_windows_list: list[Any] = list(generated_windows)
        window_aggregates: list[Any] = list(generated_windows)
    else:
        result_windows_list = list(time_windows)
        window_aggregates = [
            {
                "start_time": start,
                "end_time": end,
                **window_counts(start, end, closed_end=index == len(result_windows_list) - 1),
            }
            for index, (start, end) in enumerate(result_windows_list)
        ]

    return {
        "time_windows": result_windows_list,
        "window_aggregates": window_aggregates,
        "aggregates": {
            "total_events": total_events,
            "total_people": total_people,
            "time_span_years": time_span,
        },
    }


def map_events_to_traits(
    event_sequence: Any,
    trait_mapping: Optional[Dict[str, List[str]]] = None,
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

    The outcome measures are simple counting indices on a 0-10 scale, not validated instruments:

    * ``career_success``: ``min(10, 2 * number of job_start events)``
    * ``health_outcomes``: ``max(0, 10 - number of health_event events)``
    * ``social_stability``: ``min(10, max(0, relationship_start - relationship_end events))``

    Args:
        sequences: List of EventSequence objects
        outcome_measures: Optional list of outcome measures to analyze; the three above are the default

    Returns:
        Dictionary with trajectory analysis results

    Raises:
        ValueError: If an outcome measure is not one of the three above (an unknown name used to be skipped
            without a word, leaving it out of the result).
    """
    supported_outcomes = ("career_success", "health_outcomes", "social_stability")
    if outcome_measures is None:
        outcome_measures = list(supported_outcomes)
    unknown = [outcome for outcome in outcome_measures if outcome not in supported_outcomes]
    if unknown:
        raise ValueError(f"Unknown outcome measure(s) {unknown}; supported: {list(supported_outcomes)}")

    results: Dict[str, Any] = {
        "total_sequences": len(sequences),
        "trajectory_patterns": {},
        "outcome_analysis": {},
        "risk_factors": {},
    }

    # Analyze trajectory patterns
    trajectory_patterns: defaultdict[str, int] = defaultdict(int)

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

    # Analyze outcomes
    for outcome in outcome_measures:
        if outcome == "career_success":
            # Simple career success based on job events
            career_scores = []
            for sequence in sequences:
                job_events = sequence.get_events_by_type("job_start")
                career_score = min(10, len(job_events) * 2)  # two points per job started, capped at 10
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
        n_health_events = len(sequence.get_events_by_type("health_event"))
        if n_health_events > 5:
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

    results: Dict[str, Any] = {
        "age_ranges": age_ranges,
        "period_importance": {},
        "event_clusters": {},
    }

    # Analyze event density in each age range
    for start_age, end_age in age_ranges:
        period_key = f"{start_age}_{end_age}"
        total_events = 0
        event_types: defaultdict[str, int] = defaultdict(int)

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


def _poisson_forecast(n_observed: int, exposure_years: float, horizon_years: float) -> Dict[str, Any]:
    """Forecast of the events in the next *horizon_years* from *n_observed* events in *exposure_years*.

    Events are a Poisson process with an unknown rate. With Jeffreys' prior the rate's posterior is
    ``Gamma(n + 1/2, exposure)``, so the count in the horizon follows a negative binomial with size ``n + 1/2``
    and success probability ``exposure / (exposure + horizon)``. ``expected`` is the maximum-likelihood
    extrapolation ``n / exposure * horizon`` and ``lower`` / ``upper`` are the 2.5% and 97.5% quantiles of the
    predictive count (an equal-tailed 95% interval).
    """
    expected = n_observed / exposure_years * horizon_years
    size = n_observed + 0.5
    probability = exposure_years / (exposure_years + horizon_years)
    quantiles = np.atleast_1d(scipy_stats.nbinom.ppf([0.025, 0.975], size, probability))
    return {"expected": expected, "lower": int(quantiles[0]), "upper": int(quantiles[1])}


def predict_life_course_outcomes(sequences: List[EventSequence], prediction_horizon: float = 5.0) -> Dict[str, Any]:
    """Predict future life course outcomes based on current trajectories.

    Event timestamps are taken to be in days. For each person the rate of job starts, health events and social
    events over the observed span (at least one year, so a short record does not extrapolate to an absurd rate)
    is extrapolated linearly over ``prediction_horizon`` years. The count is truncated to a whole number, and the
    uncertainty is a 95% Poisson-process predictive interval (see :func:`_poisson_forecast`) for the
    additional events, reported in ``confidence_intervals`` as bounds on the total count (observed plus
    additional).

    Args:
        sequences: List of EventSequence objects with historical data
        prediction_horizon: Years to predict into the future

    Returns:
        Dictionary with ``prediction_horizon``, ``outcome_predictions`` (keyed ``<person_id>_career``,
        ``_health`` and ``_social``) and ``confidence_intervals`` (the same keys, each with ``predicted``,
        ``lower_bound``, ``upper_bound`` and ``level``).
    """
    validation.validate_not_empty(sequences, "sequences")
    validation.validate_range(prediction_horizon, min_val=0.1, name="prediction_horizon")

    results: Dict[str, Any] = {
        "prediction_horizon": prediction_horizon,
        "outcome_predictions": {},
        "confidence_intervals": {},
    }

    predictions: Dict[str, Any] = {}
    confidence_intervals: Dict[str, Any] = {}

    def record(key: str, current: int, forecast: Dict[str, Any]) -> int:
        predicted = int(forecast["expected"])
        confidence_intervals[key] = {
            "predicted": current + predicted,
            "lower_bound": current + forecast["lower"],
            "upper_bound": current + forecast["upper"],
            "level": 0.95,
        }
        return predicted

    for sequence in sequences:
        person_id = sequence.person_id
        exposure_years = max(1.0, sequence.duration / 365.25)

        # Predict career trajectory
        job_events = sequence.get_events_by_type("job_start")
        key = f"{person_id}_career"
        predicted_jobs = record(
            key, len(job_events), _poisson_forecast(len(job_events), exposure_years, prediction_horizon)
        )
        predictions[key] = {
            "current_jobs": len(job_events),
            "predicted_additional_jobs": predicted_jobs,
            "observed_years": exposure_years,
        }

        # Predict health trajectory
        health_events = sequence.get_events_by_type("health_event")
        key = f"{person_id}_health"
        predicted_health_events = record(
            key, len(health_events), _poisson_forecast(len(health_events), exposure_years, prediction_horizon)
        )
        predictions[key] = {
            "current_health_events": len(health_events),
            "predicted_health_events": predicted_health_events,
            "risk_level": (
                "high" if predicted_health_events > 2 else "moderate" if predicted_health_events > 0 else "low"
            ),
            "observed_years": exposure_years,
        }

        # Predict social trajectory
        social_events = sequence.get_events_by_type("social_event")
        key = f"{person_id}_social"
        predicted_social_events = record(
            key, len(social_events), _poisson_forecast(len(social_events), exposure_years, prediction_horizon)
        )
        predictions[key] = {
            "current_social_events": len(social_events),
            "predicted_social_events": predicted_social_events,
            "observed_years": exposure_years,
        }

    results["outcome_predictions"] = predictions
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
    event_types: defaultdict[str, int] = defaultdict(int)
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


def _recorded_outcomes(sequence: EventSequence) -> Optional[set[str]]:
    """The outcome names recorded for a person in ``metadata["outcomes"]``, or None when none are recorded.

    A sequence without that entry has an unknown outcome status, which is not the same as having no outcome:
    it is left out of the outcome analysis instead of being counted as negative.
    """
    recorded = sequence.metadata.get("outcomes")
    if recorded is None:
        return None
    if isinstance(recorded, str):
        return {recorded}
    try:
        return set(recorded)
    except TypeError:
        raise errors.ValidationError(
            f"metadata['outcomes'] of sequence {sequence.person_id!r} must be a collection of outcome names, "
            f"got {type(recorded).__name__}"
        ) from None


def _rank_sum_p_value(with_outcome: List[float], without_outcome: List[float]) -> float:
    """Two-sided Mann-Whitney U p-value; 1.0 when every value is the same (scipy returns NaN there)."""
    if len(set(with_outcome) | set(without_outcome)) == 1:
        return 1.0
    return float(scipy_stats.mannwhitneyu(with_outcome, without_outcome, alternative="two-sided").pvalue)


def _compare_feature(with_outcome: List[float], without_outcome: List[float]) -> Optional[Dict[str, Any]]:
    """Mean of a trajectory feature in the groups with and without the outcome, and a rank-sum test."""
    if not with_outcome or not without_outcome:
        return None
    mean_with = statistics.fmean(with_outcome)
    mean_without = statistics.fmean(without_outcome)
    return {
        "n_with_outcome": len(with_outcome),
        "n_without_outcome": len(without_outcome),
        "mean_with_outcome": mean_with,
        "mean_without_outcome": mean_without,
        "difference": mean_with - mean_without,
        "p_value": _rank_sum_p_value(with_outcome, without_outcome),
    }


def _event_type_associations(rows: List[Tuple[set[str], bool]]) -> Dict[str, Dict[str, Any]]:
    """Association of having each event type with the outcome, as Fisher exact tests on 2x2 tables.

    ``rows`` holds, per sequence with a known outcome status, the set of its event types and whether it has the
    outcome. The result is ordered by p-value (ties by event type).
    """
    n_with = sum(1 for _, has_outcome in rows if has_outcome)
    n_without = len(rows) - n_with
    associations: Dict[str, Dict[str, Any]] = {}
    for event_type in sorted({event_type for types, _ in rows for event_type in types}):
        exposed_with = sum(1 for types, has_outcome in rows if has_outcome and event_type in types)
        exposed_without = sum(1 for types, has_outcome in rows if not has_outcome and event_type in types)
        unexposed_with = n_with - exposed_with
        unexposed_without = n_without - exposed_without
        denominator = exposed_without * unexposed_with
        associations[event_type] = {
            "n_exposed": exposed_with + exposed_without,
            "n_exposed_with_outcome": exposed_with,
            "n_unexposed_with_outcome": unexposed_with,
            # Sample odds ratio; None when a zero cell makes it infinite or undefined
            "odds_ratio": exposed_with * unexposed_without / denominator if denominator else None,
            "p_value": float(
                scipy_stats.fisher_exact([[exposed_with, exposed_without], [unexposed_with, unexposed_without]]).pvalue
            ),
        }
    # p-values that agree to 12 decimals are ties (the exact sums behind them can differ in the last bits)
    return dict(sorted(associations.items(), key=lambda item: (round(item[1]["p_value"], 12), item[0])))


def _analyze_outcome(
    sequences: List[EventSequence], recorded: List[Optional[set[str]]], outcome: str
) -> Dict[str, Any]:
    """Prevalence of one outcome and how it relates to the trajectories, among sequences with a known status."""
    observed = [(seq, outcome in names) for seq, names in zip(sequences, recorded) if names is not None]
    n_with = sum(1 for _, has_outcome in observed if has_outcome)
    analysis: Dict[str, Any] = {
        "status": "ok" if observed else "no_outcome_data",
        "n_sequences": len(sequences),
        "n_with_outcome_data": len(observed),
        "n_missing_outcome_data": len(sequences) - len(observed),
        "n_with_outcome": n_with,
        "n_without_outcome": len(observed) - n_with,
        "prevalence": n_with / len(observed) if observed else None,
        "trajectory_comparison": {},
        "event_type_association": {},
    }
    if not n_with or n_with == len(observed):
        return analysis  # a comparison needs sequences both with and without the outcome

    def feature_values(has_outcome: bool, feature: str) -> List[float]:
        values = {
            "n_events": [float(len(seq.events)) for seq, has in observed if has == has_outcome],
            # duration and complexity are only defined for sequences that have events
            "duration": [float(seq.duration) for seq, has in observed if has == has_outcome and seq.events],
            "complexity": [
                len({event.event_type for event in seq.events}) / len(seq.events)
                for seq, has in observed
                if has == has_outcome and seq.events
            ],
        }
        return values[feature]

    analysis["trajectory_comparison"] = {
        feature: _compare_feature(feature_values(True, feature), feature_values(False, feature))
        for feature in ("n_events", "duration", "complexity")
    }
    analysis["event_type_association"] = _event_type_associations(
        [({event.event_type for event in seq.events}, has_outcome) for seq, has_outcome in observed]
    )
    return analysis


def analyze_life_course(sequences: List[EventSequence], outcomes: List[str] | None = None) -> Dict[str, Any]:
    """Comprehensive analysis of life course trajectories.

    Args:
        sequences: List of event sequences to analyze
        outcomes: Optional list of outcome names to analyze. A person's outcomes are recorded in
            ``sequence.metadata["outcomes"]`` as a collection of the outcome names that person experienced
            (an empty collection means none). A sequence without that entry has an unknown status and is
            left out of the analysis of every outcome, never counted as negative.

    Returns:
        Dictionary containing comprehensive life course analysis: ``n_sequences``, ``n_sequences_with_events``,
        ``mean_duration`` (mean over the sequences that have events, in timestamp units), ``duration_stats``,
        ``event_frequency``, ``complexity_stats`` and ``trajectory_patterns``. When no sequence has any events
        there is nothing to summarise: ``mean_duration``, ``duration_stats`` and ``complexity_stats`` are None
        and ``event_frequency`` is empty (``{"error": ...}`` is returned for an empty list of sequences).
        With ``outcomes`` it also holds ``outcome_analysis``, one entry per requested outcome:

        - ``status``: ``"ok"``, or ``"no_outcome_data"`` when no sequence records its outcomes (then
          ``prevalence`` is None and nothing else is computed)
        - ``n_sequences``, ``n_with_outcome_data``, ``n_missing_outcome_data``, ``n_with_outcome`` and
          ``n_without_outcome``: the counts behind the analysis (the last two among sequences with data)
        - ``prevalence``: ``n_with_outcome / n_with_outcome_data``
        - ``trajectory_comparison``: for ``n_events``, ``duration`` (in timestamp units) and ``complexity``
          (distinct event types / events), the mean among sequences with and without the outcome, their
          difference and the two-sided Mann-Whitney U p-value; empty unless both groups exist, and a
          feature is None when a group has no sequence that defines it (duration and complexity need events)
        - ``event_type_association``: for every event type, how many sequences have it, how many of those
          have the outcome, the sample odds ratio (None when a zero cell makes it infinite or undefined) and
          the two-sided Fisher exact p-value of the 2x2 table (event type present or not, outcome or not),
          ordered by p-value. The p-values are not adjusted for the number of event types tested.

    Raises:
        errors.ValidationError: If a sequence's ``metadata["outcomes"]`` is not a collection of names
    """
    if not sequences:
        return {"error": "No sequences provided"}

    # Only sequences with events have a duration, a complexity or event frequencies; with none there is nothing
    # to summarise, and the summaries below say so explicitly (None / empty) instead of being left out.
    durations = [seq.duration for seq in sequences if seq.events]
    results: Dict[str, Any] = {
        "n_sequences": len(sequences),
        "n_sequences_with_events": len(durations),
        "mean_duration": statistics.mean(durations) if durations else None,
        "components": {},
        "duration_stats": None,
        "event_frequency": {},
        "complexity_stats": None,
    }

    # Basic statistics
    if durations:
        results["duration_stats"] = {
            "mean": statistics.mean(durations),
            "median": statistics.median(durations),
            "min": min(durations),
            "max": max(durations),
            "std": statistics.stdev(durations) if len(durations) > 1 else 0,
        }

    # Event frequency analysis
    all_events: list[Any] = []
    for seq in sequences:
        all_events.extend([event.event_type for event in seq.events])

    if all_events:
        event_counts: defaultdict[str, int] = defaultdict(int)
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

    # Outcome analysis: prevalence and trajectory associations among sequences with recorded outcomes
    if outcomes:
        recorded_outcomes = [_recorded_outcomes(seq) for seq in sequences]
        results["outcome_analysis"] = {
            outcome: _analyze_outcome(sequences, recorded_outcomes, outcome) for outcome in outcomes
        }

    # Trajectory patterns
    results["trajectory_patterns"] = identify_trajectory_patterns(sequences)

    return results
