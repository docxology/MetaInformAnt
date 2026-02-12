# Life Course Analysis

Life course phenotype analysis for extracting temporal traits from event
sequences, modelling trajectories, and predicting outcomes.

## Key Concepts

**Event sequences** represent ordered life events (education, career, health,
social, mobility) for individuals. The analysis pipeline extracts phenotypic
traits from these sequences, aggregates them across time windows and
populations, identifies critical developmental periods, and predicts future
outcomes via simple linear extrapolation.

**Trajectory classification** categorises individuals by activity level
(high/moderate/low) and stability (stable/unstable) based on event counts and
transition patterns. Transition probabilities between event types are computed
to reveal common pathways through life stages.

## Data Structures

### Event

```python
from metainformant.phenotype.analysis.life_course import Event

event = Event(
    timestamp=25.0,
    event_type="education_complete",
    description="PhD earned",
    metadata={"level": "doctorate"},
    confidence=0.95,
)
```

| Field        | Type             | Description                          |
|-------------|------------------|--------------------------------------|
| timestamp   | float            | Time of the event                    |
| event_type  | str              | Category of event                    |
| description | str              | Free-text description                |
| metadata    | Dict[str, Any]   | Arbitrary key-value metadata         |
| confidence  | float            | Confidence score in [0.0, 1.0]       |

### EventSequence

```python
from metainformant.phenotype.analysis.life_course import EventSequence

seq = EventSequence(person_id="P001", events=[event])
seq.add_event(another_event)          # inserts in sorted order
seq.duration                          # total timespan
seq.event_count                       # number of events
seq.get_events_in_range(20.0, 30.0)   # filter by time window
seq.get_events_by_type("job_start")   # filter by type
```

## Function Reference

### extract_phenotypes_from_events

```python
def extract_phenotypes_from_events(
    event_sequence: EventSequence,
    trait_mapping: Optional[Dict[str, List[str]]] = None,
) -> Dict[str, Any]
```

Extracts phenotypic traits (education level, career stability, health score,
social connectivity, mobility) from an event sequence. The default trait
mapping covers five domains; provide a custom mapping to override.

### aggregate_temporal_phenotypes

```python
def aggregate_temporal_phenotypes(
    sequences: List[EventSequence],
    time_windows: List[Tuple[float, float]],
    trait_categories: List[str],
) -> Dict[str, Any]
```

Aggregates extracted phenotypes across individuals and time windows. Numeric
traits produce mean/median/min/max/count; categorical traits produce frequency
distributions with most-common values.

### analyze_life_course_trajectories

```python
def analyze_life_course_trajectories(
    sequences: List[EventSequence],
    outcome_measures: Optional[List[str]] = None,
) -> Dict[str, Any]
```

Classifies trajectories by activity level and stability, then computes outcome
scores for career success, health, and social stability. Identifies risk
factors (high health events, job instability).

### identify_critical_periods

```python
def identify_critical_periods(
    sequences: List[EventSequence],
    age_ranges: List[Tuple[float, float]],
) -> Dict[str, Any]
```

Analyses event density in each age range and flags periods with event rates
exceeding one standard deviation above the mean as high-significance clusters.

### predict_life_course_outcomes

```python
def predict_life_course_outcomes(
    sequences: List[EventSequence],
    prediction_horizon: float = 5.0,
) -> Dict[str, Any]
```

Linear extrapolation of career, health, and social trajectories over the
prediction horizon (in years). Returns predictions with simplified 20% margin
confidence intervals.

### identify_trajectory_patterns

```python
def identify_trajectory_patterns(
    sequences: List[EventSequence],
) -> Dict[str, Any]
```

Classifies sequences by activity and diversity (complex_high_activity,
repetitive_high_activity, diverse_moderate, focused_moderate, diverse_low,
simple_low). Computes transition probability matrices between event types.

### analyze_life_course

```python
def analyze_life_course(
    sequences: List[EventSequence],
    outcomes: List[str] | None = None,
) -> Dict[str, Any]
```

Comprehensive analysis combining duration stats, event frequencies, sequence
complexity metrics, outcome analysis, and trajectory pattern identification.

### create_life_course_report

```python
def create_life_course_report(
    sequences: List[EventSequence],
    output_path: Optional[str | Path] = None,
) -> str
```

Generates a formatted text report with individual counts, event type
distributions, trajectory patterns, outcome measures, and risk factors.
Optionally saves to a file.

## Usage Example

```python
from metainformant.phenotype.analysis.life_course import (
    Event, EventSequence, analyze_life_course, create_life_course_report,
)

events = [
    Event(timestamp=18.0, event_type="education_complete", metadata={"level": "bachelor"}),
    Event(timestamp=22.0, event_type="job_start"),
    Event(timestamp=30.0, event_type="health_event"),
]
seq = EventSequence(person_id="P001", events=events)

results = analyze_life_course([seq], outcomes=["career_success"])
report = create_life_course_report([seq], output_path="output/phenotype/report.txt")
```

## Related Modules

- `metainformant.phenotype.visualization` -- plot life course trajectories
- `metainformant.phenotype.behavior` -- behavioural sequence analysis
- `metainformant.phenotype.workflow` -- pipeline orchestration
- `metainformant.life_events` -- event sequence analysis
