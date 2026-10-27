# Life Events: Event Data Structures

Core data structures for representing, storing, and querying life event sequences.

## Concept Overview

The life events module models individual life trajectories as timestamped event sequences.
Three primary structures form the foundation:

- **Event** -- A single occurrence with type, timestamp, domain, and arbitrary attributes
- **EventSequence** -- An ordered collection of events for one person
- **EventDatabase** -- A queryable collection of sequences with persistence

Events are domain-tagged (education, career, health, family, finance) enabling
domain-specific filtering and cross-domain analysis.

## Data Structures

### Event

```python
from datetime import datetime

from metainformant.life_events.core.events import Event

event = Event(
    event_type="job_start",
    timestamp=datetime(2024, 3, 15),
    domain="career",
    attributes={"company": "Acme", "role": "engineer"}
)
```

**Fields:**

| Field | Type | Description |
|-------|------|-------------|
| `event_type` | `str` | Category identifier (e.g., `"graduation"`, `"diagnosis"`); must not be empty |
| `timestamp` | `datetime` | When the event occurred; POSIX seconds and ISO-format strings are converted to `datetime` |
| `domain` | `str` | Life domain such as education, career, health, family, finance; must not be empty |
| `attributes` | `Dict[str, Any]` | Arbitrary key-value metadata (default: empty) |
| `person_id` | `Optional[str]` | Optional owner identifier (default: `None`) |
| `description` | `str` | Optional free text (default: empty) |
| `confidence` | `float` | Confidence in the record (default: `1.0`) |

### EventSequence

```python
from metainformant.life_events.core.events import EventSequence

seq = EventSequence(
    person_id="P001",
    events=[event1, event2, event3],
    metadata={"cohort": "2024", "region": "west"}
)
```

The constructor sorts the events by timestamp.

**Key Methods:**

| Method | Signature | Returns |
|--------|-----------|---------|
| `filter_by_domain` | `(domain: str) -> EventSequence` | New sequence with matching domain |
| `filter_by_time_range` | `(start: datetime, end: datetime) -> EventSequence` | New sequence within range (inclusive) |
| `filter_by_time` | `(start, end) -> EventSequence` | Same as `filter_by_time_range`; accepts `datetime` or POSIX seconds |
| `get_events_by_type` | `(event_type: str) -> List[Event]` | Events of one type |
| `to_dataframe` | `() -> pd.DataFrame` | Tabular representation |
| `get_event_types` | `() -> List[str]` | Unique event types |
| `get_domains` | `() -> List[str]` | Unique domains |
| `to_dict` / `from_dict` | `() -> Dict` / `(data) -> EventSequence` | Serialization (`from_dict` is a class method) |

`len(seq)` is the number of events and `seq[i]` the i-th event in time order.

### EventDatabase

```python
from metainformant.life_events.core.events import EventDatabase

db = EventDatabase(
    sequences=[seq1, seq2, seq3],
    metadata={"study": "longitudinal_2024"}
)
```

**Key Methods:**

| Method | Signature | Returns |
|--------|-----------|---------|
| `add_sequence` | `(sequence: EventSequence) -> None` | Adds and indexes sequence |
| `get_sequence` | `(person_id: str) -> Optional[EventSequence]` | Lookup by person ID |
| `filter_by_domain` | `(domain: str) -> EventDatabase` | New DB filtered by domain |
| `filter_by_time_range` | `(start, end) -> EventDatabase` | New DB within time range |
| `get_statistics` | `() -> Dict[str, Any]` | Summary stats (`n_sequences`, `n_events`, `unique_event_types`, `unique_domains`, per-domain and per-type counts, `avg_events_per_sequence`) |
| `save_json` | `(path: str) -> None` | Persist to JSON; the parent directory must exist |
| `load_json` | `(path: str) -> EventDatabase` | Class method to restore from JSON |

## Usage Examples

### Building a Life Trajectory

```python
from datetime import datetime
from metainformant.life_events.core.events import Event, EventSequence, EventDatabase

events = [
    Event("college_start", datetime(2018, 9, 1), "education", {"school": "State U"}),
    Event("internship", datetime(2020, 6, 1), "career", {"company": "TechCo"}),
    Event("graduation", datetime(2022, 5, 15), "education", {"degree": "BS"}),
    Event("job_start", datetime(2022, 7, 1), "career", {"role": "analyst"}),
]

seq = EventSequence(person_id="P001", events=events)
career_events = seq.filter_by_domain("career")
print(f"Career events: {len(career_events.events)}")  # 2
```

### Querying a Database

```python
from pathlib import Path

db = EventDatabase(sequences=[seq])
stats = db.get_statistics()
print(f"Total sequences: {stats['n_sequences']}")
print(f"Domains: {stats['domains']}")

# Persist and reload
Path("output/life_events").mkdir(parents=True, exist_ok=True)
db.save_json("output/life_events/cohort_db.json")
restored = EventDatabase.load_json("output/life_events/cohort_db.json")
```

## Utility Functions

The `core.utils` module provides generation and validation helpers:

| Function | Signature | Purpose |
|----------|-----------|---------|
| `generate_cohort_sequences` | `(n_sequences, avg_events_per_sequence, ..., random_state=None)` | Synthetic cohort data (a list of sequences, or a dictionary of cohorts when `n_cohorts` is set) |
| `generate_synthetic_life_events` | `(n_sequences=100, ..., random_state=None) -> (sequences, outcomes)` | Life-stage sequences with optional outcomes, rare events, missing data, and timestamp noise (see below) |
| `generate_realistic_life_events` | `(n_events, start_age, end_age, seed) -> EventSequence` | One sequence of life-stage-typical event types (at most `n_events` events) |
| `validate_sequence` | `(sequence) -> (bool, List[str])` | Validity flag plus the issues and warnings found |
| `get_event_statistics` | `(sequences) -> Dict` | Aggregate statistics |
| `convert_sequences_to_tokens` | `(sequences) -> List[List[str]]` | `"domain:event_type"` tokens for ML models |

### generate_synthetic_life_events

```python
from metainformant.life_events.core.utils import generate_synthetic_life_events

sequences, outcomes = generate_synthetic_life_events(
    n_sequences=200,
    generate_outcomes=True,
    outcome_relationship="event_based",
    temporal_patterns="clustered",
    rare_event_probability=0.05,
    missing_data_probability=0.1,
    temporal_noise=0.02,
    random_state=42,
)
```

| Argument | Default | Meaning |
|----------|---------|---------|
| `n_sequences` | `100` | Number of sequences |
| `min_events_per_sequence`, `max_events_per_sequence` | `3`, `20` | Range of the number of events drawn per sequence (before missing data removes any) |
| `generate_outcomes` | `False` | Return a 0/1 outcome per sequence; `outcomes` is `None` otherwise |
| `outcome_relationship` | `"random"` | `"random"`: a fair coin flip. `"event_based"`: 1 if the complete record contains a diagnosis, divorce, or financial difficulty. `"temporal"`: 1 if the complete record has more than 12 events |
| `event_types`, `domains` | built-in lists | Event types and domains to draw from |
| `temporal_patterns` | `"realistic"` | `"realistic"`: life-stage chains with plausible offsets. `"random"`: the same events with independent uniform timestamps inside the realistic timeline's span. `"clustered"`: the same events grouped into up to three bursts |
| `rare_event_probability` | `0.0` | Probability, independently per event, that the event is replaced by a rare event (`serious_accident`, `bereavement`, `natural_disaster`, `windfall`, or `business_failure`; see `RARE_EVENT_DOMAINS`) at the same time; the number of events is unchanged |
| `missing_data_probability` | `0.0` | Probability, independently per event, that the event is missing from the observed sequence; sequences can become shorter than `min_events_per_sequence` or empty |
| `temporal_noise` | `0.0` | Standard deviation of Gaussian jitter added to every timestamp, as a fraction of the sequence's time span; events are re-sorted afterwards |
| `random_state` | `None` | Seed |

Outcomes describe the complete record, rare events included, and are computed before events go missing.

`ValueError` is raised for an unknown `outcome_relationship` or `temporal_patterns`, a probability outside [0, 1], or a negative `temporal_noise`. The function has no `**kwargs`: any other keyword argument (for example `transition_probabilities`, `seasonal_patterns`, or `cooccurrence_patterns`) raises `TypeError`.

## Configuration

Life events workflows are configured with `LifeEventsWorkflowConfig`, loaded from YAML, TOML, or JSON:

```yaml
work_dir: output/life_events/work
embedding:
  embedding_dim: 100
  window_size: 5
  epochs: 10
model:
  model_type: embedding
  task_type: classification
  random_state: 42
workflow:
  save_model: true
```

Load with `load_life_events_config(config_file)`; `work_dir` is required. The `LE_` prefix applies `LE_THREADS`, `LE_WORK_DIR`, and `LE_LOG_DIR` to the loaded values. The keys that `analyze_life_course` reads, and the full template, are described in [Workflow](workflow.md#configuration).

## See Also

- [Models](models.md) -- Prediction models using event sequences
- [Survival](survival.md) -- Time-to-event analysis
- [Visualization](visualization.md) -- Plotting event data
