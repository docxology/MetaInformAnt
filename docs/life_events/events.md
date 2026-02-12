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
| `event_type` | `str` | Category identifier (e.g., `"graduation"`, `"diagnosis"`) |
| `timestamp` | `datetime` | When the event occurred |
| `domain` | `str` | Life domain: education, career, health, family, finance |
| `attributes` | `Dict[str, Any]` | Arbitrary key-value metadata |

### EventSequence

```python
from metainformant.life_events.core.events import EventSequence

seq = EventSequence(
    person_id="P001",
    events=[event1, event2, event3],
    metadata={"cohort": "2024", "region": "west"}
)
```

**Key Methods:**

| Method | Signature | Returns |
|--------|-----------|---------|
| `add_event` | `(event: Event) -> None` | Adds event to sequence |
| `filter_by_domain` | `(domain: str) -> EventSequence` | New sequence with matching domain |
| `filter_by_time_range` | `(start: datetime, end: datetime) -> EventSequence` | New sequence within range |
| `to_dataframe` | `() -> pd.DataFrame` | Tabular representation |
| `get_event_types` | `() -> List[str]` | Unique event types |
| `get_domains` | `() -> List[str]` | Unique domains |

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
| `get_statistics` | `() -> Dict[str, Any]` | Summary stats (counts, domains, types) |
| `save_json` | `(path: str) -> None` | Persist to JSON |
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
db = EventDatabase(sequences=[seq1, seq2, seq3])
stats = db.get_statistics()
print(f"Total sequences: {stats['n_sequences']}")
print(f"Domains: {stats['domains']}")

# Persist and reload
db.save_json("output/life_events/cohort_db.json")
restored = EventDatabase.load_json("output/life_events/cohort_db.json")
```

## Utility Functions

The `core.utils` module provides generation and validation helpers:

| Function | Signature | Purpose |
|----------|-----------|---------|
| `generate_cohort_sequences` | `(n_sequences, avg_events_per_sequence, ...)` | Generate synthetic cohort data |
| `generate_synthetic_life_events` | `(n_sequences, ...)` | Random event sequences |
| `generate_realistic_life_events` | `(n_events, start_age, end_age, seed)` | Age-appropriate events |
| `validate_sequence` | `(sequence) -> bool` | Check sequence validity |
| `get_event_statistics` | `(sequences) -> Dict` | Aggregate statistics |
| `convert_sequences_to_tokens` | `(sequences) -> List[List[str]]` | Tokenize for ML models |

## Configuration

Life events configuration uses the `LE_` environment prefix:

```yaml
work_dir: output/life_events
embedding_dim: 100
model_type: embedding
learning_rate: 0.001
batch_size: 32
epochs: 100
random_seed: 42
```

Load with `load_life_events_config(config_file, prefix="LE")`.

## See Also

- [Models](models.md) -- Prediction models using event sequences
- [Survival](survival.md) -- Time-to-event analysis
- [Visualization](visualization.md) -- Plotting event data
