# Life Events Core

Core data structures and utilities for life event sequence analysis: event/sequence models, configuration, synthetic generation, and validation.

## Contents

| File | Purpose |
|------|---------|
| `events.py` | Event, EventSequence, and EventDatabase data classes |
| `config.py` | LifeEventsWorkflowConfig dataclass, YAML config loading, defaults |
| `utils.py` | Synthetic event generation, sequence statistics, validation, tokenization |

## Key Classes

| Class | Description |
|-------|-------------|
| `Event` | Single life event with type, timestamp, domain, and metadata |
| `EventSequence` | Ordered collection of Events representing an individual's life course |
| `EventDatabase` | Storage and retrieval interface for multiple EventSequences |
| `LifeEventsWorkflowConfig` | Configuration for life events analysis pipelines |

## Key Functions

| Function | Description |
|----------|-------------|
| `load_life_events_config()` | Load YAML config with environment variable overrides (prefix: LE) |
| `create_default_config()` | Generate default LifeEventsWorkflowConfig for a work directory |
| `generate_realistic_life_events()` | Create synthetic event sequences with domain-aware patterns |
| `get_event_statistics()` | Compute summary statistics across a list of EventSequences |
| `load_sequences_from_json()` | Load EventSequences from JSON file |
| `validate_sequence()` | Validate an EventSequence for temporal ordering and completeness |
| `generate_cohort_sequences()` | Generate synthetic cohort with configurable parameters |

## Usage

```python
from metainformant.life_events.core.events import Event, EventSequence
from metainformant.life_events.core.utils import generate_realistic_life_events

sequences = generate_realistic_life_events(n_individuals=100, n_events=50)
stats = get_event_statistics(sequences)
```
