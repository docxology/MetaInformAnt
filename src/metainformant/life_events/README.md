# Life Events Module

The `life_events` module provides comprehensive tools for analyzing human life courses as temporal event sequences, enabling prediction of life outcomes from event data using NLP-inspired methods.

## Overview

This module treats life courses as sequences of events (similar to how NLP treats text as sequences of words), allowing the application of sequence modeling techniques to predict outcomes like early mortality, personality traits, and health outcomes.

## Key Features

- **Event Sequence Representation**: Structured data types for temporal event sequences
- **Event Embeddings**: Learn dense vector representations of events using Word2Vec-style methods
- **Sequence Embeddings**: Aggregate event sequences into fixed-size vectors
- **Domain-Specific Embeddings**: Separate embedding spaces for different life domains
- **Temporal Analysis**: Filtering and querying by time periods
- **Integration**: Seamless integration with ML, phenotype, and visualization modules

## Core Components

### Event Data Structures (`events.py`)

#### Event
Individual life event with temporal and domain information:

```python
from metainformant.life_events import Event
from datetime import datetime

event = Event(
    event_type="diagnosis",
    timestamp=datetime(2020, 1, 15),
    domain="health",
    attributes={"condition": "diabetes", "severity": "moderate"}
)
```

#### EventSequence
Container for temporal event sequence of a single person:

```python
from metainformant.life_events import EventSequence, Event

sequence = EventSequence(
    person_id="person_001",
    events=[
        Event("degree", datetime(2010, 6, 1), "education", {"degree": "BS"}),
        Event("job_change", datetime(2015, 3, 1), "occupation", {"company": "TechCorp"}),
        Event("diagnosis", datetime(2020, 1, 15), "health", {"condition": "diabetes"}),
    ]
)

# Filter by domain
health_events = sequence.filter_by_domain("health")

# Filter by time
recent_events = sequence.filter_by_time(
    start_time=datetime(2015, 1, 1),
    end_time=datetime(2025, 1, 1)
)
```

#### EventDatabase
Collection of multiple event sequences:

```python
from metainformant.life_events import EventDatabase

database = EventDatabase(sequences=[sequence1, sequence2, sequence3])
stats = database.get_statistics()
```

### Event Embeddings (`embeddings.py`)

#### Learn Event Embeddings

```python
from metainformant.life_events import learn_event_embeddings

# Convert events to token format: "domain:event_type"
sequences = [
    ["health:diagnosis", "occupation:job_change", "income:raise"],
    ["education:degree", "occupation:job_change", "address:move"],
]

embeddings = learn_event_embeddings(
    sequences,
    embedding_dim=100,
    window_size=5,
    method="skipgram",
    epochs=10
)

# Access embedding for an event
diagnosis_embedding = embeddings["health:diagnosis"]
```

#### Sequence Embeddings

```python
from metainformant.life_events import sequence_embeddings

# Aggregate sequences into fixed-size vectors
seq_embeddings = sequence_embeddings(
    sequences,
    event_embeddings=embeddings,
    method="mean",
    temporal_weighting=True
)
```

#### Domain-Specific Embeddings

```python
from metainformant.life_events import domain_specific_embeddings

sequences = [
    ["diagnosis", "job_change", "raise"],
    ["degree", "job_change", "move"],
]
domains = [
    ["health", "occupation", "income"],
    ["education", "occupation", "address"],
]

domain_embeddings = domain_specific_embeddings(
    sequences,
    domains,
    embedding_dim=100
)
```

## Integration with Other Modules

### With ML Module

```python
from metainformant.life_events import learn_event_embeddings, sequence_embeddings
from metainformant.ml import BiologicalClassifier

# Learn embeddings
embeddings = learn_event_embeddings(sequences)

# Convert sequences to vectors
X = sequence_embeddings(sequences, embeddings)

# Train classifier
classifier = BiologicalClassifier(algorithm="random_forest")
classifier.fit(X, y)
```

### With Visualization Module

```python
from metainformant.life_events import learn_event_embeddings
from metainformant.ml.dimensionality import biological_embedding
from metainformant.visualization import plots

# Learn embeddings
embeddings = learn_event_embeddings(sequences)

# Reduce to 2D for visualization
embedding_matrix = np.array(list(embeddings.values()))
reduced = biological_embedding(embedding_matrix, method="umap", n_components=2)

# Visualize
plots.scatterplot(reduced["embedding"][:, 0], reduced["embedding"][:, 1])
```

## Data Format

### Event Sequence JSON Format

```json
{
  "person_id": "person_001",
  "events": [
    {
      "event_type": "diagnosis",
      "timestamp": "2020-01-15T00:00:00",
      "domain": "health",
      "attributes": {
        "condition": "diabetes",
        "severity": "moderate"
      }
    }
  ],
  "metadata": {}
}
```

## Usage Patterns

### Loading Event Data

```python
from metainformant.life_events import EventSequence
from metainformant.core import io

# Load from JSON
data = io.load_json("events.json")
sequence = EventSequence.from_dict(data)
```

### Analyzing Event Patterns

```python
# Get unique event types
event_types = sequence.get_event_types()

# Get domains covered
domains = sequence.get_domains()

# Filter by domain
health_sequence = sequence.filter_by_domain("health")
```

### Learning and Using Embeddings

```python
# Convert sequences to token lists
sequences_tokens = [
    [f"{e.domain}:{e.event_type}" for e in seq.events]
    for seq in database.sequences
]

# Learn embeddings
embeddings = learn_event_embeddings(sequences_tokens)

# Get sequence-level embeddings
seq_embeddings = sequence_embeddings(sequences_tokens, embeddings)
```

## Performance Considerations

- Embedding learning scales with vocabulary size and number of sequences
- Sequence embeddings use efficient NumPy operations
- Filtering operations are O(n) where n is number of events
- All operations use real implementations (no mocks)

## Testing

Comprehensive tests cover:
- Event data structure validation
- Embedding learning algorithms
- Sequence aggregation methods
- Integration with other modules

Run tests with:
```bash
pytest tests/test_life_events_*.py
```

## Dependencies

- `numpy`: Core numerical operations
- `pandas`: Data structure support
- Standard library: `datetime`, `collections`

Optional dependencies for advanced features:
- Deep learning frameworks (for future transformer models)
- SHAP (for model interpretation)

