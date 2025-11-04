# Life Events Module

The `life_events` module provides comprehensive tools for analyzing human life courses as temporal event sequences, enabling prediction of life outcomes from event data using NLP-inspired methods.

## Overview

This module treats life courses as sequences of events (similar to how NLP treats text as sequences of words), allowing the application of sequence modeling techniques to predict outcomes like early mortality, personality traits, and health outcomes. The implementation is inspired by research that represents human lives using NLP techniques to examine the evolution and predictability of life outcomes.

## Key Features

- **Event Sequence Representation**: Structured data types for temporal event sequences with domain categorization
- **Event Embeddings**: Learn dense vector representations of events using Word2Vec-style methods (Skip-gram, CBOW)
- **Sequence Embeddings**: Aggregate event sequences into fixed-size vectors for ML models
- **Domain-Specific Embeddings**: Separate embedding spaces for different life domains (health, education, occupation, etc.)
- **Temporal Analysis**: Filtering and querying by time periods
- **Prediction Models**: Classification and regression models for predicting life outcomes
- **Visualization**: Timeline plots, embedding visualizations, attention heatmaps, importance plots
- **Model Interpretation**: Event importance, temporal patterns, feature attribution
- **Workflow Integration**: End-to-end pipelines for complete analysis
- **Utility Functions**: Data loading, validation, and statistics computation
- **Integration**: Seamless integration with ML, phenotype, visualization, and multi-omics modules

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
from metainformant.life_events import EventSequence, load_sequences_from_json
from metainformant.core import io

# Load from JSON (convenience function)
sequences = load_sequences_from_json("events.json")

# Or load manually
data = io.load_json("events.json")
sequence = EventSequence.from_dict(data)
```

### Utility Functions

```python
from metainformant.life_events import (
    validate_sequence,
    convert_sequences_to_tokens,
    get_event_statistics
)

# Validate sequence
is_valid, errors = validate_sequence(sequence)
if not is_valid:
    print("Validation errors:", errors)

# Convert to tokens for ML
tokens = convert_sequences_to_tokens(sequences)

# Get statistics
stats = get_event_statistics(sequences)
print(f"Total events: {stats['total_events']}")
print(f"Domains: {stats['domains']}")
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

### Prediction Models (`models.py`)

#### EventSequencePredictor

Train models to predict outcomes from event sequences:

```python
from metainformant.life_events import EventSequencePredictor, learn_event_embeddings
import numpy as np

# Prepare sequences and outcomes
sequences = [
    ["health:diagnosis", "occupation:job_change"],
    ["education:degree", "occupation:job_change"],
]
y = np.array([0, 1])  # Binary classification labels

# Train predictor
predictor = EventSequencePredictor(
    model_type="embedding",
    task_type="classification",
    embedding_dim=100,
    random_state=42
)

# Learn embeddings first
embeddings = learn_event_embeddings(sequences, embedding_dim=100, random_state=42)

# Fit model
predictor.fit(sequences, y, event_embeddings=embeddings)

# Make predictions
predictions = predictor.predict(sequences)
probabilities = predictor.predict_proba(sequences)
```

### Visualization (`visualization.py`)

#### Plot Event Timeline

Visualize individual life course as timeline:

```python
from metainformant.life_events import EventSequence, Event, plot_event_timeline
from datetime import datetime

sequence = EventSequence(
    person_id="person_001",
    events=[
        Event("degree", datetime(2010, 6, 1), "education"),
        Event("job_change", datetime(2015, 3, 1), "occupation"),
        Event("diagnosis", datetime(2020, 1, 15), "health"),
    ]
)

fig = plot_event_timeline(sequence, output_path="output/timeline.png")
```

#### Plot Event Embeddings

Visualize event embeddings in 2D or 3D:

```python
from metainformant.life_events import learn_event_embeddings, plot_event_embeddings

sequences = [
    ["health:diagnosis", "occupation:job_change"],
    ["education:degree", "occupation:job_change"],
]

embeddings = learn_event_embeddings(sequences, embedding_dim=100)

# Plot with UMAP reduction
fig = plot_event_embeddings(
    embeddings,
    method="umap",
    n_components=2,
    output_path="output/embeddings.png"
)
```

#### Plot Prediction Importance

Visualize which events are most important for predictions:

```python
from metainformant.life_events import (
    EventSequencePredictor, learn_event_embeddings, event_importance,
    plot_prediction_importance
)
import numpy as np

# Train model and get importance
sequences = [["health:diagnosis", "occupation:job_change"]]
y = np.array([0])
predictor = EventSequencePredictor(random_state=42)
embeddings = learn_event_embeddings(sequences, random_state=42)
predictor.fit(sequences, y, event_embeddings=embeddings)

importance = event_importance(predictor, sequences, embeddings)

# Plot importance
fig = plot_prediction_importance(
    importance,
    top_n=10,
    output_path="output/importance.png"
)
```

### Model Interpretation (`interpretability.py`)

#### Event Importance

Identify which events contribute most to predictions:

```python
from metainformant.life_events import (
    EventSequencePredictor, learn_event_embeddings, event_importance
)
import numpy as np

sequences = [
    ["health:diagnosis", "occupation:job_change", "income:raise"],
    ["education:degree", "occupation:job_change", "address:move"],
]
y = np.array([0, 1])

predictor = EventSequencePredictor(random_state=42)
embeddings = learn_event_embeddings(sequences, random_state=42)
predictor.fit(sequences, y, event_embeddings=embeddings)

# Compute importance
importance = event_importance(
    predictor,
    sequences,
    embeddings,
    method="permutation"
)

# Events with highest importance
sorted_events = sorted(importance.items(), key=lambda x: x[1], reverse=True)
print("Top important events:", sorted_events[:5])
```

#### Temporal Patterns

Analyze which time periods are critical for predictions:

```python
from metainformant.life_events import temporal_patterns
import numpy as np

sequences = [
    ["health:diagnosis", "occupation:job_change"],
    ["education:degree", "occupation:job_change"],
]
predictions = np.array([0.8, 0.3])

patterns = temporal_patterns(sequences, predictions)
print("Position importance:", patterns["position_importance"])
```

#### Feature Attribution

Compute SHAP-style feature attribution:

```python
from metainformant.life_events import feature_attribution

attribution = feature_attribution(
    predictor,
    sequences,
    embeddings,
    use_shap=False  # Uses permutation method if SHAP not available
)
```

### Workflow Examples (`workflow.py`)

#### Complete Life Course Analysis

End-to-end analysis from sequences to predictions:

```python
from metainformant.life_events import EventSequence, Event, analyze_life_course
from datetime import datetime
import numpy as np

# Create sequences
sequences = [
    EventSequence(
        "person_001",
        [Event("degree", datetime(2010, 6, 1), "education")]
    ),
    EventSequence(
        "person_002",
        [Event("diagnosis", datetime(2020, 1, 15), "health")]
    ),
]

# Optional outcomes for training
outcomes = np.array([0, 1])

# Run complete analysis
results = analyze_life_course(
    sequences,
    outcomes=outcomes,
    output_dir="output/life_events"
)

print(f"Trained model: {results['model_type']}")
print(f"Accuracy: {results.get('accuracy', 'N/A')}")
```

#### Compare Populations

Compare event patterns between groups:

```python
from metainformant.life_events import compare_populations

group1 = [EventSequence("p1", [Event("degree", datetime(2010, 1, 1), "education")])]
group2 = [EventSequence("p2", [Event("job_change", datetime(2015, 1, 1), "occupation")])]

comparison = compare_populations(group1, group2)

print("Common event types:", comparison["comparison"]["common_event_types"])
print("Unique to group 1:", comparison["comparison"]["unique_to_group1"])
```

#### Intervention Analysis

Analyze effects of interventions:

```python
from metainformant.life_events import intervention_analysis
from datetime import datetime
import numpy as np

sequences = [
    EventSequence(
        "person_001",
        [
            Event("degree", datetime(2010, 1, 1), "education"),
            Event("job_change", datetime(2015, 1, 1), "occupation"),
            Event("diagnosis", datetime(2020, 1, 1), "health"),
        ]
    )
]

intervention_time = datetime(2017, 1, 1).timestamp()
pre_outcomes = np.array([0.5])
post_outcomes = np.array([0.8])

results = intervention_analysis(
    sequences,
    intervention_time=intervention_time,
    pre_intervention_outcomes=pre_outcomes,
    post_intervention_outcomes=post_outcomes
)

print("Outcome change:", results["outcome_change"]["mean"])
```

## CLI Usage

The module provides CLI commands through the main METAINFORMANT CLI:

### Learn Event Embeddings

```bash
metainformant life-events embed \
    --input data/event_sequences.json \
    --output output/life_events/embeddings \
    --embedding-dim 100 \
    --window-size 5 \
    --epochs 10
```

### Predict Outcomes

```bash
metainformant life-events predict \
    --events data/test_sequences.json \
    --model output/life_events/model.json \
    --outcome mortality \
    --output output/life_events/predictions
```

### Interpret Predictions

```bash
metainformant life-events interpret \
    --model output/life_events/model.json \
    --sequences data/test_sequences.json \
    --output output/life_events/interpretation
```

## Configuration

Workflow functions support configuration files for parameter control:

### Example Configuration (YAML)

```yaml
embedding_dim: 100
window_size: 5
epochs: 10
model_type: embedding
task_type: classification
random_state: 42
```

### Using Configuration

```python
from metainformant.life_events import analyze_life_course

results = analyze_life_course(
    sequences,
    outcomes=outcomes,
    config_path="config/life_events.yaml",
    output_dir="output/life_events"
)
```

## Integration Examples

### With Phenotype Module

Extract phenotypes from event sequences:

```python
from metainformant.life_events import EventSequence, Event
from metainformant.phenotype import extract_phenotypes_from_events
from datetime import datetime

sequence = EventSequence(
    "person_001",
    [
        Event("diabetes", datetime(2020, 1, 1), "health", {"severity": "moderate"}),
        Event("bachelors", datetime(2010, 6, 1), "education", {"degree": "BS"}),
    ]
)

phenotypes = extract_phenotypes_from_events(sequence)
print("Health events:", phenotypes["health_events"])
print("Education achievements:", phenotypes["education_achievements"])
```

### With Networks Module

Create event co-occurrence networks:

```python
from metainformant.life_events import EventDatabase
from metainformant.networks import graph

database = EventDatabase(sequences=[seq1, seq2, seq3])

# Build co-occurrence network
# (Implementation would extract event pairs and create network)
```

### With Multi-omics Module

Integrate events as one data layer:

```python
from metainformant.life_events import EventDatabase, sequence_embeddings, learn_event_embeddings
from metainformant.multiomics import MultiOmicsData, joint_pca
import pandas as pd

# Convert events to feature matrix
embeddings = learn_event_embeddings(sequences_tokens)
event_features = sequence_embeddings(sequences_tokens, embeddings)

# Create multi-omics data
omics_data = MultiOmicsData(
    genomics=genomic_df,
    transcriptomics=transcriptomic_df,
    # Add events as custom layer
    metadata=pd.DataFrame(event_features)
)

# Joint analysis
joint_embeddings, loadings, variance = joint_pca(omics_data, n_components=10)
```

## Performance Considerations

- Embedding learning scales with vocabulary size and number of sequences
- Sequence embeddings use efficient NumPy operations
- Filtering operations are O(n) where n is number of events
- All operations use real implementations (no mocks)
- For large datasets, consider using domain-specific embeddings to reduce vocabulary size

## Troubleshooting

### Empty Sequences Error

If you get "sequences list cannot be empty" error:
- Ensure your input data contains at least one sequence
- Check that sequences are properly loaded from files
- Verify that filtering operations haven't removed all events

### Embedding Dimension Issues

If embeddings are too large or too small:
- Adjust `embedding_dim` parameter (typical values: 50-200)
- Consider vocabulary size when choosing dimension
- Use domain-specific embeddings for large vocabularies

### Model Not Fitted Error

If you get "Model must be fitted" error:
- Call `fit()` before `predict()` or `predict_proba()`
- Ensure training data is not empty
- Check that outcomes match sequence length

### Missing Optional Dependencies

If visualization functions fail:
- Install matplotlib: `pip install matplotlib`
- Visualization functions are optional and will gracefully degrade
- Check that output directory is writable

## Testing

Comprehensive tests cover:
- Event data structure validation
- Embedding learning algorithms
- Sequence aggregation methods
- Prediction models (classification and regression)
- Workflow functions
- Visualization functions
- Interpretability tools
- Integration with other modules

Run tests with:
```bash
pytest tests/test_life_events_*.py
```

Run specific test suites:
```bash
pytest tests/test_life_events_events.py
pytest tests/test_life_events_embeddings.py
pytest tests/test_life_events_models.py
pytest tests/test_life_events_workflow.py
pytest tests/test_life_events_visualization.py
pytest tests/test_life_events_interpretability.py
pytest tests/test_life_events_utils.py
pytest tests/test_life_events_integration.py
```

## Dependencies

### Required
- `numpy>=1.21.0`: Core numerical operations
- `pandas>=1.3.0`: Data structure support
- Standard library: `datetime`, `collections`, `pathlib`

### Optional
- `matplotlib>=3.5.0`: For visualization functions
- `scikit-learn>=1.1.0`: For ML model evaluation metrics
- `scipy>=1.9.0`: For statistical tests in intervention analysis
- `shap>=0.40.0`: For advanced feature attribution (not yet fully integrated)
- `torch` or `tensorflow`: For future transformer/LSTM models (not yet fully implemented)

All optional dependencies are handled defensively - the module will work without them, but some features may be limited.

