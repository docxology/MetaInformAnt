# Life Events: Prediction Models

Machine learning models for predicting outcomes from life event sequences.

## Concept Overview

The models submodule provides several architectures for learning from event sequences:

- **EventSequencePredictor** -- General-purpose predictor with embedding-based features
- **LSTMSequenceModel** -- Long Short-Term Memory recurrent model
- **GRUSequenceModel** -- Gated Recurrent Unit model (lighter than LSTM)
- **MultiTaskPredictor** -- Simultaneous prediction of multiple outcomes
- **EnsemblePredictor** -- Weighted combination of multiple models
- **SurvivalPredictor** -- Time-to-event prediction with Cox or embedding methods

All models share a consistent `fit` / `predict` interface and support serialization.

## Model Reference

### EventSequencePredictor

The primary entry point for event-based prediction.

```python
from metainformant.life_events.models.predictor import EventSequencePredictor

predictor = EventSequencePredictor(
    model_type="embedding",   # "embedding", "frequency", or "sequential"
    embedding_dim=100,
    random_seed=42
)
predictor.fit(sequences, outcomes, task="classification")
predictions = predictor.predict(test_sequences)
probabilities = predictor.predict_proba(test_sequences)
```

**Methods:** `fit(sequences, outcomes, task)`, `predict(sequences)`, `predict_proba(sequences)`, `save_model(path)`, `load_model(path)`

### LSTMSequenceModel / GRUSequenceModel

Recurrent models for capturing dependencies in event sequences. GRU has fewer
parameters (faster training, similar performance).

```python
from metainformant.life_events.models.sequence_models import LSTMSequenceModel, GRUSequenceModel

model = LSTMSequenceModel(
    embedding_dim=100, hidden_dim=64, num_layers=1,
    dropout=0.1, epochs=50, batch_size=32, learning_rate=0.001
)
model.fit(sequences, targets)
predictions = model.predict(test_sequences)
```

### MultiTaskPredictor

Predicts multiple outcomes simultaneously, sharing learned representations.

```python
from metainformant.life_events.models.statistical_models import MultiTaskPredictor

predictor = MultiTaskPredictor(
    task_types={"employment": "classification", "income": "regression"},
    embedding_dim=100, hidden_dim=64
)
predictor.fit(sequences, outcomes_dict)
results = predictor.predict(test_sequences)
```

### EnsemblePredictor

Combines multiple predictors with optional weighting.

```python
from metainformant.life_events.models.predictor import EnsemblePredictor

ensemble = EnsemblePredictor(
    models=[predictor1, predictor2, predictor3],
    weights=[0.5, 0.3, 0.2]   # Optional; uniform if None
)
```

### SurvivalPredictor

Time-to-event prediction combining embeddings with survival analysis.

```python
from metainformant.life_events.models.statistical_models import SurvivalPredictor

surv = SurvivalPredictor(method="cox", embedding_dim=50)
surv.fit(sequences, event_times, event_occurred)
risk_scores = surv.predict(test_sequences)
survival_curves = surv.predict_survival_function(test_sequences, times=[1, 5, 10])
```

## Embedding Functions

Standalone functions for learning event vector representations.

| Function | Signature | Purpose |
|----------|-----------|---------|
| `learn_event_embeddings` | `(sequences, embedding_dim=100, window_size=5, min_count=1, sg=1, epochs=5)` | Word2Vec-style event embeddings |
| `biological_embedding` | `(sequences, embedding_type="event_type")` | Domain-aware embeddings (event_type, domain, temporal) |
| `domain_specific_embeddings` | `(sequences, domains, embedding_dim=50)` | Separate embeddings per life domain |

```python
from metainformant.life_events.models.embeddings import learn_event_embeddings

result = learn_event_embeddings(sequences, embedding_dim=100, window_size=5)
# result["embeddings"]   -> np.ndarray (vocab_size, embedding_dim)
# result["vocabulary"]   -> Dict[str, int] event counts
# result["event_to_idx"] -> Dict[str, int] index mapping
```

## Typical Workflow

```python
from metainformant.life_events.core.utils import generate_cohort_sequences
from metainformant.life_events.models.predictor import EventSequencePredictor

sequences = generate_cohort_sequences(n_sequences=500, avg_events_per_sequence=15)
outcomes = [1 if len(s.events) > 10 else 0 for s in sequences]

model = EventSequencePredictor(model_type="embedding", embedding_dim=64)
model.fit(sequences[:400], outcomes[:400], task="classification")
predictions = model.predict(sequences[400:])
model.save_model("output/life_events/trained_model.json")
```

## Configuration

| Field | Default | Env Override |
|-------|---------|-------------|
| `embedding_dim` | `100` | `LE_EMBEDDING_DIM` |
| `model_type` | `"embedding"` | `LE_MODEL_TYPE` |
| `learning_rate` | `0.001` | `LE_LEARNING_RATE` |
| `batch_size` | `32` | `LE_BATCH_SIZE` |
| `epochs` | `100` | `LE_EPOCHS` |

## See Also

- [Events](events.md) -- Input data structures
- [Survival](survival.md) -- Dedicated survival analysis
- [Workflow](workflow.md) -- End-to-end analysis pipelines
