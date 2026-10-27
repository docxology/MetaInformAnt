# Life Events Models

Models for life event sequences: vector representations of events, outcome and survival predictors, and next-event models.

## Contents

| File | Purpose |
|------|---------|
| `embeddings.py` | Event embedding learning: word2vec, type-based, domain-based, and temporal embeddings |
| `predictor.py` | `EventSequencePredictor` (`model_type` `"embedding"` or `"simple"`) and `EnsemblePredictor` |
| `statistical_models.py` | `MultiTaskPredictor` (ridge and logistic task heads) and `SurvivalPredictor` (Cox, exponential, Kaplan-Meier), both NumPy |
| `sequence_models.py` | `LSTMSequenceModel` and `GRUSequenceModel`: next-event models that require PyTorch |

## Key Functions

| Function | Description |
|----------|-------------|
| `learn_event_embeddings()` | Train word2vec (skip-gram or CBOW, negative sampling) event embeddings from co-occurrence within a window; options `method`, `min_count`, `negative`, `learning_rate`, `combine`; returns each event's input plus output vector by default (`combine="input"`: input vector only); seeded via `random_state` |
| `biological_embedding()` | Create embeddings using biological event structure (type, domain, temporal) |
| `domain_specific_embeddings()` | Train one `learn_event_embeddings` model per life domain; returns a dictionary keyed by domain |

## Internal Helpers

| Function | Description |
|----------|-------------|
| `_create_event_type_embeddings()` | Build deterministic bit-pattern embeddings from a stable hash of each event type name |
| `_create_domain_embeddings()` | Build embeddings that capture domain-level relationships |
| `_create_temporal_embeddings()` | Build embeddings encoding temporal proximity patterns |

## Usage

```python
from metainformant.life_events.models.embeddings import learn_event_embeddings, biological_embedding

embeddings = learn_event_embeddings(sequences, embedding_dim=64, window_size=5, epochs=10, method="skipgram", random_state=0)
bio_emb = biological_embedding(sequences, embedding_type="event_type")
```
