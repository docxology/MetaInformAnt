# Life Events Models

Embedding models for life event sequences: learning vector representations from event types, domains, and temporal patterns for downstream prediction and clustering.

## Contents

| File | Purpose |
|------|---------|
| `embeddings.py` | Event embedding learning: type-based, domain-based, and temporal embeddings |

## Key Functions

| Function | Description |
|----------|-------------|
| `learn_event_embeddings()` | Train event embeddings from co-occurrence in event sequences |
| `biological_embedding()` | Create embeddings using biological event structure (type, domain, temporal) |
| `domain_specific_embeddings()` | Generate embeddings specialized to a particular event domain |

## Internal Helpers

| Function | Description |
|----------|-------------|
| `_create_event_type_embeddings()` | Build embeddings based on event type co-occurrence |
| `_create_domain_embeddings()` | Build embeddings that capture domain-level relationships |
| `_create_temporal_embeddings()` | Build embeddings encoding temporal proximity patterns |

## Usage

```python
from metainformant.life_events.models.embeddings import learn_event_embeddings, biological_embedding

embeddings = learn_event_embeddings(sequences, dim=64, window=5)
bio_emb = biological_embedding(sequences, embedding_type="event_type")
```
