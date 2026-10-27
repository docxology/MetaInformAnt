# Life Events Analysis

Interpretability and explainability methods for life event sequence models: attention weights, feature attribution, temporal pattern extraction, and Shapley value computation.

## Contents

| File | Purpose |
|------|---------|
| `interpretability.py` | Model interpretability: attention, importance, attribution, Shapley values |

## Key Functions

| Function | Description |
|----------|-------------|
| `attention_weights()` | Extract and normalize attention weights from sequence models |
| `event_importance()` | Rank events by their predictive importance for outcomes |
| `feature_attribution()` | Attribute model predictions to input features via gradient methods |
| `temporal_patterns()` | Discover recurring temporal patterns within event sequences |
| `learn_event_embeddings()` | Skip-gram word2vec embeddings from event co-occurrence: a wrapper over `models.embeddings.learn_event_embeddings` with this module's defaults (`embedding_dim=50`, `window_size=2`, `epochs=10`) that returns `{}` instead of raising when no event reaches `min_count` |
| `shapley_values()` | Shapley values of a sequence's events for a predictor's output: exact up to 10 events, seeded permutation sampling above; the values sum to `v(all) - v(empty sequence)` |

## Usage

```python
from metainformant.life_events.analysis.interpretability import (
    attention_weights, feature_attribution, shapley_values
)

weights = attention_weights(model, sequence)
attributions = feature_attribution(model, sequence, target="outcome")
shap = shapley_values(model, sequence, n_samples=100)
```
