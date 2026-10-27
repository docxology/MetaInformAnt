# Life Events Analysis

Interpretability and explainability methods for life event sequence models: attention weights, feature attribution, temporal pattern extraction, and Shapley value computation.

## Contents

| File | Purpose |
|------|---------|
| `interpretability.py` | Model interpretability: attention, importance, attribution, Shapley values |

## Key Functions

| Function | Description |
|----------|-------------|
| `attention_weights()` | Attention weights of a model that exposes them; for a scikit-learn model without attention, its `feature_importances_` or `coef_` |
| `event_importance()` | Score events by frequency (default) or, for a fitted predictor with `method="permutation"`, by the total absolute change in its output when the event is removed |
| `feature_attribution()` | Attribute a fitted predictor's outputs to events by ablation (each event removed in turn) or, with `use_shap=True` and `shap` installed, by SHAP |
| `temporal_patterns()` | Discover recurring temporal patterns within event sequences |
| `learn_event_embeddings()` | Skip-gram word2vec embeddings from event co-occurrence: a wrapper over `models.embeddings.learn_event_embeddings` with this module's defaults (`embedding_dim=50`, `window_size=2`, `epochs=10`) that returns `{}` instead of raising when no event reaches `min_count` |
| `shapley_values()` | Shapley values of a sequence's events for a predictor's output: exact up to 10 events, seeded permutation sampling above; the values sum to `v(all) - v(empty sequence)`. Absent events are removed, not replaced by background values, so `background_sequences` is deprecated and ignored: a non-empty value emits a `DeprecationWarning` |

## Usage

```python
from metainformant.life_events.analysis.interpretability import (
    attention_weights, feature_attribution, shapley_values
)

# token_sequences: lists of "domain:event_type" tokens; predictor: a fitted EventSequencePredictor;
# embeddings: the dictionary returned by learn_event_embeddings
weights = attention_weights(predictor, token_sequences)
attributions = feature_attribution(predictor, token_sequences, embeddings)
shap = shapley_values(predictor, token_sequences[0], n_samples=100)
```

`shapley_values()` removes absent events from the sequence rather than replacing them with values drawn from a background dataset. Its `background_sequences` argument is therefore deprecated and ignored: it is kept so that existing positional calls keep working, and passing a non-empty value emits a `DeprecationWarning`.
