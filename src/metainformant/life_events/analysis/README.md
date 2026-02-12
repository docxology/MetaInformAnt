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
| `learn_event_embeddings()` | Learn vector embeddings from event co-occurrence patterns |
| `shapley_values()` | Compute Shapley values for event contribution to predictions |

## Usage

```python
from metainformant.life_events.analysis.interpretability import (
    attention_weights, feature_attribution, shapley_values
)

weights = attention_weights(model, sequence)
attributions = feature_attribution(model, sequence, target="outcome")
shap = shapley_values(model, sequence, n_samples=100)
```
