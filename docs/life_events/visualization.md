# Life Events: Visualization

Comprehensive plotting library for life event sequences and analysis results.

## Concept Overview

The visualization submodule provides 18+ specialized plot types for exploring
event data, embeddings, predictions, and temporal patterns. All functions follow
a consistent interface:

- Accept data as the first argument
- Optional `output_path` for saving (displays if `None`)
- Optional `figsize`, `title`, and domain-specific parameters
- Return a matplotlib `Figure` object

Requires `matplotlib` and optionally `seaborn` for enhanced styling.

## Function Reference

All functions are imported from `metainformant.life_events.visualization` submodules (`timeline`, `network`, `statistical`).

### Event Timelines

| Function | Key Parameters | Purpose |
|----------|---------------|---------|
| `plot_event_timeline` | `(sequence, output_path, figsize, show_domains, title)` | Single sequence timeline, color-coded by domain |
| `plot_domain_timeline` | `(sequences, output_path, figsize, domains)` | Multiple sequences aligned by domain |

```python
from metainformant.life_events.visualization.timeline import plot_event_timeline

fig = plot_event_timeline(
    sequence,
    output_path="output/timeline.png",
    figsize=(14, 4),
    show_domains=True,
    title="Life Trajectory: P001"
)
```

### Embedding Visualizations

| Function | Key Parameters | Purpose |
|----------|---------------|---------|
| `plot_event_embeddings` | `(embeddings, output_path, method, figsize, title)` | 2D PCA or t-SNE projection |
| `plot_embedding_clusters` | `(embeddings, clusters, output_path)` | Embeddings with cluster labels |

```python
fig = plot_event_embeddings(
    embeddings,
    output_path="output/embeddings.png",
    method="tsne"   # "pca" or "tsne"
)
```

### Prediction & Importance

| Function | Key Parameters | Purpose |
|----------|---------------|---------|
| `plot_prediction_importance` | `(feature_importance, output_path, top_n)` | Bar chart of feature importance |
| `plot_prediction_accuracy` | `(y_true, y_pred, output_path)` | Predicted vs actual scatter/confusion |
| `plot_attention_heatmap` | `(attention_weights, output_path)` | Attention weight matrix |

### Distribution & Frequency

| Function | Key Parameters | Purpose |
|----------|---------------|---------|
| `plot_domain_distribution` | `(sequences, output_path)` | Event counts per domain |
| `plot_event_frequency_heatmap` | `(sequences, output_path)` | Event type frequencies across sequences |
| `plot_outcome_distribution` | `(outcomes, output_path)` | Outcome value distribution |
| `plot_sequence_length_distribution` | `(sequences, output_path)` | Histogram of sequence lengths |

### Temporal & Network

| Function | Key Parameters | Purpose |
|----------|---------------|---------|
| `plot_temporal_density` | `(sequences, output_path)` | KDE of event occurrence over time |
| `plot_temporal_patterns` | `(sequences, importance_scores, output_path)` | Temporal patterns with importance overlay |
| `plot_event_cooccurrence` | `(sequences, output_path)` | Co-occurrence matrix |
| `plot_transition_network` | `(sequences, output_path, min_transitions)` | Directed event-to-event transition graph |

### Comparison

| Function | Key Parameters | Purpose |
|----------|---------------|---------|
| `plot_population_comparison` | `(group1_sequences, group2_sequences, output_path)` | Side-by-side population comparison |
| `plot_intervention_effects` | `(pre_sequences, post_sequences, output_path)` | Before/after intervention comparison |
| `plot_sequence_similarity` | `(sequences, output_path)` | Pairwise similarity heatmap |

## Usage Examples

### Complete Visualization Pipeline

```python
from metainformant.life_events.visualization.timeline import (
    plot_event_timeline,
    plot_domain_distribution,
    plot_temporal_density,
)
from metainformant.life_events.visualization.statistical import plot_event_frequency_heatmap
from metainformant.life_events.visualization.network import plot_transition_network

# Individual trajectory
plot_event_timeline(seq, "output/timeline.png", show_domains=True)

# Cohort-level plots
plot_domain_distribution(sequences, "output/domains.png")
plot_event_frequency_heatmap(sequences, "output/freq_heatmap.png")
plot_transition_network(sequences, "output/network.png", min_transitions=2)
plot_temporal_density(sequences, "output/density.png")
```

### Embedding Visualization

```python
from metainformant.life_events.models.embeddings import learn_event_embeddings
from metainformant.life_events.visualization.statistical import (
    plot_event_embeddings,
    plot_embedding_clusters,
)

emb = learn_event_embeddings(sequences, embedding_dim=50)
plot_event_embeddings(emb, "output/embed_pca.png", method="pca")
plot_event_embeddings(emb, "output/embed_tsne.png", method="tsne")
```

## See Also

- [Events](events.md) -- Data structures used as plot inputs
- [Models](models.md) -- Generating predictions and embeddings to plot
- [Workflow](workflow.md) -- Automated visualization in pipelines
