# Life Events Visualization

Plotting functions for life event data: event timelines, transition networks, embedding scatter plots, attention heatmaps, population comparisons, and statistical summaries.

## Contents

| File | Purpose |
|------|---------|
| `timeline.py` | Event timeline plots, domain distribution, temporal density charts |
| `network.py` | Event transition network graph visualization |
| `statistical.py` | Embedding scatter, attention heatmaps, importance plots, population comparison |

## Key Functions

| Function | Description |
|----------|-------------|
| `plot_event_timeline()` | Timeline visualization of events for an individual |
| `plot_domain_distribution()` | Bar chart of event counts by domain |
| `plot_domain_timeline()` | Domain-colored timeline across multiple individuals |
| `plot_temporal_density()` | Kernel density estimate of event timing |
| `plot_temporal_patterns()` | Recurring pattern visualization across sequences |
| `plot_transition_network()` | Directed graph of event-to-event transitions |
| `plot_event_embeddings()` | 2D scatter of learned event embeddings |
| `plot_attention_heatmap()` | Attention weight heatmap for sequence models |
| `plot_prediction_importance()` | Feature importance bar chart for predictions |
| `plot_intervention_effects()` | Before/after intervention effect visualization |
| `plot_population_comparison()` | Side-by-side comparison of two population cohorts |
| `plot_event_frequency_heatmap()` | Event type frequency heatmap across cohorts |

## Usage

```python
from metainformant.life_events.visualization.timeline import plot_event_timeline
from metainformant.life_events.visualization.statistical import plot_event_embeddings

fig = plot_event_timeline(sequence, output_path="output/timeline.png")
fig = plot_event_embeddings(embeddings, labels=clusters, output_path="output/embed.png")
```
