# Animations

Animation functions for dynamic visualizations including time series, evolutionary processes, clustering iterations, network evolution, and trajectory inference.

## Functions

### `animate_time_series(series, *, interval_ms=100, repeat=False, init_points=1, ax=None)`

Animate a time series as a growing line.

**Example:**
```python
from metainformant.visualization import animate_time_series

fig, anim = animate_time_series([0, 1, 3, 2, 5], interval_ms=150)
# anim.save("series.mp4")  # optional - requires matplotlib writer
```

### `animate_evolution(generations, *, interval_ms=200, ax=None, **kwargs)`

Animate evolutionary process over generations.

### `animate_clustering(iterations, *, interval_ms=300, ax=None, **kwargs)`

Animate clustering iterations.

### `animate_network(network_states, *, interval_ms=200, ax=None, **kwargs)`

Animate network evolution over time.

### `animate_trajectory(trajectory_data, *, interval_ms=200, ax=None, **kwargs)`

Animate trajectory inference over time.
