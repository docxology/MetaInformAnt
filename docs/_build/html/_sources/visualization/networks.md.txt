# Network Plots

Network visualization functions for network graphs including basic network plots, circular layouts, hierarchical layouts, force-directed layouts, and community-based visualizations.

## Functions

### `network_plot(nodes, edges, *, node_sizes=None, node_colors=None, ax=None, **kwargs)`

Create a network graph visualization.

**Example:**
```python
from metainformant.visualization import network_plot

nodes = ['A', 'B', 'C']
edges = [('A', 'B'), ('B', 'C')]
ax = network_plot(nodes, edges)
```

### `circular_network_plot(nodes, edges, *, node_sizes=None, node_colors=None, ax=None, **kwargs)`

Create a circular network plot.

### `hierarchical_network_plot(nodes, edges, root=None, *, node_sizes=None, node_colors=None, ax=None, **kwargs)`

Create a hierarchical network plot.

### `force_directed_plot(nodes, edges, *, node_sizes=None, node_colors=None, iterations=50, ax=None, **kwargs)`

Create a force-directed network plot.

### `community_network_plot(nodes, edges, communities=None, *, node_sizes=None, ax=None, **kwargs)`

Create a network plot colored by community membership.

