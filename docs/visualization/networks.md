# Network Plots

Network visualization functions for network graphs including basic network plots, circular layouts, hierarchical layouts, force-directed layouts, and community-based visualizations. They live in `metainformant.visualization.genomics.networks`, take a NetworkX graph `G`, return the `matplotlib.axes.Axes`, and accept `ax=` (draw on existing axes) and `output_path=` (save the figure at 300 dpi). Extra keyword arguments such as `node_color`, `node_size`, `edge_color`, `width`, `alpha` and `with_labels` are passed to `networkx.draw`; `figsize` sizes a figure the function creates.

## Functions

### `plot_network_basic(G, *, ax=None, output_path=None, **kwargs)`

Create a network graph visualization (spring layout) from a NetworkX graph. Extra keyword
arguments such as `node_color` and `node_size` are passed to `networkx.draw`, and
`layout_kwargs={...}` is passed to `networkx.spring_layout` (for example `{'seed': 3}` for a reproducible layout).

**Example:**
```python
import networkx as nx
from metainformant.visualization.genomics.networks import plot_network_basic

G = nx.Graph()
G.add_nodes_from(['A', 'B', 'C'])
G.add_edges_from([('A', 'B'), ('B', 'C')])
ax = plot_network_basic(G)
```

### `plot_network_circular(G, *, ax=None, output_path=None, **kwargs)`

Draw the graph on a circular layout (`networkx.circular_layout`). Nodes are green by default; there is no `layout_kwargs` option.

**Example:**
```python
import networkx as nx
from metainformant.visualization.genomics.networks import plot_network_circular

G = nx.karate_club_graph()
ax = plot_network_circular(G, with_labels=False)
```

### `plot_network_hierarchical(G, *, ax=None, output_path=None, **kwargs)`

Draw the graph in a hierarchical style. An undirected graph (or a directed graph with no root, meaning no node without incoming edges) is drawn with a shell layout, and arrows are drawn only for directed graphs. For a directed graph with at least one root, every node is placed in a column by its distance from the nearest root: every root is in the first column, so a graph with several roots is laid out from all of them. Nodes that no root reaches, such as a cycle with no root of its own, go in a final column. If the layered layout cannot be built, the function logs a warning and falls back to a spring layout.

**Example:**
```python
import networkx as nx
from metainformant.visualization.genomics.networks import plot_network_hierarchical

# Undirected graph: shell layout
G = nx.Graph([('root', 'a'), ('root', 'b'), ('a', 'c'), ('a', 'd'), ('b', 'e')])
ax = plot_network_hierarchical(G)

# Directed graph with two roots (r1 and r2): columns by distance from the nearest root
D = nx.DiGraph([('r1', 'a'), ('r2', 'b'), ('a', 'c'), ('b', 'c'), ('c', 'd')])
ax = plot_network_hierarchical(D)
```

### `plot_network_force_directed(G, *, ax=None, output_path=None, **kwargs)`

Draw the graph with the Fruchterman-Reingold force-directed layout (`networkx.fruchterman_reingold_layout`). Node labels are off by default (`with_labels=True` turns them on), and `layout_kwargs={...}` is passed to the layout function.

**Example:**
```python
import networkx as nx
from metainformant.visualization.genomics.networks import plot_network_force_directed

G = nx.karate_club_graph()
ax = plot_network_force_directed(G, layout_kwargs={'seed': 1})
```

### `plot_community_network(G, communities, *, ax=None, output_path=None, **kwargs)`

Draw the graph on a spring layout with nodes colored by community. `communities` is a non-empty dict mapping each node to an integer community id; nodes missing from the dict are drawn gray and listed as "Unassigned" in the legend, which otherwise names each community. `layout_kwargs={...}` is passed to `networkx.spring_layout`. The node colors come from the communities, so `node_color` is not accepted.

**Example:**
```python
import networkx as nx
from metainformant.visualization.genomics.networks import plot_community_network

G = nx.karate_club_graph()
communities = {node: int(G.nodes[node]['club'] != 'Mr. Hi') for node in G}
ax = plot_community_network(G, communities, layout_kwargs={'seed': 7})
```
