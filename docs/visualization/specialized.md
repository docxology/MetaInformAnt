# Specialized Visualization

This document provides comprehensive documentation for specialized visualization types in METAINFORMANT, including Venn diagrams, Sankey diagrams, chord diagrams, alluvial plots, and other advanced chart types used in bioinformatics.

## Overview

Specialized visualization includes advanced plotting types commonly used in systems biology, set analysis, flow visualization, and complex data relationships. These tools complement the core plotting modules with publication-quality specialized charts.

All functions live in `metainformant.visualization.plots.specialized`; the examples import that module as `spec_viz`. Each function returns the `matplotlib.axes.Axes`, except `plot_sankey_diagram` and `plot_upset_plot`, which return a figure and take no `ax` argument. Most accept `output_path=` to save the figure and `figsize=` to size a figure they create. The circular plots (chord diagram, circular bar plot, circular network layout) draw on polar axes: they create one for you, and an `ax` you pass in must be a polar axes (for example one created with `projection="polar"`).

## Module Functions

### Set Analysis

#### Venn Diagrams

`plot_venn_diagram(sets, *, ax=None, output_path=None, figsize=(8, 8), **kwargs)` takes a dict that maps each name to a Python `set`. Exactly two or three sets are supported. With the optional `matplotlib-venn` package installed it draws `venn2` or `venn3` and passes extra keyword arguments to them; otherwise a built-in fallback draws overlapping circles labelled with element counts (without the set names for three sets) and ignores extra keyword arguments. Any other number of sets raises an error (`ValueError` with matplotlib-venn, `ImportError` without it).

```python
from metainformant.visualization.plots import specialized as spec_viz

# Create Venn diagram for set intersections
sets = {
    'Dataset_A': {'gene1', 'gene2', 'gene3', 'gene4'},
    'Dataset_B': {'gene2', 'gene3', 'gene5', 'gene6'},
    'Dataset_C': {'gene1', 'gene3', 'gene5', 'gene7'}
}

ax = spec_viz.plot_venn_diagram(sets, figsize=(8, 8))
```

#### UpSet Plots

`plot_upset_plot(data, *, output_path=None, figsize=(10, 6), **kwargs)` draws an UpSet plot of any number of sets and returns the matplotlib `Figure`. It requires the optional `upsetplot` package (an `ImportError` is raised without it) and passes extra keyword arguments to `upsetplot.plot`.

```python
from metainformant.visualization.plots import specialized as spec_viz

# Create UpSet plot for complex set intersections (requires upsetplot)
data = {
    'Condition_1': {'A', 'B', 'C', 'D'},
    'Condition_2': {'B', 'C', 'E', 'F'},
    'Condition_3': {'A', 'C', 'E', 'G'},
    'Condition_4': {'B', 'D', 'F', 'H'}
}

fig = spec_viz.plot_upset_plot(data, figsize=(10, 6))
```

### Flow Visualization

#### Sankey Diagrams

`plot_sankey_diagram(flows, *, output_path=None, figsize=(10, 8), **kwargs)` takes a list of `(source, target, value)` tuples and returns a figure. There is no `ax` argument and no `interactive` switch; which kind of figure you get depends on whether Plotly can be imported:

- **With Plotly**: an interactive Plotly Sankey figure. Extra keyword arguments go to `fig.update_layout`, and `output_path` writes an HTML file (its suffix is replaced by `.html`). The node indices of the links come from the helper `_build_sankey_links`: each distinct name in `flows` is one node (a name that is both a source and a target is a single node), numbered in the order it first appears, so an index is never reassigned once it exists.
- **Without Plotly**: a simple matplotlib flow diagram (sources on the left, targets on the right, line width and marker size from the value) returned as a matplotlib `Figure`. Sources and targets keep the order in which they first appear in `flows`, so the same input always gives the same row order. Extra keyword arguments are ignored and `output_path` saves an image, so use an image extension.

```python
from metainformant.visualization.plots import specialized as spec_viz

# Create Sankey diagram for flow data
flows = [
    ('Source_A', 'Target_X', 10),
    ('Source_A', 'Target_Y', 5),
    ('Source_B', 'Target_X', 8),
    ('Source_B', 'Target_Z', 12),
    ('Source_C', 'Target_Y', 6)
]

fig = spec_viz.plot_sankey_diagram(flows, figsize=(10, 8))
```

#### Alluvial Diagrams

`plot_alluvial_diagram(data, stages, *, ax=None, output_path=None, figsize=(12, 8), **kwargs)` takes a DataFrame and the list of column names that are the stages. If the optional `alluvial` package can be imported the function delegates to `alluvial.plot(data, stages, ax=ax, **kwargs)`. Otherwise it draws a simplified fallback with one line per DataFrame row connecting that row's values in the stage columns; it does not draw flow widths and does not use a count column.

```python
# Create alluvial diagram for transitions
import pandas as pd
from metainformant.visualization.plots import specialized as spec_viz

data = pd.DataFrame({
    'Stage_1': ['A', 'A', 'B', 'B', 'C', 'C'],
    'Stage_2': ['X', 'Y', 'X', 'Y', 'X', 'Z'],
    'Stage_3': ['P', 'P', 'Q', 'Q', 'R', 'R'],
    'Count': [5, 3, 4, 6, 2, 8]
})

ax = spec_viz.plot_alluvial_diagram(data, ['Stage_1', 'Stage_2', 'Stage_3'])
```

### Network Analysis

#### Chord Diagrams

`plot_chord_diagram(matrix, labels=None, *, ax=None, output_path=None, figsize=(8, 8), **kwargs)` takes a square NumPy array of pairwise values and an optional list with one label per row. Nodes sit evenly around a circle and one chord is drawn for each pair `i < j` (the upper triangle), with its height scaled by the value relative to the largest value. Extra keyword arguments are ignored.

```python
import numpy as np
from metainformant.visualization.plots import specialized as spec_viz

# Create chord diagram for relationships
matrix = np.array([
    [0, 5, 3, 2],
    [5, 0, 4, 1],
    [3, 4, 0, 6],
    [2, 1, 6, 0]
])
labels = ['Group_A', 'Group_B', 'Group_C', 'Group_D']

ax = spec_viz.plot_chord_diagram(matrix, labels, figsize=(8, 8))
```

#### Circular Network Layout

`plot_network_circular_layout(graph, *, ax=None, output_path=None, figsize=(10, 10), **kwargs)` draws a NetworkX graph with its nodes evenly spaced on a circle of polar axes, with node names written around it. Extra keyword arguments are ignored; for a plain circular layout that accepts `networkx.draw` options see `plot_network_circular` in [Network Plots](networks.md).

```python
import networkx as nx
from metainformant.visualization.plots import specialized as spec_viz

# Plot network in circular layout
G = nx.erdos_renyi_graph(20, 0.1, seed=1)
ax = spec_viz.plot_network_circular_layout(G, figsize=(10, 10))
```

### Multi-dimensional Data

#### Radar Charts

Radar charts are not part of this module: `plot_radar_chart` lives in `metainformant.visualization.plots.multidim` (see [Multi-dimensional Plots](multidim.md)). Each numeric column of the DataFrame is one axis (at least three are needed) and each row is one polygon, labelled by its index.

```python
import pandas as pd
from metainformant.visualization.plots.multidim import plot_radar_chart

# Create radar chart for multi-dimensional data
data = pd.DataFrame({
    'Trait_A': [0.8, 0.6, 0.7],
    'Trait_B': [0.6, 0.8, 0.7],
    'Trait_C': [0.9, 0.5, 0.7],
    'Trait_D': [0.4, 0.6, 0.7],
    'Trait_E': [0.7, 0.8, 0.7]
}, index=['Sample_1', 'Sample_2', 'Sample_3'])

ax = plot_radar_chart(data, figsize=(8, 8))
```

#### Circular Bar Plots

`plot_circular_barplot(values, labels=None, *, ax=None, output_path=None, figsize=(8, 8), **kwargs)` draws one bar per value at evenly spaced angles; `values` must be a NumPy array and `labels` an optional list written around the circle. Extra keyword arguments go to `Axes.bar`; `width` (0.4), `bottom` (0.0) and `alpha` (0.7) are defaults that your keyword arguments override.

```python
import numpy as np
from metainformant.visualization.plots import specialized as spec_viz

# Create circular bar plot
values = np.random.default_rng(0).random(12) * 100
labels = [f'Category_{i+1}' for i in range(12)]

ax = spec_viz.plot_circular_barplot(values, labels, figsize=(8, 8), color='tab:blue')
```

## Integration Examples

### With Multi-Omics Data
```python
import numpy as np
import pandas as pd
from metainformant.multiomics.analysis import integration
from metainformant.visualization.plots import specialized as spec_viz

# Integrate multi-omics data: one DataFrame per layer, samples as rows
rng = np.random.default_rng(0)
samples = [f'S{i}' for i in range(20)]
omics_data = {
    'rna': pd.DataFrame(rng.random((20, 50)), index=samples, columns=[f'gene{i}' for i in range(50)]),
    'protein': pd.DataFrame(rng.random((20, 30)), index=samples, columns=[f'gene{i}' for i in range(20, 50)]),
    'metabolomics': pd.DataFrame(rng.random((20, 20)), index=samples, columns=[f'met{i}' for i in range(20)]),
}

integrated = integration.integrate_omics_data(data=omics_data)

# Venn diagram of feature overlaps between the layers
# (use plot_upset_plot instead for more than three layers)
feature_sets = {layer: set(df.columns) for layer, df in integrated.data.items()}
ax = spec_viz.plot_venn_diagram(feature_sets)
```

### With GWAS Results
```python
from metainformant.visualization.plots import specialized as spec_viz

# GWAS analysis results
# Create Sankey diagram for workflow steps (illustrative values: % of the raw variants)
workflow_flows = [
    ('Raw_Variants', 'QC_Filtered', 95),
    ('Raw_Variants', 'Removed_LowQual', 5),
    ('QC_Filtered', 'Associated_SNPs', 10),
    ('QC_Filtered', 'Non_Associated', 85),
    ('Associated_SNPs', 'Significant_Hits', 3),
    ('Associated_SNPs', 'Suggestive_Hits', 7)
]

fig = spec_viz.plot_sankey_diagram(workflow_flows)
```

### With Ontology Analysis
```python
import numpy as np
from metainformant.visualization.plots import specialized as spec_viz

# GO enrichment results
# Create chord diagram for GO term relationships
go_terms = ['GO:0008150', 'GO:0003674', 'GO:0005575', 'GO:0009987']
relationship_matrix = np.random.default_rng(0).random((4, 4))  # stand-in for a term-term similarity matrix
relationship_matrix = (relationship_matrix + relationship_matrix.T) / 2
np.fill_diagonal(relationship_matrix, 0)

ax = spec_viz.plot_chord_diagram(relationship_matrix, go_terms)
```

## Output Options

Static plots are saved at 300 dpi and missing parent directories are created:
```python
from metainformant.visualization.plots import specialized as spec_viz

sets = {'A': {'gene1', 'gene2', 'gene3'}, 'B': {'gene2', 'gene3', 'gene4'}}
flows = [('Source_A', 'Target_X', 10), ('Source_B', 'Target_X', 8)]

# Save static plots
ax = spec_viz.plot_venn_diagram(sets, output_path="venn_diagram.png")

# Sankey: an image without Plotly; with Plotly the suffix is replaced and "sankey.html" is written
fig = spec_viz.plot_sankey_diagram(flows, output_path="sankey.png")
```

UpSet plots save an image the same way (this needs `upsetplot`):
```python
from metainformant.visualization.plots import specialized as spec_viz

data = {'Set1': {'A', 'B', 'C'}, 'Set2': {'B', 'C', 'D'}, 'Set3': {'A', 'D', 'E'}}
fig = spec_viz.plot_upset_plot(data, output_path="upset.png")
```

## Dependencies

- **Required**: matplotlib, numpy, pandas
- **Venn Diagrams**: matplotlib-venn (optional; a built-in fallback handles two or three sets)
- **UpSet Plots**: upsetplot (optional, but `plot_upset_plot` raises `ImportError` without it)
- **Alluvial Diagrams**: alluvial (optional; a simplified fallback is used without it)
- **Sankey Diagrams**: plotly (optional, in the `visualization` extra; a matplotlib fallback is used without it)
- **Network Analysis**: networkx (a core dependency)

## Performance Considerations

- **Large Datasets**: Venn diagrams support only two or three sets; use UpSet plots for more
- **Complex Networks**: Circular layouts scale to ~100 nodes
- **Interactive Plots**: The interactive Sankey view needs Plotly (the `visualization` extra); without it a simple static diagram is drawn
- **Memory Usage**: Chord diagrams require O(n²) memory for dense matrices

## Color Schemes

Recommended color schemes for your own specialized figures (most functions here use fixed colors; only `plot_circular_barplot` takes a `color=` argument directly):
```python
import matplotlib.pyplot as plt
import numpy as np

# Venn diagrams
venn_colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FFEAA7']

# Sankey diagrams
sankey_colors = {
    'source': '#FF6B6B',
    'target': '#4ECDC4',
    'intermediate': '#45B7D1'
}

# Chord diagrams
chord_cmap = 'coolwarm'  # Diverging colormap for relationships

# Radar charts
radar_colors = plt.cm.tab10(np.linspace(0, 1, 10))  # Up to 10 series
```

## Examples

### Complete Specialized Visualization Workflow
```python
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from metainformant.visualization.plots import specialized as spec_viz
from metainformant.visualization.plots.multidim import plot_radar_chart

# Create comprehensive specialized visualization suite
fig = plt.figure(figsize=(18, 12))

# Venn diagram (regular axes)
sets = {
    'Group_A': set(range(50)),
    'Group_B': set(range(25, 75)),
    'Group_C': set(range(40, 90))
}
spec_viz.plot_venn_diagram(sets, ax=fig.add_subplot(2, 3, 1))

# Chord diagram (needs polar axes)
rng = np.random.default_rng(0)
matrix = rng.random((5, 5))
matrix = (matrix + matrix.T) / 2
np.fill_diagonal(matrix, 0)
labels = ['A', 'B', 'C', 'D', 'E']
spec_viz.plot_chord_diagram(matrix, labels, ax=fig.add_subplot(2, 3, 2, projection='polar'))

# Radar chart: one polygon per row, one axis per column
radar_data = pd.DataFrame(
    {
        'Trait_A': [0.8, 0.6, 0.7],
        'Trait_B': [0.6, 0.8, 0.7],
        'Trait_C': [0.7, 0.9, 0.8],
        'Trait_D': [0.9, 0.6, 0.6],
        'Trait_E': [0.5, 0.8, 0.9],
    },
    index=['Sample1', 'Sample2', 'Sample3'],
)
plot_radar_chart(radar_data, ax=fig.add_subplot(2, 3, 3, projection='polar'))

# Circular bar plot (needs polar axes)
values = rng.random(8) * 100
labels = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug']
spec_viz.plot_circular_barplot(values, labels, ax=fig.add_subplot(2, 3, 4, projection='polar'))

# Simplified alluvial fallback: one line per row across numeric stage columns
stages = pd.DataFrame({
    'Stage1': [1.0, 2.0, 3.0],
    'Stage2': [2.0, 2.5, 1.0],
    'Stage3': [1.5, 3.0, 2.0],
})
spec_viz.plot_alluvial_diagram(stages, ['Stage1', 'Stage2', 'Stage3'], ax=fig.add_subplot(2, 3, 5))

plt.tight_layout()
plt.savefig("specialized_visualizations.png", dpi=300, bbox_inches='tight')

# Sankey creates its own figure, so it cannot be placed in the grid above
flows = [
    ('Input', 'Process_A', 40),
    ('Input', 'Process_B', 60),
    ('Process_A', 'Output_1', 30),
    ('Process_A', 'Output_2', 10),
    ('Process_B', 'Output_2', 40),
    ('Process_B', 'Output_3', 20)
]
sankey_fig = spec_viz.plot_sankey_diagram(flows)
```

UpSet plots (which need the optional `upsetplot` package) also create their own figure:
```python
from metainformant.visualization.plots import specialized as spec_viz

upset_data = {
    'Set1': {'A', 'B', 'C', 'D'},
    'Set2': {'B', 'C', 'E', 'F'},
    'Set3': {'A', 'D', 'E', 'G'}
}
upset_fig = spec_viz.plot_upset_plot(upset_data)
```

## Troubleshooting

### Common Issues

1. **Missing Dependencies**: Install optional packages for full functionality (`plot_upset_plot` raises `ImportError` without `upsetplot`)
2. **Large Data**: Reduce complexity for better performance
3. **Polar Axes Errors**: Chord diagrams, circular bar plots and circular network layouts raise an `AttributeError` when given regular (non-polar) axes
4. **Layout Issues**: Adjust figure sizes for complex diagrams

### Data Requirements

- **Venn Diagrams**: Dictionary of set objects (exactly two or three sets)
- **Chord Diagrams**: Square relationship matrix (NumPy array)
- **Sankey Diagrams**: List of (source, target, value) tuples
- **Alluvial Diagrams**: Pandas DataFrame with one column per stage (numeric values for the built-in fallback)
- **Radar Charts**: DataFrame with variables as columns and observations as rows

## Related Documentation

- **[Core Visualization](basic.md)**: Basic plotting functions
- **[Statistical Plots](statistical.md)**: Statistical visualizations
- **[Network Visualization](networks.md)**: Network analysis plots
- **[Multi-dimensional Plots](multidim.md)**: Complex data visualization

This module provides advanced specialized visualization capabilities for complex bioinformatics data analysis and presentation.
