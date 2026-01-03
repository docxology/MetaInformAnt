# Specialized Visualization

This document provides comprehensive documentation for specialized visualization types in METAINFORMANT, including Venn diagrams, Sankey diagrams, chord diagrams, alluvial plots, and other advanced chart types used in bioinformatics.

## Overview

Specialized visualization includes advanced plotting types commonly used in systems biology, set analysis, flow visualization, and complex data relationships. These tools complement the core plotting modules with publication-quality specialized charts.

## Module Functions

### Set Analysis

#### Venn Diagrams
```python
from metainformant.visualization import specialized as spec_viz

# Create Venn diagram for set intersections
sets = {
    'Dataset_A': {'gene1', 'gene2', 'gene3', 'gene4'},
    'Dataset_B': {'gene2', 'gene3', 'gene5', 'gene6'},
    'Dataset_C': {'gene1', 'gene3', 'gene5', 'gene7'}
}

ax = spec_viz.plot_venn_diagram(sets, figsize=(8, 8))
```

#### UpSet Plots
```python
# Create UpSet plot for complex set intersections
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
```python
# Create Sankey diagram for flow data
flows = [
    ('Source_A', 'Target_X', 10),
    ('Source_A', 'Target_Y', 5),
    ('Source_B', 'Target_X', 8),
    ('Source_B', 'Target_Z', 12),
    ('Source_C', 'Target_Y', 6)
]

# Static matplotlib version
ax = spec_viz.plot_sankey_diagram(flows, figsize=(10, 8))

# Interactive Plotly version
fig = spec_viz.plot_sankey_diagram(flows, figsize=(10, 8), interactive=True)
```

#### Alluvial Diagrams
```python
# Create alluvial diagram for transitions
import pandas as pd

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
```python
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
```python
# Plot network in circular layout
import networkx as nx

G = nx.erdos_renyi_graph(20, 0.1)
ax = spec_viz.plot_network_circular_layout(G, figsize=(10, 10))
```

### Multi-dimensional Data

#### Radar Charts
```python
# Create radar chart for multi-dimensional data
import pandas as pd

data = pd.DataFrame({
    'Sample_1': [0.8, 0.6, 0.9, 0.4, 0.7],
    'Sample_2': [0.6, 0.8, 0.5, 0.6, 0.8],
    'Sample_3': [0.7, 0.7, 0.7, 0.7, 0.7]
}, index=['Trait_A', 'Trait_B', 'Trait_C', 'Trait_D', 'Trait_E'])

ax = spec_viz.plot_radar_chart(data, figsize=(8, 8))
```

#### Circular Bar Plots
```python
# Create circular bar plot
values = np.random.rand(12) * 100
labels = [f'Category_{i+1}' for i in range(12)]

ax = spec_viz.plot_circular_barplot(values, labels, figsize=(8, 8))
```

## Integration Examples

### With Multi-Omics Data
```python
from metainformant.multiomics import integration, visualization as multi_viz
from metainformant.visualization import specialized as spec_viz

# Integrate multi-omics data
omics_data = {
    'rna': pd.DataFrame(np.random.rand(100, 50)),
    'protein': pd.DataFrame(np.random.rand(100, 30)),
    'metabolomics': pd.DataFrame(np.random.rand(100, 20))
}

integrated = integration.integrate_omics_data(**omics_data)

# Create specialized visualizations
# UpSet plot for feature overlaps
feature_sets = {
    'RNA': set(range(50)),
    'Protein': set(range(50, 80)),
    'Metabolites': set(range(80, 100))
}
fig = spec_viz.plot_upset_plot(feature_sets)
```

### With GWAS Results
```python
from metainformant.gwas import visualization as gwas_viz
from metainformant.visualization import specialized as spec_viz

# GWAS analysis results
# Create Sankey diagram for workflow steps
workflow_flows = [
    ('Raw_Variants', 'QC_Filtered', 950000),
    ('Raw_Variants', 'Removed_LowQual', 50000),
    ('QC_Filtered', 'Associated_SNPs', 5000),
    ('QC_Filtered', 'Non_Associated', 945000),
    ('Associated_SNPs', 'Significant_Hits', 500),
    ('Associated_SNPs', 'Suggestive_Hits', 4500)
]

fig = spec_viz.plot_sankey_diagram(workflow_flows)
```

### With Ontology Analysis
```python
from metainformant.ontology import go, visualization as ont_viz
from metainformant.visualization import specialized as spec_viz

# GO enrichment results
# Create chord diagram for GO term relationships
go_terms = ['GO:0008150', 'GO:0003674', 'GO:0005575', 'GO:0009987']
relationship_matrix = np.random.rand(4, 4)
relationship_matrix = (relationship_matrix + relationship_matrix.T) / 2
np.fill_diagonal(relationship_matrix, 0)

ax = spec_viz.plot_chord_diagram(relationship_matrix, go_terms)
```

## Output Options

All specialized visualization functions support:
```python
# Save static plots
ax = spec_viz.plot_venn_diagram(sets, output_path="venn_diagram.png")

# Save interactive plots
fig = spec_viz.plot_sankey_diagram(flows, output_path="sankey.html")

# UpSet plots
fig = spec_viz.plot_upset_plot(data, output_path="upset.png")
```

## Dependencies

- **Required**: matplotlib, numpy
- **Venn Diagrams**: matplotlib-venn (optional, fallback provided)
- **UpSet Plots**: upsetplot (optional)
- **Alluvial Diagrams**: alluvial (optional)
- **Interactive Plots**: plotly (optional)
- **Network Analysis**: networkx (optional)

## Performance Considerations

- **Large Datasets**: Venn diagrams work best with <5 sets
- **Complex Networks**: Circular layouts scale to ~100 nodes
- **Interactive Plots**: Plotly Sankey diagrams handle thousands of flows
- **Memory Usage**: Chord diagrams require O(n²) memory for dense matrices

## Color Schemes

Recommended color schemes for specialized plots:
```python
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
from metainformant.visualization import specialized as spec_viz
import numpy as np
import pandas as pd

# Create comprehensive specialized visualization suite
fig, axes = plt.subplots(2, 3, figsize=(18, 12))

# Venn diagram
sets = {
    'Group_A': set(range(50)),
    'Group_B': set(range(25, 75)),
    'Group_C': set(range(40, 90))
}
spec_viz.plot_venn_diagram(sets, ax=axes[0,0])

# Chord diagram
matrix = np.random.rand(5, 5)
matrix = (matrix + matrix.T) / 2
np.fill_diagonal(matrix, 0)
labels = ['A', 'B', 'C', 'D', 'E']
spec_viz.plot_chord_diagram(matrix, labels, ax=axes[0,1])

# Radar chart
radar_data = pd.DataFrame({
    'Sample1': [0.8, 0.6, 0.9, 0.7, 0.5],
    'Sample2': [0.6, 0.8, 0.7, 0.6, 0.8],
    'Sample3': [0.7, 0.7, 0.8, 0.9, 0.6]
})
spec_viz.plot_radar_chart(radar_data, ax=axes[0,2])

# Circular bar plot
values = np.random.rand(8) * 100
labels = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug']
spec_viz.plot_circular_barplot(values, labels, ax=axes[1,0])

# Alluvial diagram (placeholder - requires alluvial package)
# alluvial_data = pd.DataFrame({
#     'Stage1': ['A', 'B', 'C'] * 10,
#     'Stage2': ['X', 'Y', 'Z'] * 10,
#     'Count': np.random.randint(1, 10, 30)
# })
# spec_viz.plot_alluvial_diagram(alluvial_data, ['Stage1', 'Stage2'], ax=axes[1,1])

# Sankey diagram (static version)
flows = [
    ('Input', 'Process_A', 40),
    ('Input', 'Process_B', 60),
    ('Process_A', 'Output_1', 30),
    ('Process_A', 'Output_2', 10),
    ('Process_B', 'Output_2', 40),
    ('Process_B', 'Output_3', 20)
]
spec_viz.plot_sankey_diagram(flows, ax=axes[1,2])

plt.tight_layout()
plt.savefig("specialized_visualizations.png", dpi=300, bbox_inches='tight')

# Create separate interactive plots
# UpSet plot
upset_data = {
    'Set1': {'A', 'B', 'C', 'D'},
    'Set2': {'B', 'C', 'E', 'F'},
    'Set3': {'A', 'D', 'E', 'G'}
}
upset_fig = spec_viz.plot_upset_plot(upset_data)

# Interactive Sankey
interactive_sankey = spec_viz.plot_sankey_diagram(flows, interactive=True)
```

## Troubleshooting

### Common Issues

1. **Missing Dependencies**: Install optional packages for full functionality
2. **Large Data**: Reduce complexity for better performance
3. **Color Conflicts**: Use custom color schemes for clarity
4. **Layout Issues**: Adjust figure sizes for complex diagrams

### Data Requirements

- **Venn Diagrams**: Dictionary of set objects (up to 5 sets recommended)
- **Chord Diagrams**: Square relationship matrix
- **Sankey Diagrams**: List of (source, target, value) tuples
- **Alluvial Diagrams**: Pandas DataFrame with stage columns
- **Radar Charts**: DataFrame with variables as rows, observations as columns

## Related Documentation

- **[Core Visualization](basic.md)**: Basic plotting functions
- **[Statistical Plots](statistical.md)**: Statistical visualizations
- **[Network Visualization](networks.md)**: Network analysis plots
- **[Multi-dimensional Plots](multidim.md)**: Complex data visualization

This module provides advanced specialized visualization capabilities for complex bioinformatics data analysis and presentation.
