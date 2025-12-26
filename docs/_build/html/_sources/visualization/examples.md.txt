# Visualization Examples

Comprehensive usage examples for the visualization module.

## Basic Plotting

```python
from metainformant.visualization import lineplot, scatter_plot, bar_plot, heatmap
import numpy as np

# Line plot
ax = lineplot(None, [1, 4, 2, 8, 5], label="Data")
ax.set_xlabel("Index")
ax.set_ylabel("Value")

# Scatter plot
x = np.random.random(100)
y = np.random.random(100)
ax = scatter_plot(x, y, xlabel="X", ylabel="Y")

# Bar plot
ax = bar_plot(["A", "B", "C"], [10, 20, 15])

# Heatmap
data = np.random.random((10, 10))
ax = heatmap(data, annot=True)
```

## Statistical Visualization

```python
from metainformant.visualization import histogram, box_plot, qq_plot, correlation_heatmap
import pandas as pd
import numpy as np

# Histogram
data = np.random.normal(0, 1, 1000)
ax = histogram(data, bins=30)

# Box plot
data = [np.random.normal(0, 1, 100) for _ in range(3)]
ax = box_plot(data, labels=["A", "B", "C"])

# Q-Q plot
pvals = np.random.uniform(0, 1, 1000)
ax = qq_plot(pvals.tolist())

# Correlation heatmap
df = pd.DataFrame(np.random.random((10, 5)))
ax = correlation_heatmap(df)
```

## Genomic Visualization

```python
from metainformant.visualization import manhattan_plot, volcano_plot, regional_plot
import pandas as pd
import numpy as np

# Manhattan plot
data = pd.DataFrame({
    'position': range(1000, 10000, 100),
    'pvalue': np.random.uniform(1e-9, 1, 90),
    'chromosome': ['chr1'] * 45 + ['chr2'] * 45
})
data['neg_log10_p'] = -np.log10(data['pvalue'])
ax = manhattan_plot(data, 'position', 'neg_log10_p', 'chromosome')

# Volcano plot
data = pd.DataFrame({
    'log2fc': np.random.normal(0, 1, 100),
    'pvalue': np.random.uniform(0, 1, 100)
})
data['neg_log10_p'] = -np.log10(data['pvalue'])
ax = volcano_plot(data, 'log2fc', 'neg_log10_p')
```

## Dimensionality Reduction

```python
from metainformant.visualization import pca_plot, umap_plot, pca_scree_plot
import pandas as pd
import numpy as np

# PCA plot
data = pd.DataFrame({
    'PC1': np.random.normal(0, 1, 100),
    'PC2': np.random.normal(0, 1, 100),
    'group': ['A'] * 50 + ['B'] * 50
})
ax = pca_plot(data, hue='group')

# PCA scree plot
variance = np.array([0.4, 0.3, 0.2, 0.1])
ax = pca_scree_plot(variance)
```

## Multi-panel Figures

```python
from metainformant.visualization.layout import create_multi_panel, add_shared_axis_labels
import matplotlib.pyplot as plt

fig, axes = create_multi_panel(4, layout='grid')
# Plot on each axis...
add_shared_axis_labels(fig, xlabel='Time', ylabel='Value')
plt.savefig('output/multi_panel.png')
```

## Publication-quality Export

```python
from metainformant.visualization.export import save_figure_multiformat
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
# Create plot...
save_figure_multiformat(fig, 'output/plot', formats=['png', 'pdf', 'svg'], dpi=300)
```

