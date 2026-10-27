# Visualization Examples

Comprehensive usage examples for the visualization module.

## Basic Plotting

```python
from metainformant.visualization.plots.basic import lineplot, scatter_plot, bar_plot, heatmap
import numpy as np

# Line plot (y omitted: x is plotted against its index)
ax = lineplot(np.array([1, 4, 2, 8, 5]), label="Data")
ax.set_xlabel("Index")
ax.set_ylabel("Value")

# Scatter plot
x = np.random.random(100)
y = np.random.random(100)
ax = scatter_plot(x, y)
ax.set_xlabel("X")
ax.set_ylabel("Y")

# Bar plot
ax = bar_plot(["A", "B", "C"], np.array([10, 20, 15]))

# Heatmap
data = np.random.random((10, 10))
ax = heatmap(data)
```

## Statistical Visualization

```python
from metainformant.visualization.analysis.statistical import histogram, box_plot, qq_plot, correlation_heatmap
import pandas as pd
import numpy as np

# Histogram
data = np.random.normal(0, 1, 1000)
ax = histogram(data, bins=30)

# Box plot
data = [np.random.normal(0, 1, 100) for _ in range(3)]
ax = box_plot(data)
ax.set_xticklabels(["A", "B", "C"])

# Q-Q plot (p-values against the uniform distribution)
pvals = np.random.uniform(0, 1, 1000)
ax = qq_plot(pvals, distribution="uniform")

# Correlation heatmap
df = pd.DataFrame(np.random.random((10, 5)))
ax = correlation_heatmap(df)
```

## Genomic Visualization

```python
from metainformant.visualization.genomics.genomics import manhattan_plot, volcano_plot, regional_plot
import pandas as pd
import numpy as np

# Manhattan plot (raw p-values; -log10 is applied by the function)
data = pd.DataFrame({
    'position': range(1000, 10000, 100),
    'pvalue': np.random.uniform(1e-9, 1, 90),
    'chromosome': ['chr1'] * 45 + ['chr2'] * 45
})
ax = manhattan_plot(data, chr_col='chromosome', pos_col='position', pval_col='pvalue')

# Volcano plot
data = pd.DataFrame({
    'log2fc': np.random.normal(0, 1, 100),
    'pvalue': np.random.uniform(0, 1, 100)
})
ax = volcano_plot(data, log2fc_col='log2fc', pval_col='pvalue')
```

## Dimensionality Reduction

```python-snippet
from metainformant.visualization.plots.general import pca_plot
from metainformant.visualization.analysis.dimred import plot_umap
from metainformant.gwas.visualization.population import pca_scree_plot
import pandas as pd
import numpy as np

# PCA plot (pre-computed components plus the variance explained, in percent)
data = pd.DataFrame({
    'PC1': np.random.normal(0, 1, 100),
    'PC2': np.random.normal(0, 1, 100),
    'PC1_variance': 40.0,
    'PC2_variance': 30.0,
    'group': ['A'] * 50 + ['B'] * 50
})
ax = pca_plot(data, hue='group')

# PCA scree plot (GWAS module; saves the figure and returns a status dict)
variance = [0.4, 0.3, 0.2, 0.1]
result = pca_scree_plot(variance, output_file="output/pca_scree.png")
```

## Multi-panel Figures

```python-snippet
from metainformant.visualization.dashboards.composite import multi_panel
from metainformant.visualization.plots.basic import lineplot, scatter_plot
import numpy as np

# Each panel is a dict with a plotting callable, its keyword arguments and an optional title
x = np.linspace(0, 10, 100)
fig = multi_panel(
    [
        {"func": lineplot, "kwargs": {"x": x, "y": np.sin(x)}, "title": "Sine"},
        {"func": scatter_plot, "kwargs": {"x": x, "y": np.cos(x)}, "title": "Cosine"},
    ],
    ncols=2,
    output_path='output/multi_panel.png',
)
```

## Publication-quality Export

```python-snippet
from metainformant.visualization.config.conventions import save_figure_deterministic
from pathlib import Path
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
# Create plot...
Path('output').mkdir(exist_ok=True)
for fmt in ['png', 'pdf']:
    save_figure_deterministic(fig, f'output/plot.{fmt}', dpi=300)
```

