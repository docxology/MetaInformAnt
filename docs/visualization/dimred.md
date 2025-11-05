# Dimensionality Reduction Plots

Visualization functions for dimensionality reduction techniques including PCA, UMAP, t-SNE, and related diagnostic plots.

## Functions

### `pca_plot(data, *, pc_x=1, pc_y=2, hue=None, ax=None, **kwargs)`

Create a PCA scatter plot.

**Example:**
```python
from metainformant.visualization import pca_plot
import pandas as pd
import numpy as np

data = pd.DataFrame({
    'PC1': np.random.normal(0, 1, 100),
    'PC2': np.random.normal(0, 1, 100),
    'group': ['A'] * 50 + ['B'] * 50
})
ax = pca_plot(data, hue='group')
```

### `umap_plot(data, *, x_col=0, y_col=1, hue=None, ax=None, title='UMAP Plot', **kwargs)`

Create a UMAP visualization plot.

### `tsne_plot(data, *, x_col=0, y_col=1, hue=None, ax=None, title='t-SNE Plot', **kwargs)`

Create a t-SNE visualization plot.

### `pca_scree_plot(explained_variance, *, n_components=None, ax=None, title='PCA Scree Plot', **kwargs)`

Create a PCA scree plot showing variance explained.

### `pca_loadings_plot(loadings, pc_x=0, pc_y=1, *, feature_names=None, ax=None, title='PCA Loadings Plot', **kwargs)`

Create a PCA loadings plot.

### `biplot(scores, loadings, *, pc_x=0, pc_y=1, feature_names=None, sample_names=None, ax=None, title='PCA Biplot', **kwargs)`

Create a PCA biplot showing both samples and loadings.

