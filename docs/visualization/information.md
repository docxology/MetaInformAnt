# Information Theory Plots

Information theory visualization functions for information-theoretic analysis including entropy plots, mutual information plots, information profiles, Rényi entropy spectra, and information networks.

## Functions

### `entropy_plot(entropies, *, positions=None, ax=None, title='Entropy Plot', **kwargs)`

Plot entropy values across positions or sequences.

**Example:**
```python
from metainformant.visualization import entropy_plot
import numpy as np

entropies = np.random.uniform(0, 2, 100)
ax = entropy_plot(entropies)
```

### `mutual_information_plot(mi_matrix, labels=None, *, ax=None, title='Mutual Information Matrix', **kwargs)`

Plot mutual information matrix as heatmap.

### `information_profile_plot(profile_data, *, figsize=None, **kwargs)`

Plot information profile visualization.

### `renyi_spectrum_plot(alpha_range, entropies, *, ax=None, title='Rényi Entropy Spectrum', **kwargs)`

Plot Rényi entropy as a function of order α.

### `information_network_plot(mi_matrix, labels=None, *, threshold=0.1, ax=None, title='Information Network', **kwargs)`

Plot network visualization based on mutual information matrix.

