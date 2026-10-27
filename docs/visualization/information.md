# Information Theory Plots

Information theory visualization functions for information-theoretic analysis including entropy profiles, mutual information matrices, Rényi entropy spectra, information landscapes, and information networks. They live in `metainformant.visualization.analysis.information`, return the `matplotlib.axes.Axes`, and accept `ax=` (draw on existing axes) and `output_path=` (save the figure at 300 dpi). Inputs are type-checked, so pass a `dict`, `list` or `numpy.ndarray` exactly as described below.

## Functions

### `plot_entropy_profile(entropy_data, *, ax=None, output_path=None, **kwargs)`

Plot Shannon entropy as a function of k-mer size. `entropy_data` maps a label (for example a
sequence or sample name) to a list of entropy values, one per k-mer size (k = 1, 2, ...). Extra keyword
arguments go to `Axes.plot` (`label` and `marker` are set by the function).

**Example:**
```python
from metainformant.visualization.analysis.information import plot_entropy_profile

entropy_data = {
    "seq_1": [1.98, 3.90, 5.71, 7.35],
    "seq_2": [1.85, 3.52, 5.01, 6.28],
}
ax = plot_entropy_profile(entropy_data)
```

### `plot_mutual_information_matrix(mi_matrix, labels=None, *, ax=None, output_path=None, **kwargs)`

Plot a mutual information matrix as a heatmap with a colorbar. `mi_matrix` must be a square NumPy array; `labels`, if given, must be a `list` with one label per row and column.

**Example:**
```python
import numpy as np
from metainformant.visualization.analysis.information import plot_mutual_information_matrix

mi_matrix = np.array([
    [1.00, 0.35, 0.10],
    [0.35, 1.00, 0.22],
    [0.10, 0.22, 1.00],
])
ax = plot_mutual_information_matrix(mi_matrix, labels=["gene_a", "gene_b", "gene_c"])
```

### `plot_renyi_spectra(renyi_data, alpha_values, *, ax=None, output_path=None, **kwargs)`

Plot Rényi entropy as a function of the order α. `renyi_data` maps a label to a list of entropy values, one per entry of `alpha_values`, and `alpha_values` is a `list` of the α values used. When 1.0 is among the α values a dashed vertical line marks Shannon entropy (α = 1). Extra keyword arguments go to `Axes.plot` (`label` and `marker` are set by the function).

**Example:**
```python
from metainformant.visualization.analysis.information import plot_renyi_spectra

alpha_values = [0.5, 1.0, 2.0, 4.0]
renyi_data = {
    "seq_1": [1.95, 1.90, 1.78, 1.61],
    "seq_2": [1.70, 1.62, 1.45, 1.20],
}
ax = plot_renyi_spectra(renyi_data, alpha_values)
```

### `plot_information_landscape(landscape_data, x_coords=None, y_coords=None, *, ax=None, output_path=None, **kwargs)`

Plot a 2-D array of information values as a filled contour map with a colorbar. Rows of `landscape_data` run along y and columns along x; `x_coords` and `y_coords` (NumPy arrays of matching length) default to the array indices.

**Example:**
```python
import numpy as np
from metainformant.visualization.analysis.information import plot_information_landscape

x = np.linspace(0, 1, 30)
y = np.linspace(0, 1, 30)
landscape = np.sin(3 * x)[None, :] * np.cos(3 * y)[:, None]
ax = plot_information_landscape(landscape, x_coords=x, y_coords=y)
```

### `plot_information_network(nodes, edges, *, ax=None, output_path=None, **kwargs)`

Plot a directed information flow network with NetworkX. `nodes` is a list of node names and `edges` a list of `(source, target, weight)` tuples whose endpoints must all appear in `nodes`. Node size follows the node's degree and edge width follows the weight, so use non-negative weights.

**Example:**
```python
from metainformant.visualization.analysis.information import plot_information_network

nodes = ["A", "B", "C", "D"]
edges = [("A", "B", 0.8), ("B", "C", 0.4), ("A", "D", 0.2)]
ax = plot_information_network(nodes, edges)
```
