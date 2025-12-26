# Genomics Plots

Genomic visualization functions for GWAS and sequence analysis including Manhattan plots, volcano plots, regional plots, circular Manhattan plots, chromosome ideograms, coverage plots, and variant visualizations.

## Functions

### `manhattan_plot(data, x_col, y_col, chromosome_col='chromosome', *, p_threshold=5e-8, highlight_color='red', ax=None, **kwargs)`

Create a Manhattan plot for genome-wide association studies.

**Example:**
```python
from metainformant.visualization import manhattan_plot
import pandas as pd
import numpy as np

data = pd.DataFrame({
    'position': [1000, 2000, 3000],
    'pvalue': [1e-6, 1e-9, 0.01],
    'chromosome': ['chr1', 'chr1', 'chr2']
})
data['neg_log10_p'] = -np.log10(data['pvalue'])
ax = manhattan_plot(data, 'position', 'neg_log10_p', 'chromosome')
```

### `volcano_plot(data, x_col, y_col, *, p_threshold=0.05, fc_threshold=1.0, ax=None, **kwargs)`

Create a volcano plot for differential expression analysis.

**Example:**
```python
from metainformant.visualization import volcano_plot
import pandas as pd
import numpy as np

data = pd.DataFrame({
    'log2fc': [-2, 1, 0.5, -1.5],
    'pvalue': [0.001, 0.01, 0.5, 0.001]
})
data['neg_log10_p'] = -np.log10(data['pvalue'])
ax = volcano_plot(data, 'log2fc', 'neg_log10_p')
```

### `regional_plot(data, chromosome, start, end, x_col='position', y_col='neg_log10_p', *, p_threshold=5e-8, ax=None, **kwargs)`

Create a regional plot for a specific genomic region.

**Example:**
```python
from metainformant.visualization import regional_plot
import pandas as pd
import numpy as np

data = pd.DataFrame({
    'position': range(1000000, 1001000, 100),
    'neg_log10_p': np.random.uniform(0, 5, 10)
})
ax = regional_plot(data, 'chr1', 1000000, 1001000)
```

### `circular_manhattan_plot(data, x_col, y_col, chromosome_col='chromosome', *, p_threshold=5e-8, ax=None, **kwargs)`

Create a circular Manhattan plot for genome-wide visualization.

**Example:**
```python
from metainformant.visualization import circular_manhattan_plot
import pandas as pd
import numpy as np

data = pd.DataFrame({
    'position': [1000, 2000, 3000],
    'pvalue': [1e-6, 1e-9, 0.01],
    'chromosome': ['chr1', 'chr1', 'chr2']
})
data['neg_log10_p'] = -np.log10(data['pvalue'])
ax = circular_manhattan_plot(data, 'position', 'neg_log10_p', 'chromosome')
```

### `chromosome_ideogram(chromosomes, positions, values, *, threshold=None, ax=None, **kwargs)`

Create a chromosome ideogram with marked positions.

**Example:**
```python
from metainformant.visualization import chromosome_ideogram

ax = chromosome_ideogram(
    ['chr1', 'chr2'],
    [1000000, 2000000],
    [5.0, 6.0],
    threshold=5.0
)
```

### `coverage_plot(positions, coverage, *, ax=None, title='Coverage Plot', **kwargs)`

Create a coverage plot for sequencing data.

**Example:**
```python
from metainformant.visualization import coverage_plot
import numpy as np

positions = np.arange(1000, 2000, 10)
coverage = np.random.poisson(50, len(positions))
ax = coverage_plot(positions, coverage)
```

### `variant_plot(positions, ref_alleles, alt_alleles, frequencies=None, *, ax=None, title='Variant Plot', **kwargs)`

Create a variant visualization plot.

**Example:**
```python
from metainformant.visualization import variant_plot

ax = variant_plot(
    [1000, 2000, 3000],
    ['A', 'T', 'G'],
    ['G', 'C', 'A'],
    [0.1, 0.2, 0.05]
)
```

