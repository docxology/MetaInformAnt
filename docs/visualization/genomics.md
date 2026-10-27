# Genomics Plots

Genomic visualization functions for GWAS and sequence analysis including Manhattan plots, volcano plots, regional plots, circular Manhattan plots, a human chromosome overview, coverage plots, and variant position plots. They live in `metainformant.visualization.genomics.genomics`, return the `matplotlib.axes.Axes`, and accept `ax=` (draw on existing axes) and `output_path=` (save the figure at 300 dpi).

## Functions

### `manhattan_plot(data, *, chr_col='CHR', pos_col='BP', pval_col='P', ax=None, output_path=None, **kwargs)`

Create a Manhattan plot for genome-wide association studies. `pval_col` holds raw p-values (the
function applies -log10 itself) and the genome-wide significance line is drawn at 5e-8.

**Example:**
```python
from metainformant.visualization.genomics.genomics import manhattan_plot
import pandas as pd

data = pd.DataFrame({
    'position': [1000, 2000, 3000],
    'pvalue': [1e-6, 1e-9, 0.01],
    'chromosome': ['chr1', 'chr1', 'chr2']
})
ax = manhattan_plot(data, chr_col='chromosome', pos_col='position', pval_col='pvalue')
```

### `volcano_plot(data, *, log2fc_col='log2FoldChange', pval_col='padj', ax=None, output_path=None, **kwargs)`

Create a volcano plot for differential expression analysis. `pval_col` holds raw (or adjusted)
p-values, which the function converts with -log10; points are colored using fixed thresholds of
p = 0.05 and |log2 fold change| = 1.

**Example:**
```python
from metainformant.visualization.genomics.genomics import volcano_plot
import pandas as pd

data = pd.DataFrame({
    'log2fc': [-2, 1, 0.5, -1.5],
    'pvalue': [0.001, 0.01, 0.5, 0.001]
})
ax = volcano_plot(data, log2fc_col='log2fc', pval_col='pvalue')
```

### `regional_plot(data, *, chr, start, end, ax=None, output_path=None, **kwargs)`

Create a regional association plot for the window `start`-`end` (base pairs) on chromosome `chr`. `chr`, `start` and `end` are keyword-only, and `start` must be smaller than `end`. `data` is a DataFrame with the fixed column names `CHR`, `BP` and `P` (raw p-values, converted with -log10); there are no column-name options. A dashed line marks genome-wide significance (5e-8), and a window with no rows draws a "No data in region" message instead of points. Extra keyword arguments (`s`, `alpha`, ...) are passed to `Axes.scatter`.

**Example:**
```python
from metainformant.visualization.genomics.genomics import regional_plot
import numpy as np
import pandas as pd

rng = np.random.default_rng(1)
data = pd.DataFrame({
    'CHR': ['chr1'] * 200,
    'BP': np.sort(rng.integers(1_000_000, 1_100_000, 200)),
    'P': rng.uniform(1e-8, 1, 200),
})
ax = regional_plot(data, chr='chr1', start=1_000_000, end=1_100_000)
```

### `circular_manhattan_plot(data, *, ax=None, output_path=None, **kwargs)`

Create a circular Manhattan plot. Chromosomes (sorted by name) share the circle equally, each point is placed by its position within its chromosome and plotted at a radius of -log10(p). `data` needs the fixed column names `CHR`, `BP` and `P` (raw p-values). The function creates polar axes; an `ax` you pass in must be a polar axes. Extra keyword arguments are passed to `Axes.scatter`.

**Example:**
```python
from metainformant.visualization.genomics.genomics import circular_manhattan_plot
import numpy as np
import pandas as pd

rng = np.random.default_rng(1)
data = pd.DataFrame({
    'CHR': np.repeat(['chr1', 'chr2', 'chr3'], 100),
    'BP': np.tile(np.arange(100) * 10_000, 3),
    'P': rng.uniform(1e-8, 1, 300),
})
ax = circular_manhattan_plot(data)
```

### `chromosome_ideogram(*, ax=None, output_path=None, **kwargs)`

Draw an overview of the 24 human chromosomes (1-22, X, Y) as horizontal bars scaled to the longest chromosome, using built-in GRCh38 lengths. The function takes no data: it cannot mark positions or values, and `figsize` is the only keyword argument it reads from `**kwargs`. To show per-position results use `manhattan_plot` or `circular_manhattan_plot`.

**Example:**
```python
from metainformant.visualization.genomics.genomics import chromosome_ideogram

ax = chromosome_ideogram(figsize=(10, 6))
```

### `coverage_plot(coverage, positions, *, ax=None, output_path=None, **kwargs)`

Create a coverage plot for sequencing data: a line with the area under it filled. The argument order is coverage depth first, then genomic positions; both must be NumPy arrays of equal length. Extra keyword arguments are passed to both `Axes.plot` and `Axes.fill_between` (the fill defaults to `alpha=0.3`).

**Example:**
```python
from metainformant.visualization.genomics.genomics import coverage_plot
import numpy as np

positions = np.arange(1000, 2000, 10)
coverage = np.random.default_rng(0).poisson(50, len(positions))
ax = coverage_plot(coverage, positions)
```

### `variant_plot(variants, *, ax=None, output_path=None, **kwargs)`

Plot variant positions as points on a single track with a small fixed vertical jitter so nearby variants stay visible. `variants` is a DataFrame that must contain a `POS` column; other columns (alleles, frequencies) are not used. Extra keyword arguments (`s`, `alpha`, ...) are passed to `Axes.scatter`.

**Example:**
```python
from metainformant.visualization.genomics.genomics import variant_plot
import pandas as pd

variants = pd.DataFrame({
    'POS': [1000, 2000, 3000],
    'REF': ['A', 'T', 'G'],
    'ALT': ['G', 'C', 'A'],
})
ax = variant_plot(variants)
```
