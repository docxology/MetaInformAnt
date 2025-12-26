# Expression Plots

Expression analysis visualization functions including expression heatmaps, enrichment plots, gene expression plots, differential expression plots, and log fold change visualizations.

## Functions

### `expression_heatmap(data, *, row_cluster=True, col_cluster=True, cmap='RdYlBu_r', ax=None, **kwargs)`

Create an expression heatmap with clustering.

**Example:**
```python
from metainformant.visualization import expression_heatmap
import pandas as pd
import numpy as np

data = pd.DataFrame(np.random.random((10, 5)))
ax = expression_heatmap(data)
```

### `enrichment_plot(data, x_col, y_col, *, p_threshold=0.05, ax=None, **kwargs)`

Create an enrichment plot for pathway/gene set analysis.

### `gene_expression_plot(gene_name, expression_data, sample_groups=None, *, ax=None, title=None, **kwargs)`

Plot expression levels for a single gene across samples.

### `differential_expression_plot(data, gene_col='gene', log2fc_col='log2fc', pvalue_col='pvalue', *, top_n=20, ax=None, **kwargs)`

Plot top differentially expressed genes.

### `log_fold_change_plot(data, log2fc_col='log2fc', group_col=None, *, ax=None, title='Log2 Fold Change Distribution', **kwargs)`

Plot distribution of log fold changes.

