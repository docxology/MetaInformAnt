# Quality Control Plots

Quality control visualization functions for quality control metrics including QC metrics plots, quality score plots, per-base quality plots, adapter content plots, and sequence length distributions.

## Functions

### `qc_metrics_plot(metrics, *, ncols=3, figsize=None, **kwargs)`

Plot multiple quality control metrics.

**Example:**
```python
from metainformant.visualization import qc_metrics_plot
import numpy as np

metrics = {
    'total_counts': np.random.poisson(1000, 100),
    'n_genes': np.random.poisson(2000, 100),
    'pct_mt': np.random.uniform(0, 10, 100)
}
fig = qc_metrics_plot(metrics)
```

### `quality_score_plot(quality_scores, *, ax=None, title='Quality Score Distribution', **kwargs)`

Plot distribution of quality scores.

### `per_base_quality_plot(positions, quality_scores, *, ax=None, title='Per-Base Quality Scores', **kwargs)`

Plot quality scores across read positions.

### `adapter_content_plot(positions, adapter_content, *, threshold=0.1, ax=None, title='Adapter Content', **kwargs)`

Plot adapter content across read positions.

### `sequence_length_distribution(lengths, *, ax=None, title='Sequence Length Distribution', **kwargs)`

Plot distribution of sequence lengths.

