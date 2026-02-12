# Visualization Analysis

Analytical visualization functions for dimensionality reduction, quality control, information theory, statistical diagnostics, and time series.

## Contents

| File | Purpose |
|------|---------|
| `dimred.py` | PCA, UMAP, t-SNE scatter plots and biplots |
| `information.py` | Entropy profiles, mutual information matrices, Renyi spectra |
| `quality.py` | FASTQ quality metrics, GC distribution, adapter content plots |
| `quality_assessment.py` | Coverage uniformity, error profiles, batch effect QC |
| `quality_omics.py` | VCF, single-cell, protein structure, multi-omics quality plots |
| `quality_sequencing.py` | Per-base quality, duplication levels, k-mer profiles |
| `statistical.py` | Histogram, boxplot, violin, QQ, ROC, correlation heatmap |
| `timeseries.py` | Time series, autocorrelation, seasonal decomposition, forecast |

## Key Functions

| Function | Description |
|----------|-------------|
| `plot_pca()` | PCA scatter plot with optional grouping |
| `plot_umap()` | UMAP embedding visualization |
| `plot_entropy_profile()` | Positional entropy across sequence or features |
| `plot_mutual_information_matrix()` | Pairwise MI heatmap |
| `plot_quality_metrics()` | Multi-panel FASTQ quality summary |
| `histogram()` | Statistical histogram with optional density overlay |
| `violin_plot()` | Violin plot for distribution comparison |
| `plot_time_series()` | Time series line plot with annotations |

## Usage

```python
from metainformant.visualization.analysis.dimred import plot_pca
from metainformant.visualization.analysis.statistical import histogram

plot_pca(data, color_by="group", output_path="output/pca.png")
histogram(values, bins=50, output_path="output/hist.png")
```
