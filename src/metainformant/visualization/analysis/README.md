# Visualization Analysis

Analytical visualization functions for dimensionality reduction, quality control, information theory, statistical diagnostics, and time series.

## Contents

| File | Purpose |
|------|---------|
| `dimred.py` | PCA, UMAP, t-SNE scatter plots and biplots |
| `information.py` | Entropy profiles, mutual information matrices, Renyi spectra |
| `math_plots.py` | Population genetics and mathematical biology plots (allele frequency spectra, coalescent trees, epidemic models, evolutionary trajectories) |
| `quality.py` | FASTQ quality metrics, GC distribution, adapter content plots |
| `quality_assessment.py` | Coverage uniformity, error profiles, batch effect QC |
| `quality_omics.py` | VCF, single-cell, protein structure, multi-omics quality plots |
| `quality_sequencing.py` | Per-base quality, duplication levels, k-mer profiles |
| `simulation_plots.py` | Simulation result plots and animations (sequence evolution, population dynamics, agent-based models, sensitivity analysis) |
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
| `plot_allele_frequency_spectrum()` | Allele frequency spectrum histogram |
| `plot_coalescent_tree()` | Coalescent tree visualization |
| `plot_epidemic_model_simulation()` | SIR/SEIR epidemic model curves |
| `plot_sequence_evolution()` | Sequence evolution over generations |
| `plot_population_dynamics_simulation()` | Population dynamics simulation results |
| `animate_population_dynamics()` | Population dynamics animation |

## Usage

```python
from metainformant.visualization.analysis.dimred import plot_pca
from metainformant.visualization.analysis.statistical import histogram
from metainformant.visualization.analysis.math_plots import plot_allele_frequency_spectrum
from metainformant.visualization.analysis.simulation_plots import plot_sequence_evolution

plot_pca(data, color_by="group", output_path="output/pca.png")
histogram(values, bins=50, output_path="output/hist.png")
plot_allele_frequency_spectrum(freqs, output_path="output/afs.png")
plot_sequence_evolution(sequences, output_path="output/seq_evolution.png")
```
