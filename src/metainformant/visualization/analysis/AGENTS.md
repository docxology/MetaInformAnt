# AI Agents: Statistical Analysis Visualization

Function index and technical details for the `visualization.analysis` module.

## Function Index

### `metainformant.visualization.analysis.statistical`
- **`histogram(data: np.ndarray, bins: int = 30) -> Axes`**: Plots a frequency distribution.
- **`box_plot(data: list[np.ndarray], labels: list[str]) -> Axes`**: Standard box-and-whisker plot.
- **`correlation_heatmap(corr_matrix: pd.DataFrame) -> Axes`**: Visualizes correlation matrices.

### `metainformant.visualization.analysis.dimred`
- **`plot_pca(pca_results: np.ndarray, labels: list[int]) -> Axes`**: Plots PCA results.
- **`plot_umap(umap_results: np.ndarray, labels: list[int]) -> Axes`**: Plots UMAP results.

### `metainformant.visualization.analysis.quality`
- **`plot_quality_metrics(metrics: dict[str, Any]) -> Axes`**: Summarizes quality control results.
- **`plot_gc_distribution(seqs: list[str]) -> Axes`**: Plots GC content distribution.
