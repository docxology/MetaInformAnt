# SPEC: Genomics Visualization

## Standards

1. **Output Format**: All plots return a `matplotlib.axes.Axes` object to allow for further customization by the caller.
2. **Data structures**: Functions are optimized to handle `pandas.DataFrame` and `numpy.ndarray` as primary inputs.
3. **Scientific Rigor**: Plot scales (e.g., -log10 p-value) and thresholds conform to standard bioinformatics practices.

## Design Decisions

- **Modularity**: Visualizations are categorized into `genomics.py`, `expression.py`, `networks.py`, and `trees.py` to keep the codebase maintainable.
- **Interactivity**: While primarily returning static plots, the functions are designed to be compatible with interactive backends where possible.
