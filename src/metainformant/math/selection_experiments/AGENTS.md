# AI Agents in Selection Experiment Development

This document captures AI assistance for the natural selection simulation suite.

## AI Contributions

### Model Implementation
**Code Assistant Agent** translated the theoretical framework into code by:
- Implementing `model.py` simulation kernels for structural/quality trait dynamics
- Encoding Price equation components and signal fidelity transformations
- Building stochastic experiment runners used by the CLI entry points

### CLI and Visualization
**Code Assistant Agent** authored:
- `cli.py` commands (`replay`, `abstract`) with uv-friendly argument parsing
- Plotting utilities in `plotting.py` that generate publication-quality PNGs
- Reproducibility hooks (deterministic seeds, output path controls)

### Documentation and Context
**Documentation Agent** assisted with:
- The theoretical background and references in `README.md`
- Experiment-by-experiment commentary and usage instructions
- Cross-links to related math module resources

## Maintenance
- Validate plots and metrics against the cited literature when parameters change
- Synchronize README tables with available CLI flags and outputs
- Run real simulations after significant refactors to confirm reproducibility guarantees


