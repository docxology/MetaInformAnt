# VISUALIZATION

## Overview
Command-line helpers for GWAS Visualization workflows. Scripts should remain thin wrappers around `src/metainformant/` implementations and be run from the repository root with `uv`.

## Contents
- [generate_missing_plots.py](generate_missing_plots.py) — regenerates the PCA scree plot from saved PCA results

The former `visualizations.py` module (kinship heatmap, PCA scatter/scree
plots) moved to `src/metainformant/gwas/visualization/structure_plots.py`.

## Usage
```bash
uv run python scripts/gwas/visualization/generate_missing_plots.py --help
```
