# AGENTS.md — `MetaInformAnt/tests/visualization/`

Tests for the visualization module (`src/metainformant/visualization/`): 27 test modules plus `__init__.py`.

## Layout

| Area | Test modules |
|------|--------------|
| `plots/` | `test_visualization_basic.py`, `test_visualization_animations.py`, `test_visualization_multidim.py`, `test_visualization_specialized.py`, `test_cross_species.py`, `test_sensitivity_endpoints.py` |
| `analysis/` | `test_visualization_quality.py`, `test_visualization_statistical.py`, `test_visualization_dimred.py`, `test_visualization_information.py`, `test_visualization_timeseries.py` |
| `genomics/` | `test_visualization_genomics.py`, `test_visualization_expression.py`, `test_visualization_networks.py`, `test_visualization_trees.py`, `test_visualization_phylo.py`, `test_circular_tree_layout.py` |
| `config/` and export | `test_conventions.py`, `test_palette_single_source.py`, `test_deterministic_export_formats.py`, `test_subfigure_axes_export.py` |
| `dashboards/` | `test_interactive_html_fallback.py`, `test_interactive_scatter_fallback.py` |
| Whole module | `test_visualization.py`, `test_visualization_comprehensive.py` |
| Defect regressions | `test_visualization_defect_fixes.py`, `test_visualization_followup_fixes.py` |

The two defect-regression modules read the data that landed on the artists (line data, collection offsets and sizes, colormaps, patch geometry, axes counts) rather than only the return value. They cover:

- keyword arguments that override a function's defaults (`cmap`, `marker`, `levels` and the `networkx.draw` options of the information plots; `color`, `linestyle` and `alpha` of the forecast and trend plots; `width`, `bottom` and `alpha` of the circular bar plot);
- the node indices of the Sankey links and the first-seen order of the matplotlib Sankey fallback;
- the circular Manhattan plot (`s` and `alpha` on every chromosome, single-position chromosomes at the middle of their sector);
- the hierarchical network layout (columns by distance from the nearest root, unreachable nodes in a final column);
- the optional recombination-rate track of `genomics.regional_plot`;
- `multi_panel` with a single panel.

Plots from other domains are tested beside their domain: GWAS plots in `tests/gwas/test_gwas_visualization_*.py`, DNA population plots in `tests/dna/test_dna_population_viz*.py` and `tests/dna/test_dna_population_summary_plot.py`, and the protein sequence logo in `tests/protein/test_protein_sequence_logo.py`.

## Conventions

- Real implementations with small deterministic data; the lexical no-mocks gate applies (no `MagicMock`/`unittest.mock`).
- Plotly is optional: tests of the Plotly-free fallbacks skip when it is installed.
- Run: `env -u VIRTUAL_ENV .venv/bin/python -m pytest -q tests/visualization -p no:cacheprovider`.

Repo-wide policy: see the repository-root `AGENTS.md`.
