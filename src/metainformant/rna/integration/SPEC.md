# RNA Integration Adapters Specification

## Purpose

Hold the RNA domain's cross-domain imports in one place so the package-boundary
check (`scripts/quality/check_module_boundaries.py`, enforced by
`tests/core/test_module_boundaries.py`) can admit them explicitly.

## Contract

- Modules here may import other `metainformant` domains; RNA engine, retrieval
  and analysis modules import those dependencies from here instead.
- Adapters re-export or wrap the target domain's public API without changing
  its behavior, signatures or outputs.
- `campaign_figures.plot_acquisition_observation` is the visualization
  function itself, so observation bundles render identically through the
  adapter.
