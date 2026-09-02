# Agent Directives: scripts/popgen

## Role
Thin orchestrator scripts for the population genetics workflow
(generate -> analyze -> report -> visualize).

## Contents
- `analysis.py` - end-to-end CLI entry point
- `analyze.py` - analysis orchestration helper
- `generate_dataset.py` - synthetic dataset generation
- `report.py` - markdown report writer
- `visualize.py` - plot generation

## Rules
- Scripts are thin wrappers: analysis business logic lives in
  `metainformant.popgen`; generators in
  `metainformant.simulation.models.popgen`
- Keep import paths in sync with the actual module homes (this directory
  previously imported from defunct `metainformant.dna.population_analysis`
  and `metainformant.simulation.popgen` paths)
- Follow REAL IMPLEMENTATION policy; use `uv` for Python dependencies
