# AGENTS.md — `MetaInformAnt/tests/visualization/`

tests for the visualization module.
Files (verified 2026-09-05): 19 test modules + `__init__.py`, covering `plots/` (basic, animations, multidim, cross_species, specialized), `analysis/` (quality, statistical, dimred, information, timeseries), `genomics/` (genomics, expression, networks, trees/phylo), `config/` (conventions, palette single-source), and top-level (`test_visualization.py`, `test_visualization_comprehensive.py`).


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV .venv/bin/python -m pytest -q tests/visualization -p no:cacheprovider` (verified 2026-09-05: 312 passed, 2 skipped pre-review).
Repo-wide policy: see the repository-root `AGENTS.md`.