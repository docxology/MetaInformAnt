# AGENTS.md — `MetaInformAnt/projects/template_bioinformatics_project/scripts`

Verified against disk 2026-08-30 (doc-realization fleet pass). Thin orchestrator scripts: `01_process_data.py`, `02_analyze_results.py`, `03_visualize.py`, plus `99_create_synthetic_data.py` for deterministic synthetic input.

## Layout

- `01_process_data.py`
- `02_analyze_results.py`
- `03_visualize.py`
- `99_create_synthetic_data.py`

## Gotchas

- Run in numeric order; 99 exists so the pipeline works without real data.
