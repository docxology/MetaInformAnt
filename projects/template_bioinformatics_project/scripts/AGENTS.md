# AGENTS.md — `scripts/` — Thin Orchestrator Scripts

Technical specification for this project's pipeline scripts (verified 2026-09-06:
4 files, 0 subdirs). This is a standalone template project nested inside
MetaInformAnt/projects/ — see its root `AGENTS.md` and `SPEC.md` for the
project's own rules, and `../doc/architecture.md` for the pattern rationale.

## Script Inventory

| Script | Pattern | Delegates to | Outputs |
|--------|---------|--------------|---------|
| `01_process_data.py` | Stage 01 Thin Orchestrator | `src/template_bioinformatics_project/processing.py` | `data/processed/processed_data.csv`, `logs/01_process_data.log` |
| `02_analyze_results.py` | Stage 02 Thin Orchestrator | `src/template_bioinformatics_project/analysis.py` | `results/tables/summary_statistics.csv`, `results/tables/correlation_matrix.csv` (summary/correlation), `results/tables/pca_summary.json` (pca), `results/tables/analysis_metadata.json`, `logs/02_analyze_results.log` |
| `03_visualize.py` | Stage 03 Thin Orchestrator | `src/template_bioinformatics_project/visualization.py` | `results/figures/distribution_grid.{fmt}`, `results/figures/correlation_heatmap.{fmt}` (if Stage 2 wrote a correlation matrix), `logs/03_visualize.log` |
| `99_create_synthetic_data.py` | Stage 99 Thin Orchestrator | `src/template_bioinformatics_project/synthetic.py` | `data/raw/samples_A.csv`, `data/raw/samples_B.csv`, `data/raw/metadata.yaml` |

## Design Contract

- Scripts are orchestration only: argparse, `sys.path` bootstrap, logging, and
  one delegated call into `src/template_bioinformatics_project/`.
- No domain algorithms, data wrangling, statistics, or plotting inline. No
  reimplementing `metainformant.*` or library logic — if business logic is
  needed, it goes into the appropriate `src/template_bioinformatics_project/`
  module and is imported.
- No hardcoded paths: every path and threshold routes through
  `config/default.yaml` (loaded via `common.load_config` /
  `common.load_optional_config`).
- Every stage script: `--config` + `--force` flags, idempotency guard on
  existing outputs, structured log to `logs/`, graceful `sys.exit(1)` on fatal
  errors.
- All reusable script behavior must be covered by tests in `../tests/`:
  `test_library.py` unit-tests the src modules directly;
  `test_pipeline.py` exercises the scripts end-to-end via subprocess
  (Real-Implementation policy, no mocks).
