# scripts

Thin orchestrator scripts for the template bioinformatics pipeline. All pipeline
logic lives in `src/template_bioinformatics_project/` — scripts only parse
arguments, bootstrap `sys.path`, configure logging, and make one delegated call.
See `../AGENTS.md`, `../doc/architecture.md`, and `AGENTS.md` (this directory).

## Inventory

| Script | Purpose | Delegates to | Command |
|--------|---------|--------------|---------|
| `01_process_data.py` | Stage 1 — ingest raw CSVs, filter, z-score normalise | `template_bioinformatics_project.processing.process_data` | `uv run scripts/01_process_data.py --config config/default.yaml [--force]` |
| `02_analyze_results.py` | Stage 2 — summary statistics, correlation, optional PCA | `template_bioinformatics_project.analysis.run_analysis` | `uv run scripts/02_analyze_results.py --config config/default.yaml [--force]` |
| `03_visualize.py` | Stage 3 — distribution grid + correlation heatmap figures | `template_bioinformatics_project.visualization.run_visualize` | `uv run scripts/03_visualize.py --config config/default.yaml [--force]` |
| `99_create_synthetic_data.py` | Generate synthetic raw CSVs + provenance metadata for testing | `template_bioinformatics_project.synthetic.run` | `uv run scripts/99_create_synthetic_data.py [--n-samples 200] [--seed 2026] [--force]` |

All scripts accept `--config` (default `config/default.yaml`); stage scripts
additionally accept `--force` to bypass the idempotency guard. Stage 99 also
accepts `--n-samples`, `--n-features`, and `--seed`.

## Module logic (src/template_bioinformatics_project/)

- `processing.py` — Stage 1: CSV discovery, concatenation, missing-value and
  row-count filters, z-score normalisation
- `analysis.py` — Stage 2: `compute_summary_statistics`, `compute_correlation_matrix`,
  `compute_pca_summary`, orchestration via `run_analysis`
- `visualization.py` — Stage 3: `plot_distribution_grid`, `plot_correlation_heatmap`,
  orchestration via `run_visualize`
- `synthetic.py` — Stage 99: `generate_sample_dataframe`, `write_csv`,
  `write_metadata`, orchestration via `run`
- `common.py` — shared `load_config` / `load_optional_config` YAML helpers

## Outputs

All paths resolve from `config/default.yaml`: Stage 1 → `data/processed/`,
Stage 2 → `results/tables/` (`summary_statistics.csv`, `correlation_matrix.csv`
for `summary`/`correlation` methods, `pca_summary.json` for `pca`,
`analysis_metadata.json`), Stage 3 → `results/figures/`, all stages → `logs/`.
