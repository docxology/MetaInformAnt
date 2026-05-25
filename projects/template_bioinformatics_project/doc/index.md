# Template Bioinformatics Project — Documentation

Welcome to the documentation for the **MetaInformAnt Template Bioinformatics Project** — the canonical starting point for new standalone bioinformatics projects within the MetaInformAnt ecosystem.

> [!NOTE]
> This project follows the **Thin Orchestration Pattern**. All heavy computational logic lives in `metainformant.*` library modules. Scripts in `scripts/` are configurable, idempotent wrappers that delegate to these modules.

## Documentation Map

| Document | Description |
| :--- | :--- |
| [Architecture Guide](architecture.md) | Thin Orchestration Pattern, MetaInformAnt module mapping, extension guidelines |
| [Configuration Reference](configuration.md) | Every key in `config/default.yaml`, with descriptions and defaults |
| [Data Management](data.md) | Directory conventions, file formats, immutability rules, provenance |

### Per-Stage References

| Stage | Script | Doc |
| :--- | :--- | :--- |
| 1 — Process Data | `scripts/01_process_data.py` | [stages/01_process_data.md](stages/01_process_data.md) |
| 2 — Analyze Results | `scripts/02_analyze_results.py` | [stages/02_analyze_results.md](stages/02_analyze_results.md) |
| 3 — Visualise | `scripts/03_visualize.py` | [stages/03_visualize.md](stages/03_visualize.md) |
| 99 — Synthetic Data | `scripts/99_create_synthetic_data.py` | [stages/99_synthetic_data.md](stages/99_synthetic_data.md) |

## Project Layout

```text
template_bioinformatics_project/
 AGENTS.md # AI agent responsibilities & standards
 SPEC.md # Living project specification
 README.md # User-facing quickstart
 pyproject.toml # uv-compatible project manifest
 main.py # Top-level CLI pipeline runner
 run.sh # Shell pipeline runner (stage-selectable)
 config/
 default.yaml # Centralized configuration
 data/
 raw/ # Immutable raw inputs
 processed/ # Stage 1 outputs
 scripts/ # Numbered thin-orchestrator scripts
 results/
 figures/ # Stage 3 plots
 tables/ # Stage 2 summary tables
 logs/ # Per-stage execution logs
 tests/ # Real-Implementation end-to-end tests
 doc/ # This documentation suite
 stages/ # Per-stage technical references
```

## Prerequisites

| Tool | Version | Purpose |
| :--- | :--- | :--- |
| `uv` | ≥ 0.4 | Python environment and dependency management |
| Python | ≥ 3.11 | Runtime |
| `metainformant` | latest | Core library (when delegating beyond stdlib) |

Install `uv`:
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

## Quick Start

```bash
# 1. Install project dependencies
cd projects/template_bioinformatics_project
uv sync

# 2. Generate synthetic test data
uv run scripts/99_create_synthetic_data.py

# 3. Run the full pipeline
./run.sh
# or equivalently:
uv run main.py --stage all

# 4. Run individual stages
./run.sh --stage process
./run.sh --stage analyze
./run.sh --stage visualize

# 5. Use a custom config
./run.sh --config config/my_experiment.yaml

# 6. Clean all generated outputs
./run.sh --clean

# 7. Run tests
uv run pytest tests/ -v
```

## Pipeline Execution Flow

```text
 Starting pipeline (stage=all)

==> Stage: Data Processing   (01)   → data/processed/processed_data.csv
==> Stage: Results Analysis  (02)   → results/tables/summary_statistics.csv
==> Stage: Visualization     (03)   → results/figures/distribution_grid.png
                                       results/figures/correlation_heatmap.png

 Pipeline finished successfully.
 Raw data: data/raw/
 Processed: data/processed/
 Results: results/
 Logs: logs/
```

Each stage is **idempotent**: re-running any stage after completion safely skips redundant work.  Use `--force` to unconditionally reprocess.

## Output Routing

All output paths are defined in `config/default.yaml` under `paths`. No files should ever be written to the project root.

| Config Key | Default | Contents |
| :--- | :--- | :--- |
| `paths.data_raw` | `data/raw/` | Immutable raw inputs (CSV, YAML metadata) |
| `paths.data_processed` | `data/processed/` | Cleaned, normalised data from Stage 1 |
| `paths.results_tables` | `results/tables/` | Summary statistics, correlation matrix (Stage 2) |
| `paths.results_figures` | `results/figures/` | Matplotlib figures (Stage 3) |
| `paths.logs` | `logs/` | Per-stage execution logs |

> [!IMPORTANT]
> If you observe stray `.csv` or `.json` files in the project root, a script is writing outside the configured paths. Check the `paths` section in your config.

## Troubleshooting

| Symptom | Cause | Fix |
| :--- | :--- | :--- |
| `No CSV files found in data/raw/` | Raw data not generated | Run `uv run scripts/99_create_synthetic_data.py` first |
| Stage 1 output already exists | Idempotency guard activated | Use `--force` flag to reprocess |
| Stage 2 exits non-zero | Stage 1 not completed | Run Stage 1 first |
| `uv: command not found` | uv not installed | `curl -LsSf https://astral.sh/uv/install.sh \| sh` |
| Tests fail with import errors | Dependencies not installed | `uv sync --extra dev` |
| Figures look wrong | Wrong `color_palette` or `dpi` | Adjust `visualization` section in config |
