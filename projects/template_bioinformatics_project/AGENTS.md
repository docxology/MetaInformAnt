# AI Agents Documentation — Template Bioinformatics Project

This project follows the **MetaInformAnt Thin Orchestration Pattern**. All heavy computational logic resides in `metainformant.*` library modules. Scripts in `scripts/` are configurable, idempotent wrappers that delegate to these modules.

## Agent Responsibilities

### Code Assistant Agent

- Enforces the **Thin Orchestration Pattern**: scripts call `metainformant.core`, `metainformant.dna`, `metainformant.rna`, etc. — no reimplementing algorithms inline.
- Validates correct API usage: function signatures, parameter names, and return types.
- Ensures all file outputs route to `data/`, `results/`, or `logs/` — **never** the project root.
- Maintains **Zero-Mock** compliance: all tests exercise real code against real or synthetic data.
- Ensures every script has: argparse CLI with `--config`, idempotency guard, timed logging, and graceful error handling.

### Documentation Agent

- Maintains `README.md` and `AGENTS.md` at every directory level.
- Keeps per-stage docs in `doc/stages/` synchronized with script implementations.
- Ensures `doc/index.md`, `doc/architecture.md`, `doc/configuration.md`, and `doc/data.md` reflect the current pipeline.

### Workflow Orchestrator Agent

- Manages stage ordering and dependencies in `run.sh` and `main.py`.
- Ensures idempotency (skip-if-output-exists) logic in every script.
- Validates that `config/default.yaml` covers all parameters used by scripts.

## Pipeline Stages

| Stage | Script | Key MetaInformAnt Module |
| :--- | :--- | :--- |
| 1 | `scripts/01_process_data.py` | `metainformant.core.io`, `metainformant.core.paths` |
| 2 | `scripts/02_analyze_results.py` | `metainformant.core.io`, `metainformant.core.validation` |
| 3 | `scripts/03_visualize.py` | `metainformant.core.io` |
| 99 | `scripts/99_create_synthetic_data.py` | stdlib + pyyaml |

## Key Standards

- **No hardcoded paths**: all paths come from `config/default.yaml`.
- **Immutability**: `data/raw/` is read-only; all mutations go to `data/processed/`.
- **Idempotency**: every stage checks for existing outputs and skips if present.
- **Logging**: every script writes structured logs to `logs/`.
- **Testing**: `tests/` contains real end-to-end pipeline tests (Zero-Mock).

## Documentation

- [Complete Documentation](doc/index.md)
- [Architecture Guide](doc/architecture.md)
- [Configuration Reference](doc/configuration.md)
- [Data Management](doc/data.md)
- [Per-Stage Docs](doc/stages/)
