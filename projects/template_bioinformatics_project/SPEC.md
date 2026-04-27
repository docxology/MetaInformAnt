# Project Specification — Template Bioinformatics Project

## Purpose

This repository is the **canonical template** for new bioinformatics projects within the MetaInformAnt ecosystem. When forking this template:

1. Replace project-specific config values in `config/default.yaml`.
2. Rename and adapt scripts in `scripts/` to your domain.
3. Update `README.md`, `AGENTS.md`, and `doc/` to reflect your project.
4. Add phenotype/metadata files to `data/raw/` or `data/phenotypes/`.

## Design Constraints

| Constraint | Rationale |
| :--- | :--- |
| Scripts are thin orchestrators | Algorithmic logic belongs in `metainformant.*` modules — not duplicated in every project |
| All paths from config | Prevents environment-specific breakage; ensures reproducibility |
| `data/raw/` is immutable | Raw data is the source of truth; never overwrite it |
| Idempotent stages | Re-running any stage should be safe and fast (skip existing outputs) |
| Zero-Mock tests | Tests must exercise real code; no patching or stubbing of logic |
| `uv` for environments | Reproducible, fast, PEP 723—compatible dependency management |

## Directory Layout

```text
template_bioinformatics_project/
 AGENTS.md # AI agent responsibilities and standards
 SPEC.md # This specification
 README.md # User-facing quickstart and overview
 pyproject.toml # uv-compatible project manifest
 main.py # Top-level CLI entry point
 run.sh # Master shell pipeline runner
 config/
 default.yaml # Centralized configuration (all parameters)
 data/
 raw/ # Immutable raw inputs (FASTQ, CSV, YAML metadata)
 processed/ # Intermediate outputs from Stage 1
 scripts/
 01_process_data.py # Stage 1: ingest and process raw data
 02_analyze_results.py # Stage 2: downstream statistical analysis
 03_visualize.py # Stage 3: plot generation
 99_create_synthetic_data.py # Test data generator
 results/
 figures/ # Matplotlib output figures
 logs/ # Per-stage execution logs
 tests/
 conftest.py # Shared pytest fixtures
 test_pipeline.py # End-to-end pipeline tests
 doc/
 index.md # Documentation hub
 architecture.md # Thin orchestration pattern guide
 configuration.md # Config reference
 data.md # Data management guide
 stages/ # Per-stage technical references
 01_process_data.md
 02_analyze_results.md
 03_visualize.md
 99_synthetic_data.md
```

## Extension Patterns

### Adding a New Stage

1. Create `scripts/NN_stage_name.py`.
2. Load config with `--config` argparse argument.
3. Implement idempotency: check if output exists, skip if so.
4. Delegate computation to `metainformant.*`.
5. Write structured logs to `logs/NN_stage_name.log`.
6. Register the stage in `run.sh` and `main.py`.
7. Add a doc file at `doc/stages/NN_stage_name.md`.

### Adding New Config Parameters

1. Add the key to `config/default.yaml` with a sensible default and inline comment.
2. Document it in `doc/configuration.md`.
3. Access it in scripts via `config['section']['key']`.

### Swapping in a Domain-Specific MetaInformAnt Module

Replace calls to `metainformant.core.*` with the appropriate domain module (e.g., `metainformant.dna.*`, `metainformant.rna.*`, `metainformant.gwas.*`) once your project has a concrete biological domain.

## Versioning

This template uses CalVer (`YYYY.MM.DD`) in `pyproject.toml`. When forking, reset to `0.1.0` and adopt SemVer for your project-specific releases.
