# Specification: quality

## Scope

`config/quality/` holds YAML/JSON configuration for the quality module (QC metrics, contamination detection — see `src/metainformant/quality/`). Currently the directory contains no YAML config file; it exists as a documentation slot with these files: `AGENTS.md`, `README.md`, `PAI.md`, `SPEC.md`, and `mypy_error_budget.txt` (a 4-byte mypy error counter used by typed-ci workflows).

## Architecture

- **Component**: Configuration under the `config/` domain layer (see root `SPEC.md`: "config/ — YAML configuration templates").
- **Config loading**: Shared helper `metainformant.core.utils.config` (`load_mapping_from_file` and friends) with environment-variable overrides using domain prefixes.
- **Consumer**: The quality module reads config through `metainformant.core`, never parses YAML directly.

## Data Structures

- **Format**: YAML (plain mapping). No schema model specific to quality exists yet; when a config file is added it should mirror the section style of sibling templates (e.g. `config/singlecell/singlecell_template.yaml`: `input`, `qc`, `normalization`, `output`, `performance`).
- `mypy_error_budget.txt`: plain text integer.

## API Definition

This directory exposes no code API. Contract:

- New quality configs must be loadable via `metainformant.core.utils.config.load_mapping_from_file("config/quality/<name>.yaml")` and validate against the quality module's schema before committing (rule from `AGENTS.md`: "Validate with schema before committing new configs").
- Environment overrides follow the repo-wide prefix convention documented in root `SPEC.md` § "Configuration with Environment Overrides".
- Per `AGENTS.md`: dependency management via `uv` only; REAL IMPLEMENTATION policy applies — tests use real config files, never mocks.
