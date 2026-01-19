# SPEC: Script Orchestration

The `scripts/` directory contains production-ready entry points for executing complex biological workflows.

## Orchestration Standards

### 1. Thin Wrappers
Scripts must not contain complex business logic. They should be responsible for:
- Argument parsing
- Environment setup
- Calling core methods from `src/metainformant/`
- Logging and progress reporting

### 2. Modular Structure
Scripts are organized by biological domain, mirroring the `src/` directory structure.

## Discovery Logic

Automated tools (like `scripts/run_all_scripts.py`) discover orchestrators by looking for:
- Python files starting with `run_` or named `orchestrate.py`.
- Directories containing a "Triple Play" set of documentation.

## Technical Patterns

### Config-Driven Execution
Most scripts accept a `--config` argument pointing to a YAML or TOML file. This configuration is merged with defaults using `metainformant.core.config`.

### Environment Compliance
Scripts that depend on external tools (e.g., R, Amalgkit, GATK) must perform a pre-flight check using `metainformant.core.execution` to ensure the environment is correctly configured.

## TUI vs CLI
- **Production Scripts**: Should use the TUI (Terminal UI) for progress monitoring in long-running jobs.
- **Utility Scripts**: Simple command-line output for one-off tasks.
