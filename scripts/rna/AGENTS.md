# Agent Directives: scripts/rna

## Role

RNA-seq analysis and amalgkit workflow scripts.

## Key Scripts

- `run_all_species.sh` - **Primary pipeline** — runs all species sequentially with per-sample concurrency
- `run_workflow.py` - Single-species amalgkit workflow orchestrator
- `run_workflow_tui.py` - TUI-based workflow runner
- `check_environment.py` - Verify environment setup
- `check_pipeline_status.py` - Per-species progress dashboard
- `verify_rna.py` - Validate RNA module functionality
- `setup_genome.py` - Genome preparation for quantification
- `discover_species.py` - Discover available species
- `filter_valid_samples.py` - Filter samples for processing
- `validate_all_species_workflow.py` - Cross-species validation
- `validate_configs.py` - Config YAML validation
- `_setup_utils.py` - Shared setup utilities
- `install_r_deps.R` - R dependency installation
- `install_r_packages.sh` - R package installation script
- `run_rna_tests.sh` - RNA test runner

## Usage

```bash
# Run full pipeline (all species)
nohup bash scripts/rna/run_all_species.sh > output/amalgkit/run_all_species_incremental.log 2>&1 &

# Run single species
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pbarbatus.yaml --stream --chunk-size 6

# Check progress
.venv/bin/python scripts/package/generate_custom_summary.py

# Check environment
python3 scripts/rna/check_environment.py
```

## 📋 Code Quality Policy

All scripts **MUST** be:

1. **Functional** — Real implementations only. No stubs or placeholders.
2. **Modular** — Use `argparse` for all CLI arguments. No hardcoded paths.
3. **Tested** — Integration tests covering all scripts. No mocks.
4. **Documented** — Docstrings, `--help` output, and README entries.

## 🚫 NO_MOCKING_POLICY

> **NEVER use `unittest.mock`, `pytest-mock`, `MagicMock`, or `patch` in tests.**

Use real filesystem operations, real configs, and real function calls. See `tests/NO_MOCKING_POLICY.md`.
