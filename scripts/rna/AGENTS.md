# Agent Directives: scripts/rna

## Role

RNA-seq analysis and amalgkit workflow scripts.

## Key Scripts

- `run_all_species.py` - Configuration-derived multi-species launcher
- `process_species.py` - Single-species producer using the shared streaming engine
- `check_environment.py` - Verify environment setup
- `check_pipeline_status.py` - Run-state and downstream evidence checker
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

## ENA Contributions
- **ENA-First Strategy**: Robust ENA download logic (`download_ena_robust`) integrated into orchestrators.
- **Failover**: Automatic failover to `fasterq-dump` from NCBI SRA.

## Usage

```bash
# Inspect the multi-species plan
uv run python scripts/rna/run_all_species.py \
  --config-dir projects/hymenoptera_amalgkit/config/amalgkit \
  --data-root "$AMALGKIT_DATA_ROOT" --dry-run

# Run single species
uv run python scripts/rna/process_species.py \
  --species pogonomyrmex_barbatus \
  --config-dir projects/hymenoptera_amalgkit/config/amalgkit \
  --data-root "$AMALGKIT_DATA_ROOT"

# Generate current project evidence
bash projects/hymenoptera_amalgkit/scripts/verify_setup.sh \
  --data-root "$AMALGKIT_DATA_ROOT" --require-data --report

# Check environment
uv run python scripts/rna/check_environment.py
```

## Code Quality Policy

All scripts **MUST** be:

1. **Functional** — Real implementations only. No inert placeholders.
2. **Modular** — Use `argparse` for all CLI arguments. No hardcoded paths.
3. **Tested** — Integration tests covering all scripts. Real implementations.
4. **Documented** — Docstrings, `--help` output, and README entries.

## real-implementation policy

> **NEVER use `test-double modules`, `test-double plugins`, `test-double classes`, or `patch` in tests.**

Use real filesystem operations, real configs, and real function calls. See `tests/REAL_IMPLEMENTATION_TESTING_POLICY.md`.
