# Agent Directives: scripts/rna

## Role

RNA-seq analysis and amalgkit workflow scripts.

## Key Scripts

- `run_all_species.py` - Multi-species streaming producer (the campaign launcher is `projects/hymenoptera_amalgkit/scripts/run_full_campaign.sh`)
- `process_species.py` - Single-species producer using the shared streaming engine
- `check_environment.py` - Verify environment setup
- `check_pipeline_status.py` - Run-state and downstream evidence checker
- `verify_rna.py` - Validate RNA module functionality
- `filter_valid_samples.py` - Filter samples for processing
- `validate_all_species_workflow.py` - Cross-species validation
- `validate_configs.py` - Config YAML validation
- `_setup_utils.py` - Shared setup utilities
- `install_r_packages.R` - R dependency installation
- `install_r_packages.sh` - R package installation script
- `run_rna_tests.sh` - RNA test runner
- `prepare_cohort_host.py` - Prepare a fresh data root for frozen-cohort quantification on a new host
- `import_quantifications.py` - Import verified quantifications from other data roots
- `export_quant_snapshot.py` - Export verified, reusable quantification outputs to a portable snapshot
- `reclaim_quantified_raw.py` - Reclaim raw inputs only for verified, reusable quantification outputs
- `reconcile_quantification_compatibility.py` - Classify and safely restore existing quantification outputs
- `campaign_watcher.py` - Detached, read-only campaign watcher
- `benchmark_campaign_scale.py` - Campaign-scale benchmark on disposable fixtures
- `new_clade_scaffold.py` - Scaffold a sibling clade project from the Hymenoptera template pattern

## ENA Contributions
- **ENA-First Strategy**: `retrieval.ena_downloader.ENADownloader` and the streaming producer download from ENA first.
- **Failover**: Automatic failover to `fasterq-dump` from NCBI SRA.

## Usage

```bash
# Inspect the multi-species plan
uv run python scripts/rna/run_all_species.py \
  --config-dir projects/hymenoptera_amalgkit/config/amalgkit \
  --data-root "$AMALGKIT_DATA_ROOT" --dry-run

# Run one cohort species (diagnostic; production uses the project launcher)
env -u VIRTUAL_ENV uv run --frozen --no-sync python scripts/rna/process_species.py \
  --species pogonomyrmex_barbatus \
  --config-dir projects/hymenoptera_amalgkit/config/amalgkit \
  --data-root "$AMALGKIT_DATA_ROOT" \
  --cohort "$AMALGKIT_DATA_ROOT/cohort/hymenoptera27-20261005"

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
