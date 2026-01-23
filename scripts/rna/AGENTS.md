# Agent Directives: scripts/rna

## Role
RNA-seq analysis and amalgkit workflow scripts.

## Key Scripts
- `run_workflow.py` - Main amalgkit workflow runner
- `run_workflow_tui.py` - TUI-based workflow runner
- `check_environment.py` - Verify environment setup
- `verify_rna.py` - Validate RNA module functionality
- `verify_fallback.py` - Test fallback mechanisms
- `setup_genome.py` - Genome preparation for quantification
- `discover_species.py` - Discover available species
- `filter_valid_samples.py` - Filter samples for processing
- `process_samples_sequential.py` - Sequential sample processing
- `recover_missing_batch.py` - Batch recovery for failed samples
- `recover_missing_parallel.py` - Parallel recovery processing
- `validate_all_species_workflow.py` - Cross-species validation
- `_setup_utils.py` - Shared setup utilities
- `install_r_deps.R` - R dependency installation
- `install_r_packages.sh` - R package installation script
- `run_rna_tests.sh` - RNA test runner

## Usage
```bash
# Run workflow
uv run python scripts/rna/run_workflow.py --config config/amalgkit/species.yaml

# Check environment
uv run python scripts/rna/check_environment.py
```
