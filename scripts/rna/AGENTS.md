# AI Agents in RNA Script Development

This document records AI involvement in the RNA processing scripts housed here.

## Script Consolidation (December 2025)

**Code Assistant Agent** consolidated 20+ scripts into 4-5 thin orchestrators:

### Phase 1: Move Logic to src/

Created new modules in `src/metainformant/rna/`:
- **`orchestration.py`**: Workflow orchestration functions
  - `run_workflow_for_species()` - Run workflow steps for a species
  - `cleanup_unquantified_samples()` - Quantify and cleanup samples
  - `monitor_workflows()` - Real-time monitoring
  - `discover_species_configs()` - Discover species from configs

- **`cleanup.py`**: Cleanup functions
  - `cleanup_partial_downloads()` - Clean up partial downloads
  - `fix_abundance_naming_for_species()` - Fix abundance file naming

- **`discovery.py`**: Species discovery and config generation
  - `search_species_with_rnaseq()` - Search NCBI SRA for species
  - `get_genome_info()` - Get genome assembly information
  - `generate_config_yaml()` - Generate amalgkit YAML configs

Extended existing modules:
- **`genome_prep.py`**: Added orchestration functions
  - `verify_genome_status()` - Verify genome and index status
  - `orchestrate_genome_setup()` - Complete genome setup pipeline

- **`monitoring.py`**: Added assessment functions
  - `assess_all_species_progress()` - Assess progress for all species
  - `initialize_progress_tracking()` - Initialize progress tracking

### Phase 2: Create Thin Orchestrators

Created thin wrapper scripts that call src methods:
- **`run_workflow.py`**: Main workflow orchestrator
  - Calls `metainformant.rna.orchestration.run_workflow_for_species()`
  - Single-species sequential execution
  - Status checking, cleanup operations

- **`setup_genome.py`**: Genome setup orchestrator
  - Calls `metainformant.rna.genome_prep.orchestrate_genome_setup()`
  - Complete genome setup pipeline

- **`discover_species.py`**: Discovery/config generator
  - Calls `metainformant.rna.discovery` functions
  - Discover species and generate configs

- **`check_environment.py`**: Made thinner
  - Only calls `metainformant.rna.environment.validate_environment()`

### Phase 3: Remove Obsolete Scripts

Deleted 20+ obsolete scripts:
- Workflow orchestration: `orchestrate_workflows.py`, `run_multi_species.py`, `run_all_species_parallel.py`, `workflow_ena_integrated.py`, `batch_download_species.py`, `assess_progress.py`, `initialize_progress_tracking.py`, `cleanup_progress_state.py`, `cleanup_partial_downloads.py`, `emergency_cleanup.py`, `fix_abundance_naming.py`
- Genome setup: `orchestrate_genome_setup.py`, `verify_genomes_and_indexes.py`, `download_missing_genomes.py`, `prepare_transcriptomes.py`, `build_kallisto_indexes.py`, `run_genome_setup.sh`
- Discovery: `discover_ant_species_with_rnaseq.py`, `discover_ant_rnaseq_by_genus.py`, `generate_ant_configs_with_genomes.py`
- Environment: `check_r_dependencies.py` (merged into `check_environment.py`)

### Phase 4: Update Documentation

Updated `README.md` and `AGENTS.md` to reflect new consolidated structure.

## Design Principles

1. **Thin Scripts**: Scripts only parse args and call src methods
2. **Single-Species Sequential**: Focus on one species at a time (not parallel multi-species)
3. **Hybrid Approach**: Use metainformant for complex logic, direct CLI for simple steps
4. **All Logic in src**: No business logic in scripts directory
5. **Backward Compatibility**: New orchestrators support same use cases as old scripts

## Previous AI Contributions

### Pipeline Automation
**Code Assistant Agent** implemented:
- `batch_ena.py` for resilient ENA downloads, concurrent quantification, and log management
- Kallisto execution orchestration with dynamic thread allocation
- Multi-species workflow orchestration

### Monitoring and Testing
**Code Assistant Agent** developed:
- Monitoring functions in `src/metainformant/rna/monitoring.py` for programmatic status checking
- `run_workflow.py --status` for comprehensive workflow status
- Test scripts in `scripts/rna/` for integration testing (test_*.py, run_e2e_pbarbatus.py)
- Comprehensive test suites in `tests/rna/`

### Data Management
**Code Assistant Agent** created:
- Cleanup functions for safe deletion of quantified files
- Progress tracking and assessment functions

### Download Optimization
**Code Assistant Agent** implemented:
- ENA API integration for direct FASTQ URL discovery
- wget-based downloads with --continue for resume capability
- Parallel download orchestration
- Automatic retry logic

### Operational Guidance
**Documentation Agent** authored:
- Usage instructions in `README.md`
- Environment and dependency reminders
- Comprehensive workflow documentation

## Maintenance Practices

- All business logic in `src/metainformant/rna/`
- Scripts are thin wrappers that call src methods
- Single-species sequential execution focus
- Hybrid approach: metainformant for complex logic, direct CLI for simple steps
- Update documentation when adding new functionality
- Keep tests in `tests/rna` aligned with automation behavior

## Script Organization

- **Main orchestrator**: `run_workflow.py` (single-species end-to-end workflows) ⭐ **RECOMMENDED**
- **Genome setup**: `setup_genome.py` (genome preparation pipeline)
- **Discovery**: `discover_species.py` (species discovery and config generation)
- **Environment**: `check_environment.py` (environment validation)
- **Amalgkit workflows**: Centralized in `scripts/rna/amalgkit/` (run_amalgkit.sh, verify_workflow.sh)
- **Utility scripts**: `convert_existing_sra.py` (SRA to FASTQ conversion), `fix_tmp_space.sh` (temp cleanup)
- **Test scripts**: `test_*.py`, `run_e2e_pbarbatus.py` (integration testing, not for production)

**For end-to-end workflows**: Use `run_workflow.py` which calls `execute_workflow()` to run all 11 amalgkit steps with automatic genome setup and per-sample processing.

All scripts follow repository standards with no hardcoded paths and proper integration with `output/` directory structure.
