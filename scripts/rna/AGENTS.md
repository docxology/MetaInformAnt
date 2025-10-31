# AI Agents in RNA Script Development

This document records AI involvement in the high-throughput RNA processing scripts housed here.

## AI Contributions

### Pipeline Automation
**Code Assistant Agent** implemented:
- `batch_ena.py` for resilient ENA downloads, concurrent quantification, and log management
- `run_multi_species_amalgkit.py` for multi-species orchestration with cross-species analysis
- `run_multi_species_sequential.py` for disk-space-friendly sequential processing
- Kallisto execution orchestration with dynamic thread allocation
- Multi-species workflow orchestration

### Monitoring and Testing
**Code Assistant Agent** developed:
- `monitor_workflow.py` for real-time workflow monitoring dashboard
- `monitor_amalgkit_progress.sh` for multi-species progress tracking
- `test_pbarbatus_workflow.py` for workflow testing
- `test_single_species.py` for single species validation
- `test_skip_logic.py` and `verify_skip_logic.sh` for skip logic verification

### Data Management
**Code Assistant Agent** created:
- `cleanup_quantified_sra.sh` for safe deletion of quantified SRA files
- `list_unquantified.sh` for reporting unquantified samples
- `quant_downloaded_samples.py` for quantifying already-downloaded samples

### Download Optimization
**Code Assistant Agent** implemented:
- `force_fasterq.sh` for SRA FASTQ processing with fallbacks
- `force_fasterq_parallel.sh` for parallel processing
- `process_one_srr.sh` for single SRR processing with multiple sources (ENA, AWS ODP, SRA)

### Operational Guidance
**Documentation Agent** authored:
- Usage instructions in `README.md`
- Environment and dependency reminders for wget, kallisto, and uv workflows
- Integration notes tying script output paths to `docs/rna/examples`
- Comprehensive workflow documentation

### Validation Support
**Code Assistant Agent** cross-checked:
- CLI argument handling with amalgkit metadata expectations
- File placement within `output/` per repository rules
- Error handling branches used during large cohort reprocessing
- Skip logic for already-quantified samples

## Maintenance Practices
- Capture real run logs when updating retry semantics or concurrency parameters
- Reflect new script flags in both this AGENTS file and the accompanying README
- Keep tests in `tests/rna` aligned with the automation behavior documented here
- Remove obsolete scripts with hardcoded user-specific paths
- Maintain consistency across multi-species workflows

## Script Organization
- **Amalgkit workflows**: Centralized in `scripts/rna/amalgkit/`
- **Multi-species orchestration**: `run_multi_species_*.py` scripts
- **Download utilities**: `batch_ena.py`, `force_fasterq*.sh`, `process_one_srr.sh`
- **Monitoring**: `monitor_*.py/sh` scripts
- **Testing**: `test_*.py`, `verify_*.sh` scripts
- **Data management**: `cleanup_*.sh`, `list_*.sh`, `quant_*.py` scripts

All scripts follow repository standards with no hardcoded paths and proper integration with `output/` directory structure.


