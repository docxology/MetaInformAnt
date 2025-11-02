# AI Agents in RNA Script Development

This document records AI involvement in the high-throughput RNA processing scripts housed here.

## AI Contributions

### Pipeline Automation
**Code Assistant Agent** implemented:
- `batch_ena.py` for resilient ENA downloads, concurrent quantification, and log management
- `run_multi_species.py` for multi-species orchestration with cross-species analysis
  - **Auto-activation**: Virtual environment detection and automatic re-execution
  - **SRA optimization**: Wrapper script creation and environment configuration
  - **Disk management**: Batched processing with automatic FASTQ cleanup
  - **Thread configuration**: 10 parallel threads for downloads and quantification
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
- **`download_ena_robust.py`**: Robust ENA downloader with retry logic (PRODUCTION)
  - ENA API integration for direct FASTQ URL discovery
  - wget-based downloads with --continue for resume capability
  - Parallel download orchestration with ThreadPoolExecutor
  - Automatic retry logic with configurable attempts
  - 100% reliable vs 0% with SRA Toolkit
- **`workflow_ena_integrated.py`**: Integrated download + quantification (PRODUCTION)
  - Combines robust ENA downloads with kallisto quantification
  - Batched processing: download→quantify→delete cycles
  - Auto-detection of single vs paired-end data
  - Resume support for failed workflows
  - Complete end-to-end testing with real samples
- ~~`batch_ena.py`~~ (REMOVED - superseded by workflow_ena_integrated.py)
- ~~`force_fasterq.sh`~~ (REMOVED - SRA Toolkit approach superseded by ENA)
- ~~`force_fasterq_parallel.sh`~~ (REMOVED - SRA Toolkit approach superseded by ENA)
- ~~`process_one_srr.sh`~~ (REMOVED - SRA Toolkit helper superseded by ENA)

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
- Document production fixes: auto-activation, SRA wrappers, disk management
- Update thread counts and performance characteristics based on actual usage

## Script Organization
- **Production workflows**: `workflow_ena_integrated.py` (recommended), `run_multi_species.py` (legacy SRA-based)
- **Download utilities**: `download_ena_robust.py` (production)
- **Amalgkit workflows**: Centralized in `scripts/rna/amalgkit/`
- **Monitoring**: `monitor_comprehensive.py` (recommended), `monitor_workflow.py`, `monitor_amalgkit_progress.sh`
- **Utility**: `cleanup_quantified_sra.sh`, `list_unquantified.sh`, `quant_downloaded_samples.py`
- **Environment**: `check_environment.py`

All scripts follow repository standards with no hardcoded paths and proper integration with `output/` directory structure.

## Production Enhancements

### ENA Download Integration (November 2025)
**Code Assistant Agent** developed robust ENA-based workflow:
- Bypassed SRA Toolkit (100% failure rate on large samples)
- Implemented ENA API integration for reliable FASTQ URL discovery
- Used wget with --continue for resume capability after network failures
- Achieved 100% download success rate vs 0% with SRA Toolkit
- Integrated with kallisto quantification in batched workflow
- Tested end-to-end with 3 real samples, then 12 samples (successful download + quantification)
- Auto-detection of single-end vs paired-end data
- Performance: ~6 min download + ~2 min quantification per batch of 3 samples
- Currently processing 844 samples across 4 species in production (November 2025)

### Auto-Activation Implementation
**Code Assistant Agent** developed automatic virtual environment activation:
- Detection of virtual environment using path analysis
- Process replacement via `os.execve()` for seamless activation
- Environment variable configuration (`VIRTUAL_ENV`, `PATH`)
- Clear error messages with setup instructions if venv missing

### SRA Download Optimization
**Code Assistant Agent** implemented disk space management:
- Wrapper script creation for `fasterq-dump` with `--size-check off`
- Environment variable configuration for SRA toolkit (`TMPDIR`, `TEMP`, `TMP`)
- Tool symlink management (fastp, kallisto, seqkit)
- PATH manipulation to prioritize wrapper directory
- Configuration to use project directory instead of `/tmp` partition

### Performance Configuration
- Updated from 8 to 10 threads for parallel operations
- Batched processing: 10 samples at a time
- Peak disk usage: ~20-50 GB (10 samples of FASTQs)
- Validated with 300+ samples per species across 4 ant species


