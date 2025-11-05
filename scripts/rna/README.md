# RNA-seq Processing Scripts

Scripts for RNA-seq data processing, including SRA download, quantification, and multi-species workflows.

## Directory Structure

```
scripts/rna/
├── amalgkit/                       # Amalgkit-specific workflow scripts
│   ├── run_amalgkit.sh            # Comprehensive pipeline orchestrator
│   ├── verify_workflow.sh         # Workflow validation
│   └── README.md                  # Detailed amalgkit documentation
├── orchestrate_workflows.py       # **UNIFIED ORCHESTRATOR** ⭐⭐⭐ (START HERE)
│                                  # Status, monitoring, cleanup, workflow execution
├── workflow_ena_integrated.py     # Integrated ENA download + quantification (PRODUCTION) ⭐
├── run_multi_species.py           # Multi-species with cross-species analysis (legacy SRA-based)
├── run_all_species_parallel.py    # Parallel execution of all species workflows
├── batch_download_species.py      # Configurable batch download for multiple species
├── check_environment.py           # Environment validation
├── cleanup_partial_downloads.py   # Clean up partial/failed downloads
├── fix_abundance_naming.py        # Fix abundance file naming for merge compatibility
├── _setup_utils.py                # Shared setup utilities (venv, dependencies)
├── README.md                      # This file
└── AGENTS.md                      # AI agent documentation
```

## Quick Start: Unified Orchestrator ⭐

**The new `orchestrate_workflows.py` consolidates all workflow management into one intelligent script.**

### Common Tasks

```bash
# 1. Status and monitoring (replaces unified_status.py, monitor_comprehensive.py, etc.)
python3 scripts/rna/orchestrate_workflows.py --status                    # Brief status
python3 scripts/rna/orchestrate_workflows.py --status --detailed         # Detailed status with categories
python3 scripts/rna/orchestrate_workflows.py --monitor                   # Real-time monitoring
python3 scripts/rna/orchestrate_workflows.py --monitor --watch 60        # Watch mode with interval

# 2. Full assessment
python3 scripts/rna/orchestrate_workflows.py --assess

# 3. Cleanup: quantify downloaded samples and delete FASTQs
python3 scripts/rna/orchestrate_workflows.py --cleanup-unquantified

# 4. Run merge+curate+sanity for species ready
python3 scripts/rna/orchestrate_workflows.py --steps merge curate sanity --auto-species

# 5. Resume incomplete downloads
python3 scripts/rna/orchestrate_workflows.py --resume-downloads

# 6. Run specific steps for specific species
python3 scripts/rna/orchestrate_workflows.py --species cfloridanus --steps merge curate sanity

# 7. Full pipeline for one species
python3 scripts/rna/orchestrate_workflows.py --species pbarbatus --steps metadata select getfastq quant merge curate sanity
```

### Key Features

- **Status & Monitoring**: Comprehensive status checks and real-time monitoring (replaces 11+ status/monitoring scripts)
- **Smart Assessment**: Comprehensive status check showing progress, readiness, and recommendations
- **Flexible Execution**: Run any combination of steps for any combination of species
- **Auto-Selection**: Automatically choose species ready for specific steps (--auto-species)
- **Cleanup**: Quantify downloaded samples and cleanup FASTQs in one command
- **Resume Support**: Restart incomplete downloads intelligently
- **Auto-Discovery**: Automatically discovers all species from config files

## Script Consolidation

All scripts in `scripts/rna/` have been consolidated to use thin orchestration patterns that call tested, documented methods in `metainformant.rna`. Core functionality (quantification, deletion, status checking, etc.) is now in the metainformant module, making scripts maintainable and consistent.

**Key Functions in `metainformant.rna`:**
- `metainformant.rna.steps.quantify_sample()` - Quantify a single sample
- `metainformant.rna.steps.delete_sample_fastqs()` - Delete FASTQ files for a sample
- `metainformant.rna.steps.process_sample_pipeline()` - Complete pipeline: SRA→FASTQ→Quant→Delete
- `metainformant.rna.steps.convert_sra_to_fastq()` - Convert SRA files to FASTQ
- `metainformant.rna.count_quantified_samples()` - Count quantified samples
- `metainformant.rna.find_unquantified_samples()` - Find unquantified samples
- `metainformant.rna.check_workflow_progress()` - Get workflow progress
- `metainformant.rna.analyze_species_status()` - Comprehensive status analysis

**Consolidated Functionality:**

All status, monitoring, and workflow management functionality has been consolidated into `orchestrate_workflows.py`:

- **Status scripts** (11 removed): `check_status.py`, `comprehensive_status.py`, `detailed_progress.py`, `full_assessment.py`, `get_current_status.py`, `quick_status.py`, `unified_status.py`, `monitor_workflow.py`, `monitor_comprehensive.py`, `check_progress_now.py`, `run_assessment.py`
  → Use `orchestrate_workflows.py --status` or `--status --detailed`

- **Monitoring scripts** (2 removed): `monitor_comprehensive.py`, `monitor_workflow.py`
  → Use `orchestrate_workflows.py --monitor`

- **Quant/cleanup scripts** (3 removed): `quant_and_cleanup.py`, `quant_downloaded_samples.py`, `manual_quant_cleanup.py`
  → Use `orchestrate_workflows.py --cleanup-unquantified` or `--steps quant`

- **Workflow management scripts** (4 removed): `restart_all_workflows.py`, `check_and_restart_workflows.py`, `process_sra_samples.py`, `analyze_sample_status.py`
  → Use `orchestrate_workflows.py --resume-downloads` or `--steps <step>`

- **Download scripts** (1 removed): `download_ena_robust.py`
  → Functionality integrated into `workflow_ena_integrated.py`

- **Shell scripts** (6 removed): `run_all_ant_species.sh`, `run_top10_ant_species.sh`, `run_batch2_ant_species.sh`, `cleanup_quantified_sra.sh`, `list_unquantified.sh`, `monitor_amalgkit_progress.sh`
  → Use `orchestrate_workflows.py` or `run_all_species_parallel.py`

## Available Scripts

### `workflow_ena_integrated.py` ⭐ **PRODUCTION**
**Robust integrated download + quantification workflow**

Production-ready workflow bypassing SRA Toolkit issues:

**Features:**
- **Direct ENA downloads**: Fetches FASTQs directly from European Nucleotide Archive API
- **Robust retry logic**: Uses wget with --continue for resume capability and automatic retries
- **Batched processing**: Download N samples → Quantify → Delete FASTQs → Repeat
- **Disk-friendly**: Only one batch of FASTQs on disk at a time
- **Auto-detection**: Handles both single-end and paired-end data
- **Resume support**: Skips already-quantified samples automatically

**Why ENA over SRA Toolkit:**
- SRA Toolkit downloads fail frequently (~100% failure rate on large samples)
- ENA provides direct FASTQ files (no SRA→FASTQ conversion needed)
- wget --continue allows proper resumption after network interruptions
- Much more reliable for large-scale downloads

**Usage:**
```bash
# Full workflow with default settings (12 samples/batch, 12 threads)
python3 scripts/rna/workflow_ena_integrated.py \
  --config config/amalgkit/amalgkit_cfloridanus.yaml \
  --batch-size 12 \
  --threads 12

# Test with 3 samples
python3 scripts/rna/workflow_ena_integrated.py \
  --config config/amalgkit/amalgkit_cfloridanus.yaml \
  --batch-size 3 \
  --threads 8 \
  --max-samples 3

# Resume (skip download, only quantify existing FASTQs)
python3 scripts/rna/workflow_ena_integrated.py \
  --config config/amalgkit/amalgkit_cfloridanus.yaml \
  --skip-download
```

**Performance:**
- Download: ~6 minutes per 3 samples (varies by sample size)
- Quantification: ~36 seconds per sample (single-end, 8 threads)
- Cleanup: Instant FASTQ deletion after quantification
- Peak disk: ~1.5 GB per sample (batch size × 1.5 GB)

**Requirements:**
- wget (for robust downloads)
- kallisto (for quantification)
- Kallisto index must exist (built once per species)

### `download_ena_robust.py` ⭐
**Standalone robust ENA downloader**

Used internally by `workflow_ena_integrated.py`, but can be used standalone:

**Usage:**
```bash
python3 scripts/rna/download_ena_robust.py \
  --metadata output/amalgkit/cfloridanus/work/metadata/metadata.tsv \
  --out-dir output/amalgkit/cfloridanus/fastq \
  --threads 12 \
  --max-retries 3
```

### `run_multi_species.py`
**Legacy multi-species workflow (SRA Toolkit-based)**

Alternative workflow using SRA Toolkit instead of ENA direct downloads:
- Auto-activation and environment management
- Batched processing (10 samples at a time)
- Cross-species analysis (CSTMM, CSCA)
- Complete pipeline automation

**Note:** Consider using `workflow_ena_integrated.py` for better reliability.

**Usage:**
```bash
# Process all discovered species
python3 scripts/rna/run_multi_species.py
```

### Status and Monitoring

All status and monitoring functionality is now in `orchestrate_workflows.py`:

**Status:**
```bash
# Brief status (replaces check_status.py, quick_status.py, etc.)
python3 scripts/rna/orchestrate_workflows.py --status

# Detailed status with sample categories (replaces unified_status.py --detailed)
python3 scripts/rna/orchestrate_workflows.py --status --detailed
```

**Monitoring:**
```bash
# Real-time monitoring (replaces monitor_comprehensive.py, monitor_workflow.py)
python3 scripts/rna/orchestrate_workflows.py --monitor

# Watch mode with custom interval
python3 scripts/rna/orchestrate_workflows.py --monitor --watch 60
```

### Utility Scripts

#### `cleanup_partial_downloads.py`
**Clean up partial and failed downloads**

Removes samples with partial FASTQ/SRA files that aren't quantified:
- Safe deletion with dry-run option
- Frees disk space for retrying downloads
- Processes all species automatically

**Usage:**
```bash
# Dry run (see what would be deleted)
python3 scripts/rna/cleanup_partial_downloads.py --dry-run

# Actually delete partial downloads
python3 scripts/rna/cleanup_partial_downloads.py --execute
```

**Note:** For quantifying downloaded samples and cleaning up FASTQs, use `orchestrate_workflows.py --cleanup-unquantified`.

#### `batch_download_species.py` ⭐ **NEW**
**Configurable batch download for multiple species in parallel**

Downloads samples for multiple species simultaneously with configurable parallelism:
- Configurable species count and threads per species
- Automatic virtual environment activation
- Cloud acceleration (AWS, GCP, NCBI)
- Processes all species in batches

**Usage:**
```bash
# Default: 3 species × 10 threads = 30 total downloads
python3 scripts/rna/batch_download_species.py

# Custom: 4 species × 12 threads = 48 total downloads
python3 scripts/rna/batch_download_species.py --species-count 4 --threads-per-species 12

# Limited: 2 species × 8 threads = 16 total downloads
python3 scripts/rna/batch_download_species.py --species-count 2 --threads-per-species 8
```

See `docs/rna/BATCH_DOWNLOAD_CONFIGURATION.md` for complete configuration guide.

#### `fix_abundance_naming.py`
**Fix abundance file naming for amalgkit merge compatibility**

Creates symlinks from `abundance.tsv` to `{SRR}_abundance.tsv` for amalgkit merge compatibility:

**Usage:**
```bash
python3 scripts/rna/fix_abundance_naming.py
```

**Note:** For quantifying downloaded samples and cleaning up FASTQs, use `orchestrate_workflows.py --cleanup-unquantified`.

### Testing Scripts

All test scripts have been moved to `tests/test_rna_ena_workflow.py`.
See that file for comprehensive integration tests of the ENA workflow.

## Disk Space Management

**Automatic Disk Space Management:**

All workflows and scripts now automatically delete FASTQ/SRA files after quantification to manage disk space. This is handled by `metainformant.rna.steps.delete_sample_fastqs()` which is called automatically by:

- `batch_download_species.py` - Deletes immediately after per-sample quant
- `workflow_ena_integrated.py` - Deletes after batch quant
- `orchestrate_workflows.py --cleanup-unquantified` - Deletes after quant during cleanup

**Batched Processing:**
- Downloads N samples from ENA (parallel, robust)
- Quantifies all downloaded samples with kallisto
- Deletes FASTQ files immediately after quantification (via metainformant functions)
- Repeats with next batch
- Peak usage: ~1.5 GB per sample × batch size

**FASTQ Cleanup:**
- Automatic cleanup via `metainformant.rna.steps.delete_sample_fastqs()`
- Manual cleanup available via `orchestrate_workflows.py --cleanup-unquantified`
- Partial/failed downloads cleanup via `cleanup_partial_downloads.py`
- Quantification files retained permanently (~2 MB per sample)
- Low disk space warnings indicate the delete pipeline is working correctly

## Output Structure

All outputs go to `output/amalgkit/{species}/`:
- `quant/` - Quantification results (abundance.tsv, ~2 MB per sample)
- `work/metadata/` - Filtered metadata
- `fastq/` - FASTQ files (flat structure, sample directories directly in fastq/, automatically cleaned after quantification)
- `work/index/` - Kallisto index files
- `logs/` - Processing logs from workflow_ena_integrated.py

**Structure Standardization (Nov 2025):**
- ✅ All species use flat `fastq/` structure (no nested `getfastq/` subdirectories)
- ✅ Sample directories: `fastq/{SAMPLE_ID}/` containing `.fastq.gz` or `.sra` files
- ✅ All utility scripts updated to use standardized paths
- ✅ Ensures consistent behavior and prevents false positives/negatives in workflow logic

See `docs/rna/examples/` for complete documentation.

## Troubleshooting

**ENA Downloads:**
- The workflow uses wget with automatic retry and resume
- If downloads fail, re-run the workflow - it will resume automatically
- Check network connectivity if many downloads fail

**Virtual Environment:**
- `run_multi_species.py` auto-activates venv if available
- If venv missing, script provides setup instructions
- Manual activation: `source .venv/bin/activate`

## Examples

See `docs/rna/examples/pbarbatus_analysis.md` for a complete workflow example.
