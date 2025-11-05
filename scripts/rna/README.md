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
├── workflow_ena_integrated.py     # Integrated ENA download + quantification (PRODUCTION) ⭐
├── download_ena_robust.py         # Robust ENA downloader with retry logic ⭐
├── run_multi_species.py           # Multi-species with cross-species analysis (legacy SRA-based)
├── restart_all_workflows.py       # Restart all workflows in parallel
├── check_environment.py           # Environment validation
├── cleanup_quantified_sra.sh      # Safe deletion of FASTQ files after quantification
├── monitor_comprehensive.py       # Comprehensive real-time monitoring
├── check_and_restart_workflows.py # Smart status check and conditional restart
├── [Legacy monitoring scripts]    # Various status/progress scripts
├── README.md                      # This file
└── AGENTS.md                      # AI agent documentation
```

## Quick Start: Unified Orchestrator ⭐

**The new `orchestrate_workflows.py` consolidates all workflow management into one intelligent script.**

### Common Tasks

```bash
# 1. Check status of all workflows
python3 scripts/rna/orchestrate_workflows.py --assess

# 2. Cleanup: quantify downloaded samples and delete FASTQs
python3 scripts/rna/orchestrate_workflows.py --cleanup-unquantified

# 3. Run merge+curate+sanity for species ready
python3 scripts/rna/orchestrate_workflows.py --steps merge curate sanity --auto-species

# 4. Resume incomplete downloads
python3 scripts/rna/orchestrate_workflows.py --resume-downloads

# 5. Run specific steps for specific species
python3 scripts/rna/orchestrate_workflows.py --species cfloridanus --steps merge curate sanity

# 6. Full pipeline for one species
python3 scripts/rna/orchestrate_workflows.py --species pbarbatus --steps metadata select getfastq quant merge curate sanity
```

### Key Features

- **Smart Assessment**: Comprehensive status check showing progress, readiness, and recommendations
- **Flexible Execution**: Run any combination of steps for any combination of species
- **Auto-Selection**: Automatically choose species ready for specific steps (--auto-species)
- **Cleanup**: Quantify downloaded samples and cleanup FASTQs in one command
- **Resume Support**: Restart incomplete downloads intelligently

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

**Obsolete Scripts (Replaced by `unified_status.py`):**
- `check_status.py` → Use `unified_status.py`
- `comprehensive_status.py` → Use `unified_status.py --detailed`
- `detailed_progress.py` → Use `unified_status.py`
- `full_assessment.py` → Use `unified_status.py --detailed`
- `get_current_status.py` → Use `unified_status.py`
- `quick_status.py` → Use `unified_status.py`

These scripts are kept for reference but should not be used. All functionality is available in `unified_status.py` which uses metainformant functions.

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

### Monitoring Scripts

#### `monitor_comprehensive.py` ⭐
**Comprehensive real-time workflow monitor**

Tracks all 4 species simultaneously with detailed progress:
- Sample counts (quantified/total)
- Current batch numbers
- Downloading sample counts
- FASTQ directory sizes
- Time estimates

**Usage:**
```bash
python3 scripts/rna/monitor_comprehensive.py
```

#### `monitor_workflow.py`
**Alternative monitoring dashboard**

Real-time monitoring dashboard with:
- Species progress tracking
- Disk usage monitoring
- Running process detection
- Batch activity logs

**Usage:**
```bash
python3 scripts/rna/monitor_workflow.py
```

#### `monitor_amalgkit_progress.sh`
**Simple progress monitor**

Lightweight bash-based monitoring (no Python dependencies):

**Usage:**
```bash
bash scripts/rna/monitor_amalgkit_progress.sh

# Or watch mode
watch -n 60 bash scripts/rna/monitor_amalgkit_progress.sh
```

#### `unified_status.py` ⭐ **NEW**
**Unified status checking script (replaces multiple status scripts)**

Replaces: `check_status.py`, `comprehensive_status.py`, `detailed_progress.py`, `full_assessment.py`, `get_current_status.py`, `quick_status.py`

**Usage:**
```bash
# Brief status
python3 scripts/rna/unified_status.py

# Detailed status with sample categories
python3 scripts/rna/unified_status.py --detailed

# Filter by species
python3 scripts/rna/unified_status.py --species cfloridanus --detailed
```

**Output:**
- Per-species progress and quantified counts
- Sample categories (quantified and deleted, downloading, failed, etc.)
- Overall summary statistics

#### `restart_all_workflows.py` ⭐
**Restart all workflows in parallel**

Convenient script to restart all 4 species workflows simultaneously:

**Usage:**
```bash
python3 scripts/rna/restart_all_workflows.py
```

This will start all workflows with:
- 12 samples per batch
- 12 parallel threads
- Background execution (nohup)
- Timestamped log files

#### `check_and_restart_workflows.py`
**Smart status check and conditional restart**

Checks workflow status and only restarts inactive workflows:

**Usage:**
```bash
python3 scripts/rna/check_and_restart_workflows.py
```

### Utility Scripts

#### `list_unquantified.sh`
**Generate reports of samples needing quantification**

Identifies samples with downloaded FASTQ or SRA files but no quantification output:
- Works with standardized flat `fastq/` structure
- Detects both `.fastq.gz` (ENA downloads) and `.sra` (SRA Toolkit) files
- Scans all species directories
- Reports size and sample count
- Creates sample lists in output/amalgkit/

**Usage:**
```bash
bash scripts/rna/list_unquantified.sh
```

**Output:**
- `output/amalgkit/{species}_unquantified.txt` - Sample lists
- Console report with sizes and counts

#### `cleanup_quantified_sra.sh`
**Safe deletion of FASTQ files after quantification**

Reclaims disk space by removing FASTQ files after successful quantification.
Works with both SRA Toolkit downloads (.sra files) and ENA downloads (.fastq.gz files):
- Verifies quantification completion before deletion
- Detailed logging of operations
- Safe: skips unquantified samples
- Uses standardized flat `fastq/` structure

**Usage:**
```bash
# Preview what will be deleted (dry run)
bash scripts/rna/cleanup_quantified_sra.sh

# Execute cleanup
bash scripts/rna/cleanup_quantified_sra.sh --execute
```

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

#### `quant_downloaded_samples.py`
**Quantify already-downloaded samples**

Automatically discovers all species and quantifies downloaded but unquantified samples:
- Auto-activates virtual environment
- Discovers all species configs automatically
- Handles both FASTQ and SRA files
- Deletes FASTQs after successful quantification

**Usage:**
```bash
# Quantify all downloaded samples across all species
python3 scripts/rna/quant_downloaded_samples.py
```

#### `analyze_sample_status.py` ⭐ **NEW**
**Comprehensive sample status analysis**

Categorizes samples into:
- Quantified and deleted
- Quantified but not deleted
- Currently downloading
- Failed download
- Undownloaded

**Usage:**
```bash
# Analyze all samples across all species
python3 scripts/rna/analyze_sample_status.py
```

#### `cleanup_partial_downloads.py` ⭐ **NEW**
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

#### `manual_quant_cleanup.py`
**Manual quantification and cleanup utility**

Sequential processing for manual control:

**Usage:**
```bash
python3 scripts/rna/manual_quant_cleanup.py \
  --config config/amalgkit/amalgkit_cfloridanus.yaml
```

#### `quant_and_cleanup.py`
**Batch quantification and cleanup**

Parallel quantification followed by cleanup:

**Usage:**
```bash
python3 scripts/rna/quant_and_cleanup.py \
  --config config/amalgkit/amalgkit_cfloridanus.yaml \
  --threads 12
```

### Testing Scripts

All test scripts have been moved to `tests/test_rna_ena_workflow.py`.
See that file for comprehensive integration tests of the ENA workflow.

## Disk Space Management

**Automatic Disk Space Management:**

All workflows and scripts now automatically delete FASTQ/SRA files after quantification to manage disk space. This is handled by `metainformant.rna.steps.delete_sample_fastqs()` which is called automatically by:

- `batch_download_species.py` - Deletes immediately after per-sample quant
- `process_sample_pipeline()` - Automatically deletes after quant
- `quant_downloaded_samples.py` - Deletes after quant
- `quant_and_cleanup.py` - Deletes after quant
- `orchestrate_workflows.py` - Deletes after quant during cleanup
- `manual_quant_cleanup.py` - Deletes after quant

**Batched Processing:**
- Downloads N samples from ENA (parallel, robust)
- Quantifies all downloaded samples with kallisto
- Deletes FASTQ files immediately after quantification (via metainformant functions)
- Repeats with next batch
- Peak usage: ~1.5 GB per sample × batch size

**FASTQ Cleanup:**
- Automatic cleanup via `metainformant.rna.steps.delete_sample_fastqs()`
- Manual cleanup available via `cleanup_quantified_sra.sh`
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
