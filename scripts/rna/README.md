# RNA-seq Processing Scripts

Thin orchestrator scripts for RNA-seq data processing. All business logic is in `src/metainformant/rna/`.

## Directory Structure

```
scripts/rna/
├── amalgkit/                       # Amalgkit-specific workflow scripts
│   ├── run_amalgkit.sh            # Comprehensive pipeline orchestrator
│   ├── verify_workflow.sh         # Workflow validation
│   └── README.md                  # Detailed amalgkit documentation
├── run_workflow.py                # Main workflow orchestrator ⭐
├── setup_genome.py                # Genome setup orchestrator
├── discover_species.py             # Species discovery and config generator
├── check_environment.py           # Environment validation
├── _setup_utils.py                # Shared setup utilities (venv, dependencies)
├── README.md                      # This file
└── AGENTS.md                      # AI agent documentation
```

## Quick Start

### Environment Setup

**All setup uses `uv` for package management** (required):

```bash
# Install uv if not already installed
curl -LsSf https://astral.sh/uv/install.sh | sh

# Create virtual environment
uv venv

# Install metainformant
uv pip install -e . --python .venv/bin/python3

# Install amalgkit
uv pip install git+https://github.com/kfuku52/amalgkit --python .venv/bin/python3
```

**Note**: Scripts automatically detect and use the virtual environment (`.venv` or `/tmp/metainformant_venv` on filesystems without symlink support). No manual activation needed when running scripts.

**External Drive Support**: On filesystems with limitations (e.g., ext6 without symlink support), the system automatically:
- Uses `/tmp/metainformant_venv` for virtual environment
- Uses `/tmp/uv-cache` for UV package cache (avoids symlink errors)
- Handles all filesystem limitations transparently

See [docs/rna/EXTERNAL_DRIVE_SETUP.md](../../docs/rna/EXTERNAL_DRIVE_SETUP.md) for complete documentation.

### Run End-to-End Workflow for a Species

**Recommended**: Use `run_workflow.py` for complete end-to-end execution:

```bash
# Full end-to-end workflow (all steps: metadata → integrate → config → select → getfastq → quant → merge → cstmm → curate → csca → sanity)
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml

# Specific steps only
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps getfastq quant merge

# Check status
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --status

# Cleanup unquantified samples
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --cleanup-unquantified
```

**What happens automatically**:
- Genome download and kallisto index preparation (if genome config exists)
- Metadata retrieval from NCBI SRA
- Per-sample processing: download → quantify → delete FASTQ (as configured with `keep_fastq: no`)
- All downstream steps: integrate, merge, cstmm, curate, csca, sanity

### Setup Genome for a Species

```bash
# Full genome setup
python3 scripts/rna/setup_genome.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml

# Verify status only
python3 scripts/rna/setup_genome.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --verify-only

# Skip specific steps
python3 scripts/rna/setup_genome.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --skip-download --skip-prepare
```

### Discover Species and Generate Configs

```bash
# Discover ant species with RNA-seq data
python3 scripts/rna/discover_species.py --output config/amalgkit/

# Generate config for specific species
python3 scripts/rna/discover_species.py --species "Camponotus floridanus" --output config/amalgkit/
```

### Check Environment

```bash
python3 scripts/rna/check_environment.py
```

## Available Scripts

### `run_workflow.py` ⭐ **Main Orchestrator (Recommended for End-to-End)**

Thin wrapper that calls `metainformant.rna.orchestration.run_workflow_for_species()` which uses `execute_workflow()` to run complete amalgkit workflows for single species.

**Features:**
- **Complete end-to-end execution**: All 11 amalgkit steps in correct order
- **Automatic genome setup**: Downloads genome, prepares transcriptome, builds kallisto index (if genome config exists)
- **Per-sample processing**: Download → quantify → delete FASTQ (as configured)
- **Status checking and progress monitoring**: Check workflow status at any time
- **Cleanup operations**: Partial downloads, unquantified samples, abundance file naming
- **Resume support**: Automatically skips completed steps

**Usage:**
```bash
# Full workflow
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml

# Specific steps
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps getfastq quant merge

# Status check
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --status

# Cleanup operations
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --cleanup-unquantified
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --cleanup-partial
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --fix-abundance-naming
```

### `setup_genome.py` **Genome Setup Orchestrator**

Thin wrapper that calls `metainformant.rna.genome_prep.orchestrate_genome_setup()` to orchestrate complete genome setup pipelines.

**Features:**
- Verify genome/index status
- Download missing genomes
- Prepare transcriptomes
- Build kallisto indexes
- Complete genome setup pipeline

**Usage:**
```bash
# Full setup
python3 scripts/rna/setup_genome.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml

# Verify only
python3 scripts/rna/setup_genome.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --verify-only

# Skip steps
python3 scripts/rna/setup_genome.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --skip-download --skip-prepare

# Custom k-mer size
python3 scripts/rna/setup_genome.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --kmer-size 23
```

### `discover_species.py` **Discovery/Config Generator**

Thin wrapper that calls `metainformant.rna.discovery` functions to discover species with RNA-seq data and generate amalgkit YAML configurations.

**Features:**
- Discover species with RNA-seq data
- Generate amalgkit YAML configurations
- Get genome assembly information

**Usage:**
```bash
# Discover ant species
python3 scripts/rna/discover_species.py --output config/amalgkit/

# Specific species
python3 scripts/rna/discover_species.py --species "Camponotus floridanus" --output config/amalgkit/

# Custom search query
python3 scripts/rna/discover_species.py --search-query 'txid1234[Organism] AND RNA-Seq[Strategy]' --output config/amalgkit/

# Minimum samples filter
python3 scripts/rna/discover_species.py --output config/amalgkit/ --min-samples 10
```

### `check_environment.py` **Environment Checker**

Thin wrapper that calls `metainformant.rna.environment.validate_environment()` to verify all required tools and dependencies.

**Usage:**
```bash
python3 scripts/rna/check_environment.py
```

### Status Checking

Use `run_workflow.py --status` for comprehensive workflow status:

```bash
# Quick status check
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --status

# Detailed status with recommendations
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --status --detailed
```

**Status output includes:**
- Step completion status
- Sample counts (downloaded, quantified, remaining)
- Disk usage information
- Error summaries
- Recommended next steps

## Script Consolidation

All scripts have been consolidated into 4-5 thin orchestrators that call methods from `src/metainformant/rna/`. All business logic has been moved to the source code, making scripts maintainable and consistent.

**Key Functions in `metainformant.rna`:**
- `metainformant.rna.orchestration.run_workflow_for_species()` - Run workflow for a species
- `metainformant.rna.orchestration.cleanup_unquantified_samples()` - Cleanup unquantified samples
- `metainformant.rna.genome_prep.orchestrate_genome_setup()` - Orchestrate genome setup
- `metainformant.rna.genome_prep.verify_genome_status()` - Verify genome status
- `metainformant.rna.discovery.search_species_with_rnaseq()` - Search for species
- `metainformant.rna.discovery.generate_config_yaml()` - Generate config files
- `metainformant.rna.cleanup.cleanup_partial_downloads()` - Cleanup partial downloads
- `metainformant.rna.monitoring.check_workflow_progress()` - Check workflow progress
- `metainformant.rna.orchestration.check_workflow_status()` - Check workflow status (unified interface)

**Consolidated Scripts:**

The following scripts have been consolidated:

- **Workflow orchestration** (12 scripts removed): `orchestrate_workflows.py`, `run_multi_species.py`, `run_all_species_parallel.py`, `workflow_ena_integrated.py`, `batch_download_species.py`, `assess_progress.py`, `initialize_progress_tracking.py`, `cleanup_progress_state.py`, `cleanup_partial_downloads.py`, `emergency_cleanup.py`, `fix_abundance_naming.py`
  → Use `run_workflow.py`

- **Genome setup** (6 scripts removed): `orchestrate_genome_setup.py`, `verify_genomes_and_indexes.py`, `download_missing_genomes.py`, `prepare_transcriptomes.py`, `build_kallisto_indexes.py`, `run_genome_setup.sh`
  → Use `setup_genome.py`

- **Discovery/config generation** (3 scripts removed): `discover_ant_species_with_rnaseq.py`, `discover_ant_rnaseq_by_genus.py`, `generate_ant_configs_with_genomes.py`
  → Use `discover_species.py`

- **Environment checks** (1 script removed): `check_r_dependencies.py`
  → Merged into `check_environment.py`

## Output Structure

All outputs go to `output/amalgkit/{species}/`:
- `quant/` - Quantification results (abundance.tsv, ~2 MB per sample)
- `work/metadata/` - Filtered metadata
- `fastq/` - FASTQ files (automatically cleaned after quantification)
- `work/index/` - Kallisto index files
- `logs/` - Processing logs

## Troubleshooting

**Environment:**
- Run `check_environment.py` to verify all dependencies
- Virtual environment is auto-activated by scripts
- **All setup uses `uv`**: `uv venv` and `uv pip install` (see Environment Setup above)
- If venv creation fails (symlink issues), scripts automatically use `/tmp/metainformant_venv`
- UV cache automatically uses `/tmp/uv-cache` to avoid symlink errors on ext6 filesystems

**Workflow:**
- Check status with `--status` flag
- Use `--cleanup-unquantified` to quantify downloaded samples
- Use `--cleanup-partial` to remove partial downloads

**Genome Setup:**
- Use `--verify-only` to check current status
- Use `--skip-*` flags to skip specific steps
- Check logs in `output/amalgkit/{species}/logs/`

## Test Scripts

These scripts are for integration testing and debugging specific functionality. They use real implementations (no mocks) and are meant for manual testing.

### `test_heartbeat.py`
Tests heartbeat logging on a single sample during getfastq step.

**Usage:**
```bash
python3 scripts/rna/test_heartbeat.py
```

### `test_enhanced_heartbeat.py`
Tests enhanced heartbeat logging with detailed progress tracking.

**Usage:**
```bash
python3 scripts/rna/test_enhanced_heartbeat.py
```

### `test_getfastq_fix.py`
Tests fixes for getfastq step issues (disk space, LITE SRA files, etc.).

**Usage:**
```bash
python3 scripts/rna/test_getfastq_fix.py
```

### `test_genome_prep.py`
Tests genome preparation functions (download, transcriptome extraction, index building).

**Usage:**
```bash
python3 scripts/rna/test_genome_prep.py
```

### `run_e2e_pbarbatus.py`
**TESTING ONLY**: Runs end-to-end workflow for Pogonomyrmex barbatus with a single sample to verify the entire system.

**Note**: This script limits to 1 sample for testing. For production workflows, use `run_workflow.py`.

**Usage:**
```bash
python3 scripts/rna/run_e2e_pbarbatus.py
```

**For automated unit tests**, see `tests/rna/` directory.

## Utility Scripts

### `convert_existing_sra.py`
Converts existing SRA files to FASTQ format. Finds all SRA files without corresponding FASTQ files and converts them using parallel processing.

**Usage:**
```bash
# Convert SRA files for a specific species
python3 scripts/rna/convert_existing_sra.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml

# Use specific number of threads
python3 scripts/rna/convert_existing_sra.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --threads 8
```

**Features:**
- Automatically finds SRA files without FASTQ counterparts
- Parallel conversion for faster processing
- Progress tracking and error handling
- Respects workflow configuration

### `fix_tmp_space.sh`
Quick utility to clean up `/tmp` space by removing pytest temp files and uv cache.

**Usage:**
```bash
bash scripts/rna/fix_tmp_space.sh
```

**What it does:**
- Checks `/tmp` disk usage
- Removes pytest temp directories (`/tmp/pytest-of-*`)
- Removes uv cache (`/tmp/uv-cache`)
- Shows large temp files for manual review

**Note**: This is a general-purpose utility. Use when `/tmp` space is low and affecting workflow execution.

## Examples

See `docs/rna/EXAMPLES.md` for complete workflow examples.
