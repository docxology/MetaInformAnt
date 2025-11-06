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

### Run Workflow for a Species

```bash
# Full workflow (all steps)
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pbarbatus.yaml

# Specific steps
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pbarbatus.yaml --steps getfastq quant merge

# Check status
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pbarbatus.yaml --status

# Cleanup unquantified samples
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pbarbatus.yaml --cleanup-unquantified
```

### Setup Genome for a Species

```bash
# Full genome setup
python3 scripts/rna/setup_genome.py --config config/amalgkit/amalgkit_pbarbatus.yaml

# Verify status only
python3 scripts/rna/setup_genome.py --config config/amalgkit/amalgkit_pbarbatus.yaml --verify-only

# Skip specific steps
python3 scripts/rna/setup_genome.py --config config/amalgkit/amalgkit_pbarbatus.yaml --skip-download --skip-prepare
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

### `run_workflow.py` ⭐ **Main Orchestrator**

Thin wrapper that calls `metainformant.rna.orchestration.run_workflow_for_species()` to run complete amalgkit workflows for single species.

**Features:**
- Single-species sequential execution
- Run all amalgkit steps: metadata → integrate → config → select → getfastq → quant → merge → cstmm → curate → csca → sanity
- Status checking and progress monitoring
- Cleanup operations (partial downloads, unquantified samples)
- Progress tracking initialization

**Usage:**
```bash
# Full workflow
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pbarbatus.yaml

# Specific steps
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pbarbatus.yaml --steps getfastq quant merge

# Status check
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pbarbatus.yaml --status

# Cleanup operations
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pbarbatus.yaml --cleanup-unquantified
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pbarbatus.yaml --cleanup-partial
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pbarbatus.yaml --fix-abundance-naming
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
python3 scripts/rna/setup_genome.py --config config/amalgkit/amalgkit_pbarbatus.yaml

# Verify only
python3 scripts/rna/setup_genome.py --config config/amalgkit/amalgkit_pbarbatus.yaml --verify-only

# Skip steps
python3 scripts/rna/setup_genome.py --config config/amalgkit/amalgkit_pbarbatus.yaml --skip-download --skip-prepare

# Custom k-mer size
python3 scripts/rna/setup_genome.py --config config/amalgkit/amalgkit_pbarbatus.yaml --kmer-size 23
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
- `metainformant.rna.monitoring.check_workflow_status()` - Check workflow status

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

**Workflow:**
- Check status with `--status` flag
- Use `--cleanup-unquantified` to quantify downloaded samples
- Use `--cleanup-partial` to remove partial downloads

**Genome Setup:**
- Use `--verify-only` to check current status
- Use `--skip-*` flags to skip specific steps
- Check logs in `output/amalgkit/{species}/logs/`

## Examples

See `docs/rna/examples/` for complete workflow examples.
