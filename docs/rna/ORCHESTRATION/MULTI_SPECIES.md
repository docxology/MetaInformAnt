# Multi-Species SRA Workflow Orchestrator

**Script**: `scripts/rna/run_multi_species.py`

## Overview

Multi-species SRA-based workflow orchestrator with auto-discovery and cross-species analysis. Best for testing, legacy environments, and multi-species coordination.

## Configuration Block

```python
# ============================================================================
# CONFIGURATION
# ============================================================================
# Scope: Multi-species SRA-based workflow with auto-discovery
# Steps: metadata → select → getfastq → quant → merge → curate → cstmm → csca → sanity
# Config: Auto-discovers all config/amalgkit/amalgkit_*.yaml files
# Threads: Per-config (default 10) or config override
# Batch Size: 10 samples per species (configurable via threads in config)
# Output: output/amalgkit/{species}/work/ per species
# Dependencies: SRA Toolkit, kallisto, fastp, seqkit, amalgkit
# Virtual Env: Auto-activates if .venv exists
# Reliability: ~0% for large samples (SRA Toolkit limitation)
# ============================================================================
```

## When to Use

Use this orchestrator when:
- ✅ Testing workflows
- ✅ Processing multiple species sequentially
- ✅ Legacy SRA Toolkit environment
- ✅ Cross-species analysis needed (CSTMM, CSCA)
- ✅ Auto-activation of venv required
- ⚠️ **NOT recommended for large samples** (>50GB) - use ENA workflow instead

## Features

- **Auto-activation**: Automatically detects and activates virtual environment
- **Auto-discovery**: Finds all species configs in `config/amalgkit/`
- **Batched processing**: 10 samples at a time (~20-50 GB peak)
- **SRA optimization**: Automatic wrapper and environment setup
- **Cross-species analysis**: CSTMM and CSCA for multi-species comparisons
- **Complete pipeline**: All 11 amalgkit steps

## Usage

### Basic Usage

```bash
# Prerequisites: Ensure .venv exists with amalgkit installed
# If not set up, run:
#   python3 -m venv .venv
#   source .venv/bin/activate
#   pip install -e .
#   pip install git+https://github.com/kfuku52/amalgkit

# Processes all discovered species (automatically activates venv)
python3 scripts/rna/run_multi_species.py

# With configurable threads:
export AK_THREADS=12
python3 scripts/rna/run_multi_species.py
```

### What It Does

1. **Auto-activates virtual environment** if `.venv` exists (must be set up first)
2. **Checks environment** for required tools (amalgkit, fasterq-dump, kallisto, etc.)
3. **Discovers configs** in `config/amalgkit/amalgkit_*.yaml`
4. **Runs workflows** for each species sequentially
5. **Cross-species analysis** if ≥2 species completed

## Method Signatures

### `ensure_venv_activated()`
Automatically activate virtual environment if needed.

**Side effects**:
- Re-executes script using venv Python if needed
- Sets VIRTUAL_ENV and PATH environment variables
- Exits with instructions if venv not found

### `check_environment_or_exit()`
Check that all required tools are available.

**Verifies**:
- Virtual environment is activated
- amalgkit is installed and accessible
- SRA Toolkit (fasterq-dump) is installed
- kallisto, fastp, seqkit are installed

**Exits** with clear error message if any critical dependency is missing.

### `discover_species_configs(config_dir: Path = Path("config/amalgkit")) -> list[tuple[str, Path]]`
Discover all species configuration files.

**Parameters**:
- `config_dir`: Directory to search for config files (default: config/amalgkit)

**Returns**: List of (species_name, config_path) tuples, sorted by filename

**Notes**: Excludes template files, extracts species name from filename

### `run_species_workflow(config_path: Path, species_name: str) -> tuple[bool, Path]`
Run full workflow for a single species using batched processing.

**Parameters**:
- `config_path`: Path to species YAML configuration file
- `species_name`: Human-readable species name for logging

**Returns**: Tuple of (success, work_dir)

**Side effects**:
- Loads workflow configuration
- Executes complete amalgkit workflow (metadata → sanity)
- Creates output directories and files
- Writes logs to work_dir/logs/

### `run_cross_species_analysis(work_dirs: list[Path], output_dir: Path = Path("output/amalgkit/cross_species"))`
Run cross-species TMM normalization and correlation analysis.

**Parameters**:
- `work_dirs`: List of work directories for all species (must have completed workflows)
- `output_dir`: Output directory for cross-species results

**Side effects**:
- Runs amalgkit cstmm (Cross-Species TMM Normalization)
- Runs amalgkit csca (Cross-Species Correlation Analysis)
- Writes results to output_dir/cstmm/ and output_dir/csca/

### `main()`
Main entry point for multi-species workflow.

**Orchestrates**:
1. Auto-activate virtual environment if needed
2. Check environment and create SRA wrapper/symlinks
3. Discover all species configs in config/amalgkit/
4. Run individual species workflows (batched processing)
5. Run cross-species analysis (CSTMM, CSCA) if ≥2 species

## Environment Setup

The script automatically configures:

- **SRA temp directory**: `output/sra_temp` (instead of `/tmp`)
- **Environment variables**: `TMPDIR`, `TEMP`, `TMP`, `NCBI_VDB_QUALITY`
- **Wrapper script**: `fasterq-dump` with `--size-check off`
- **Tool symlinks**: fastp, kallisto, seqkit in wrapper directory

## Workflow Steps

For each species:
1. **metadata** - Retrieve and filter SRA metadata
2. **config** - Generate configuration files
3. **select** - Select qualified samples
4. **getfastq** - Download FASTQ files (batched, 10 samples)
5. **quant** - Quantify expression (batched with getfastq)
6. **merge** - Merge quantification results
7. **curate** - Quality control and curation
8. **cstmm** - Cross-species TMM normalization (if ≥2 species)
9. **csca** - Cross-species correlation analysis (if ≥2 species)
10. **sanity** - Validate workflow outputs

## Output Structure

```
output/amalgkit/
├── {species1}/
│   ├── work/
│   ├── quant/
│   └── logs/
├── {species2}/
│   └── ...
└── cross_species/
    ├── cstmm/
    └── csca/
```

## Performance

- **Download speed**: Slow (SRA Toolkit conversion)
- **Success rate**: ~0% for large samples (>50GB)
- **Disk usage**: Medium (~20-50 GB peak per species)
- **Processing time**: 6-13 days for 300+ samples (without optimization)

## Limitations

⚠️ **SRA Toolkit has known issues with large samples**:
- Fails with "disk-limit exceeded" errors
- Conservative size checks cause false failures
- Slow conversion from SRA to FASTQ

**Recommendation**: Use `workflow_ena_integrated.py` for production workflows.

## Troubleshooting

### "amalgkit: command not found"
Script auto-activates venv. If still failing:
```bash
source .venv/bin/activate
uv pip install git+https://github.com/kfuku52/amalgkit
```

### "fasterq-dump: command not found"
Install SRA Toolkit:
```bash
sudo apt-get install -y sra-toolkit
```

### Large sample failures
Use ENA workflow instead:
```bash
python3 scripts/rna/workflow_ena_integrated.py --config ...
```

## Related Documentation

- **[ORCHESTRATION/README.md](README.md)**: Orchestrator comparison
- **[ENA_WORKFLOW.md](ENA_WORKFLOW.md)**: Recommended ENA-based workflow
- **[WORKFLOW.md](../WORKFLOW.md)**: Workflow planning and execution

## See Also

- **Source Script**: `scripts/rna/run_multi_species.py`
- **Configuration Guide**: [CONFIGURATION.md](../CONFIGURATION.md)
- **Setup Guide**: [GETTING_STARTED.md](../GETTING_STARTED.md)

