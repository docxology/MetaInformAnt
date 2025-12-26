# RNA Workflow Orchestration

Overview of orchestrator scripts for RNA-seq workflows, including when to use each and how to configure them.

## Quick Links

- **[API Reference](API.md#orchestration-functions)** - Orchestration function documentation
- **[Workflow Guide](workflow.md)** - Workflow planning and execution
- **[Configuration Guide](CONFIGURATION.md)** - Configuration management
- **[Main Index](README.md)** - RNA domain master index

## Overview

METAINFORMANT provides a unified orchestrator script for RNA-seq workflows:

1. **`run_workflow.py`** - **MAIN ORCHESTRATOR** ⭐⭐⭐ (RECOMMENDED)
   - Complete end-to-end workflow execution
   - Status checking and monitoring
   - Cleanup operations
   - Per-sample processing: download → quantify → delete FASTQ
   - Parallel downloads via `num_download_workers` configuration
   - All 11 amalgkit steps in correct order

**DEPRECATED**: The following scripts have been removed:
- `workflow_ena_integrated.py` - Use `run_workflow.py` instead
- `batch_download_species.py` - Use `run_workflow.py` with `num_download_workers` in config
- `run_multi_species.py` - Use `run_workflow.py` separately for each species config

## Main Orchestrator: run_workflow.py

| Feature | run_workflow.py |
|---------|------------------|
| **Primary Use** | Complete end-to-end workflow execution |
| **Download Method** | ENA direct (wget) - 100% reliability |
| **Reliability** | 100% (ENA-based downloads) |
| **Species Support** | Single species per execution (run separately for multiple species) |
| **Status/Monitoring** | ✅ Yes (--status, --detailed flags) |
| **Processing Mode** | Immediate per-sample (download → quantify → delete FASTQ) |
| **Parallel Downloads** | Configurable via `num_download_workers` in config (default: 16) |
| **Auto-cleanup** | Yes (automatic FASTQ deletion after quantification) |
| **Virtual Env** | Auto-setup and activation |
| **Best For** | Production workflows, single-species processing, maximum reliability |

## Main Orchestrator: run_workflow.py ⭐⭐⭐

**Script**: `scripts/rna/run_workflow.py`

**Best for**: All end-to-end RNA-seq workflows

**Features**:
- Complete end-to-end execution via `execute_workflow()`
- Automatic genome setup (if genome config exists)
- Per-sample processing: download → quantify → delete FASTQ
- All 11 amalgkit steps in correct order
- Status checking and cleanup operations
- Parallel downloads via `num_download_workers` configuration
- ENA-based downloads (100% reliability)

**Usage**:
```bash
# Full end-to-end workflow (all steps)
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml

# Specific steps only
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps getfastq quant merge

# Check status
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --status

# Detailed status
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --status --detailed

# Cleanup unquantified samples
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --cleanup-unquantified
```

**Configuration for parallel downloads**:
```yaml
steps:
  getfastq:
    num_download_workers: 16  # Number of parallel download processes
    threads: 24
```

**See**: [scripts/rna/README.md](../../scripts/rna/README.md) for complete documentation.

## ENA Workflow (Recommended)

**DEPRECATED**: This section describes functionality that has been consolidated into `run_workflow.py`.

**Current Recommendation**: Use `run_workflow.py` for all end-to-end workflows, which provides ENA-based downloads by default.

### Migration

**Old approach** (removed):
```bash
python3 scripts/rna/workflow_ena_integrated.py --config config.yaml
```

**New approach** (recommended):
```bash
python3 scripts/rna/run_workflow.py config.yaml
```

### Features

The `run_workflow.py` script provides all the functionality previously in `workflow_ena_integrated.py`:

- **Direct ENA Downloads**: 100% reliability using ENA API via wget
- **Parallel Downloads**: Configurable via `num_download_workers` in config file
- **Automatic Resume**: wget --continue for interrupted downloads
- **Auto-cleanup**: FASTQs deleted after quantification
- **Disk Management**: Immediate per-sample processing (only one sample's FASTQs at a time)

### Configuration

Configure parallel downloads in your config file:

```yaml
steps:
  getfastq:
    out_dir: output/amalgkit/{species}/fastq
    threads: 24
    num_download_workers: 16  # Number of parallel download processes
    aws: yes     # Use AWS Open Data Program
    gcp: yes     # Use Google Cloud
    ncbi: yes    # Use NCBI directly
    fastp: yes   # ENABLED: fastp QC tool
```

### Performance

- **Download speed**: Fast (ENA direct, no SRA conversion)
- **Success rate**: 100% (vs ~0% for SRA Toolkit on large samples)
- **Disk usage**: Low (auto-cleanup, immediate per-sample processing)
- **Processing time**: ~7.5 minutes per sample average

## Multi-Species Processing

**DEPRECATED**: This section describes functionality that has been consolidated into `run_workflow.py`.

**Current Recommendation**: Use `run_workflow.py` separately for each species config. For multiple species, run it separately for each species, either sequentially or in parallel.

### Migration

**Old approach** (removed):
```bash
python3 scripts/rna/run_multi_species.py
```

**New approach** (recommended):
```bash
# Run separately for each species
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_species1.yaml
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_species2.yaml
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_species3.yaml
```

**For parallel execution** (multiple species at once):
```bash
# Run in background for each species
nohup python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_species1.yaml > logs/species1.log 2>&1 &
nohup python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_species2.yaml > logs/species2.log 2>&1 &
nohup python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_species3.yaml > logs/species3.log 2>&1 &
```

### Features

The `run_workflow.py` script provides all the functionality previously in `run_multi_species.py`:

- **Complete Pipeline**: All 11 amalgkit steps (metadata → sanity)
- **Auto-activation**: Virtual environment detection and activation
- **Per-Sample Processing**: Immediate download → quantify → delete FASTQ
- **Parallel Downloads**: Configurable via `num_download_workers` in config file
- **Resume Support**: Automatically skips completed steps

### Configuration

Each species has its own config file with parallel download settings:

```yaml
steps:
  getfastq:
    num_download_workers: 16  # Parallel downloads for this species
    threads: 24
```

### Cross-Species Analysis

For cross-species analysis (CSTMM, CSCA), run the workflow for all species first, then run the cross-species steps:

```bash
# Process all species
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_species1.yaml
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_species2.yaml

# Run cross-species analysis (if needed)
# The workflow automatically handles cross-species steps when multiple species are complete
```

## Batch Download Configuration

**DEPRECATED**: This section describes functionality that has been consolidated into `run_workflow.py`.

**Current Recommendation**: Use `run_workflow.py` with `num_download_workers` configured in each species config file for parallel downloads.

### Migration

**Old approach** (removed):
```bash
python3 scripts/rna/batch_download_species.py --total-threads 24
```

**New approach** (recommended):
```yaml
# Configure parallel downloads in each species config file:
steps:
  getfastq:
    num_download_workers: 16  # Number of parallel download processes
    threads: 24
```

Then run `run_workflow.py` separately for each species:
```bash
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_species1.yaml
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_species2.yaml
```

For multiple species, run workflows in parallel using `nohup` or `screen`:
```bash
# Run in background for each species
nohup python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_species1.yaml > logs/species1.log 2>&1 &
nohup python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_species2.yaml > logs/species2.log 2>&1 &
```

### Features

The `run_workflow.py` script provides all the functionality previously in `batch_download_species.py`:

- **Parallel Downloads**: Configurable via `num_download_workers` in config file
- **Immediate Quantification**: Each sample quantified and FASTQs deleted as soon as download completes
- **Disk Space Management**: Immediate per-sample processing (only one sample's FASTQs at a time)
- **Auto-cleanup**: Automatic FASTQ deletion after quantification

### Configuration

Configure parallel downloads per species:

```yaml
steps:
  getfastq:
    num_download_workers: 16  # Parallel downloads for this species
    threads: 24
```

### Thread Allocation

- Each species config has its own `num_download_workers` setting
- For multiple species, run `run_workflow.py` separately for each
- To control total system load, adjust `num_download_workers` in each config file
- Run workflows in parallel (using `nohup` or `screen`) or sequentially

## Configuration

All orchestrators support configuration through:

1. **YAML config files**: Primary configuration source
2. **Command-line arguments**: Runtime overrides
3. **Environment variables**: System-wide defaults

See [CONFIGURATION.md](CONFIGURATION.md) for detailed configuration options.

## Per-Sample Processing Workflow

The orchestrator implements an immediate per-sample processing workflow:

1. **Download**: FASTQ files are downloaded for each sample
2. **Quantify**: Sample is immediately quantified against the kallisto index
3. **Delete**: FASTQ files are automatically deleted after successful quantification

This workflow ensures maximum disk efficiency - only one sample's FASTQ files exist at any time.

### Manual Per-Sample Processing

For processing individual samples manually (e.g., testing or recovery):

```python
from metainformant.rna.steps.quant import quantify_sample
from metainformant.rna.steps.getfastq import delete_sample_fastqs
from metainformant.rna.workflow import load_workflow_config
from metainformant.core.io import read_delimited
from pathlib import Path

# Load config and metadata
cfg = load_workflow_config("config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml")
metadata_file = cfg.work_dir / "metadata" / "metadata.tsv"
rows = list(read_delimited(metadata_file, delimiter="\t"))
sample_rows = [row for row in rows if row.get("run") == "SRR14740514"]

# Prepare quant params
quant_params = dict(cfg.per_step.get("quant", {}))
quant_params["out_dir"] = str(cfg.per_step.get("quant", {}).get("out_dir", cfg.work_dir / "quant"))
quant_params["threads"] = cfg.threads or 12
quant_params["work_dir"] = str(cfg.work_dir.absolute())

# Quantify
success, message, abundance_path = quantify_sample(
    sample_id="SRR14740514",
    metadata_rows=sample_rows,
    quant_params=quant_params,
    log_dir=cfg.log_dir,
)

# Delete FASTQ files after successful quantification
if success and abundance_path and abundance_path.exists():
    fastq_dir = Path(cfg.per_step.get("getfastq", {}).get("out_dir", cfg.work_dir / "fastq"))
    delete_sample_fastqs("SRR14740514", fastq_dir)
```

**Test Script**: A complete test script is available:

```bash
python3 scripts/rna/test_quantify_sample.py \
    --sample SRR14740514 \
    --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml
```

### Batch Cleanup

The `cleanup_unquantified_samples()` function processes all downloaded but unquantified samples:

```python
from metainformant.rna.orchestration import cleanup_unquantified_samples
from pathlib import Path

config_path = Path("config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml")
quantified, failed = cleanup_unquantified_samples(config_path)
```

This function:
- Finds all samples with FASTQ files but no quantification results
- Quantifies each sample using `quantify_sample()`
- Deletes FASTQ files after successful quantification using `delete_sample_fastqs()`
- Returns counts of quantified and failed samples

## Method Signatures

### run_workflow.py

The `run_workflow.py` script calls methods from `metainformant.rna.orchestration`:

```python
from metainformant.rna.orchestration import run_workflow_for_species

# Main function called by run_workflow.py
def run_workflow_for_species(
    config_path: Path,
    steps: list[str] | None = None,
    check: bool = False,
) -> dict[str, Any]:
    """Run complete workflow for a single species."""
    
# Status checking
from metainformant.rna.monitoring import check_workflow_progress, analyze_species_status

def check_workflow_progress(config_path: Path) -> dict[str, Any]:
    """Get workflow progress summary."""
    
# Cleanup operations
from metainformant.rna.orchestration import cleanup_unquantified_samples
from metainformant.rna.cleanup import cleanup_partial_downloads

def cleanup_unquantified_samples(config_path: Path) -> tuple[int, int]:
    """Quantify downloaded samples and cleanup FASTQs.
    
    Finds all samples with FASTQ files but no quantification results,
    quantifies each sample, and deletes FASTQ files after successful quantification.
    """
```

See [API.md](API.md#orchestration-functions) for complete function documentation.

## When to Use run_workflow.py

Use `run_workflow.py` for:
- ✅ All end-to-end RNA-seq workflows
- ✅ Single-species processing (recommended)
- ✅ Production workflows requiring maximum reliability
- ✅ Large samples (>50GB)
- ✅ Network interruptions are common
- ✅ Parallel downloads (configure via `num_download_workers`)

**For multiple species**: Run `run_workflow.py` separately for each species config, either sequentially or in parallel (using `nohup` or `screen`).

## Performance

| Metric | run_workflow.py |
|--------|-----------------|
| **Download Speed** | Fast (ENA direct, wget-based) |
| **Success Rate** | 100% (ENA-based downloads) |
| **Disk Usage** | Low (immediate per-sample cleanup) |
| **CPU Usage** | Configurable (via threads and num_download_workers) |
| **Parallelism** | Configurable per species (num_download_workers in config) |

## Troubleshooting

### Common Issues

**Downloads failing**: 
- Check network connectivity
- Verify NCBI_EMAIL is set
- Check that `num_download_workers` is configured in config file

**Disk space issues**:
- Use `--cleanup-unquantified` to cleanup downloaded but unquantified samples
- Use `--cleanup-partial` to remove partial downloads
- Reduce `num_download_workers` if disk space is limited

**Virtual environment issues**:
- Scripts automatically detect and activate venv (`.venv` or `/tmp/metainformant_venv`)
- If issues persist, manually activate: `source .venv/bin/activate`

**Performance issues**:
- Adjust `num_download_workers` in config file based on bandwidth and CPU
- Monitor resource usage: `python3 scripts/rna/run_workflow.py <config> --status`
- Check active processes: `ps aux | grep "run_workflow\|amalgkit" | grep -v grep`

## Related Documentation

- **[CONFIGURATION.md](CONFIGURATION.md)**: Configuration management
- **[workflow.md](workflow.md)**: Workflow planning and execution
- **[GETTING_STARTED.md](GETTING_STARTED.md)**: Setup and installation

## See Also

### Documentation
- **[API Reference](API.md#orchestration-functions)** - Complete function documentation
- **[Function Index](amalgkit/FUNCTIONS.md)** - Quick function lookup
- **[Workflow Guide](workflow.md)** - Workflow planning and execution
- **[Configuration Guide](CONFIGURATION.md)** - Configuration management
- **[Main Index](README.md)** - RNA domain master index

### Related
- **Source Scripts**: `scripts/rna/` - Implementation details
- **Module Documentation**: `src/metainformant/rna/README.md` - API reference
- **Examples**: [EXAMPLES.md](EXAMPLES.md) - Real-world usage examples

