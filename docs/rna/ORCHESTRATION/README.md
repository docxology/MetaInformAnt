# RNA Workflow Orchestration

Overview of all orchestrator scripts for RNA-seq workflows, including when to use each and how to configure them.

## Overview

METAINFORMANT provides several orchestrator scripts for RNA-seq workflows:

1. **`orchestrate_workflows.py`** - **UNIFIED ORCHESTRATOR** ⭐⭐⭐ (START HERE)
   - Status, monitoring, cleanup, and workflow execution
   - Consolidates all status/monitoring/management functionality
2. **`workflow_ena_integrated.py`** - Single-species ENA-based workflow (recommended)
3. **`run_multi_species.py`** - Multi-species SRA-based workflow (legacy)
4. **`run_all_species_parallel.py`** - Parallel execution of all species workflows
5. **`batch_download_species.py`** - Configurable batch download for multiple species

## Comparison Table

| Feature | orchestrate_workflows.py | workflow_ena_integrated.py | run_multi_species.py | batch_download_species.py | run_all_species_parallel.py |
|---------|-------------------------|---------------------------|---------------------|---------------------------|---------------------------|
| **Primary Use** | Status, monitoring, cleanup, step execution | Single-species ENA workflow | Multi-species SRA workflow | Configurable batch downloads | Parallel all-species execution |
| **Download Method** | N/A (orchestrates) | ENA direct (wget) | SRA Toolkit | Configurable (ENA/SRA) | SRA Toolkit |
| **Reliability** | N/A | 100% | ~0% for large samples | Depends on method | ~0% for large samples |
| **Species Support** | All (auto-discovered) | Single species | Multiple species | Multiple species | All (parallel) |
| **Status/Monitoring** | ✅ Yes (consolidates 11+ scripts) | No | No | No | No |
| **Processing Mode** | N/A | Immediate per-sample | Immediate per-sample | Immediate per-sample | Per-species |
| **Auto-cleanup** | Yes | Yes | Yes | Yes | Yes |
| **Virtual Env** | Auto-setup | Manual | Auto-activation | Auto-activation | Auto-setup |
| **Best For** | Status, monitoring, cleanup, workflow management | Production, reliability | Testing, legacy | Multi-species parallel | All-species parallel |

## Detailed Documentation

### 0. Unified Orchestrator (Start Here) ⭐⭐⭐

**Script**: `scripts/rna/orchestrate_workflows.py`

**Best for**: Status checking, monitoring, cleanup, and workflow execution

**Features**:
- Auto-discovers all species from config files
- Status and monitoring (replaces 11+ status/monitoring scripts)
- Cleanup and quantification management
- Flexible step execution for any species
- Resume incomplete downloads

**Usage**:
```bash
# Status and monitoring
python3 scripts/rna/orchestrate_workflows.py --status                    # Brief status
python3 scripts/rna/orchestrate_workflows.py --status --detailed         # Detailed status
python3 scripts/rna/orchestrate_workflows.py --monitor                   # Real-time monitoring

# Cleanup and workflow execution
python3 scripts/rna/orchestrate_workflows.py --cleanup-unquantified      # Quantify and cleanup
python3 scripts/rna/orchestrate_workflows.py --resume-downloads          # Resume downloads
python3 scripts/rna/orchestrate_workflows.py --steps merge curate sanity --auto-species
```

**See**: [scripts/rna/README.md](../../scripts/rna/README.md) for complete documentation.

### 1. ENA Workflow (Recommended)

**Script**: `scripts/rna/workflow_ena_integrated.py`

**Best for**: Production workflows requiring maximum reliability

**Features**:
- Direct ENA downloads with 100% reliability
- Batched processing (configurable batch size)
- Automatic resume with `wget --continue`
- Auto-cleanup (FASTQs deleted after quantification)

**Documentation**: [ENA_WORKFLOW.md](ENA_WORKFLOW.md)

**Usage**:
```bash
python3 scripts/rna/workflow_ena_integrated.py \
  --config config/amalgkit/amalgkit_cfloridanus.yaml \
  --batch-size 12 \
  --threads 12
```

### 2. Multi-Species SRA Workflow

**Script**: `scripts/rna/run_multi_species.py`

**Best for**: Testing, legacy environments, multi-species coordination

**Features**:
- Auto-activation of virtual environment
- Auto-discovery of all species configs
- Batched processing (10 samples at a time)
- SRA Toolkit optimization
- Cross-species analysis (CSTMM, CSCA)

**Documentation**: [MULTI_SPECIES.md](MULTI_SPECIES.md)

**Usage**:
```bash
# Prerequisites: Ensure .venv exists with amalgkit installed
# See GETTING_STARTED.md or RUN_ALL_SPECIES.md for setup

# All species with default threads (10 per species)
python3 scripts/rna/run_multi_species.py

# All species with configurable threads (via environment variable)
export AK_THREADS=12
python3 scripts/rna/run_multi_species.py
```

**See [RUN_ALL_SPECIES.md](../RUN_ALL_SPECIES.md) for complete guide including prerequisites and setup instructions.**

### 3. Batch Download

**Script**: `scripts/rna/batch_download_species.py`

**Best for**: Multi-species parallel downloads with configurable throughput

**Features**:
- Configurable parallelism (species count × threads per species)
- Dynamic thread allocation
- Immediate per-sample quantification
- Disk space monitoring and cleanup

**Documentation**: [BATCH_DOWNLOAD.md](BATCH_DOWNLOAD.md)

**Usage**:
```bash
# Default: 3 species × 10 threads = 30 total downloads
python3 scripts/rna/batch_download_species.py

# Custom: 4 species × 12 threads = 48 total downloads
python3 scripts/rna/batch_download_species.py \
  --species-count 4 \
  --threads-per-species 12
```

## Configuration

All orchestrators support configuration through:

1. **YAML config files**: Primary configuration source
2. **Command-line arguments**: Runtime overrides
3. **Environment variables**: System-wide defaults

See [CONFIGURATION.md](../CONFIGURATION.md) for detailed configuration options.

## Method Signatures

### workflow_ena_integrated.py

```python
def load_config(config_path: Path) -> dict:
    """Load amalgkit YAML config."""
    
def get_sample_list(metadata_path: Path) -> list[str]:
    """Get list of run IDs from metadata file."""
    
def sample_already_quantified(run_id: str, quant_dir: Path) -> bool:
    """Check if sample has abundance.tsv."""
    
def download_batch_ena(run_ids: list[str], metadata_path: Path, 
                       fastq_dir: Path, threads: int, batch_num: int) -> tuple[list[str], list[str]]:
    """Download a batch of samples using the robust ENA downloader."""
    
def quantify_batch_kallisto(run_ids: list[str], fastq_dir: Path, 
                            quant_dir: Path, index_path: Path, 
                            threads: int, batch_num: int) -> tuple[list[str], list[str]]:
    """Quantify samples using kallisto."""
```

### run_multi_species.py

```python
def ensure_venv_activated():
    """Automatically activate virtual environment if needed."""
    
def setup_sra_environment(repo_root: Path):
    """Configure SRA Toolkit environment for large downloads."""
    
def find_species_configs(config_dir: Path) -> list[Path]:
    """Find all amalgkit_*.yaml config files."""
    
def run_species_workflow(config_path: Path, repo_root: Path):
    """Run complete workflow for a single species."""
```

### batch_download_species.py

```python
def ensure_venv_activated():
    """Automatically activate virtual environment if needed."""
    
def check_disk_space(min_free_gb: float) -> bool:
    """Check if sufficient disk space is available."""
    
def cleanup_partial_downloads(output_dir: Path, threshold_gb: float):
    """Clean up partial/failed downloads when space is low."""
    
def allocate_threads(species_configs: list[Path], total_threads: int) -> dict[Path, int]:
    """Allocate threads evenly across species."""
```

## Choosing the Right Orchestrator

### Use `workflow_ena_integrated.py` when:
- ✅ You need maximum reliability (100% success rate)
- ✅ Processing single species
- ✅ Production workflows
- ✅ Large samples (>50GB)
- ✅ Network interruptions are common

### Use `run_multi_species.py` when:
- ✅ Testing workflows
- ✅ Processing multiple species sequentially
- ✅ Legacy SRA Toolkit environment
- ✅ Cross-species analysis needed
- ✅ Auto-activation of venv required

### Use `batch_download_species.py` when:
- ✅ Processing multiple species in parallel with immediate per-sample processing
- ✅ Need total thread allocation (not per-species)
- ✅ Disk space is limited (maximum efficiency: only one sample's FASTQs at a time)
- ✅ Dynamic resource allocation needed (threads redistribute as species complete)
- ✅ Maximum throughput with controlled resource usage

## Performance Comparison

| Metric | workflow_ena_integrated.py | run_multi_species.py | batch_download_species.py |
|--------|---------------------------|---------------------|---------------------------|
| **Download Speed** | Fast (ENA direct) | Slow (SRA Toolkit) | Fast (configurable) |
| **Success Rate** | 100% | ~0% (large samples) | Depends on method |
| **Disk Usage** | Low (auto-cleanup) | Low (immediate delete) | Minimal (one sample at a time) |
| **CPU Usage** | Moderate | High (SRA conversion) | Configurable |
| **Parallelism** | Per-species | Per-species | Multi-species |

## Troubleshooting

### Common Issues

**Downloads failing**: 
- Use `workflow_ena_integrated.py` for ENA-based downloads
- Check network connectivity
- Verify NCBI_EMAIL is set

**Disk space issues**:
- Use `batch_download_species.py` with automatic cleanup
- Reduce batch size or threads
- Enable auto-cleanup thresholds

**Virtual environment issues**:
- Use `run_multi_species.py` or `batch_download_species.py` (auto-activation)
- Or manually activate: `source .venv/bin/activate`

**Performance issues**:
- Adjust thread counts based on bandwidth and CPU
- Use `batch_download_species.py` for fine-grained control
- Monitor resource usage and adjust accordingly

## Related Documentation

- **[CONFIGURATION.md](../CONFIGURATION.md)**: Configuration management
- **[WORKFLOW.md](../WORKFLOW.md)**: Workflow planning and execution
- **[MULTI_SPECIES_QUICK_START.md](../MULTI_SPECIES_QUICK_START.md)**: Production workflows
- **[GETTING_STARTED.md](../GETTING_STARTED.md)**: Setup and installation

## See Also

- **Source Scripts**: `scripts/rna/` - Implementation details
- **Module Documentation**: `src/metainformant/rna/README.md` - API reference
- **Examples**: `docs/rna/examples/` - Real-world usage examples

