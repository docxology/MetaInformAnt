# RNA Configuration Guide

Complete guide to configuring RNA-seq workflows, including species profiles, workflow parameters, and batch download settings.

## Quick Links

- **[API Reference](API.md#workflow-functions)** - Workflow configuration functions
- **[Workflow Guide](workflow.md)** - Workflow planning and execution
- **[Orchestration Guide](ORCHESTRATION.md)** - Orchestrator configuration
- **[Main Index](README.md)** - RNA domain master index

## Overview

Configuration in the RNA module is managed through:
- **YAML config files**: Complete workflow configuration per species
- **Species profiles**: High-level species metadata and tissue filters
- **Step parameters**: Per-step configuration overrides
- **Command-line arguments**: Runtime parameter overrides

### Path Resolution

**All paths in configuration files are resolved relative to the repository root.** This allows configs to work whether the repository is located on `/home/q/...` or `/media/q/ext6/...` (external drives).

- **Relative paths**: Resolved relative to repository root (e.g., `output/amalgkit/...`)
- **Absolute paths**: Used as-is if provided
- **Automatic detection**: Repository root is detected by looking for `.git`, `pyproject.toml`, or `.cursorrules` markers

Example:
```yaml
# These paths work on any drive location
work_dir: output/amalgkit/my_species/work
log_dir: output/amalgkit/my_species/logs
genome:
  dest_dir: output/amalgkit/my_species/genome
```

## Quick Start

### Basic Configuration

```python
from pathlib import Path
from metainformant.rna.configs import SpeciesProfile, AmalgkitRunLayout, build_step_params

spec = SpeciesProfile(name="Apis mellifera", taxon_id=7460, tissues=["head", "abdomen"])
layout = AmalgkitRunLayout(base_dir=Path("./work/Apis_mellifera"))
params = build_step_params(spec, layout)
```

### Loading from YAML

```python
from metainformant.rna.workflow import load_workflow_config

cfg = load_workflow_config("config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml")
```

## Configuration Components

### Species Profiles

High-level helpers for building per-step parameter maps.

**Classes**: `SpeciesProfile`, `AmalgkitRunLayout`

**Function**: `build_step_params(species, layout)`

```python
from pathlib import Path
from metainformant.rna.configs import SpeciesProfile, AmalgkitRunLayout, build_step_params

spec = SpeciesProfile(name="Apis mellifera", taxon_id=7460, tissues=["head", "abdomen"])
layout = AmalgkitRunLayout(base_dir=Path("./work/Apis_mellifera"))
params = build_step_params(spec, layout)
```

**Convention**: prefer `output/` for working directories, e.g. `base_dir=Path("output/amalgkit/Apis_mellifera")`.

### YAML Configuration Files

Example YAML configuration with automatic genome setup:

```yaml
# Paths are resolved relative to repository root
work_dir: output/amalgkit/camponotus_floridanus/work
log_dir: output/amalgkit/camponotus_floridanus/logs
threads: 24

auto_install_amalgkit: true

filters:
  require_tissue: false

species_list:
  - Camponotus_floridanus

# Genome configuration - automatically downloaded and indexed if missing
genome:
  accession: GCF_003227725.1
  dest_dir: output/amalgkit/camponotus_floridanus/genome
  include:
    - genome
    - gff3
    - rna
    - cds
    - protein
  ftp_url: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/227/725/...

steps:
  metadata:
    out_dir: output/amalgkit/camponotus_floridanus/work
    search_string: '"Camponotus floridanus"[Organism] AND RNA-Seq[Strategy] AND Illumina[Platform]'
    redo: yes
  getfastq:
    out_dir: output/amalgkit/camponotus_floridanus/fastq
    threads: 24
    aws: yes
    gcp: yes
    ncbi: yes
  quant:
    out_dir: output/amalgkit/camponotus_floridanus/quant
    threads: 24
    keep_fastq: no
    redo: no
    build_index: yes  # Automatically build kallisto index if missing
  merge:
    out: output/amalgkit/camponotus_floridanus/merged/merged_abundance.tsv
    out_dir: output/amalgkit/camponotus_floridanus/merged
```

**Automatic Genome Setup**: If `genome` configuration is provided, the workflow automatically:
1. Validates genome and kallisto index status at startup
2. Downloads genome package if missing
3. Prepares transcriptome FASTA if missing
4. Builds kallisto index if `build_index: yes` and index is missing
5. Only proceeds to metadata/download steps after genome/index is confirmed ready

### Environment Variable Overrides

Configuration values can be overridden using environment variables with the `AK_` prefix:

```bash
export AK_THREADS=8
export AK_WORK_DIR=/path/to/work
```

## Immediate Processing Configuration

The immediate processing system uses total thread allocation across all species. This provides fine-grained control over download throughput and system resource usage while ensuring maximum disk efficiency.

### Quick Start

**DEPRECATED**: `batch_download_species.py` has been removed. Use `run_workflow.py` with `num_download_workers` configured in each species config file.

#### Default Configuration (Recommended)
```yaml
# Configure in each species config file:
steps:
  getfastq:
    num_download_workers: 16  # Number of parallel download processes
    threads: 24
```

```bash
# Run workflow for each species
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_species1.yaml
```

#### Custom Configuration
```yaml
# More parallelism (high-end systems)
steps:
  getfastq:
    num_download_workers: 24  # More parallel downloads
    threads: 24

# Less parallelism (limited resources)
steps:
  getfastq:
    num_download_workers: 8   # Fewer parallel downloads
    threads: 12
```

### Configuration Parameters

#### `num_download_workers` (in config file)
- **Default**: 16
- **Description**: Number of parallel download processes per species
- **Range**: 4-32 (recommended: 16 for standard systems)
- **Impact**: More workers = faster downloads but higher resource usage
- **Location**: Configure in each species config file under `steps.getfastq.num_download_workers`

#### `threads` (in config file)
- **Default**: 24
- **Description**: Threads for quantification operations
- **Range**: 8-48 (recommended: 24)
- **Impact**: More threads = faster quantification but higher CPU usage

### Thread Count Recommendations

#### System Resource Considerations

| CPU Cores | Recommended Total Threads | Distribution (20 species) |
|-----------|---------------------------|---------------------------|
| 4-8 cores | 16-24 total | 1 thread per species |
| 16 cores | 24-32 total | 1-2 threads per species |
| 32 cores | 48-64 total | 2-3 threads per species |
| 64+ cores | 64-96 total | 3-4 threads per species |

### Examples

#### Standard System (Recommended)
```yaml
# Configure in config file:
steps:
  getfastq:
    num_download_workers: 16
    threads: 24
```
```bash
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_species1.yaml
```

#### High-End System
```yaml
# More parallelism:
steps:
  getfastq:
    num_download_workers: 24
    threads: 32
```

#### Limited Resources
```yaml
# Conservative settings:
steps:
  getfastq:
    num_download_workers: 8
    threads: 12
```

### Monitoring Download Progress

#### Check Active Downloads
```bash
# View active processes
ps aux | grep "amalgkit getfastq"

# Count active threads
ps aux | grep "amalgkit getfastq" | grep -c "amalgkit"
```

#### Analyze Sample Status
```bash
# Comprehensive status analysis
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_species1.yaml --status --detailed
```

### Performance Tuning

#### Signs of Optimal Configuration
- ✅ Network bandwidth is well-utilized (check with `nload` or `iftop`)
- ✅ CPU usage is moderate (not pegged at 100%)
- ✅ Downloads complete without excessive retries
- ✅ System remains responsive

#### Signs of Too Many Threads
- ⚠️ High latency/packet loss
- ⚠️ Frequent download failures
- ⚠️ System becomes unresponsive
- ⚠️ Network saturation warnings

#### Signs of Too Few Threads
- ⚠️ Downloads are slow
- ⚠️ Network bandwidth underutilized
- ⚠️ CPU idle time high

### Troubleshooting

#### Downloads Stuck or Slow
1. Check network connectivity: `ping -c 4 8.8.8.8`
2. Reduce total threads: `--total-threads 16`
3. Check disk space: `df -h /`
4. Verify immediate processing is working (only one sample's FASTQs should exist)

#### Out of Disk Space
1. Clean up partial downloads: `python3 scripts/rna/cleanup_partial_downloads.py --execute`
2. Clean up quantified samples: `python3 scripts/rna/orchestrate_workflows.py --cleanup-unquantified` (replaces cleanup_quantified_sra.sh)
3. Check for large temp files: `du -sh output/amalgkit/*/fastq`

#### High CPU Usage
1. Reduce total threads: `--total-threads 16`
2. Check other processes: `top` or `htop`
3. Monitor thread distribution: Should see threads distributed across species

## Best Practices

1. **Start Conservative**: Begin with 24 total threads (distributed across all species)
2. **Monitor Progress**: Watch resource usage and download speeds
3. **Adjust Gradually**: Increase total threads if bandwidth and CPU allow
4. **Verify Immediate Processing**: Only one sample's FASTQs should exist at a time
5. **Thread Distribution**: Threads automatically redistribute as species complete

## Related Documentation

- **[workflow.md](workflow.md)**: Workflow planning and execution
- **[Step Documentation](amalgkit/steps/README.md)**: Individual step configuration
- **[ORCHESTRATION.md](ORCHESTRATION.md)**: Orchestrator configuration
- **[GETTING_STARTED.md](GETTING_STARTED.md)**: Setup and installation

## Disk Space Configuration

### Drive Size Detection

The system automatically detects drive size categories and adjusts defaults accordingly:

- **Large drives** (> 1TB free): Optimized for high-throughput processing
- **Medium drives** (500GB-1TB free): Balanced settings
- **Small drives** (< 500GB free): Conservative settings to prevent disk exhaustion

### Batch Size Configuration

Batch sizes are automatically calculated based on available disk space, but can be overridden:

**DEPRECATED**: Batch size configuration has been removed. The workflow uses immediate per-sample processing (download → quantify → delete FASTQ) for maximum disk efficiency.

**Current approach** (recommended):
```bash
# Parallel downloads configured via num_download_workers in config file
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml
```

**Configuration**:
```yaml
steps:
  getfastq:
    num_download_workers: 16  # Number of parallel download processes
    threads: 24
```

### Batch Size Recommendations by Drive Size

| Drive Size | Default Batch Size | Max Recommended | Notes |
|------------|-------------------|-----------------|-------|
| Large (> 1TB free) | 50 | 100 | Aggressive processing, high throughput |
| Medium (500GB-1TB) | 20-30 | 50 | Balanced processing |
| Small (< 500GB) | 8-12 | 20 | Conservative, prevents disk exhaustion |

### Temporary Directory Configuration

The system automatically selects the best temporary directory location:

**Preference order:**
1. `TMPDIR` environment variable (if set)
2. External drive `output/.tmp/` (if drive is large/medium)
3. System temp directory (`/tmp` or system default)

**Manual override:**
```bash
export TMPDIR=/path/to/temp/directory
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml
```

**Check current temp directory:**
```python
from metainformant.core.disk import get_recommended_temp_dir
from pathlib import Path

repo_root = Path(".")
temp_dir = get_recommended_temp_dir(repo_root)
print(f"Temporary directory: {temp_dir}")
```

### Disk Space Thresholds

Minimum free space and auto-cleanup thresholds adapt to drive size:

**Auto-detection (recommended):**
```bash
# Workflow automatically manages disk space with immediate per-sample cleanup
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml
```

**Note**: The workflow uses immediate per-sample processing (download → quantify → delete FASTQ), so disk space is automatically managed. Only one sample's FASTQs exist at a time.

**Default thresholds by drive size:**

| Drive Size | Min Free GB | Auto-Cleanup Threshold GB |
|------------|-------------|---------------------------|
| Large (> 1TB free) | 50 | 20 |
| Medium (500GB-1TB) | 20 | 10 |
| Small (< 500GB) | 10 | 5 |

### Maximum Batch Disk Usage

**Note**: The workflow uses immediate per-sample processing (download → quantify → delete FASTQ), so batch disk usage is not applicable. Only one sample's FASTQs exist at a time, ensuring maximum disk efficiency.

### Examples

**Large drive (6TB):**
```yaml
# Configure in config file:
steps:
  getfastq:
    num_download_workers: 24  # More parallelism for large drives
    threads: 32
```
```bash
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml
```

**Medium drive (1TB):**
```yaml
# Configure in config file:
steps:
  getfastq:
    num_download_workers: 16  # Standard parallelism
    threads: 24
```
```bash
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml
```

**Small drive (500GB):**
```yaml
# Configure in config file:
steps:
  getfastq:
    num_download_workers: 8   # Conservative parallelism
    threads: 12
```
```bash
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml
```

## See Also

### Documentation
- **[API Reference](API.md#workflow-functions)** - Workflow configuration functions
- **[Workflow Guide](workflow.md)** - Workflow planning and execution
- **[Orchestration Guide](ORCHESTRATION.md)** - Orchestrator configuration
- **[Main Index](README.md)** - RNA domain master index

### Related
- **Scripts**: `scripts/rna/run_workflow.py` - Main workflow orchestrator
- **Source Code**: `src/metainformant/rna/configs.py` - Configuration module
- **Examples**: See [EXAMPLES.md](EXAMPLES.md) for real-world configurations

