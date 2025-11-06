# RNA Configuration Guide

Complete guide to configuring RNA-seq workflows, including species profiles, workflow parameters, and batch download settings.

## Overview

Configuration in the RNA module is managed through:
- **YAML config files**: Complete workflow configuration per species
- **Species profiles**: High-level species metadata and tissue filters
- **Step parameters**: Per-step configuration overrides
- **Command-line arguments**: Runtime parameter overrides

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

cfg = load_workflow_config("config/amalgkit/amalgkit_cfloridanus.yaml")
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

Example YAML configuration:

```yaml
work_dir: output/amalgkit/camponotus_floridanus/work
log_dir: output/amalgkit/camponotus_floridanus/logs
threads: 24

auto_install_amalgkit: true

filters:
  require_tissue: false

species_list:
  - Camponotus_floridanus

steps:
  metadata:
    out_dir: output/amalgkit/camponotus_floridanus/work
    search_string: '"Camponotus floridanus"[Organism] AND RNA-Seq[Strategy] AND Illumina[Platform]'
    redo: yes
  getfastq:
    threads: 24
    aws: yes
    gcp: yes
    ncbi: yes
  quant:
    threads: 24
    keep_fastq: no
    redo: no
```

### Environment Variable Overrides

Configuration values can be overridden using environment variables with the `AK_` prefix:

```bash
export AK_THREADS=8
export AK_WORK_DIR=/path/to/work
```

## Immediate Processing Configuration

The immediate processing system uses total thread allocation across all species. This provides fine-grained control over download throughput and system resource usage while ensuring maximum disk efficiency.

### Quick Start

#### Default Configuration (Recommended)
```bash
# 24 threads TOTAL distributed across all species (minimum 1 per species)
python3 scripts/rna/batch_download_species.py --total-threads 24
```

#### Custom Configuration
```bash
# 48 threads total (more throughput on high-end systems)
python3 scripts/rna/batch_download_species.py --total-threads 48

# 16 threads total (conservative for limited resources)
python3 scripts/rna/batch_download_species.py --total-threads 16
```

### Configuration Parameters

#### `--total-threads`
- **Default**: 24
- **Description**: Total number of threads to distribute across all species
- **Range**: 8-64 (recommended: 24 for standard systems)
- **Impact**: Threads distributed evenly (minimum 1 per species). More threads = faster overall but higher resource usage
- **Distribution**: For 20 species with 24 threads: 4 species get 2 threads, 16 get 1 thread

#### `--max-species`
- **Default**: None (all species)
- **Description**: Maximum number of species to process
- **Use case**: Testing or processing subset of species

#### `--quant-threads`
- **Default**: 10
- **Description**: Separate thread pool for quantification operations
- **Range**: 4-20 (recommended: 10)
- **Impact**: Quantification runs in parallel with downloads

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
```bash
# Default: 24 threads total distributed across all species
python3 scripts/rna/batch_download_species.py
```

#### High-End System
```bash
# 48 threads total for faster processing
python3 scripts/rna/batch_download_species.py --total-threads 48
```

#### Limited Resources
```bash
# 16 threads total for conservative resource usage
python3 scripts/rna/batch_download_species.py --total-threads 16
```

#### Testing (Limited Species)
```bash
# Process only first 5 species with 24 threads total
python3 scripts/rna/batch_download_species.py --total-threads 24 --max-species 5
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
# Comprehensive status analysis (replaces analyze_sample_status.py)
python3 scripts/rna/orchestrate_workflows.py --status --detailed
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

- **[WORKFLOW.md](WORKFLOW.md)**: Workflow planning and execution
- **[STEPS.md](STEPS.md)**: Individual step configuration
- **[ORCHESTRATION/README.md](ORCHESTRATION/README.md)**: Orchestrator configuration
- **[GETTING_STARTED.md](GETTING_STARTED.md)**: Setup and installation

## Disk Space Configuration

### Drive Size Detection

The system automatically detects drive size categories and adjusts defaults accordingly:

- **Large drives** (> 1TB free): Optimized for high-throughput processing
- **Medium drives** (500GB-1TB free): Balanced settings
- **Small drives** (< 500GB free): Conservative settings to prevent disk exhaustion

### Batch Size Configuration

Batch sizes are automatically calculated based on available disk space, but can be overridden:

**Auto-detection (recommended):**
```bash
# Batch size automatically calculated from available space
python3 scripts/rna/workflow_ena_integrated.py --config config/amalgkit/amalgkit_cfloridanus.yaml
```

**Manual override:**
```bash
# Specify batch size explicitly
python3 scripts/rna/workflow_ena_integrated.py \
  --config config/amalgkit/amalgkit_cfloridanus.yaml \
  --batch-size 50 \
  --max-batch-size 100
```

**Environment variable:**
```bash
export AK_BATCH_SIZE=50
python3 scripts/rna/workflow_ena_integrated.py --config config/amalgkit/amalgkit_cfloridanus.yaml
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
python3 scripts/rna/workflow_ena_integrated.py --config ...
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
# Thresholds automatically set based on drive size
python3 scripts/rna/batch_download_species.py
```

**Manual override:**
```bash
python3 scripts/rna/batch_download_species.py \
  --min-free-gb 50.0 \
  --auto-cleanup-threshold 20.0
```

**Default thresholds by drive size:**

| Drive Size | Min Free GB | Auto-Cleanup Threshold GB |
|------------|-------------|---------------------------|
| Large (> 1TB free) | 50 | 20 |
| Medium (500GB-1TB) | 20 | 10 |
| Small (< 500GB) | 10 | 5 |

### Maximum Batch Disk Usage

Limit disk space usage per batch:

```bash
python3 scripts/rna/workflow_ena_integrated.py \
  --config config/amalgkit/amalgkit_cfloridanus.yaml \
  --max-batch-disk-gb 150.0
```

If batch size would exceed this limit, it's automatically adjusted downward.

### Examples

**Large drive (6TB):**
```bash
# Auto-detected settings (recommended)
export AK_THREADS=24
python3 scripts/rna/run_all_species_parallel.py --threads-per-species 24
# Batch size: 50 (auto-detected)
# Min free: 50GB (auto-detected)
# Temp dir: output/.tmp/ (auto-detected)
```

**Medium drive (1TB):**
```bash
# Auto-detected settings
python3 scripts/rna/workflow_ena_integrated.py \
  --config config/amalgkit/amalgkit_cfloridanus.yaml
# Batch size: 25 (auto-detected)
# Min free: 20GB (auto-detected)
```

**Small drive (500GB):**
```bash
# Conservative settings
python3 scripts/rna/workflow_ena_integrated.py \
  --config config/amalgkit/amalgkit_cfloridanus.yaml \
  --batch-size 10 \
  --min-free-gb 10
```

## See Also

- **Scripts**: `scripts/rna/batch_download_species.py` - Main batch download script
- **Source Code**: `src/metainformant/rna/configs.py` - Configuration module
- **Examples**: See `docs/rna/examples/` for real-world configurations

