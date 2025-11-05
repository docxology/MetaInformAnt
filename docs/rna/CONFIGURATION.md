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
threads: 12

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
    threads: 12
    aws: yes
    gcp: yes
    ncbi: yes
  quant:
    threads: 12
    keep_fastq: no
    redo: no
```

### Environment Variable Overrides

Configuration values can be overridden using environment variables with the `AK_` prefix:

```bash
export AK_THREADS=8
export AK_WORK_DIR=/path/to/work
```

## Batch Download Configuration

The batch download system allows you to configure how many species download in parallel and how many threads each species uses. This provides fine-grained control over download throughput and system resource usage.

### Quick Start

#### Default Configuration (Recommended)
```bash
# 3 species in parallel, 10 threads each (30 total downloads)
python3 scripts/rna/batch_download_species.py
```

#### Custom Configuration
```bash
# 4 species, 12 threads each (48 total downloads)
python3 scripts/rna/batch_download_species.py --species-count 4 --threads-per-species 12

# 2 species, 8 threads each (16 total downloads)
python3 scripts/rna/batch_download_species.py --species-count 2 --threads-per-species 8
```

### Configuration Parameters

#### `--species-count`
- **Default**: 3
- **Description**: Number of species to download in parallel
- **Range**: 1-10 (recommended: 2-5)
- **Impact**: More species = faster overall but higher resource usage

#### `--threads-per-species`
- **Default**: 10
- **Description**: Number of download threads per species
- **Range**: 4-20 (recommended: 8-12)
- **Impact**: More threads = faster per-species but diminishing returns

#### `--max-species`
- **Default**: None (all species)
- **Description**: Maximum number of species to process
- **Use case**: Testing or processing subset of species

### Thread Count Recommendations

#### Network Bandwidth Considerations

| Bandwidth | Recommended Threads | Total (3 species) |
|-----------|---------------------|------------------|
| < 100 Mbps | 8 per species | 24 total |
| 100-400 Mbps | 10 per species | 30 total |
| 400-500 Mbps | 12 per species | 36 total |
| > 500 Mbps | 16 per species | 48 total |

#### System Resource Considerations

| CPU Cores | Recommended Config | Total Threads |
|-----------|-------------------|---------------|
| 4 cores | 2 species × 8 threads | 16 total |
| 8 cores | 3 species × 10 threads | 30 total |
| 16 cores | 4 species × 12 threads | 48 total |
| 32+ cores | 5 species × 16 threads | 80 total |

### Examples

#### High-Speed Network (500+ Mbps)
```bash
python3 scripts/rna/batch_download_species.py \
    --species-count 4 \
    --threads-per-species 16
```

#### Limited Bandwidth (< 100 Mbps)
```bash
python3 scripts/rna/batch_download_species.py \
    --species-count 2 \
    --threads-per-species 8
```

#### Balanced Configuration (Default)
```bash
python3 scripts/rna/batch_download_species.py \
    --species-count 3 \
    --threads-per-species 10
```

#### Testing (Single Species)
```bash
python3 scripts/rna/batch_download_species.py \
    --species-count 1 \
    --threads-per-species 10 \
    --max-species 1
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
python3 scripts/rna/analyze_sample_status.py
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
2. Reduce threads: `--threads-per-species 8`
3. Reduce parallel species: `--species-count 2`
4. Check disk space: `df -h /`

#### Out of Disk Space
1. Clean up partial downloads: `python3 scripts/rna/cleanup_partial_downloads.py --execute`
2. Clean up quantified samples: `bash scripts/rna/cleanup_quantified_sra.sh --execute`
3. Check for large temp files: `du -sh output/amalgkit/*/fastq`

#### High CPU Usage
1. Reduce threads: `--threads-per-species 8`
2. Reduce parallel species: `--species-count 2`
3. Check other processes: `top` or `htop`

## Best Practices

1. **Start Conservative**: Begin with default (3 species × 10 threads)
2. **Monitor First Batch**: Watch resource usage and download speeds
3. **Adjust Gradually**: Increase threads/species if bandwidth allows
4. **Clean Regularly**: Clean up partial downloads to free space
5. **Use Dry Runs**: Always test cleanup with `--dry-run` first

## Related Documentation

- **[WORKFLOW.md](WORKFLOW.md)**: Workflow planning and execution
- **[STEPS.md](STEPS.md)**: Individual step configuration
- **[ORCHESTRATION/README.md](ORCHESTRATION/README.md)**: Orchestrator configuration
- **[GETTING_STARTED.md](GETTING_STARTED.md)**: Setup and installation

## See Also

- **Scripts**: `scripts/rna/batch_download_species.py` - Main batch download script
- **Source Code**: `src/metainformant/rna/configs.py` - Configuration module
- **Examples**: See `docs/rna/examples/` for real-world configurations

