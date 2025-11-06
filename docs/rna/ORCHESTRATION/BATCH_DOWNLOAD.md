# Batch Download Orchestrator

**Script**: `scripts/rna/batch_download_species.py`

## Overview

Configurable batch download orchestrator for multiple species with dynamic thread allocation and immediate per-sample quantification. Best for multi-species parallel downloads with fine-grained control.

## Configuration Block

```python
# ============================================================================
# CONFIGURATION
# ============================================================================
# Scope: Multi-species parallel batch download with configurable throughput
# Steps: getfastq → quant (per-sample, immediate)
# Config: Auto-discovers all config/amalgkit/amalgkit_*.yaml files
# Threads: Configurable total (default 24) distributed across species
# Batch Size: Per-sample processing (download → quant → delete immediately)
# Output: output/amalgkit/{species}/fastq/ and quant/ per species
# Dependencies: amalgkit, kallisto, wget (for ENA) or SRA Toolkit
# Virtual Env: Auto-activates if .venv exists
# Disk Management: Pre-flight checks, auto-cleanup when space low
# Reliability: Depends on download method (ENA recommended)
# ============================================================================
```

## When to Use

Use this orchestrator when:
- ✅ Processing multiple species in parallel
- ✅ Need fine-grained control over parallelism
- ✅ Disk space is limited
- ✅ Dynamic resource allocation needed
- ✅ Immediate quantification per sample required
- ✅ Want to monitor and adjust throughput

## Features

- **Configurable parallelism**: Adjust total threads and distribution
- **Dynamic allocation**: Threads redistribute as species complete
- **Immediate quantification**: Each sample quantified and FASTQs deleted as soon as download completes
- **Disk space monitoring**: Pre-flight checks and automatic cleanup
- **Process health checks**: Automatic recovery from failures
- **Auto-discovery**: Finds all species configs automatically

## Usage

### Basic Usage

```bash
# Default: 24 threads total, distributed across all species
python3 scripts/rna/batch_download_species.py

# Larger drive: 48 threads
python3 scripts/rna/batch_download_species.py --total-threads 48

# Custom disk space thresholds
python3 scripts/rna/batch_download_species.py \
  --min-free-gb 15.0 \
  --auto-cleanup-threshold 10.0
```

### Command-Line Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `--total-threads` | int | 24 | Total threads across all species |
| `--max-species` | int | None | Maximum species to process |
| `--check-interval` | int | 300 | Seconds between completion checks |
| `--monitor-interval` | int | 60 | Seconds between sample monitoring |
| `--quant-threads` | int | 10 | Threads for quantification pool |
| `--max-retries` | int | 3 | Maximum retries for failed processes |
| `--process-timeout` | int | 3600 | Process timeout in seconds |
| `--health-check-interval` | int | 600 | Health check interval in seconds |
| `--min-free-gb` | float | 10.0 | Minimum free disk space GB |
| `--auto-cleanup-threshold` | float | 5.0 | GB threshold for auto-cleanup |

## Method Signatures

### `is_sample_download_complete(sample_dir: Path) -> bool`
Check if a sample directory has completed download with valid FASTQ files.

**Parameters**:
- `sample_dir`: Path to sample directory

**Returns**: True if sample has valid FASTQ files (non-empty .fastq.gz files, minimum 1KB each)

### `is_sra_download_complete(sample_dir: Path) -> bool`
Check if a sample has completed SRA download (ready for conversion to FASTQ).

**Parameters**:
- `sample_dir`: Path to sample directory

**Returns**: True if sample has SRA file that appears complete (size > 100 MB, < 100 GB, not modified in last 5 minutes)

### `detect_completed_downloads(species_configs: list[tuple[str, Path]], processed_samples: set[str], quant_dirs: dict[str, Path]) -> list[tuple[str, str, Path, str]]`
Detect samples that have completed download but not yet been quantified.

**Parameters**:
- `species_configs`: List of (species_name, config_path) tuples
- `processed_samples`: Set of sample IDs already being processed
- `quant_dirs`: Dictionary mapping species_name -> quant_dir Path

**Returns**: List of (species_name, sample_id, sample_dir, status) tuples

**Side effects**: Updates processed_samples set

### `main()`
Main entry point for batch download orchestrator.

**Orchestrates**:
1. Auto-discovery of species configs
2. Dynamic thread allocation across species
3. Per-sample immediate quantification and cleanup
4. Disk space monitoring and automatic cleanup
5. Process health checks and recovery

## Thread Allocation

### Initial Distribution

- Threads distributed evenly across species (minimum 1 per species)
- Example: 24 threads, 20 species → 4 species get 2 threads, 16 get 1 thread

### Dynamic Reallocation

- As species complete, threads redistribute to remaining species
- Threads concentrate on fewer species as workflow progresses
- Example: 24 threads, 2 species remaining → [12, 12] threads

## Disk Space Management

### Pre-flight Checks

- Verifies sufficient disk space before starting (default: 10GB minimum)
- Exits with error if insufficient space

### Automatic Cleanup

- Removes partial/failed downloads when space drops below threshold (default: 5GB)
- Monitors disk space every 10 minutes during execution
- Automatically attempts cleanup and retry on disk space errors

## Workflow

1. **Downloads** run in parallel across species with dynamic thread allocation
2. **Every 60 seconds** (configurable), monitor detects completed downloads:
   - Completed SRA files → SRA→FASTQ conversion
   - Completed FASTQ files → Direct quantification
3. **Pipeline per sample**: SRA→FASTQ (if needed) → Quant → Delete
4. **Processes ALL species** (including backlogged completed downloads)
5. **No waiting** for batch completion - maximum disk efficiency

## Performance

- **Download speed**: Fast (configurable, ENA recommended)
- **Success rate**: Depends on download method (ENA: 100%, SRA: ~0% for large)
- **Disk usage**: Low (per-sample cleanup, immediate deletion)
- **CPU usage**: Configurable (adjust total-threads)

## Thread Count Recommendations

### Network Bandwidth

| Bandwidth | Recommended Total Threads |
|-----------|--------------------------|
| < 100 Mbps | 8-16 |
| 100-400 Mbps | 16-30 |
| 400-500 Mbps | 30-48 |
| > 500 Mbps | 48+ |

### System Resources

| CPU Cores | Recommended Total Threads |
|-----------|--------------------------|
| 4 cores | 8-12 |
| 8 cores | 16-24 |
| 16 cores | 30-48 |
| 32+ cores | 48+ |

## Troubleshooting

### Downloads Stuck or Slow
1. Check network connectivity
2. Reduce threads: `--total-threads 16`
3. Check disk space: `df -h /`

### Out of Disk Space
- Automatic cleanup enabled (handles this)
- Or manually: `python3 scripts/rna/cleanup_partial_downloads.py --execute`

### High CPU Usage
- Reduce threads: `--total-threads 16`
- Check other processes: `top` or `htop`

## Related Documentation

- **[ORCHESTRATION/README.md](README.md)**: Orchestrator comparison
- **[CONFIGURATION.md](../CONFIGURATION.md)**: Detailed configuration guide
- **[ENA_WORKFLOW.md](ENA_WORKFLOW.md)**: Single-species ENA workflow

## See Also

- **Source Script**: `scripts/rna/batch_download_species.py`
- **Configuration Guide**: [CONFIGURATION.md](../CONFIGURATION.md#batch-download-configuration)
- **Performance Tuning**: [CONFIGURATION.md](../CONFIGURATION.md#performance-tuning)

