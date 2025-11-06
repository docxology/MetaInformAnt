# Batch Download Configuration Guide

## Overview

The batch download system uses **total thread allocation** across all species, providing fine-grained control over download throughput and system resource usage while ensuring maximum disk efficiency through immediate per-sample processing.

## Quick Start

### Default Configuration (Recommended)
```bash
# 24 threads TOTAL distributed across all species (minimum 1 per species)
python3 scripts/rna/batch_download_species.py
```

### Custom Configuration
```bash
# 48 threads total (more throughput on high-end systems)
python3 scripts/rna/batch_download_species.py --total-threads 48

# 16 threads total (conservative for limited resources)
python3 scripts/rna/batch_download_species.py --total-threads 16
```

## Configuration Parameters

### `--total-threads`
- **Default**: 24
- **Description**: Total number of threads to distribute across all species
- **Range**: 8-64 (recommended: 24 for standard systems)
- **Impact**: Threads distributed evenly (minimum 1 per species). More threads = faster overall but higher resource usage
- **Distribution**: For 20 species with 24 threads: 4 species get 2 threads, 16 get 1 thread

### `--max-species`
- **Default**: None (all species)
- **Description**: Maximum number of species to process
- **Use case**: Testing or processing subset of species

### `--quant-threads`
- **Default**: 10
- **Description**: Separate thread pool for quantification operations
- **Range**: 4-20 (recommended: 10)
- **Impact**: Quantification runs in parallel with downloads

## Thread Allocation Strategy

The system uses **total thread allocation** where threads are distributed across all species:

1. **Initial Distribution**: Threads distributed evenly (minimum 1 per species)
2. **Dynamic Redistribution**: As species complete, threads redistribute to remaining species
3. **Concentration**: Threads concentrate on fewer species as workflow progresses
4. **Immediate Processing**: Each sample is downloaded, quantified, and FASTQs deleted immediately

### Example Distribution

For 20 species with 24 total threads:
- 4 species get 2 threads each (8 threads)
- 16 species get 1 thread each (16 threads)
- Total: 24 threads

As species complete, remaining species receive more threads automatically.

## Thread Count Recommendations

### System Resource Considerations

| CPU Cores | Recommended Total Threads | Distribution (20 species) |
|-----------|---------------------------|---------------------------|
| 4-8 cores | 16-24 total | 1 thread per species |
| 16 cores | 24-32 total | 1-2 threads per species |
| 32 cores | 48-64 total | 2-3 threads per species |
| 64+ cores | 64-96 total | 3-4 threads per species |

## Examples

### Standard System (Recommended)
```bash
# 24 threads total distributed across all species
python3 scripts/rna/batch_download_species.py --total-threads 24
```

### High-End System
```bash
# 48 threads total for faster processing
python3 scripts/rna/batch_download_species.py --total-threads 48
```

### Limited Resources
```bash
# 16 threads total for conservative resource usage
python3 scripts/rna/batch_download_species.py --total-threads 16
```

### Testing (Limited Species)
```bash
# Process only first 5 species with 24 threads total
python3 scripts/rna/batch_download_species.py --total-threads 24 --max-species 5
```

## Immediate Processing

The system processes samples immediately (not in batches):

1. **Download**: Sample FASTQ files are downloaded
2. **Immediate Quantification**: Sample is quantified as soon as download completes
3. **Immediate Deletion**: FASTQ files are deleted immediately after quantification
4. **Next Sample**: Process moves to next sample without waiting for batch completion

This ensures maximum disk efficiency - only one sample's FASTQs exist at any time.

## Monitoring Download Progress

### Check Active Downloads
```bash
# View active processes
ps aux | grep "amalgkit getfastq"

# Count active threads
ps aux | grep "amalgkit getfastq" | grep -c "amalgkit"
```

### Analyze Sample Status
```bash
# Comprehensive status analysis
python3 scripts/rna/assess_progress.py
```

### Check Disk Space
```bash
# Overall usage
df -h /
du -sh output/amalgkit/
```

## Performance Tuning

### Signs of Optimal Configuration
- ✅ Network bandwidth is well-utilized (check with `nload` or `iftop`)
- ✅ CPU usage is moderate (not pegged at 100%)
- ✅ Downloads complete without excessive retries
- ✅ System remains responsive
- ✅ Only one sample's FASTQs exist at a time (immediate processing working)

### Signs of Too Many Threads
- ⚠️ High latency/packet loss
- ⚠️ Frequent download failures
- ⚠️ System becomes unresponsive
- ⚠️ Network saturation warnings

### Signs of Too Few Threads
- ⚠️ Downloads are slow
- ⚠️ Network bandwidth underutilized
- ⚠️ CPU idle time high

## Troubleshooting

### Downloads Stuck or Slow
1. Check network connectivity: `ping -c 4 8.8.8.8`
2. Reduce total threads: `--total-threads 16`
3. Check disk space: `df -h /`
4. Verify immediate processing is working (only one sample's FASTQs should exist)

### Out of Disk Space
1. Clean up partial downloads: `python3 scripts/rna/cleanup_partial_downloads.py --execute`
2. Check for large temp files: `du -sh output/amalgkit/*/fastq`
3. Verify immediate processing: FASTQs should be deleted after quantification

### High CPU Usage
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

- **[CONFIGURATION.md](CONFIGURATION.md)**: Complete configuration guide
- **[ORCHESTRATION/BATCH_DOWNLOAD.md](ORCHESTRATION/BATCH_DOWNLOAD.md)**: Orchestration details
- **[WORKFLOW.md](WORKFLOW.md)**: Workflow planning and execution

## Related Scripts

- `scripts/rna/batch_download_species.py` - Main batch download script
- `scripts/rna/cleanup_partial_downloads.py` - Clean up partial/failed downloads
- `scripts/rna/assess_progress.py` - Analyze sample status and progress
