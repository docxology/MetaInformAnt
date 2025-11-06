# RNA-seq Immediate Processing Strategy

**Last Updated**: November 5, 2025  
**Status**: Active with immediate per-sample processing

## Overview

Processing 20 ant species with immediate per-sample processing (download → immediately quantify → immediately delete FASTQs) using 24 total threads distributed across all species.

## Rationale

✅ **Maximum Disk Efficiency**: Only one sample's FASTQs exist at any time  
✅ **Immediate Processing**: Download → immediately quantify → immediately delete  
✅ **Total Thread Allocation**: 24 threads distributed across all species (not per species)  
✅ **Dynamic Redistribution**: Threads redistribute as species complete  
✅ **Scalability**: Can process all species simultaneously without disk space issues  

## Processing Strategy

### Immediate Per-Sample Processing
- **All Species**: 20 species processed simultaneously
- **Total Samples**: 1,578 samples across all species
- **Thread Allocation**: 24 threads TOTAL distributed evenly (minimum 1 per species)
- **Processing Mode**: Each sample: download → immediately quantify → immediately delete FASTQs
- **Script**: `scripts/rna/batch_download_species.py --total-threads 24`
- **Status**: Active

## Disk Management

**Per Sample** (immediate processing):
1. Download FASTQ (~2-10 GB)
2. **Immediately** quantify with Kallisto
3. ✅ **Immediately delete FASTQ** after quantification
4. Keep results only (~10-50 MB)
5. Move to next sample

**Peak Usage**: ~2-10 GB (only one sample's FASTQs exist at a time)  
**Final Usage**: ~40-55 GB (all results across all species)

## Timeline

| Phase | Duration | Batch 1 | Batch 2 |
|-------|----------|---------|---------|
| Metadata | 1-2 hours | ✅ | ⏳ |
| Downloads | 12-36 hours | 🔄 | ⏳ |
| Quantification | 6-18 hours | ⏳ | ⏳ |
| Merge + QC | 2-4 hours | ⏳ | ⏳ |
| **Total** | **24-48 hours** | | **12-24 hours** |

**Complete Project**: 36-72 hours for all 20 species

## Monitoring

```bash
# Check status
python3 scripts/rna/get_current_status.py

# View logs
tail -f output/top10_*.log

# Check disk
df -h output/

# Active processes
ps aux | grep run_amalgkit | wc -l
```

## Launch Commands

```bash
# Immediate processing with 24 threads TOTAL distributed across all species (recommended)
python3 scripts/rna/batch_download_species.py --total-threads 24

# Or run full workflow sequentially (one species at a time)
export AK_THREADS=24
python3 scripts/rna/run_multi_species.py
```

---

See `discovery/README.md` for species details and `../amalgkit/README.md` for workflow documentation.




