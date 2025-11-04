# RNA-seq Batched Processing Strategy

**Last Updated**: November 3, 2025  
**Status**: Batch 1 Running

## Overview

Processing 20 ant species in 2 batches of 10 to manage disk space efficiently.

## Rationale

✅ **Disk Space Management**: FASTQs auto-deleted after quantification  
✅ **Resource Optimization**: System not overwhelmed  
✅ **Progress Monitoring**: Assess batch 1 before starting batch 2  
✅ **Troubleshooting**: Easier to identify and fix issues  
✅ **Flexibility**: Adjust strategy based on results  

## Batch Structure

### Batch 1 (Running)
- **Species**: 10 (top by sample count)
- **Samples**: 3,820 total
- **Status**: Active (started Nov 3, 2025 16:18 PST)
- **ETA**: 24-48 hours
- **Logs**: `output/top10_*.log`

### Batch 2 (Queued)
- **Species**: 10 (remaining)
- **Samples**: 728 total  
- **Launch**: After Batch 1 completes
- **Script**: `scripts/rna/run_batch2_ant_species.sh`
- **ETA**: 12-24 hours

## Disk Management

**Per Sample**:
1. Download FASTQ (~2-10 GB)
2. Quantify with Kallisto
3. ✅ **Delete FASTQ** immediately
4. Keep results only (~10-50 MB)

**Peak Usage**: ~50-100 GB (temporary)  
**Final Usage**: ~40-55 GB (all results)

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
# Batch 1 (already running)
bash scripts/rna/run_top10_ant_species.sh

# Batch 2 (run after Batch 1)
bash scripts/rna/run_batch2_ant_species.sh
```

---

See `discovery/README.md` for species details and `../amalgkit/README.md` for workflow documentation.



