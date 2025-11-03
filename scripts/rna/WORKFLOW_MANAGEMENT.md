# Workflow Management Guide

Complete guide for managing amalgkit RNA-seq workflows.

## Quick Commands

### Check Status
```bash
# Comprehensive status with recommendations
python3 scripts/rna/comprehensive_status.py

# Real-time monitoring
python3 scripts/rna/monitor_comprehensive.py

# Simple progress
bash scripts/rna/monitor_amalgkit_progress.sh
```

### Restart Workflows
```bash
# Restart all workflows
python3 scripts/rna/restart_all_workflows.py

# Check status and restart only if needed
python3 scripts/rna/check_and_restart_workflows.py
```

### Monitor Progress
```bash
# Watch mode (updates every 60 seconds)
watch -n 60 python3 scripts/rna/monitor_comprehensive.py

# Check individual logs
tail -f output/workflow_*_restarted_*.log

# Check running processes
ps aux | grep workflow_ena | grep -v grep
```

## Workflow Status

### Current Configuration
- **Workflow**: `workflow_ena_integrated.py` (ENA-based, robust)
- **Batch Size**: 12 samples per batch
- **Threads**: 12 parallel downloads/quantifications
- **Auto-skip**: Already quantified samples are automatically skipped
- **Cleanup**: FASTQ files deleted immediately after quantification

### Expected Progress Rates
- **Download**: ~6 minutes per 3 samples (varies by sample size)
- **Quantification**: ~36 seconds per sample (single-end, 12 threads)
- **Batch time**: ~15-20 minutes per batch of 12 samples
- **Peak disk**: ~18 GB per batch (1.5 GB × 12 samples)

## Species Status

### C. floridanus (cfloridanus)
- **Total**: 307 samples
- **Progress**: ~250/307 (81%+)
- **Remaining**: ~57 samples

### P. barbatus (pbarbatus)
- **Total**: 83 samples
- **Progress**: ~58/83 (70%+)
- **Remaining**: ~25 samples

### M. pharaonis (mpharaonis)
- **Total**: 100 samples
- **Progress**: ~65/100 (65%+)
- **Remaining**: ~35 samples

### S. invicta (sinvicta)
- **Total**: 354 samples
- **Progress**: ~159/354 (45%+)
- **Remaining**: ~195 samples

### Overall
- **Total**: 844 samples
- **Quantified**: ~532/844 (63%+)
- **Remaining**: ~312 samples

## Troubleshooting

### Workflows Not Running
1. Check status: `python3 scripts/rna/comprehensive_status.py`
2. Check logs: `ls -lht output/workflow_*.log | head -5`
3. Restart if needed: `python3 scripts/rna/restart_all_workflows.py`

### Download Failures
- Normal: Some samples fail due to network/timeout issues
- Auto-retry: Workflows retry failed downloads automatically
- Manual retry: Restart workflow to retry failed samples

### Disk Space Issues
- Batched processing keeps disk usage low (~18 GB per batch)
- FASTQ files automatically deleted after quantification
- Quantification files are small (~2 MB per sample)

### Log Files
- Located in: `output/workflow_{species}_*.log`
- Old logs can be deleted after workflow completion
- Keep only the most recent log per species

## Best Practices

1. **Regular Monitoring**: Check status daily with `comprehensive_status.py`
2. **Automatic Restart**: Use `check_and_restart_workflows.py` to auto-restart inactive workflows
3. **Log Management**: Clean up old logs periodically to save space
4. **Progress Tracking**: Use `monitor_comprehensive.py` for real-time updates

## File Locations

- **Workflow Scripts**: `scripts/rna/workflow_ena_integrated.py`
- **Configuration**: `config/amalgkit/amalgkit_*.yaml`
- **Output**: `output/amalgkit/{species}/`
- **Logs**: `output/workflow_*.log`
- **Quantification**: `output/amalgkit/{species}/quant/`
- **FASTQ**: `output/amalgkit/{species}/fastq/` (temporary, auto-deleted)

