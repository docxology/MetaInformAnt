# Workflow Management Guide

Complete guide for managing amalgkit RNA-seq workflows.

## Quick Commands

### Check Status
```bash
# Single species status (recommended)
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --status

# Detailed status with recommendations
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --status --detailed

# Check all species (using monitoring functions)
python3 -c "from metainformant.rna.monitoring import assess_all_species_progress; from pathlib import Path; print(assess_all_species_progress(Path('config/amalgkit')))"
```

### Monitor Progress
```bash
# Watch mode (updates every 60 seconds)
watch -n 60 'python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --status'

# Check individual logs
tail -f output/amalgkit/pogonomyrmex_barbatus/logs/*.log

# Check running processes
ps aux | grep amalgkit | grep -v grep
ps aux | grep fasterq-dump | grep -v grep
```

### Resume Workflows
```bash
# Workflows automatically skip completed steps
# Just re-run the workflow command to resume
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml

# Or run specific steps
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps getfastq,quant
```

## Workflow Status

### Current Configuration
- **Workflow**: `run_workflow.py` (recommended orchestrator)
- **Parallel Downloads**: Configured via `num_download_workers` in YAML
- **Auto-skip**: Already completed steps are automatically skipped
- **Cleanup**: FASTQ files deleted immediately after quantification (when using unified processing)

### Expected Progress Rates
- **Download**: Varies by sample size and network (typically 5-15 minutes per sample)
- **Quantification**: ~30-60 seconds per sample (depends on transcriptome size and threads)
- **Disk usage**: Managed automatically with immediate FASTQ deletion after quantification

## Monitoring Functions

The `src/metainformant/rna/monitoring.py` module provides programmatic access to workflow status:

```python
from metainformant.rna.monitoring import (
    analyze_species_status,      # Comprehensive status for one species
    assess_all_species_progress,  # Status for all species in config directory
    check_workflow_progress,      # Step-by-step progress
    count_quantified_samples,     # Sample counts
    find_unquantified_samples,    # List of samples needing quantification
    check_active_downloads,       # Currently downloading samples
    get_sample_status,            # Status for a specific sample
)
```

## Troubleshooting

### Workflows Not Running
1. Check status: `python3 scripts/rna/run_workflow.py --config <config> --status`
2. Check logs: `ls -lht output/amalgkit/<species>/logs/ | head -5`
3. Check environment: `python3 scripts/rna/check_environment.py`
4. Resume workflow: Re-run the workflow command (auto-skips completed steps)

### Download Failures
- Normal: Some samples fail due to network/timeout issues
- Auto-retry: Workflows retry failed downloads automatically
- Manual retry: Re-run workflow to retry failed samples
- Cleanup: Use `--cleanup-partial` to remove stuck partial downloads

### Disk Space Issues
- Unified processing keeps disk usage low (FASTQ deleted immediately after quant)
- Use `--cleanup-unquantified` to quantify downloaded samples and free space
- Quantification files are small (~2 MB per sample)

### Log Files
- Located in: `output/amalgkit/<species>/logs/`
- Old logs can be deleted after workflow completion
- Keep only the most recent logs per species

## Best Practices

1. **Regular Monitoring**: Check status frequently with `--status`
2. **Resume After Interruption**: Workflows automatically skip completed steps
3. **Log Management**: Clean up old logs periodically to save space
4. **Progress Tracking**: Use `--status --detailed` for comprehensive information
5. **Environment Checks**: Run `check_environment.py` before starting workflows

## File Locations

- **Workflow Scripts**: `scripts/rna/run_workflow.py`
- **Configuration**: `config/amalgkit/amalgkit_*.yaml`
- **Output**: `output/amalgkit/<species>/`
- **Logs**: `output/amalgkit/<species>/logs/`
- **Quantification**: `output/amalgkit/<species>/quant/`
- **FASTQ**: `output/amalgkit/<species>/fastq/` (temporary, auto-deleted after quantification)

## Related Documentation

- **Getting Started**: `docs/rna/GETTING_STARTED.md` - Complete setup and workflow guide
- **Examples**: `docs/rna/EXAMPLES.md` - Real-world workflow examples
- **Orchestration**: `docs/rna/ORCHESTRATION.md` - Multi-species and advanced workflows
- **Monitoring API**: `src/metainformant/rna/monitoring.py` - Programmatic monitoring functions
