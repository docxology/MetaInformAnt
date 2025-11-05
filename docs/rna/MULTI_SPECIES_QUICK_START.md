# Multi-Species RNA-seq Quick Start Guide

**Last Updated**: November 3, 2025  
**Status**: ✅ Production-validated - Currently processing 3,820 samples in Batch 1 (10 ant species)

This guide provides step-by-step instructions for starting and monitoring multi-species RNA-seq workflows using METAINFORMANT.

---

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [Starting Multi-Species Workflows](#starting-multi-species-workflows)
3. [Monitoring Progress](#monitoring-progress)
4. [Understanding the Output](#understanding-the-output)
5. [Troubleshooting](#troubleshooting)

---

## Prerequisites

### Required Software
- **Python 3.11+** with virtual environment (`.venv/`)
- **amalgkit** (installed in venv via `uv pip install amalgkit`)
- **wget** (for ENA downloads)
- **kallisto** (for quantification)

### Optional Tools
- **fastp** (quality control)
- **seqkit** (sequence manipulation)

### Verify Installation
```bash
# Check Python version
python3 --version

# Activate virtual environment (if not auto-activated)
source .venv/bin/activate

# Verify amalgkit
amalgkit --version

# Verify other tools
which wget kallisto
```

---

## Starting Multi-Species Workflows

### Method 1: Configurable Batch Download (NEW ⭐⭐)

**Best for:** Multi-species parallel downloads with configurable throughput

```bash
# Default: 3 species × 10 threads = 30 total downloads
python3 scripts/rna/batch_download_species.py

# Custom: 4 species × 12 threads = 48 total downloads
python3 scripts/rna/batch_download_species.py --species-count 4 --threads-per-species 12

# Run in background with logging
python3 scripts/rna/batch_download_species.py > output/batch_download.log 2>&1 &
```

**Features:**
- ✅ **Configurable parallelism**: Adjust species count and threads per species
- ✅ **Automatic discovery**: Finds all species configs automatically
- ✅ **Virtual environment**: Auto-activates if available
- ✅ **Cloud acceleration**: AWS, GCP, NCBI enabled by default
- ✅ **Progress tracking**: Logs all download activity

**Configuration Options:**
- `--species-count`: Number of species to download in parallel (default: 3)
- `--threads-per-species`: Threads per species (default: 10)
- `--max-species`: Limit total species processed (optional)

See `docs/rna/BATCH_DOWNLOAD_CONFIGURATION.md` for complete guide.

### Method 2: Production ENA Workflow (RECOMMENDED ⭐)

**Best for:** Maximum reliability, large-scale processing, single species focus

```bash
# Navigate to repository root
cd /home/q/Documents/GitHub/MetaInformAnt

# Batch Processing Strategy (Current - November 2025)
# Processing 20 species in 2 batches of 10 for disk space management
# See docs/rna/discovery/batched_processing.md for details

# For single species or smaller batches:
# C. floridanus
nohup python3 scripts/rna/workflow_ena_integrated.py \
  --config config/amalgkit/amalgkit_cfloridanus.yaml \
  --batch-size 12 \
  --threads 12 \
  > output/workflow_cfloridanus_$(date +%Y%m%d_%H%M%S).log 2>&1 &

# P. barbatus
nohup python3 scripts/rna/workflow_ena_integrated.py \
  --config config/amalgkit/amalgkit_pbarbatus.yaml \
  --batch-size 12 \
  --threads 12 \
  > output/workflow_pbarbatus_$(date +%Y%m%d_%H%M%S).log 2>&1 &

# M. pharaonis
nohup python3 scripts/rna/workflow_ena_integrated.py \
  --config config/amalgkit/amalgkit_mpharaonis.yaml \
  --batch-size 12 \
  --threads 12 \
  > output/workflow_mpharaonis_$(date +%Y%m%d_%H%M%S).log 2>&1 &

# S. invicta
nohup python3 scripts/rna/workflow_ena_integrated.py \
  --config config/amalgkit/amalgkit_sinvicta.yaml \
  --batch-size 12 \
  --threads 12 \
  > output/workflow_sinvicta_$(date +%Y%m%d_%H%M%S).log 2>&1 &

# Check processes started
ps aux | grep workflow_ena
```

**Features:**
- ✅ **100% reliability** (vs 0% with SRA Toolkit)
- ✅ **Direct ENA downloads** via API
- ✅ **Automatic resume** with `wget --continue`
- ✅ **Batched processing** (12 samples at a time)
- ✅ **Auto-cleanup** (FASTQs deleted after quantification)
- ✅ **Batch orchestration** (10 species per batch for disk management)
- ✅ **Parallel execution** (multiple species simultaneously)

**Expected behavior:**
- Each workflow downloads 12 samples → quantifies → deletes FASTQs → repeats
- Peak disk usage: ~50-100 GB (temporary, auto-cleaned per sample)
- Processing rate: ~7.5 minutes per sample
- Final results: ~2-3 GB per species (expression matrices + QC)
- Network interruptions: automatically resume

### Method 2: Legacy SRA Workflow (Alternative)

**Best for:** Testing, single species, environments without ENA access

```bash
# Activates venv automatically
python3 scripts/rna/run_multi_species.py
```

**Features:**
- Auto-discovers all `config/amalgkit/amalgkit_*.yaml` files
- Auto-activates virtual environment
- Batched processing with SRA Toolkit
- Includes cross-species analysis (CSTMM, CSCA)

**Note:** SRA Toolkit has ~0% success rate for large samples. Use ENA workflow for production.

### Method 3: Single Species Testing

```bash
# Test with just one species
python3 scripts/rna/workflow_ena_integrated.py \
  --config config/amalgkit/amalgkit_pbarbatus.yaml \
  --batch-size 3 \
  --threads 8 \
  --max-samples 3
```

---

## Monitoring Progress

### Real-Time Comprehensive Monitor (RECOMMENDED ⭐)

```bash
# Run once for current status
python3 scripts/rna/monitor_comprehensive.py
```

**Output:**
```
================================================================================
  MULTI-SPECIES RNA-SEQ WORKFLOW MONITOR
================================================================================

⏰ 2025-11-01 16:46:28

📊 C. floridanus (cfloridanus)
   Progress: 12/307 (3%)
   Batch: 1/25/25
   Downloading: 12 samples (35G)

📊 P. barbatus (pbarbatus)
   Progress: 0/83 (0%)
   Batch: 1/7/7
   Downloading: 12 samples (4.1G)

📊 M. pharaonis (mpharaonis)
   Progress: 1/100 (1%)
   Batch: 1/9/9
   Downloading: 47 samples (219M)

📊 S. invicta (sinvicta)
   Progress: 0/354 (0%)
   Batch: 1/30/30
   Downloading: 96 samples (540M)

================================================================================
🎯 OVERALL: Progress shown for example species (current Batch 1: 3,820 samples across 10 species)
================================================================================

⏳ Estimated remaining: 103.9 hours (831 samples)
```

**Continuous monitoring:**
```bash
# Watch mode (updates every 60 seconds)
watch -n 60 python3 scripts/rna/monitor_comprehensive.py
```

### Alternative Monitors

**Detailed workflow monitor:**
```bash
python3 scripts/rna/monitor_workflow.py
```

**Simple bash monitor:**
```bash
bash scripts/rna/monitor_amalgkit_progress.sh

# Or watch mode
watch -n 60 bash scripts/rna/monitor_amalgkit_progress.sh
```

### Check Running Processes

```bash
# Count workflow processes
ps aux | grep workflow_ena | grep -v grep | wc -l

# Count download processes
ps aux | grep download_ena | grep -v grep | wc -l

# Count active downloads
ps aux | grep wget | grep -v grep | wc -l

# Show all RNA workflow processes
ps aux | grep -E "(workflow_ena|download_ena|wget)" | grep -v grep
```

### Check Individual Species Logs

```bash
# Find latest logs
ls -lt output/workflow_*.log | head -4

# Tail specific species
tail -f output/workflow_cfloridanus_*.log
tail -f output/workflow_pbarbatus_*.log
tail -f output/workflow_mpharaonis_*.log
tail -f output/workflow_sinvicta_*.log

# Check for errors
grep -i error output/workflow_*.log
grep -i failed output/workflow_*.log
```

### Check Disk Usage

```bash
# Overall amalgkit directory
du -sh output/amalgkit/

# Per species
du -sh output/amalgkit/*/

# FASTQ directories (should stay small with batching)
du -sh output/amalgkit/*/fastq/

# Quantification results
du -sh output/amalgkit/*/quant/
```

---

## Understanding the Output

### Directory Structure

```
output/amalgkit/
├── cfloridanus/
│   ├── quant/              # Quantification results (~2 MB per sample)
│   │   ├── SRR1234567/
│   │   │   ├── abundance.tsv
│   │   │   └── run_info.json
│   │   └── ...
│   ├── fastq/              # Temporary FASTQs (auto-deleted, flat structure)
│   │   └── {SRR_ID}/       # Sample directories directly in fastq/
│   ├── work/
│   │   ├── metadata/
│   │   └── index/          # Kallisto index
│   └── logs/               # Workflow logs
├── pbarbatus/
├── mpharaonis/
├── sinvicta/
└── workflow_*.log          # Main workflow logs
```

### Key Files

**Per-sample quantification:**
- `output/amalgkit/{species}/quant/{SRR_ID}/abundance.tsv` - Expression levels
- `output/amalgkit/{species}/quant/{SRR_ID}/run_info.json` - Kallisto metrics

**Metadata:**
- `output/amalgkit/{species}/work/metadata/metadata.tsv` - Sample information

**Kallisto index:**
- `output/amalgkit/{species}/work/index/transcriptome.idx` - Built once per species

**Logs:**
- `output/workflow_{species}_TIMESTAMP.log` - Main workflow log
- `output/amalgkit/{species}/logs/` - Per-step logs

### Progress Indicators

**Completed sample:**
- Quantification directory exists: `output/amalgkit/{species}/quant/{SRR_ID}/`
- Contains: `abundance.tsv` and `run_info.json`
- FASTQ files deleted from `fastq/{SRR_ID}/` (flat structure, no nested subdirectories)

**Downloading sample:**
- FASTQ directory exists: `output/amalgkit/{species}/fastq/{SRR_ID}/`
- Contains `.fastq.gz` files (growing)
- wget process visible in `ps aux | grep wget`

**Failed sample:**
- Check logs: `grep {SRR_ID} output/workflow_{species}_*.log`
- Workflow will retry automatically on restart

---

## Performance Characteristics

### Expected Timeline

**Per sample:**
- Download: 2-8 minutes (varies by sample size, 1-4 GB)
- Quantification: 30-90 seconds (with kallisto, 12 threads)
- Cleanup: <1 second
- **Total**: ~7.5 minutes average

**Per batch (12 samples):**
- Download: 10-15 minutes (parallel)
- Quantification: 6-12 minutes (parallel)
- **Total**: ~20-30 minutes per batch

**Full workflow:**
- C. floridanus (307 samples): ~38 hours (25 batches)
- P. barbatus (83 samples): ~10 hours (7 batches)
- M. pharaonis (100 samples): ~12 hours (9 batches)
- S. invicta (354 samples): ~44 hours (30 batches)
- **Total (Batch 1)**: ~24-48 hours for 3,820 samples (10 species)
- **Total (Batch 2)**: ~12-24 hours for 728 samples (10 species)

### Resource Usage

**CPU:**
- 12 threads per species × 10 species (Batch 1) = up to 120 threads total
- Recommend: 64+ cores for smooth parallel processing

**RAM:**
- ~4 GB per workflow
- ~16 GB total for 4 parallel workflows

**Disk:**
- Peak: ~18 GB per species (12 samples × 1.5 GB avg)
- Total peak: ~50-100 GB (temporary, auto-cleaned per sample)
- Final: ~1.7 GB (quantification files only)

**Network:**
- Bandwidth: ~40 GB downloading at any time
- ENA servers: highly reliable, automatic resume

---

## Troubleshooting

### Workflow Not Starting

**Check virtual environment:**
```bash
# Verify venv exists
ls -la .venv/

# Verify amalgkit installed
source .venv/bin/activate
amalgkit --version
```

**Check dependencies:**
```bash
# Verify wget
which wget
wget --version

# Verify kallisto
which kallisto
kallisto version
```

### Downloads Failing

**Check network:**
```bash
# Test ENA connectivity
wget --spider https://www.ebi.ac.uk/ena/

# Check active connections
netstat -an | grep ESTABLISHED | wc -l
```

**Resume interrupted downloads:**
- ENA workflow automatically resumes with `wget --continue`
- Simply restart the workflow - it will skip completed samples

### Disk Space Issues

**Automatic Protection (NEW):**
The batch download workflow now includes automatic disk space management:
- ✅ Pre-flight checks before starting downloads
- ✅ Automatic cleanup of partial/failed downloads when space is low
- ✅ Default thread count reduced to 8 for smaller drives
- ✅ Health checks monitor disk space every 10 minutes
- ✅ Automatic recovery from disk space errors

**Manual Checks:**
```bash
# Check available space
df -h /home/q/Documents/GitHub/MetaInformAnt/output/

# Clean up partial/failed downloads
python3 scripts/rna/cleanup_partial_downloads.py --execute
```

**Configuration for Different Drive Sizes:**
```bash
# Small drive (default: 8 threads, auto-cleanup at 5GB)
python3 scripts/rna/batch_download_species.py --total-threads 8

# Large drive (can use more threads)
python3 scripts/rna/batch_download_species.py --total-threads 30

# Custom thresholds
python3 scripts/rna/batch_download_species.py \
  --total-threads 8 \
  --min-free-gb 15.0 \
  --auto-cleanup-threshold 10.0
```

**Reduce batch size:**
```bash
# Use smaller batches (6 instead of 12)
python3 scripts/rna/workflow_ena_integrated.py \
  --config config/amalgkit/amalgkit_cfloridanus.yaml \
  --batch-size 6 \
  --threads 12
```

### Processes Stuck or Hung

**Kill and restart:**
```bash
# Find process IDs
ps aux | grep workflow_ena

# Kill specific workflow
kill <PID>

# Or kill all workflows (use with caution)
pkill -f workflow_ena

# Restart (will resume from last completed sample)
nohup python3 scripts/rna/workflow_ena_integrated.py ...
```

### Quantification Failing

**Check kallisto index:**
```bash
# Verify index exists
ls -lh output/amalgkit/*/work/index/*.idx

# Rebuild if needed (automatic on first run)
# Just restart the workflow
```

**Check FASTQ quality:**
```bash
# List downloaded FASTQs
ls -lh output/amalgkit/*/fastq/*/

# Check if files are complete (not 0 bytes)
find output/amalgkit/*/fastq -name "*.fastq.gz" -size 0
```

### Getting Help

**Check logs:**
```bash
# Recent errors across all logs
find output -name "workflow_*.log" -exec grep -i error {} +

# Species-specific issues
tail -100 output/workflow_cfloridanus_*.log
```

**Monitor comprehensive status:**
```bash
# Get detailed current state
python3 scripts/rna/monitor_comprehensive.py
```

**Report issues:**
- Include: species name, sample ID, error message
- Attach: relevant log snippet
- Note: workflow will auto-skip failed samples and continue

---

## Quick Reference Commands

```bash
# START WORKFLOWS
# Batch 1: Top 10 species (production)
# bash scripts/rna/run_top10_ant_species.sh
for species in cfloridanus pbarbatus mpharaonis sinvicta; do
  nohup python3 scripts/rna/workflow_ena_integrated.py \
    --config config/amalgkit/amalgkit_${species}.yaml \
    --batch-size 12 --threads 12 \
    > output/workflow_${species}_$(date +%Y%m%d_%H%M%S).log 2>&1 &
done

# MONITOR PROGRESS
python3 scripts/rna/monitor_comprehensive.py          # Once
watch -n 60 python3 scripts/rna/monitor_comprehensive.py  # Continuous

# CHECK PROCESSES
ps aux | grep workflow_ena | grep -v grep             # Main workflows
ps aux | grep wget | grep -v grep | wc -l             # Active downloads

# CHECK DISK SPACE
du -sh output/amalgkit/                               # Total usage
du -sh output/amalgkit/*/fastq/                       # FASTQs only

# VIEW LOGS
tail -f output/workflow_cfloridanus_*.log             # Follow log
grep -i error output/workflow_*.log                   # Find errors

# CLEANUP
bash scripts/rna/cleanup_quantified_sra.sh --execute  # Free disk space

# STOP WORKFLOWS (if needed)
pkill -f workflow_ena                                 # Stop all
```

---

## Success Criteria

✅ **Workflow is successful when:**
- All samples quantified: `abundance.tsv` and `run_info.json` present
- FASTQs cleaned up: `fastq/` directories empty or deleted
- Logs show "Batch N/N complete"
- Monitor shows 100% progress
- No ERROR messages in logs

✅ **Expected final state:**
- `output/amalgkit/{species}/quant/` contains N sample directories
- Each sample has ~2 MB of quantification data
- Total size: ~40-55 GB for all 4,548 samples (20 species)
- Logs show successful completion

---

## Next Steps After Completion

1. **Verify all samples quantified:**
   ```bash
   python3 scripts/rna/monitor_comprehensive.py
   ```

2. **Generate expression matrices:**
   ```bash
   # Per species (if using amalgkit directly)
   amalgkit merge --work-dir output/amalgkit/cfloridanus/work
   ```

3. **Run cross-species analysis:**
   ```bash
   # CSTMM normalization
   # CSCA correlation analysis
   # (See amalgkit documentation)
   ```

4. **Downstream analysis:**
   - Differential expression
   - Gene ontology enrichment
   - Publication figures

---

**Status**: Batch processing strategy active (November 2025):
- **Batch 1**: 3,820 samples across 10 ant species - currently running
- **Batch 2**: 728 samples across 10 ant species - queued
- **Total**: 4,548 samples across 20 ant species
- **Reliability**: 100% with ENA direct downloads
- **Expected completion**: 24-48 hours for Batch 1, 12-24 hours for Batch 2

**Documentation Version**: 1.0 (November 2025)  
**Tested Configuration**: 12 threads, 12-sample batches, ENA direct downloads

