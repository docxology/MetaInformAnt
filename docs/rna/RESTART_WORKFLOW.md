# Restarting End-to-End Amalgkit Workflow

**Complete guide for restarting the full end-to-end amalgkit workflow for all species with configurable threads.**

## Quick Start

### Recommended: Parallel Execution (All Species Simultaneously)

```bash
# 1. Ensure /tmp space is available (if needed)
bash scripts/rna/fix_tmp_space.sh

# 2. Set threads per species (default: 1, or from AK_THREADS)
export AK_THREADS=24  # Or any number you prefer

# 3. Start all species workflows in parallel
python3 scripts/rna/run_all_species_parallel.py --threads-per-species 24

# Or use environment variable:
export AK_THREADS=24
python3 scripts/rna/run_all_species_parallel.py
```

**What this does:**
- ✅ Runs all 24 species workflows **simultaneously** (one process per species)
- ✅ Each species uses 24 threads (configurable via `--threads-per-species` or `AK_THREADS`)
- ✅ Total concurrent threads: 24 species × 24 threads = 576 threads
- ✅ Automatically handles virtual environment setup
- ✅ Downloads → Quantifies → Deletes FASTQs → Moves to next batch
- ✅ Progress tracking and dashboard updates automatically

### Alternative: Sequential Execution (One Species at a Time)

```bash
export AK_THREADS=24
python3 scripts/rna/run_multi_species.py
```

**What this does:**
- ✅ Runs species **sequentially** (one completes before next starts)
- ✅ Uses 24 threads per species (from `AK_THREADS`)
- ✅ Total concurrent threads: 24 threads (one species at a time)
- ✅ More conservative resource usage

## Prerequisites

### 1. Virtual Environment Setup

The scripts automatically discover and use virtual environments, but you need to set it up first:

```bash
# Install uv if not already installed
curl -LsSf https://astral.sh/uv/install.sh | sh

# Create virtual environment
uv venv .venv
# Or on ext6 filesystems, it will automatically use /tmp/metainformant_venv if needed

# Install dependencies
source .venv/bin/activate  # Or: source /tmp/metainformant_venv/bin/activate
uv pip install -e .
uv pip install git+https://github.com/kfuku52/amalgkit
```

**Note**: Scripts automatically discover venv location (`.venv` or `/tmp/metainformant_venv`) - no manual activation needed when running scripts.

### 2. Check Disk Space

Ensure `/tmp` has space (workflows may create temporary files):

```bash
# Check /tmp space
df -h /tmp

# If /tmp is full, clean it:
bash scripts/rna/fix_tmp_space.sh
```

### 3. Verify Configuration

All species configs should be in `config/amalgkit/amalgkit_*.yaml`. The scripts automatically discover them.

## Thread Configuration

### Option 1: Total Threads (Recommended for Immediate Processing)

```bash
# 24 threads TOTAL distributed across all species (recommended)
python3 scripts/rna/batch_download_species.py --total-threads 24
```

### Option 2: Environment Variables (For Sequential Workflows)

```bash
# Set threads per species (for sequential execution)
export AK_THREADS=24
python3 scripts/rna/run_multi_species.py
```

### Option 3: Auto-Detection (For Immediate Processing)

```bash
# Total threads distributed automatically (default: 24 total threads)
python3 scripts/rna/batch_download_species.py --total-threads 24
# Threads are distributed evenly across all species (minimum 1 per species)
# Dynamic redistribution as species complete
```

### Option 4: Edit Config Files

Edit `threads:` in each config file:
```yaml
# config/amalgkit/amalgkit_cfloridanus.yaml
threads: 24  # Change from default
```

### Thread Recommendations

| System Resources | Recommended Total Threads | Distribution (20 species) |
|------------------|---------------------------|---------------------------|
| High-end (64+ cores, 128GB+ RAM) | 48-64 total | 2-3 threads per species |
| Mid-range (32 cores, 64GB RAM) | 24-32 total | 1-2 threads per species |
| Standard (16 cores, 32GB RAM) | 16-24 total | 1 thread per species |
| Low-end (8 cores, 16GB RAM) | 8-16 total | 1 thread per species (shared) |

**Note**: Immediate processing (`batch_download_species.py`) uses **total threads distributed across all species** (minimum 1 per species). Threads redistribute dynamically as species complete. Sequential execution (`run_multi_species.py`) uses threads per species (default: 24 per species).

## Monitoring Progress

### Real-Time Dashboard

The workflow automatically generates a progress dashboard:

```bash
# View the dashboard
cat output/amalgkit/progress_dashboard.txt

# Or watch it update in real-time
watch -n 30 cat output/amalgkit/progress_dashboard.txt
```

### Status Commands

```bash
# Brief status summary
python3 scripts/rna/orchestrate_workflows.py --status

# Detailed status with categories
python3 scripts/rna/orchestrate_workflows.py --status --detailed

# Real-time monitoring (updates every 60s)
python3 scripts/rna/orchestrate_workflows.py --monitor

# Custom watch interval
python3 scripts/rna/orchestrate_workflows.py --monitor --watch 30
```

### Check Running Processes

```bash
# See active workflows
ps aux | grep -E "(workflow_ena|run_all_species_parallel|amalgkit)" | grep -v grep

# Count active downloads
ps aux | grep "amalgkit getfastq" | grep -v grep | wc -l
```

## Workflow Behavior

### Automatic Resume

The workflow automatically resumes from where it left off:

- ✅ **Skips already-quantified samples** (checks for `quant/{sample_id}/abundance.tsv`)
- ✅ **Resumes downloads** (checks for existing FASTQ files)
- ✅ **Continues from last batch** (processes remaining samples)

### Cleanup Strategy

The workflow automatically manages disk space with immediate per-sample processing:

- ✅ **Downloads** one sample
- ✅ **Immediately quantifies** the sample
- ✅ **Immediately deletes FASTQs** after quantification
- ✅ **Repeats** with next sample

This ensures maximum disk efficiency: only one sample's FASTQs exist at any time (typically ~1-5 GB peak usage per sample).

### Progress Tracking

Every sample completion (download, quant, delete) is tracked:

- ✅ **State file**: `output/amalgkit/progress_state.json`
- ✅ **Dashboard**: `output/amalgkit/progress_dashboard.txt`
- ✅ **Log file**: `output/amalgkit/progress_tracker.log`

## Troubleshooting

### Workflow Stuck or Not Progressing

```bash
# Check if processes are running
ps aux | grep -E "(workflow_ena|run_all_species)" | grep -v grep

# Check disk space
df -h /media/q/ext6 /tmp

# Check logs
tail -50 output/parallel_all_species_*.log
tail -50 output/workflow_*.log
```

### Out of Disk Space

```bash
# Clean /tmp
bash scripts/rna/fix_tmp_space.sh

# Check for quantified samples with FASTQs still present
python3 scripts/rna/emergency_cleanup.py --dry-run

# Clean them up (if needed)
python3 scripts/rna/emergency_cleanup.py
```

### Restart After Interruption

The workflow automatically resumes, but if you need to manually restart:

```bash
# Stop all running workflows
pkill -f "run_all_species_parallel"
pkill -f "workflow_ena_integrated"

# Wait a few seconds
sleep 5

# Restart
export AK_THREADS=24
python3 scripts/rna/run_all_species_parallel.py --threads-per-species 24
```

## Complete Example

```bash
# 1. Setup (one-time, if not already done)
curl -LsSf https://astral.sh/uv/install.sh | sh
uv venv .venv
source .venv/bin/activate
uv pip install -e .
uv pip install git+https://github.com/kfuku52/amalgkit

# 2. Ensure /tmp has space (if needed)
bash scripts/rna/fix_tmp_space.sh

# 3. Configure threads and batch size (optional - auto-detected if not specified)
export AK_THREADS=24
export AK_BATCH_SIZE=50  # Optional: auto-detected for large drives

# 4. Start workflow (parallel execution)
# Batch size and disk thresholds auto-detected based on drive size
python3 scripts/rna/run_all_species_parallel.py \
  --threads-per-species 24 \
  --batch-size 50 \
  --max-batch-size 100

# 5. Monitor in another terminal
watch -n 30 cat output/amalgkit/progress_dashboard.txt
```

**Note**: For large drives (6TB+), batch sizes are automatically set to 50-100 samples per batch, and temporary files use the external drive (`output/.tmp/`) instead of `/tmp`.

## See Also

- **[RUN_ALL_SPECIES.md](RUN_ALL_SPECIES.md)**: Detailed guide for all execution methods
- **[ORCHESTRATION/README.md](ORCHESTRATION/README.md)**: Orchestrator selection and features
- **[CONFIGURATION.md](CONFIGURATION.md)**: Config file structure and options
- **[WORKFLOW.md](WORKFLOW.md)**: Workflow step details

