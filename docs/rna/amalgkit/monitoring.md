# Pipeline Monitoring

How to monitor the status of active amalgkit pipeline runs.

## Live Log Monitoring

The main log for the all-species run is:

```bash
# Follow in real-time
tail -f output/amalgkit/run_all_species_incremental.log

# Filter for key status lines
grep -E "Done|SzSkip|Error|species" output/amalgkit/run_all_species_incremental.log | tail -50
```

Key log tokens:
- `✓ Done` — sample quantified successfully
- `✗ SzSkip` — sample skipped (exceeds `max_bp` size limit)
- `✗ Error` — quantification failed
- `[<species>]` — pipeline position marker

## Active Process Monitor

```bash
# Count active wget download workers
ps aux | grep wget | grep -v grep | wc -l

# List active kallisto processes
ps aux | grep 'kallisto quant' | grep -v grep

# See all amalgkit-related processes
ps aux | grep -E 'wget|kallisto|amalgkit' | grep -v grep
```

## Report Completed Samples (Disk-Based)

```bash
# Count quantified samples per species from disk
python3 scripts/rna/report_completed.py

# Or the shell version
bash scripts/rna/report_completed.sh
```

These check for `abundance.tsv` files under `output/amalgkit/*/work/quant/`.

## Single-Species Status

For species managed via `run_workflow.py`:

```bash
# High-level status
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pbarbatus.yaml --status

# With detail per step
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pbarbatus.yaml --status --detailed
```

## Quick Health Checks

```bash
# How many samples quantified for amellifera?
find output/amalgkit/amellifera/work/quant -name "abundance.tsv" | wc -l

# Any samples stuck downloading (large .fastq.gz partial files)?
find output/amalgkit -name "*.fastq.gz" -size +1G 2>/dev/null | head

# Check disk space
df -h output/

# Most recently modified quant directories
ls -lt output/amalgkit/amellifera/work/quant/ | head -20
```

## Log Quick-Counts (Without Scanning Disk)

```bash
LOG=output/amalgkit/run_all_species_incremental.log
grep  "✓ Done"   "$LOG" | wc -l   # completed
grep  "✗ SzSkip" "$LOG" | wc -l   # size-skipped
grep  "✗ Error"  "$LOG" | wc -l   # failed
```

## What Was Changed in v0.2.7

In v0.2.7 (2026-03-01):
- **Heartbeat JSON files** (`*.heartbeat.json`) were deleted and added to `.gitignore`. These were transient download-progress trackers that were never read back; the `abundance.tsv` file is the authoritative proof a sample completed.
- **`.safely_removed` markers** were deleted and added to `.gitignore`. The non-empty `abundance.tsv` makes them redundant.

See [CHANGELOG.md](../../../../CHANGELOG.md) for details.
