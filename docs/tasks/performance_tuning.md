# Performance Tuning Quick Reference

Measure the actual bottleneck before changing concurrency. RNA acquisition,
SRA extraction, compression, validation, and Kallisto quantification use
different resources and must be tuned independently.

## Observe first

```bash
uv run python scripts/rna/check_pipeline_status.py \
  --data-root "$AMALGKIT_DATA_ROOT" --verbose

uv run python projects/hymenoptera_amalgkit/scripts/report_campaign_status.py \
  --data-root "$AMALGKIT_DATA_ROOT" --json

ps -axo pid=,etime=,command= | \
  rg 'run_all_species.py|curl .*fastq|fasterq-dump|kallisto quant'
```

If the process tree shows resumable ENA transfers or `fasterq-dump` without
active Kallisto processes, increasing quantification slots will not improve
throughput and may increase external-volume contention.

## Bounded producer settings

```bash
AMALGKIT_PIPELINE_WORKERS=4 \
AMALGKIT_PIPELINE_THREADS=8 \
AMALGKIT_PIPELINE_QUANT_SLOTS=2 \
AMALGKIT_PIPELINE_FASTQ_THREADS=1 \
AMALGKIT_PIPELINE_FASTQ_SLOTS=1 \
AMALGKIT_PIPELINE_VALIDATION_SLOTS=2 \
uv run python scripts/rna/run_all_species.py \
  --config-dir projects/hymenoptera_amalgkit/config/amalgkit \
  --data-root "$AMALGKIT_DATA_ROOT" \
  --workers 4 --threads 8 --quant-slots 2 \
  --fastq-threads 1 --fastq-slots 1 --validation-slots 2
```

Keep `--max-in-flight` bounded when the metadata inventory is large. Use the
campaign report and free-space checks to confirm that a change improves the
measured rate without increasing failed or non-current outputs.

## Local profiling

```bash
uv run python scripts/rna/analyze_processing_times.py
python -m cProfile -o output/processing.pstats \
  scripts/rna/process_species.py --species apis_mellifera --dry-run
```

Write profiles and plots under `output/` or the selected campaign data root.
Do not benchmark by starting a competing producer against a live data root.

## Reliability rules

- Preserve partial transfer files for resumable retries.
- Do not delete raw inputs until current quantification provenance validates.
- Keep producer and downstream locks separate.
- Recheck the progress database, process tree, disk space, and recent log after
  each tuning change.
- Report observed throughput separately from queue ETA and scientific
  completion.
