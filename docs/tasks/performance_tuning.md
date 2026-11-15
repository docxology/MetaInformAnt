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

Without `--workers` and `--threads` the producer derives every budget from
the host; `--print-profile` prints the result. A bounded diagnostic run with
explicit budgets:

```bash
env -u VIRTUAL_ENV uv run --frozen --no-sync python scripts/rna/run_all_species.py \
  --config-dir projects/hymenoptera_amalgkit/config/amalgkit \
  --data-root "$AMALGKIT_DATA_ROOT" \
  --cohort "$AMALGKIT_DATA_ROOT/cohort/hymenoptera27-20261005" \
  --workers 4 --threads 8 --quant-slots 2 \
  --fastq-threads 1 --fastq-slots 1 --validation-slots 2
```

The same budgets can be set through `PIPELINE_WORKERS`, `PIPELINE_THREADS`
and the `AMALGKIT_PIPELINE_*` variables (`AMALGKIT_PIPELINE_QUANT_SLOTS`,
`AMALGKIT_PIPELINE_FASTQ_THREADS`, `AMALGKIT_PIPELINE_FASTQ_SLOTS`,
`AMALGKIT_PIPELINE_VALIDATION_SLOTS`); a command-line option takes
precedence.

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
