# Pipeline Monitoring

How to monitor the status of active amalgkit pipeline runs.

## Live Log Monitoring

The main log for the all-species run is:

```bash
# Follow the current campaign log in real time
tail -f "$AMALGKIT_DATA_ROOT"/results/full_campaign_*.log

# Filter for key status lines
grep -E "Downloaded|Quantified|Reclaimed|Done|Failed|ERROR|Phase" \
  "$AMALGKIT_DATA_ROOT"/results/full_campaign_*.log | tail -50
```

Key log tokens:
- `Downloaded` — a valid local FASTQ was acquired
- `Quantified` or `Done` — sample quantification completed
- `Reclaimed` — validated raw inputs were removed after quantification
- `Failed` or `ERROR` — acquisition, quantification, or validation failure
- `Phase 1 Complete` — all configured sample tasks have entered the executor

The scheduling line also records the requested and effective worker counts,
total and per-worker `quant_threads`, `fasterq_threads`, `fasterq_slots`,
`compression_threads`, `max_in_flight`, and host CPU count. Treat those values
as the authoritative resource profile; do not infer concurrency from the
number of configured species.

## Active Process Monitor

```bash
# List active acquisition, quantification, and orchestration processes
pgrep -af 'curl .*fastq|fasterq-dump|kallisto quant|run_all_species.py' || true

# List active kallisto processes
pgrep -af 'kallisto quant' || true

# See the launcher and Python workers
pgrep -af 'run_full_campaign.sh|run_all_species.py|streaming_orchestrator' || true
```

If the campaign was interrupted, do not start a second copy while one is
running. Check the process list above, then resume in a persistent `tmux`
session from the current external data root:

```bash
tmux new-session -d -s hymenoptera-campaign \
  'env AMALGKIT_DATA_ROOT="$AMALGKIT_DATA_ROOT" \
    AMALGKIT_PIPELINE_PROFILE=local-throughput \
    AMALGKIT_PIPELINE_THREADS=10 \
    AMALGKIT_PIPELINE_FASTQ_SLOTS=1 \
    AMALGKIT_RECLAIM_RAW_AFTER_QUANT=yes \
    bash projects/hymenoptera_amalgkit/scripts/run_full_campaign.sh'
tmux capture-pane -pt hymenoptera-campaign:0 -S -40
# Reattach interactively when desired: tmux attach -t hymenoptera-campaign
```

Before changing the extraction lane, inspect the resolved settings without
creating a lock or touching the data root:

```bash
bash projects/hymenoptera_amalgkit/scripts/run_full_campaign.sh --print-config
```

Keep `AMALGKIT_PIPELINE_FASTQ_SLOTS=1` while host system free space is below
`32 GiB`. A value of `2` is appropriate only after a measured extraction
bottleneck and a fresh storage/I/O check; Blue free space alone is not enough
because `fasterq-dump` temporary files use the host scratch volume.

The worker loop pauses new downloads when the external root has less than
the configured external free-space floor (8 GiB by default) or the configured
host-system reserve is reached. Existing
quantification workers and provenance-gated reclamation can continue while
the acquisition workers wait.
After metadata discovery, samples with existing final raw inputs are scheduled
before shared SRA-cache databases, resumable partials, and new acquisitions,
so a resumed campaign uses mounted data immediately and does not redownload
validated inputs.

The current evidence-aware status command checks the progress database and
downstream outputs under the selected external data root. A readable
`abundance.tsv` without a verified quantification-contract sidecar is not
counted as eligible quantification. Compatible runtime drift remains eligible
and is reported separately in the audit counts.

```bash
uv run python scripts/rna/check_pipeline_status.py \
  --data-root "$AMALGKIT_DATA_ROOT" --verbose
```

## Downstream checkpoint status

The current finalization helper is the only writer for merge, wsfilter,
finalize, and sanity outputs. It is safe to rerun after the producer stops:

```bash
bash projects/hymenoptera_amalgkit/scripts/run_all_finalization.sh \
  --data-root "$AMALGKIT_DATA_ROOT" --dry-run
```

Its execution mode skips a stage only when the contract provenance sidecar and
all recorded input/output digests still match.

## Quick Health Checks

```bash
# How many current quantification outputs exist for one species?
find "$AMALGKIT_DATA_ROOT/apis_mellifera_all/work/quant" \
  -mindepth 2 -maxdepth 2 -name abundance.tsv -size +100c | wc -l

# Any samples stuck downloading (large .fastq.gz partial files)?
find "$AMALGKIT_DATA_ROOT" -path '*/work/getfastq/*' \
  \( -name '*.fastq.gz' -o -name '*.sra' -o -name '*.part' \) -size +1G 2>/dev/null | head

# Check disk space
df -h "$AMALGKIT_DATA_ROOT"

# Most recently modified quant directories
ls -lt "$AMALGKIT_DATA_ROOT/apis_mellifera_all/work/quant/" | head -20
```

## Log Quick-Counts (Without Scanning Disk)

```bash
LOG=$(ls -t "$AMALGKIT_DATA_ROOT"/results/full_campaign_*.log | head -1)
grep -E "Quantified|Done" "$LOG" | wc -l # completed
grep -E "Reclaimed" "$LOG" | wc -l # raw-input reclamations
grep -E "Failed|ERROR" "$LOG" | wc -l # failed
```

## Compatibility-first quantification

AMALGKIT runtime version drift is recorded in the quantification audit but does
not invalidate a complete, checksum-verified quantification contract. The
status report separates current, version-drift-compatible, legacy-unverified,
and invalid outputs. Legacy-unverified outputs are quarantined rather than
automatically downloaded again. Quarantined is the fail-closed
provenance-audit state — outputs that cannot be certified current against the
recorded reference manifest — not an error bucket. Quarantined samples
re-enter as re-quantification candidates under the preserve requantification
policy; only `failed` rows are genuine errors.

Run the compatibility reconciliation as a dry run first:

    uv run python scripts/rna/reconcile_quantification_compatibility.py --data-root "$AMALGKIT_DATA_ROOT"

Apply the audit only after reviewing its manifest:

    uv run python scripts/rna/reconcile_quantification_compatibility.py --data-root "$AMALGKIT_DATA_ROOT" --apply

The producer defaults to preserving compatible outputs. Explicit rebuilds use
the preserve, version-drift, or all requantification policy.

## Current cleanup contract

After a non-empty abundance table and exact current quantification sidecar are
written, the streaming orchestrator removes only that sample's validated raw
FASTQ/SRA inputs. Set `AMALGKIT_RECLAIM_RAW_AFTER_QUANT=no` to retain them.
Use `scripts/rna/reclaim_quantified_raw.py` for a separate dry-run or audited
reclamation pass. Stale logs, heartbeats, and removed-stage directories are
archived by `scripts/rna/clean_external_artifacts.py`; they are not evidence
of a current run.

See the project commit history for details.
