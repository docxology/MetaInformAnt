# Amalgkit configuration FAQ

The canonical Hymenoptera configuration inventory is
`projects/hymenoptera_amalgkit/config/amalgkit/`. The selected
`AMALGKIT_DATA_ROOT`, current SQLite progress, provenance receipts, and
regenerated reports determine what actually completed.

## What is the command chain?

```text
metadata → select → getfastq → integrate → quant → merge
         → wsfilter → finalize → sanity
```

Inspect it without starting a run:

```bash
export AMALGKIT_DATA_ROOT=/Volumes/blue/data/amalgkit
uv run python scripts/rna/run_all_species.py \
  --config-dir projects/hymenoptera_amalgkit/config/amalgkit \
  --data-root "$AMALGKIT_DATA_ROOT" --dry-run
```

Run the complete campaign through the project lock owner:

```bash
bash projects/hymenoptera_amalgkit/scripts/run_full_campaign.sh
```

## How are resources selected?

Use the bounded project profile or explicit script options:

```bash
AMALGKIT_PIPELINE_PROFILE=local-throughput \
AMALGKIT_PIPELINE_THREADS=10 \
AMALGKIT_PIPELINE_QUANT_SLOTS=10 \
AMALGKIT_PIPELINE_FASTQ_SLOTS=1 \
bash projects/hymenoptera_amalgkit/scripts/run_full_campaign.sh
```

The local-throughput profile uses `2 x host CPUs` of sample workers bounded to
`16..24`, ten total quant threads, six validation slots, and twice the actual
worker count in flight. Check both external and system free-space floors before
raising `AMALGKIT_PIPELINE_FASTQ_SLOTS`; acquisition and mounted-volume
contention can dominate wall time even when CPU is idle. Preview the resolved
budget without starting a producer with:

```bash
bash projects/hymenoptera_amalgkit/scripts/run_full_campaign.sh --print-config
```

## How do I verify a run?

```bash
uv run python projects/hymenoptera_amalgkit/scripts/report_campaign_status.py \
  --data-root "$AMALGKIT_DATA_ROOT"
uv run python projects/hymenoptera_amalgkit/scripts/generate_pipeline_report.py \
  --data-root "$AMALGKIT_DATA_ROOT"
```

The report distinguishes configured, materialized, selected, quantified,
finalized, and analysis-ready species. A non-empty directory or process exit
code is not sufficient.

## What does a timeout mean?

A timeout is unavailable evidence within a bounded window. Preserve the
failure state, inspect the species log and receipt columns, and retry at a
resume boundary. Do not convert a timeout into a zero or a completed sample.

## How are incomplete species handled?

Keep their metadata, selected rows, logs, partial outputs, and failure records.
Repair only the missing current stage, regenerate the report, and admit the
species to cross-species analysis only after all matrix and provenance gates
pass.

## How are tissue labels handled?

Shared mappings and explicit sample/project patches live in the canonical
project config directory:

```text
projects/hymenoptera_amalgkit/config/amalgkit/tissue_mapping.yaml
projects/hymenoptera_amalgkit/config/amalgkit/tissue_patches.yaml
```

Source labels and recoding decisions remain in the evidence record.

## Related resources

- [Project configuration README](../../projects/hymenoptera_amalgkit/config/amalgkit/README.md)
- [Project running guide](../../projects/hymenoptera_amalgkit/doc/00_setup/04_running_the_pipeline.md)
- [Analysis readiness](../../projects/hymenoptera_amalgkit/doc/manuscript/analysis_readiness.md)
