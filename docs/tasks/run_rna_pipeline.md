# RNA-seq Pipeline Quick Reference

This guide documents the current ENA-first, resumable Amalgkit producer. It
does not claim cohort completion or biological inference from a partial data
root.

## Prepare and inspect

```bash
export AMALGKIT_DATA_ROOT=/path/to/amalgkit-data
uv run python scripts/rna/validate_all_species_workflow.py
uv run python scripts/rna/run_all_species.py \
  --config-dir projects/hymenoptera_amalgkit/config/amalgkit \
  --data-root "$AMALGKIT_DATA_ROOT" --dry-run
```

Reference setup is configuration-driven:

```bash
uv run python scripts/rna/setup_genome.py \
  --config projects/hymenoptera_amalgkit/config/amalgkit/amalgkit_apis_mellifera.yaml \
  --verify-only
```

Use `--dry-run` first for a real setup. Reference and index outputs belong in
the selected data root or another documented `output/` location.

## Run the producer

Run one configured species:

```bash
uv run python scripts/rna/process_species.py \
  --species apis_mellifera \
  --config-dir projects/hymenoptera_amalgkit/config/amalgkit \
  --data-root "$AMALGKIT_DATA_ROOT" \
  --dry-run
```

Run the selected cohort:

```bash
uv run python scripts/rna/run_all_species.py \
  --config-dir projects/hymenoptera_amalgkit/config/amalgkit \
  --data-root "$AMALGKIT_DATA_ROOT" \
  --workers 4 --threads 8 --quant-slots 2
```

The producer owns metadata, selection, acquisition, validation, integration,
and quantification. It resumes from `pipeline_progress.db`, retains failed
and partial states, and writes provenance before an output becomes reusable.

## Monitor and finalize

```bash
uv run python scripts/rna/check_pipeline_status.py \
  --data-root "$AMALGKIT_DATA_ROOT" --verbose

env -u VIRTUAL_ENV uv run python \
  projects/hymenoptera_amalgkit/scripts/report_campaign_status.py \
  --data-root "$AMALGKIT_DATA_ROOT" --json
```

Do not start downstream finalization while the producer lock or producer
process is active. After a durable producer stop, use the project checkpoint
runner for `merge -> wsfilter -> finalize -> sanity`; cross-species analysis
requires finalized, provenance-qualified inputs.

## Evidence boundaries

- Executable readiness means the code, dependencies, and commands validate.
- Cohort readiness requires a current manifest and reconciled progress database.
- Descriptive analysis requires current finalized matrices and receipts.
- Biological inference requires an explicitly approved analysis design and
  evidence bundle.

Plain or stale `abundance.tsv` files remain comparator/audit evidence unless
their current provenance sidecars validate the exact source, configuration,
reference, and table hash.
