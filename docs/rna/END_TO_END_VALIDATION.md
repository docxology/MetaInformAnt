# End-to-end validation

Validation is evidence production, not a process-exit check. Use the current
project validators in the following order.

## Static and configuration checks

```bash
uv run python scripts/rna/validate_configs.py
uv run python scripts/rna/validate_all_species_workflow.py
uv run python projects/hymenoptera_amalgkit/scripts/validate_project_docs.py
```

## Runtime status checks

```bash
export AMALGKIT_DATA_ROOT=/Volumes/blue/data/amalgkit
uv run python projects/hymenoptera_amalgkit/scripts/report_campaign_status.py \
  --data-root "$AMALGKIT_DATA_ROOT"
uv run python projects/hymenoptera_amalgkit/scripts/generate_pipeline_report.py \
  --data-root "$AMALGKIT_DATA_ROOT"
```

The report must distinguish configured, materialized, selected, quantified,
finalized, and analysis-ready species. It must retain failed, incomplete, and
unavailable states.

## Completion gates

For each species, verify:

1. metadata and selected tables are readable and hash-recorded;
2. selected run IDs reconcile to quantified samples;
3. every abundance table is non-empty, structurally valid, and sidecar-bound;
4. merge, within-species filtering, finalization, and sanity receipts validate;
5. the final matrix contains the expected sample and feature axes;
6. the evidence manifest is regenerated from the current data root;
7. only species passing all gates enter cross-species analysis.

## Cross-species gate

```bash
uv run python projects/hymenoptera_amalgkit/scripts/prepare_cross_species_inputs.py \
  --data-root "$AMALGKIT_DATA_ROOT" --check
uv run python projects/hymenoptera_amalgkit/scripts/run_cross_species_analysis.py \
  --data-root "$AMALGKIT_DATA_ROOT" --check
```

Do not convert missing evidence or timeouts into zeros. Preserve the status and
the run record, then repair or explicitly exclude the affected species.
