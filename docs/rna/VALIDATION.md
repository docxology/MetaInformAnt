# RNA validation protocol

The current protocol validates configuration, production state, downstream
matrices, and provenance in separate read-only or lock-owned phases.

## Configuration and plan

```bash
uv run python scripts/rna/validate_configs.py
uv run python scripts/rna/validate_all_species_workflow.py
uv run python scripts/rna/run_all_species.py \
  --config-dir projects/hymenoptera_amalgkit/config/amalgkit \
  --data-root "$AMALGKIT_DATA_ROOT" --dry-run
```

These commands validate YAML shape, the current command registry, selection
rules, and the fixed nine-stage plan without contacting public archives.

## Production state

```bash
uv run python projects/hymenoptera_amalgkit/scripts/report_campaign_status.py \
  --data-root "$AMALGKIT_DATA_ROOT"
uv run python projects/hymenoptera_amalgkit/scripts/generate_pipeline_report.py \
  --data-root "$AMALGKIT_DATA_ROOT"
```

Check that the report records configured species, metadata rows, selected rows,
quantified rows, valid abundance files, stage receipts, and failures.

## Matrix and provenance gates

For every species admitted to downstream analysis:

- metadata and selected run IDs reconcile;
- abundance tables have the expected columns and non-zero content;
- current quantification sidecars hash the exact abundance files;
- merge, `wsfilter`, `finalize`, and `sanity` receipts validate their inputs;
- the final matrix has stable sample and feature axes;
- the evidence manifest and manuscript status are regenerated.

## Cross-species analysis

```bash
uv run python projects/hymenoptera_amalgkit/scripts/prepare_cross_species_inputs.py \
  --data-root "$AMALGKIT_DATA_ROOT" --check
uv run python projects/hymenoptera_amalgkit/scripts/run_cross_species_analysis.py \
  --data-root "$AMALGKIT_DATA_ROOT" --check
```

The cross-species lane is conditional. Incomplete species remain visible in
the status ledger and are not silently imputed into a comparison.
