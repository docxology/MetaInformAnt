# Configuration maintenance record — Amalgkit

## Scope

`config/amalgkit/` contains the checked-in 27-species Hymenoptera YAML
inventory, templates, and shared tissue maps. The inventory defines planned
scope; it does not assert that external data or downstream matrices exist.

## Current contract

Every species configuration uses the current chain:

```text
metadata → select → getfastq → integrate → quant → merge
         → wsfilter → finalize → sanity
```

Provider flags are explicit Amalgkit data-source controls. Large inputs and
outputs resolve from `AMALGKIT_DATA_ROOT` when the project runners are used.
Do not replace a current report with a status copied from another data root.

## Validation

From the repository root:

```bash
uv run python scripts/rna/validate_configs.py
uv run python scripts/rna/validate_all_species_workflow.py
```

For the full project evidence boundary:

```bash
bash projects/hymenoptera_amalgkit/scripts/verify_setup.sh \
  --data-root "$AMALGKIT_DATA_ROOT" --require-data --report
```

Configuration validation confirms syntax and required fields only. The report
and evidence manifest confirm materialized outputs, sample identity, matrix
availability, and cross-species cohort integrity.

## Maintenance rules

- Update the configuration walkthrough, tests, and 27-species inventory
  together when adding or removing a species.
- Keep tissue recoding in `tissue_mapping.yaml` and
  `tissue_patches.yaml`; preserve source metadata.
- Use `redo: no` unless a deliberate rerun has been recorded.
- Regenerate the report, manuscript status, manifest, and evidence ledger after
  any configuration or data-root change.
