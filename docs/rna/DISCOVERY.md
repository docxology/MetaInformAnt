# Species and sample discovery

Discovery is configuration-derived and data-root-specific. The current
launcher excludes template, test, and cross-species YAML files, then resolves
the selected root for the SQLite progress database and per-species outputs.

## Inspect configured species

```bash
uv run python scripts/rna/validate_all_species_workflow.py
uv run python scripts/rna/run_all_species.py \
  --config-dir projects/hymenoptera_amalgkit/config/amalgkit \
  --data-root "$AMALGKIT_DATA_ROOT" --dry-run
```

The dry-run inventory is deterministic for a fixed configuration directory. It
does not query NCBI or ENA and does not create progress state.

## Metadata and selection

In frozen-cohort mode (`--cohort`, the launcher's only mode) the producer does
not run `metadata` or `select`: its tasks come from the cohort's frozen
`inputs/<species>/work/metadata/metadata_selected.tsv`, which is never
rewritten. Without `--cohort`, the producer runs `metadata` and `select` only
with `--refresh-metadata`; a missing or stale metadata table is otherwise an
error. Selection is defined by the shared `select_rules.tsv` and the species
query in its YAML. A refresh run writes the resulting tables under:

```text
$AMALGKIT_DATA_ROOT/<species>/work/metadata/metadata.tsv
$AMALGKIT_DATA_ROOT/<species>/work/metadata/metadata_selected.tsv
```

Keep the source metadata, selected rows, rule-file hash, and query record
together. A discovered accession is not a selected sample and a selected
sample is not a quantified sample.

## Single-species diagnostic

```bash
uv run python scripts/rna/process_species.py \
  --species camponotus_floridanus \
  --config-dir projects/hymenoptera_amalgkit/config/amalgkit \
  --data-root "$AMALGKIT_DATA_ROOT" --dry-run
```

For real processing, use the same command without `--dry-run` only when the
full campaign lock is not held. The single-species command uses the exact
same streaming producer and progress schema as the cohort command.

## Status interpretation

```bash
uv run python projects/hymenoptera_amalgkit/scripts/report_campaign_status.py \
  --data-root "$AMALGKIT_DATA_ROOT"
```

Treat `pending`, `downloading`, `quantifying`, `failed`, unavailable metadata,
and missing outputs as distinct states. Do not infer completion from a folder,
an accession count, or a stale report.
