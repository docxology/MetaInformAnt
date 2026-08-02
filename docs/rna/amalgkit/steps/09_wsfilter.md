# `amalgkit wsfilter`

`wsfilter` is the current within-species filtering step after merge. It
creates filtered expression inputs for finalization and records the filtering
parameters used for the run.

## Example

```bash
amalgkit wsfilter \
  --out_dir /Volumes/blue/data/amalgkit/pogonomyrmex_barbatus/work \
  --input_dir /Volumes/blue/data/amalgkit/pogonomyrmex_barbatus/work/merge \
  --metadata /Volumes/blue/data/amalgkit/pogonomyrmex_barbatus/work/metadata/metadata_selected.tsv \
  --batch no \
  --threads 8
```

Use only options shown by `amalgkit wsfilter --help` for the installed
version. The MetaInformAnt wrapper filters its own orchestration settings
before invoking the CLI.

## Validation

Before execution, verify that the input directory contains expression tables
and that metadata run IDs reconcile with the table columns. After execution,
record the output table, excluded samples/features, parameters, and row/column
counts in the project evidence manifest.
