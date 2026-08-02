# `amalgkit finalize`

`finalize` produces the analysis-ready expression tables for a species after
current filtering. It is the required downstream endpoint for the default
per-species workflow.

## Example

```bash
amalgkit finalize \
  --out_dir /Volumes/blue/data/amalgkit/pogonomyrmex_barbatus/work \
  --input_dir /Volumes/blue/data/amalgkit/pogonomyrmex_barbatus/work/wsfilter \
  --metadata /Volumes/blue/data/amalgkit/pogonomyrmex_barbatus/work/metadata/metadata_selected.tsv \
  --batch no \
  --threads 8
```

The exact `--norm`, batch, and optional effect-handling settings must be
recorded with the run. Do not describe the output as biologically normalized
or batch-free unless those properties were actually established by the
declared method and diagnostics.

## Completion evidence

Require a non-empty final expression table, a metadata table or explicit
sample map, matching feature/sample axes, a successful command record, and a
sanity check. The project report generator uses these checks rather than
directory existence alone.
