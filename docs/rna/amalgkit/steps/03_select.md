# `amalgkit select`

Amalgkit `0.16.59` reads the selection policy from a validated
`select_rules.tsv` file. Selection parameters are not emitted as unsupported
ad-hoc command-line flags by the MetaInformAnt workflow.

## Project command

```bash
amalgkit select \
  --out_dir /path/to/work \
  --metadata /path/to/work/metadata/metadata.tsv \
  --select_rules_tsv config/amalgkit/select_rules.tsv \
  --threads 8
```

When `--select_rules_tsv` is omitted, Amalgkit looks for
`<out_dir>/select_rules.tsv`. The project passes the checked-in path
explicitly so the selection policy is reproducible and independent of a
stale data directory.

## Current project rule contract

The Hymenoptera rule file records:

| Parameter | Value | Meaning |
|---|---:|---|
| `min_nspots` | `0` | No minimum spot cutoff is imposed by this cohort contract |
| `max_sample` | `999999` | Large safety bound; ordinary species are not downsampled by configuration |
| `mark_missing_rank` | `none` | Missing taxonomy rank alone does not exclude a sample |
| `mark_redundant_biosamples` | `no` | Preserve redundant records for auditability |

No tissue-group parameter is enabled. This project intentionally retains a
broad public RNA-seq cohort; biological tissue harmonization is reported as a
separate metadata/design issue rather than silently imposed by a plant-oriented
default rule set.

## Outputs

The current command writes selected metadata and audit tables below the
workspace metadata directory. Depending on the batch mode, inspect:

```text
work/metadata/
├── metadata_selected.tsv
├── metadata_qualified.tsv
├── metadata_excluded.tsv
├── select_summary.tsv
└── select_queue.tsv
```

The exact generated names are recorded in the workflow manifest and status
report. A selected table is not read evidence; it becomes analysis evidence
only after `getfastq`, quantification, merge, filtering, finalization, and
sanity checks produce validated outputs for the same run IDs.

## Validation

```bash
uv run python scripts/rna/validate_configs.py
amalgkit select --help
```

The configuration validator parses the checked-in rule file with the
installed Amalgkit selection parser and rejects unsupported direct YAML keys.
Do not infer completion from a non-empty directory or from a historical
sample count.

## References

- [Amalgkit repository and current command documentation](https://github.com/kfuku52/amalgkit)
- [`select_rules.tsv` in this repository](../../../../config/amalgkit/select_rules.tsv)
