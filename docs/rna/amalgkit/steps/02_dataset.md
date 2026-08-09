# `amalgkit dataset`

`dataset` is the current Amalgkit workspace and rule-set initialization
command. It prepares the files consumed by `metadata`, `select`, and
downstream commands; it is not a per-species completion stage.

## Inspect available assets

```bash
amalgkit dataset --list
```

The installed 0.16.38 CLI is authoritative for the available dataset names
and bundled rule sets.

## Initialize a workspace

```bash
amalgkit dataset \
  --name init \
  --out_dir /path/to/work \
  --overwrite no
```

Use a new or explicitly selected workspace. Never overwrite a mounted
analysis root merely to make a command line example succeed.

## Export a bundled rule set

```bash
amalgkit dataset \
  --name init \
  --rule_set base \
  --out_dir /path/to/work \
  --overwrite no
```

Bundled rule sets are useful starting points, but they are not automatically
appropriate for this cohort. The checked-in Hymenoptera contract uses
`config/amalgkit/select_rules.tsv`, which deliberately imposes no
plant-oriented tissue restriction and is validated against Amalgkit 0.16.38.

## Provenance

Record the exact command, Amalgkit version, rule-file hash, workspace path,
and resulting metadata hash in the run ledger. A generated workspace does not
prove that metadata was retrieved, samples were selected, or reads were
processed.

## References

- [Amalgkit repository and current command documentation](https://github.com/kfuku52/amalgkit)
- [Selection rules](03_select.md)
