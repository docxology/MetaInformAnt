# Current RNA configuration

The canonical Hymenoptera inventory is
`projects/hymenoptera_amalgkit/config/amalgkit/`. Each runnable species file
is named `amalgkit_<species>.yaml`; template, test, and cross-species files are
excluded from production discovery.

## Minimal structure

```yaml
work_dir: output/amalgkit/apis_mellifera/work
log_dir: output/amalgkit/apis_mellifera/logs
threads: 60
species_list:
  - Apis_mellifera
taxon_id: 7460

genome:
  accession: GCF_003254395.2
  dest_dir: output/amalgkit/apis_mellifera/genome

steps:
  metadata:
    search_string: '"Apis mellifera"[Organism] AND RNA-Seq[Strategy]'
    redo: no
  select:
    select_rules_tsv: config/amalgkit/select_rules.tsv
  getfastq:
    redo: no
  integrate: {}
  quant:
    redo: no
  merge: {}
  wsfilter: {}
  finalize: {}
  sanity: {}
```

The supported step options are validated by
`scripts/rna/validate_configs.py` against the pinned Amalgkit command
registry. Selection thresholds belong in `select_rules.tsv`, not in unrelated
producer settings.

## Path contract

Relative paths in project YAML files are resolved from the project/repository
root as appropriate. Large artifacts are remapped below
`AMALGKIT_DATA_ROOT` by the current workflow engine. Keep metadata, selected
tables, indexes, quantification tables, stage outputs, logs, and receipts
under the same data root for a run.

## Runtime settings

Use CLI options for direct producers and the project shell variables for the
full campaign:

```bash
export AMALGKIT_DATA_ROOT=/Volumes/blue/data/amalgkit
export AMALGKIT_PIPELINE_WORKERS=8
export AMALGKIT_PIPELINE_THREADS=8
export AMALGKIT_PIPELINE_QUANT_SLOTS=4
export AMALGKIT_PIPELINE_FASTQ_THREADS=1
export AMALGKIT_PIPELINE_FASTQ_SLOTS=1
export AMALGKIT_PIPELINE_COMPRESSION_THREADS=1
export AMALGKIT_PIPELINE_VALIDATION_SLOTS=4
export AMALGKIT_PIPELINE_MAX_IN_FLIGHT=16
```

The launcher checks both external and system free-space floors before creating
new work. Do not increase workers or quant slots based only on an idle CPU
snapshot; acquisition and mounted-volume contention may be the bottleneck.

## Validate and inspect

```bash
uv run python scripts/rna/validate_configs.py
uv run python scripts/rna/validate_all_species_workflow.py
uv run python scripts/rna/run_all_species.py \
  --config-dir projects/hymenoptera_amalgkit/config/amalgkit \
  --data-root "$AMALGKIT_DATA_ROOT" --dry-run
```

All three commands are read-only. They do not contact ENA or modify the
campaign database.

## Add a species

1. Copy the project template to `amalgkit_<species>.yaml`.
2. Set `species_list`, taxonomy, reference accession, paths, and metadata query.
3. Use the shared `select_rules.tsv` unless a reviewed rule file is required.
4. Run both validators and a single-species dry run.
5. Execute only after the reference/index and storage contract are verified.

## References

- [Project configuration README](../../projects/hymenoptera_amalgkit/config/amalgkit/README.md)
- [Project storage contract](../../projects/hymenoptera_amalgkit/doc/01_infrastructure/02_storage_contract.md)
- [Current workflow](workflow.md)
