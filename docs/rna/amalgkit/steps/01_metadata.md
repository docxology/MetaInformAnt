# `amalgkit metadata`

This page documents the metadata stage in the pinned Amalgkit `0.16.32`
contract. The installed CLI help is authoritative for defaults and option
names; MetaInformAnt passes only options in the current command registry.

## Purpose

`metadata` queries and normalizes public SRA metadata before selection. It is
an input-acquisition step, not evidence that reads were downloaded or
quantified. A reproducible run records the search string, email/contact
configuration, output path, source release, and generated metadata hashes.

## Current example

```bash
amalgkit metadata \
  --out_dir /path/to/work \
  --search_string '"Camponotus floridanus"[Organism] AND RNA-Seq[Strategy] AND Illumina[Platform]' \
  --resolve_names yes \
  --threads 8
```

For the project configurations, the wrapper supplies the species-specific
search expression and remaps `output/amalgkit` to `AMALGKIT_DATA_ROOT` when
that environment variable is set.

## Important 0.16.32 options

| Option | Role |
|---|---|
| `--out_dir` | Workspace containing metadata outputs |
| `--search_string` | Entrez search expression |
| `--species_tsv` | Batch species input with `scientific_name` and optional `species_token` |
| `--organ_terms_tsv` | Optional title-term mapping for title-based modes |
| `--mode` | Metadata query mode, including title-based modes |
| `--resolve_names yes|no` | Resolve taxonomy names while preserving the source name |
| `--download_dir` / `--download_lock_dir` | Shared cache and cross-process lock locations |
| `--internal_jobs` / `--internal_cpu_budget` | Current internal concurrency controls |
| `--ncbi_metadata_max_concurrency` | Shared Entrez request limit |

Use `amalgkit metadata --help` from the pinned environment before adding an
option to a YAML file.

## Output contract

The normal per-species workspace contains:

```text
work/
└── metadata/
    ├── metadata.tsv
    ├── metadata_original.tsv
    └── metadata_specieswise/       # batch/species-table mode when enabled
```

Exact filenames can vary with the query mode, so downstream steps should use
the generated metadata path recorded in the workflow manifest. The project
report separately records total metadata rows, selected rows, valid abundance
files, and selected-ID intersections.

## Taxonomy and integrity

Taxonomy columns are metadata fields, not a selection guarantee. Missing or
inconsistent run accessions are written to the project validation ledger and
must remain visible in the report. Do not repair a run identifier by guessing
from a directory name.

## Project hand-off

After metadata completes, initialize or verify the current selection rules in
[02_dataset.md](02_dataset.md), then run
[03_select.md](03_select.md). The default project rules are in
`config/amalgkit/select_rules.tsv` and are checked against Amalgkit 0.16.32
before a configuration is accepted.

## References

- [Amalgkit repository and current command documentation](https://github.com/kfuku52/amalgkit)
- [NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25499/)
- [SRA](https://www.ncbi.nlm.nih.gov/sra/)
