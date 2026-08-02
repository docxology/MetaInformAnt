# RNA file and storage contract

The repository stores code, configurations, tests, compact evidence, and
documentation. Raw reads, indexes, and large expression matrices live under
the explicitly configured `AMALGKIT_DATA_ROOT`.

## Canonical per-species layout

```text
$AMALGKIT_DATA_ROOT/<species>/
├── work/
│   ├── metadata/
│   ├── fastq/                 # temporary; removed after successful quantification
│   ├── quant/<run>/           # abundance.tsv or quant.sf
│   ├── merge/
│   ├── wsfilter/
│   ├── finalize/
│   ├── sanity/
│   ├── logs/
│   └── amalgkit.manifest.jsonl
├── index/                     # reference index when species-local
├── merged/                    # accepted external Amalgkit merge layout
└── evidence/                  # compact checksums and summaries
```

The planner resolves the configured `work_dir`; a report may inspect an
accepted external merge root when it is explicitly recorded, but it must not
silently combine species directories.

## Path rules

1. Use absolute paths for mounted data in real runs.
2. Keep one species and one configuration hash per work root.
3. Do not symlink a different species' quantification directory into a merge
   root.
4. Treat FASTQ as temporary only after the quantification output and manifest
   record have passed validation.
5. Save compact result checksums and counts; do not copy bulk reads into Git.

## Verification

```bash
export AMALGKIT_DATA_ROOT=/Volumes/blue/data/amalgkit
uv run python scripts/rna/check_pipeline_status.py --data-root "$AMALGKIT_DATA_ROOT"
```

For the complete Hymenoptera inventory and evidence contract, see
`projects/hymenoptera_amalgkit/doc/03_data_management/README.md`.
