# Current Amalgkit command reference

The project pins Amalgkit 0.16.33. The command registry in
`metainformant.rna.amalgkit.commands` is the source of truth for accepted
options; confirm the installed binary with `amalgkit --help`.

## Per-species stages

| Stage | Owner | Required evidence |
|---|---|---|
| `metadata` | streaming producer | source metadata table and query record |
| `select` | streaming producer | selected table and rule-file hash |
| `getfastq` | streaming producer | validated FASTQ or bounded failure record |
| `integrate` | streaming producer | integrated metadata with local paths |
| `quant` | streaming producer | abundance table and exact quantification sidecar |
| `merge` | checkpoint runner | merged matrix and input manifest |
| `wsfilter` | checkpoint runner | filtered matrix and exclusion record |
| `finalize` | checkpoint runner | final matrix and normalization record |
| `sanity` | checkpoint runner | structural and provenance checks |

## Inspect the current inventory

```bash
export AMALGKIT_DATA_ROOT=/Volumes/blue/data/amalgkit
uv run python scripts/rna/run_all_species.py \
  --config-dir projects/hymenoptera_amalgkit/config/amalgkit \
  --data-root "$AMALGKIT_DATA_ROOT" --dry-run
```

## Run the campaign and inspect state

```bash
bash projects/hymenoptera_amalgkit/scripts/run_full_campaign.sh
uv run python projects/hymenoptera_amalgkit/scripts/report_campaign_status.py \
  --data-root "$AMALGKIT_DATA_ROOT"
```

## Direct command inspection

```bash
amalgkit metadata --help
amalgkit select --help
amalgkit getfastq --help
amalgkit integrate --help
amalgkit quant --help
amalgkit merge --help
amalgkit wsfilter --help
amalgkit finalize --help
amalgkit sanity --help
```

Do not copy flags between releases. Record the exact installed version and
command line in the run evidence.
