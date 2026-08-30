# External data-root setup

Large Amalgkit inputs and outputs belong on a dedicated data volume. The
current workflow accepts that location through `AMALGKIT_DATA_ROOT` and uses
it for work directories, logs, progress state, and evidence.

## Prepare and verify the mount

```bash
export AMALGKIT_DATA_ROOT=/Volumes/external_drive/Data/amalgkit
test -d "$AMALGKIT_DATA_ROOT"
test -w "$AMALGKIT_DATA_ROOT"
df -h "$AMALGKIT_DATA_ROOT"
uv run python scripts/rna/check_environment.py
```

The host system volume also needs free space for Python, temporary archives,
and process metadata. The project launcher checks both reserves before a run.

## Inspect without writing

```bash
uv run python scripts/rna/run_all_species.py \
  --config-dir projects/hymenoptera_amalgkit/config/amalgkit \
  --data-root "$AMALGKIT_DATA_ROOT" --dry-run
uv run python projects/hymenoptera_amalgkit/scripts/report_campaign_status.py \
  --data-root "$AMALGKIT_DATA_ROOT"
```

The dry run is the safe first check after mounting a volume. It confirms the
config inventory and selected root without starting acquisition.

## Execute and resume

```bash
bash projects/hymenoptera_amalgkit/scripts/run_full_campaign.sh
```

The campaign is resumable through SQLite state, `.part` downloads, current
quantification sidecars, and stage receipts. Preserve those files when a run
is interrupted. Never start a second producer against the same data root.

## Storage safety

- Keep the campaign lock and progress database on the selected root.
- Reclaim raw reads only after contract provenance validates the exact abundance
  file and the project receipt records the decision.
- Preserve metadata, selection tables, indexes, logs, receipts, and failed
  samples for reproducibility.
- Use the project cleanup tools with an explicit data root and report output.

See the [storage contract](../../projects/hymenoptera_amalgkit/doc/01_infrastructure/02_storage_contract.md)
for the complete retention and free-space rules.
