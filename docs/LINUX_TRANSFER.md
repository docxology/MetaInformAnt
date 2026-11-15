# Moving the Amalgkit campaign to another host

A fresh host for the frozen Hymenoptera cohort follows the project's hand-off
guide, `projects/hymenoptera_amalgkit/doc/01_infrastructure/22_compute_handoff.md`
(in the private project repository that is checked out as the
`projects/hymenoptera_amalgkit` submodule). The guide's sequence is:

1. `projects/hymenoptera_amalgkit/scripts/bootstrap_compute_host.sh` prepares
   the host and its data root (options are in the guide).
2. `bash projects/hymenoptera_amalgkit/scripts/run_full_campaign.sh --species NAME`
   runs that host's species in frozen-cohort mode.

Each host quantifies only the species assigned to it (no species runs on two
hosts) and owns one data root (`AMALGKIT_DATA_ROOT`).
The campaign writer lock on `<data-root>/.campaign.flock` admits a single
writer: the launcher holds it for its whole run, and a second launcher,
producer or finalization pass against the same data root exits 73.

The library-level steps behind the bootstrap, installing each cohort species'
kallisto index and reference manifest under the new data root, are described
in [Preparing a fresh data root for cohort species](rna/amalgkit/steps/06_quant_advanced.md#preparing-a-fresh-data-root-for-cohort-species).

## Environment

```bash
uv sync --frozen --extra rna
env -u VIRTUAL_ENV uv run --frozen --no-sync python scripts/rna/check_pipeline_status.py \
  --data-root "$AMALGKIT_DATA_ROOT" --verbose
```

`uv sync --frozen --extra rna` installs the pinned Amalgkit 0.16.108 from
`uv.lock`. The status command reads `<data-root>/pipeline_progress.db` and does
not take the writer lock.
