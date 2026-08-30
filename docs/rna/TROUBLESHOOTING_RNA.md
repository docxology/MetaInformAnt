# RNA troubleshooting

## The CLI is unavailable

```bash
uv run python -c 'from metainformant.rna.amalgkit import check_cli_available; print(check_cli_available())'
amalgkit --help
```

Use the project environment and remove unrelated `VIRTUAL_ENV` contamination
when invoking `uv`:

```bash
env -u VIRTUAL_ENV uv run python scripts/rna/check_environment.py
```

## A step rejects an option

Compare the command with the installed help output. MetaInformAnt filters
wrapper-only settings, but a project configuration should still use current
option names and types. Save the rejected command and environment version in
the run record.

## No input tables found

Check the configured `input_dir`, species name, and data root. For
`wsfilter`, the input is normally the current merge output. For `finalize`, it
is normally the `wsfilter` output. Require a non-empty TSV with a valid header;
an empty directory is not a successful prerequisite.

## Quantification count disagrees with metadata

Compare run IDs using the project reconciliation tools. Inspect duplicate
accessions, failed downloads, invalid abundance files, and rows excluded by
the selection filters. Do not inflate counts by counting directories.

## The report changed unexpectedly

Regenerate the evidence manifest first, then the pipeline report. Compare
configuration hashes, data-root path, tool versions, selected rows, valid
abundance files, and downstream output counts. Preserve the previous report as
a dated snapshot rather than overwriting it while dependencies are still
being checked.

## Disk pressure

Keep raw reads and indexes on the external volume, quantify in bounded batches,
and remove temporary reads only after successful validated quantification.
Check active processes before any cleanup. Never recursively delete the active
data root.

## Incomplete species

An absent data directory is a data-availability state. Keep the configuration,
report it as not materialized, and either complete the species run or narrow
the declared analysis scope. Never mark it complete because another artifact
mentions the species.
## After the reference-manifest fix: verifying quarantine re-queue completion

Before the reference alias manifest became content-deterministic (2026-08-30),
every orchestrator restart changed the manifest hash and quarantined all prior
quantification of the affected species ("reference manifest checksum
mismatch"). A one-time re-queue sweep re-audits quarantined samples after the
fix deploys; verify it with read-only access to the progress DB:

```bash
sqlite3 "file:/Volumes/external_drive/Data/amalgkit/pipeline_progress.db?mode=ro" \
  "SELECT status, reason, COUNT(*) FROM quantification_audit GROUP BY 1, 2;"
```

Post-migration expectations:

- The quarantine count for the reason `reference manifest checksum mismatch`
  returns to its pre-campaign baseline (quarantines caused by genuinely
  missing/corrupt outputs, not by restart-varying manifest bytes). The
  baseline is per-species and per-day in `quantification_audit.checked_at`;
  do not compare a running campaign's instantaneous count against zero.
- Samples re-audited after the deploy either re-quantify once (if their
  recorded manifest hash is the old timestamped one) or return to
  `quantified`; no further `reference manifest checksum mismatch` rows should
  appear with fresh `checked_at` timestamps. If new ones do, a restart-varying
  byte still feeds the manifest digest — stop and re-audit the writer.

Distinguishing legacy-sidecar quarantines from new ones:

- `quantification_audit.checked_at` is the wall-clock time of the audit.
  Rows predating the fix deploy are legacy artifacts of the timestamped
  manifest; rows with post-deploy timestamps describe the sample's current
  state.
- `quantification_audit.contract_id` is the content-deterministic
  quantification contract digest. Legacy rows whose sidecars predate the
  contract schema carry a different contract id than samples re-quantified
  after the deploy; group on `(species, contract_id)` to see which contract
  generation a quarantine belongs to:

```bash
sqlite3 "file:/Volumes/external_drive/Data/amalgkit/pipeline_progress.db?mode=ro" \
  "SELECT species, substr(contract_id, 1, 12), COUNT(*) \
   FROM quantification_audit WHERE status = 'invalid' \
   GROUP BY 1, 2 ORDER BY 3 DESC LIMIT 20;"
```

A single dominant contract id per species after the sweep confirms all
surviving quarantines share one cause (and one fix path), rather than a mix of
legacy timestamp artifacts and fresh contract mismatches.

## The raw validation marker or SRA witness keeps re-validating

Both witnesses are content-deterministic (no wall-clock fields) and are
compared against the on-disk inputs' size/mtime fingerprint on resume. If a
sample's cached validation is replayed from scratch on every restart:

1. Confirm the input files were not touched between restarts (an mtime bump is
   enough to invalidate a witness).
2. Confirm the writer on this deploy produced a deterministic payload:

```bash
uv run python - <<'PY'
import json, subprocess, time
from pathlib import Path
# Any deterministic writer: two writes with the same inputs must be
# byte-identical. See tests/rna for the idempotency regression suite.
PY
```

3. If the payload is deterministic but the witness still misses, the input
   fingerprint (size/mtime) changed — that is expected invalidation, not a
   determinism bug.
