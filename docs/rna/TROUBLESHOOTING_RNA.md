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
