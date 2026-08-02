# Amalgkit test-contract note

The active RNA wrapper contract targets Amalgkit `0.16.32`, pinned to the
release tag `v0.16.32` in `pyproject.toml`. The test suite must exercise the
current command surface and must not describe retired command-line options or
an older release as active functionality.

## Current coverage

- `test_rna_amalgkit_current.py` checks the command registry, current wrappers,
  allow-lists, exact-version parsing, and the single-process `cstmm` contract.
- `test_rna_amalgkit_cli_args.py` checks current boolean option serialization,
  metadata options, and removal of MetaInformAnt-only values before CLI calls.
- Workflow and streaming tests verify current `clean_fastq` handling,
  metadata/run-ID integrity, output paths, and real subprocess dispatch.
- The full parent test command is run without marker filters; the resulting
  log is checked for both `SKIPPED` and `deselected` before evidence is
  reported.

## Current command model

The installed release exposes `metadata`, `dataset`, `select`, `getfastq`,
`quant`, `merge`, `busco`, `cstmm`, `wsfilter`, `csfilter`, `finalize`,
`sanity`, `rerun`, and `integrate`. The default species workflow uses the
project-specific sequence documented in `docs/rna/amalgkit/README.md`; BUSCO,
cross-species normalization/filtering, dataset initialization, and failed-run
recovery remain explicit optional branches.

Selection parameters are validated through the checked-in
`config/amalgkit/select_rules.tsv` file, which is parsed by the installed
Amalgkit selection-rule reader.

## Re-run

From the repository root:

```bash
uv run python scripts/rna/validate_configs.py
uv run python -m pytest -q --run-slow
```

The command output, interpreter, Amalgkit version, repository revision, data
root, and no-skip/no-deselection check belong in the run record for any
publication or release claim.
