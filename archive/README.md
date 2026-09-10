# Archive

One-off migration tooling retained for provenance. Nothing here is imported by
`src/` or exercised by the test suite; treat checked-out source, tests, and
maintained docs under `scripts/` as the source of truth.

## Contents

- `scripts_reorganize/` — import-map construction and mechanical import
  rewrite scripts used during the 0.4.0 module reorganization
  (`__init__.py` shortcut imports → canonical submodule paths). The generated
  `import_map_report.txt` is a historical snapshot and may not describe the
  current checkout. See `scripts_reorganize/README.md` for original usage.
