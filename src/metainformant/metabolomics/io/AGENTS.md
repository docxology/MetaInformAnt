# AGENTS.md — `MetaInformAnt/src/metainformant/metabolomics/io`

Source module under the METAINFORMANT bioinformatics toolkit (`metabolomics` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `formats.py` — Mass spectrometry file format reading and writing.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-metabolomics-io` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
