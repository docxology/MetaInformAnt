# AGENTS.md — `MetaInformAnt/src/metainformant/protein/sequence`

Source module under the METAINFORMANT bioinformatics toolkit (`protein` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `alignment.py` — Protein sequence alignment utilities.
- `proteomes.py` — Protein proteome analysis and taxonomy utilities.
- `sequences.py` — Protein sequence processing and analysis utilities.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-protein-sequence` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
