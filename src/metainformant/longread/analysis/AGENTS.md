# AGENTS.md — `MetaInformAnt/src/metainformant/longread/analysis`

Source module under the METAINFORMANT bioinformatics toolkit (`longread` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `modified_bases.py` — Modified base detection from long-read sequencing data.
- `phasing.py` — Haplotype phasing from long-read sequencing data.
- `structural.py` — Structural variant detection from long-read sequencing data.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-longread-analysis` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
