# AGENTS.md — `MetaInformAnt/src/metainformant/epigenome/assays`

Source module under the METAINFORMANT bioinformatics toolkit (`epigenome` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `atacseq.py` — ATAC-seq analysis and processing.
- `chipseq.py` — ChIP-seq analysis and processing.
- `methylation.py` — DNA methylation analysis and processing.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-epigenome-assays` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
