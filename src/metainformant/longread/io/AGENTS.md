# AGENTS.md — `MetaInformAnt/src/metainformant/longread/io`

Source module under the METAINFORMANT bioinformatics toolkit (`longread` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `bam.py` — Long-read BAM file processing with methylation and supplementary alignment support.
- `fast5.py` — FAST5/POD5 file reading for Oxford Nanopore signal data.
- `formats.py` — Format conversion utilities for long-read sequencing data.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-longread-io` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
