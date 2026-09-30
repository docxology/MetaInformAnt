# AGENTS.md — `MetaInformAnt/src/metainformant/metagenomics/shotgun`

Source module under the METAINFORMANT bioinformatics toolkit (`metagenomics` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `assembly.py` — Metagenome assembly from shotgun sequencing data.
- `binning.py` — Metagenomic binning of assembled contigs.
- `profiling.py` — Community profiling from shotgun metagenomic reads.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-metagenomics-shotgun` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
