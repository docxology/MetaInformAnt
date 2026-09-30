# AGENTS.md — `MetaInformAnt/src/metainformant/dna/variation`

Source module under the METAINFORMANT bioinformatics toolkit (`dna` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `calling.py` — Variant calling utilities for DNA sequence analysis.
- `mutations.py` — DNA mutation detection and analysis utilities.
- `variants.py` — DNA variant analysis and VCF file processing.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-dna-variation` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
