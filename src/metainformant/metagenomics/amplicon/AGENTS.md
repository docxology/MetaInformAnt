# AGENTS.md — `MetaInformAnt/src/metainformant/metagenomics/amplicon`

Source module under the METAINFORMANT bioinformatics toolkit (`metagenomics` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `asv_denoising.py` — ASV denoising for amplicon metagenomics.
- `otu_clustering.py` — OTU clustering for amplicon metagenomics.
- `taxonomy.py` — Taxonomic classification for amplicon sequences.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-metagenomics-amplicon` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
