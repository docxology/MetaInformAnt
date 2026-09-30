# AGENTS.md — `MetaInformAnt/src/metainformant/rna/splicing`

Source module under the METAINFORMANT bioinformatics toolkit (`rna` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `detection.py` — Alternative splicing detection and quantification for RNA-seq data.
- `isoforms.py` — Isoform quantification and splice graph analysis for RNA-seq data.
- `splice_analysis.py` — Splicing event classification, quantification, and differential analysis.
- `splice_sites.py` — Splice site detection and scoring for RNA-seq data.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-rna-splicing` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
