# AGENTS.md — `MetaInformAnt/src/metainformant/visualization/genomics`

Source module under the METAINFORMANT bioinformatics toolkit (`visualization` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `expression.py` — Gene expression analysis visualization functions.
- `genomics.py` — Genomic data visualization functions for GWAS and genomic analysis.
- `networks.py` — Network visualization functions for graph analysis.
- `trees.py` — Phylogenetic tree visualization functions.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-visualization-genomics` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
