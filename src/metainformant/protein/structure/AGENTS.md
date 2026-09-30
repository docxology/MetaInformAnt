# AGENTS.md — `MetaInformAnt/src/metainformant/protein/structure`

Source module under the METAINFORMANT bioinformatics toolkit (`protein` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `alphafold.py` — AlphaFold protein structure prediction integration.
- `analysis.py` — Protein structure analysis utilities.
- `contacts.py` — Protein contact analysis utilities.
- `general.py` — Protein structure analysis utilities.
- `io.py` — Protein structure file I/O utilities.
- `pdb.py` — PDB file parsing and manipulation utilities.
- `secondary.py` — Protein secondary structure prediction utilities.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-protein-structure` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
