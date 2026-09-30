# AGENTS.md — `MetaInformAnt/src/metainformant/dna/alignment`

Source module under the METAINFORMANT bioinformatics toolkit (`dna` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `distances.py` — Evolutionary distance calculations for DNA sequences.
- `msa.py` — Multiple sequence alignment for DNA sequences.
- `pairwise.py` — Sequence alignment utilities for DNA sequences.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-dna-alignment` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
