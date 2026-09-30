# AGENTS.md — `MetaInformAnt/src/metainformant/dna/annotation`

Source module under the METAINFORMANT bioinformatics toolkit (`dna` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `functional.py` — Functional annotation utilities for DNA variants and protein sequences.
- `gene_annotation.py` — Gene annotation, classification, and regulatory element analysis.
- `gene_finding.py` — Gene finding and ORF prediction utilities.
- `gene_prediction.py` — Gene annotation and prediction utilities.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-dna-annotation` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
