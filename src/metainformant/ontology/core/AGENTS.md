# AGENTS.md — `MetaInformAnt/src/metainformant/ontology/core`

Source module under the METAINFORMANT bioinformatics toolkit (`ontology` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `go.py` — Gene Ontology analysis and enrichment.
- `go_api.py` — QuickGO REST API client for live Gene Ontology annotation.
- `hpo.py` — Human Phenotype Ontology (HPO) client and GWAS trait mapper.
- `obo.py` — OBO format parsing and processing.
- `types.py` — Ontology data types and structures.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-ontology-core` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
