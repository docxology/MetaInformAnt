# AGENTS.md — `MetaInformAnt/src/metainformant/networks/analysis`

Source module under the METAINFORMANT bioinformatics toolkit (`networks` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `community.py` — Community detection algorithms for biological networks.
- `graph.py` — Graph construction and manipulation utilities for METAINFORMANT.
- `graph_algorithms.py` — Graph algorithms: metrics, similarity, filtering, centrality, paths.
- `graph_core.py` — Graph core: BiologicalNetwork class and IO/construction utilities.
- `pathway.py` — Pathway analysis and enrichment for biological networks.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-networks-analysis` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
