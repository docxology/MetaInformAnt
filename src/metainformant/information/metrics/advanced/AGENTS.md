# AGENTS.md — `MetaInformAnt/src/metainformant/information/metrics/advanced`

Source module under the METAINFORMANT bioinformatics toolkit (`information` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `channel.py` — Channel capacity and rate-distortion functions for information theory.
- `decomposition.py` — Partial Information Decomposition (PID) measures for biological data.
- `fisher_rao.py` — Fisher-Rao distance, natural gradient, and related information geometry measures.
- `geometry.py` — Information geometry measures for statistical manifolds.
- `hypothesis.py` — Information-theoretic hypothesis testing for biological data.
- `information_projection.py` — Information projection, divergences, channel capacity, and related measures.
- `semantic.py` — Semantic information theory measures for biological ontologies.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-information-metrics-advanced` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
