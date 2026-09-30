# AGENTS.md — `MetaInformAnt/src/metainformant/life_events/models`

Source module under the METAINFORMANT bioinformatics toolkit (`life_events` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `embeddings.py` — Event embedding and representation learning for life events analysis.
- `predictor.py` — Embedding-based prediction models for life event sequences.
- `sequence_models.py` — Neural sequence prediction models for life events.
- `statistical_models.py` — Statistical and multi-task prediction models for life events.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-life_events-models` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
