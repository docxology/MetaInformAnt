# AGENTS.md — `MetaInformAnt/src/metainformant/rna/engine`

Source module under the METAINFORMANT bioinformatics toolkit (`rna` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `discovery.py` — RNA species discovery and genome configuration generation.
- `pipeline.py` — RNA-seq pipeline utilities and high-level workflow orchestration.
- `progress_dashboard.py` — Pipeline progress dashboard — mosaic graphical abstract.
- `progress_db.py` — SQLite-backed progress tracking for the current RNA-seq pipeline.
- `provenance.py` — Current-method provenance for per-sample RNA quantification outputs.
- `raw_cleanup.py` — Safe per-sample reclamation of transient RNA-seq inputs.
- `species.py` — Shared species/configuration discovery helpers for RNA workflows.
- `sra_extraction.py` — SRA file extraction and fallback recovery utilities.
- `streaming_orchestrator.py` — Streaming RNA-seq Orchestrator implementation.
- `workflow.py` — RNA-seq workflow orchestration and configuration management.
- `workflow_cleanup.py` — Workflow cleanup and disk management utilities.
- `workflow_core.py` — RNA-seq workflow configuration and data classes.
- `workflow_execution.py` — RNA-seq workflow execution engine.
- `workflow_planning.py` — RNA-seq workflow planning, step defaults, and pre-execution preparation.
- `workflow_steps.py` — Workflow step helpers - prerequisite validation, post-step actions, and setup.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-rna-engine` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
