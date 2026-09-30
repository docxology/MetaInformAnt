# AGENTS.md — `MetaInformAnt/src/metainformant/gwas/workflow`

Source module under the METAINFORMANT bioinformatics toolkit (`gwas` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `rna_quantification.py` — GWAS-facing adapter for RNA quantification discovery.
- `workflow.py` — GWAS workflow orchestration and configuration management.
- `workflow_config.py` — GWAS workflow configuration and data loading utilities.
- `workflow_execution.py` — GWAS workflow execution functions.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-gwas-workflow` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
