# AGENTS.md — `MetaInformAnt/src/metainformant/rna/analysis`

Source module under the METAINFORMANT bioinformatics toolkit (`rna` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `across_species_orchestrator.py` — Across-species downstream analysis orchestrator.
- `cross_species.py` — Cross-species gene expression comparison module.
- `expression.py` — Differential expression analysis module for RNA-seq data.
- `expression_analysis.py` — Differential expression analysis, PCA, and visualization data preparation.
- `expression_core.py` — Count matrix normalization, size factor estimation, and gene filtering.
- `ortholog_mapping.py` — Ortholog mapping integration module.
- `protein_integration.py` — RNA-Protein integration and translation efficiency analysis.
- `qc.py` — RNA-seq quality control module.
- `qc_filtering.py` — RNA-seq quality control filtering, bias detection, and report generation.
- `qc_metrics.py` — RNA-seq quality control metric computation.
- `validation.py` — RNA-seq workflow sample validation utilities.
- `within_species_orchestrator.py` — Within-species downstream analysis orchestrator.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-rna-analysis` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
