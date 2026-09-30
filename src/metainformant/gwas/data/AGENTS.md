# AGENTS.md — `MetaInformAnt/src/metainformant/gwas/data`

Source module under the METAINFORMANT bioinformatics toolkit (`gwas` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `config.py` — Configuration management for GWAS workflows.
- `download.py` — Data download utilities for GWAS analysis.
- `expression.py` — Expression Data Loading Module.
- `extraction.py` — Data extraction utilities for GWAS pipeline."""
- `genome.py` — Apis mellifera (Amel_HAv3.1) genome mapping and annotation utilities.
- `metadata.py` — Sample metadata loading, validation, and merging for GWAS analysis.
- `sra_download.py` — SRA data download utilities for GWAS.
- `traits.py` — Phenotype/trait data loading utilities for GWAS analysis.
- `vcf_utils.py` — VCF file utilities for GWAS pipeline operations.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-gwas-data` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
