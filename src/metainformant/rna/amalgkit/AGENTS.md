# AGENTS.md — `MetaInformAnt/src/metainformant/rna/amalgkit`

Source module under the METAINFORMANT bioinformatics toolkit (`rna` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `__main__.py` — Run the current streaming Amalgkit producer for one configured species.
- `_amalgkit_impl.py` — Amalgkit CLI integration for RNA-seq workflow orchestration.
- `amalgkit.py` — Public facade for Amalgkit CLI integration.
- `commands.py` — Canonical command names for the supported Amalgkit release."""
- `genome_prep.py` — RNA genome and transcriptome preparation for Kallisto quantification.
- `index_prep.py` — Index Preparation and Complexity Management.
- `metadata_filter.py` — Metadata filtering utilities for RNA-seq workflows.
- `metadata_utils.py` — Utilities for Amalgkit metadata manipulation."""
- `sra_environment.py` — Campaign-local subprocess environments for SRA and FASTQ tooling."""
- `tissue_normalizer.py` — Tissue metadata normalization for RNA-seq workflows.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-rna-amalgkit` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
