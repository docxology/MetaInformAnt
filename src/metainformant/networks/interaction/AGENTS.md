# AGENTS.md — `MetaInformAnt/src/metainformant/networks/interaction`

Source module under the METAINFORMANT bioinformatics toolkit (`networks` domain).
Local-only path under `projects/ongoing/` (never committed).

## Layout

- `_ppi_impl.py` — Protein-protein interaction network analysis.
- `ppi.py` — Compatibility facade for protein-protein interaction analysis.
- `regulatory.py` — Regulatory network analysis and modeling.
- `regulatory_analysis.py` — Regulatory network analysis: inference, motifs, cascades, validation.
- `regulatory_core.py` — Regulatory network core: GeneRegulatoryNetwork class and network operations.

## Invariants & gotchas

- Real implementations only (no placeholders); file I/O via `metainformant.core.io`, logging via `metainformant.core.utils.logging`.
- All outputs go to `output/`, temp files to `.tmp/`; use `uv` only for package operations.
- Tests: `bash scripts/package/test.sh --pattern "<pattern>"` or `pytest tests/<domain>/ -v` from repo root.
- `.cursor/skills/metainformant-src-metainformant-networks-interaction` mirrors this folder as a Cursor skill; after moving/renaming AGENTS.md files, regenerate via `uv run python scripts/package/generate_cursor_skills.py`.
