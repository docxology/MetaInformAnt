---
name: metainformant-src-metainformant-phenotype
description: METAINFORMANT rules for directory src/metainformant/phenotype. Use when editing, adding tests, or reviewing code under this path. Read the linked AGENTS.md first; use uv only, write outputs to output/, real implementations.
---

# METAINFORMANT — `src/metainformant/phenotype`

Before editing files in this subtree:

- Read [`AGENTS.md`](../../../src/metainformant/phenotype/AGENTS.md) for this folder (canonical technical context).
- Optional overview: [`README.md`](../../../src/metainformant/phenotype/README.md).
- Global rules: [`CLAUDE.md`](../../../CLAUDE.md) at repo root (uv, `output/`, `.tmp/`, real implementations).
- Testing policy: [`docs/REAL_IMPLEMENTATION_POLICY.md`](../../../docs/REAL_IMPLEMENTATION_POLICY.md).
- Use `metainformant.core.io` for file I/O and `metainformant.core.utils.logging` for logs.

## Module surface (generated, validated)
- Purpose: Phenotype module for MetaInformAnt.

Comprehensive phenotypic trait analysis including:
- Morphological measurements and allometric analysis
- Behavioral sequences, ethograms, and diversity metrics
- Chemical profiles (GC-MS, CHC) with distance and marker detection
- Acoustic signals with spectral and temporal analysis
- Electronic tracking (RFID, video, GPS) with movement ecology
- Life course trajectory analysis
- Cross-omic integration (phenotype-genotype, trait-expression)
- Configurable analysis pipelines
- Public submodules: `analysis`, `behavior`, `chemical`, `data`, `electronic`, `gwas_integration`, `integration`, `morphological`, `sonic`, `visualization`, `workflow`.
- Canonical import: `import metainformant.phenotype` (submodules: `from metainformant import phenotype` then `phenotype.<submodule>`).
- Test entry point: `uv run pytest tests/phenotype -q` (one pytest directory per invocation).

Keep changes scoped; match existing patterns in this directory.
