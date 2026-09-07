---
name: metainformant-src-metainformant-metabolomics
description: METAINFORMANT rules for directory src/metainformant/metabolomics. Use when editing, adding tests, or reviewing code under this path. Read the linked AGENTS.md first; use uv only, write outputs to output/, real implementations.
---

# METAINFORMANT — `src/metainformant/metabolomics`

Before editing files in this subtree:

- Read [`AGENTS.md`](../../../src/metainformant/metabolomics/AGENTS.md) for this folder (canonical technical context).
- Optional overview: [`README.md`](../../../src/metainformant/metabolomics/README.md).
- Global rules: [`CLAUDE.md`](../../../CLAUDE.md) at repo root (uv, `output/`, `.tmp/`, real implementations).
- Testing policy: [`docs/REAL_IMPLEMENTATION_POLICY.md`](../../../docs/REAL_IMPLEMENTATION_POLICY.md).
- Use `metainformant.core.io` for file I/O and `metainformant.core.utils.logging` for logs.

## Module surface (generated, validated)
- Purpose: Metabolomics analysis module for METAINFORMANT.

This module provides tools for metabolite identification, mass spectrometry
data processing, pathway mapping, and metabolite-gene integration analysis.
- Public submodules: `analysis`, `io`, `pathways`, `visualization`.
- Canonical import: `import metainformant.metabolomics` (submodules: `from metainformant import metabolomics` then `metabolomics.<submodule>`).
- Test entry point: `uv run pytest tests/metabolomics -q` (one pytest directory per invocation).

Keep changes scoped; match existing patterns in this directory.
