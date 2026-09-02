---
name: metainformant-src-metainformant-multiomics
description: METAINFORMANT rules for directory src/metainformant/multiomics. Use when editing, adding tests, or reviewing code under this path. Read the linked AGENTS.md first; use uv only, write outputs to output/, real implementations.
---

# METAINFORMANT — `src/metainformant/multiomics`

Before editing files in this subtree:

- Read [`AGENTS.md`](../../../src/metainformant/multiomics/AGENTS.md) for this folder (canonical technical context).
- Optional overview: [`README.md`](../../../src/metainformant/multiomics/README.md).
- Global rules: [`CLAUDE.md`](../../../CLAUDE.md) at repo root (uv, `output/`, `.tmp/`, real implementations).
- Testing policy: [`docs/REAL_IMPLEMENTATION_POLICY.md`](../../../docs/REAL_IMPLEMENTATION_POLICY.md).
- Use `metainformant.core.io` for file I/O and `metainformant.core.utils.logging` for logs.

## Module surface (generated, validated)
- Purpose: Multi-omics integration module for METAINFORMANT.

This module provides comprehensive tools for integrating multiple omics data layers,
enabling cross-platform analysis of genomics, transcriptomics, epigenomics, and
proteomics data.
- Public submodules: `analysis`, `methods`, `pathways`, `survival`, `visualization`.
- Canonical import: `import metainformant.multiomics` (submodules: `from metainformant import multiomics` then `multiomics.<submodule>`).
- Test entry point: `uv run pytest tests/multiomics -q` (one pytest directory per invocation).

Keep changes scoped; match existing patterns in this directory.
