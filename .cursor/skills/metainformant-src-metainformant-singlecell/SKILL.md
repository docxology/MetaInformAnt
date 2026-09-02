---
name: metainformant-src-metainformant-singlecell
description: METAINFORMANT rules for directory src/metainformant/singlecell. Use when editing, adding tests, or reviewing code under this path. Read the linked AGENTS.md first; use uv only, write outputs to output/, real implementations.
---

# METAINFORMANT — `src/metainformant/singlecell`

Before editing files in this subtree:

- Read [`AGENTS.md`](../../../src/metainformant/singlecell/AGENTS.md) for this folder (canonical technical context).
- Optional overview: [`README.md`](../../../src/metainformant/singlecell/README.md).
- Global rules: [`CLAUDE.md`](../../../CLAUDE.md) at repo root (uv, `output/`, `.tmp/`, real implementations).
- Testing policy: [`docs/REAL_IMPLEMENTATION_POLICY.md`](../../../docs/REAL_IMPLEMENTATION_POLICY.md).
- Use `metainformant.core.io` for file I/O and `metainformant.core.utils.logging` for logs.

## Module surface (generated, validated)
- Purpose: Single-cell analysis module for METAINFORMANT.

This module provides tools for single-cell data analysis, including
dimensionality reduction, clustering, trajectory inference, multi-sample
integration, cell type annotation, differential expression, and RNA
velocity estimation.
- Public submodules: `analysis`, `celltyping`, `data`, `differential`, `doublet`, `io`, `velocity`, `visualization`.
- Canonical import: `import metainformant.singlecell` (submodules: `from metainformant import singlecell` then `singlecell.<submodule>`).
- Test entry point: `uv run pytest tests/singlecell -q` (one pytest directory per invocation).

Keep changes scoped; match existing patterns in this directory.
