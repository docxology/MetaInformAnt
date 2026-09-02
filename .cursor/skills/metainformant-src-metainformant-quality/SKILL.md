---
name: metainformant-src-metainformant-quality
description: METAINFORMANT rules for directory src/metainformant/quality. Use when editing, adding tests, or reviewing code under this path. Read the linked AGENTS.md first; use uv only, write outputs to output/, real implementations.
---

# METAINFORMANT — `src/metainformant/quality`

Before editing files in this subtree:

- Read [`AGENTS.md`](../../../src/metainformant/quality/AGENTS.md) for this folder (canonical technical context).
- Optional overview: [`README.md`](../../../src/metainformant/quality/README.md).
- Global rules: [`CLAUDE.md`](../../../CLAUDE.md) at repo root (uv, `output/`, `.tmp/`, real implementations).
- Testing policy: [`docs/REAL_IMPLEMENTATION_POLICY.md`](../../../docs/REAL_IMPLEMENTATION_POLICY.md).
- Use `metainformant.core.io` for file I/O and `metainformant.core.utils.logging` for logs.

## Module surface (generated, validated)
- Purpose: Quality control analysis module for METAINFORMANT.

This module provides tools for assessing data quality, including
contamination detection, sequencing metrics, and FastQ file handling.
- Public submodules: `analysis`, `batch`, `io`, `reporting`.
- Canonical import: `import metainformant.quality` (submodules: `from metainformant import quality` then `quality.<submodule>`).
- Test entry point: `uv run pytest tests/quality -q` (one pytest directory per invocation).

Keep changes scoped; match existing patterns in this directory.
