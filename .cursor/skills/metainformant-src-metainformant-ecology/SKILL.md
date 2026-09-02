---
name: metainformant-src-metainformant-ecology
description: METAINFORMANT rules for directory src/metainformant/ecology. Use when editing, adding tests, or reviewing code under this path. Read the linked AGENTS.md first; use uv only, write outputs to output/, real implementations.
---

# METAINFORMANT — `src/metainformant/ecology`

Before editing files in this subtree:

- Read [`AGENTS.md`](../../../src/metainformant/ecology/AGENTS.md) for this folder (canonical technical context).
- Optional overview: [`README.md`](../../../src/metainformant/ecology/README.md).
- Global rules: [`CLAUDE.md`](../../../CLAUDE.md) at repo root (uv, `output/`, `.tmp/`, real implementations).
- Testing policy: [`docs/REAL_IMPLEMENTATION_POLICY.md`](../../../docs/REAL_IMPLEMENTATION_POLICY.md).
- Use `metainformant.core.io` for file I/O and `metainformant.core.utils.logging` for logs.

## Module surface (generated, validated)
- Purpose: Ecology and community analysis module for METAINFORMANT.

This module provides tools for ecological analysis, including community
structure, diversity metrics, ordination, indicator species, functional
ecology, macroecology, and ecological visualization.
- Public submodules: `analysis`, `phylogenetic`, `traits`, `visualization`.
- Canonical import: `import metainformant.ecology` (submodules: `from metainformant import ecology` then `ecology.<submodule>`).
- Test entry point: `uv run pytest tests/ecology -q` (one pytest directory per invocation).

Keep changes scoped; match existing patterns in this directory.
