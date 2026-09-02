---
name: metainformant-src-metainformant-structural-variants
description: METAINFORMANT rules for directory src/metainformant/structural_variants. Use when editing, adding tests, or reviewing code under this path. Read the linked AGENTS.md first; use uv only, write outputs to output/, real implementations.
---

# METAINFORMANT — `src/metainformant/structural_variants`

Before editing files in this subtree:

- Read [`AGENTS.md`](../../../src/metainformant/structural_variants/AGENTS.md) for this folder (canonical technical context).
- Optional overview: [`README.md`](../../../src/metainformant/structural_variants/README.md).
- Global rules: [`CLAUDE.md`](../../../CLAUDE.md) at repo root (uv, `output/`, `.tmp/`, real implementations).
- Testing policy: [`docs/REAL_IMPLEMENTATION_POLICY.md`](../../../docs/REAL_IMPLEMENTATION_POLICY.md).
- Use `metainformant.core.io` for file I/O and `metainformant.core.utils.logging` for logs.

## Module surface (generated, validated)
- Purpose: Structural Variant (SV) analysis module for METAINFORMANT.

This module provides comprehensive tools for detecting, annotating, filtering,
and visualizing structural variants including copy number variations (CNVs),
deletions, duplications, inversions, translocations, and insertions.

Config prefix: SV_
- Public submodules: `annotation`, `detection`, `filtering`, `population`, `visualization`.
- Canonical import: `import metainformant.structural_variants` (submodules: `from metainformant import structural_variants` then `structural_variants.<submodule>`).
- Test entry point: `uv run pytest tests/structural_variants -q` (one pytest directory per invocation).

Keep changes scoped; match existing patterns in this directory.
