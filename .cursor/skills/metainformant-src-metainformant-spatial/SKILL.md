---
name: metainformant-src-metainformant-spatial
description: METAINFORMANT rules for directory src/metainformant/spatial. Use when editing, adding tests, or reviewing code under this path. Read the linked AGENTS.md first; use uv only, write outputs to output/, real implementations.
---

# METAINFORMANT — `src/metainformant/spatial`

Before editing files in this subtree:

- Read [`AGENTS.md`](../../../src/metainformant/spatial/AGENTS.md) for this folder (canonical technical context).
- Optional overview: [`README.md`](../../../src/metainformant/spatial/README.md).
- Global rules: [`CLAUDE.md`](../../../CLAUDE.md) at repo root (uv, `output/`, `.tmp/`, real implementations).
- Testing policy: [`docs/REAL_IMPLEMENTATION_POLICY.md`](../../../docs/REAL_IMPLEMENTATION_POLICY.md).
- Use `metainformant.core.io` for file I/O and `metainformant.core.utils.logging` for logs.

## Module surface (generated, validated)
- Purpose: Spatial transcriptomics module for METAINFORMANT.
- Public submodules: `analysis`, `communication`, `deconvolution`, `integration`, `io`, `niche`, `visualization`.
- Canonical import: `import metainformant.spatial` (submodules: `from metainformant import spatial` then `spatial.<submodule>`).
- Test entry point: `uv run pytest tests/spatial -q` (one pytest directory per invocation).

Keep changes scoped; match existing patterns in this directory.
