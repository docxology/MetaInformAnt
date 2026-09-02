---
name: metainformant-src-metainformant-visualization
description: METAINFORMANT rules for directory src/metainformant/visualization. Use when editing, adding tests, or reviewing code under this path. Read the linked AGENTS.md first; use uv only, write outputs to output/, real implementations.
---

# METAINFORMANT — `src/metainformant/visualization`

Before editing files in this subtree:

- Read [`AGENTS.md`](../../../src/metainformant/visualization/AGENTS.md) for this folder (canonical technical context).
- Optional overview: [`README.md`](../../../src/metainformant/visualization/README.md).
- Global rules: [`CLAUDE.md`](../../../CLAUDE.md) at repo root (uv, `output/`, `.tmp/`, real implementations).
- Testing policy: [`docs/REAL_IMPLEMENTATION_POLICY.md`](../../../docs/REAL_IMPLEMENTATION_POLICY.md).
- Use `metainformant.core.io` for file I/O and `metainformant.core.utils.logging` for logs.

## Module surface (generated, validated)
- Purpose: Visualization and plotting utilities module for METAINFORMANT.
- Public submodules: `analysis`, `config`, `dashboards`, `genomics`, `interactive_dashboards`, `plots`.
- Canonical import: `import metainformant.visualization` (submodules: `from metainformant import visualization` then `visualization.<submodule>`).
- Test entry point: `uv run pytest tests/visualization -q` (one pytest directory per invocation).

Keep changes scoped; match existing patterns in this directory.
