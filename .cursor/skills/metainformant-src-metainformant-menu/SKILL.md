---
name: metainformant-src-metainformant-menu
description: METAINFORMANT rules for directory src/metainformant/menu. Use when editing, adding tests, or reviewing code under this path. Read the linked AGENTS.md first; use uv only, write outputs to output/, real implementations.
---

# METAINFORMANT — `src/metainformant/menu`

Before editing files in this subtree:

- Read [`AGENTS.md`](../../../src/metainformant/menu/AGENTS.md) for this folder (canonical technical context).
- Optional overview: [`README.md`](../../../src/metainformant/menu/README.md).
- Global rules: [`CLAUDE.md`](../../../CLAUDE.md) at repo root (uv, `output/`, `.tmp/`, real implementations).
- Testing policy: [`docs/REAL_IMPLEMENTATION_POLICY.md`](../../../docs/REAL_IMPLEMENTATION_POLICY.md).
- Use `metainformant.core.io` for file I/O and `metainformant.core.utils.logging` for logs.

## Module surface (generated, validated)
- Purpose: Interactive menu and CLI interface for METAINFORMANT.

This module provides the interactive command-line interface and menu system,
allowing users to discover and execute analysis workflows.
- Public submodules: `core`, `ui`.
- Canonical import: `import metainformant.menu` (submodules: `from metainformant import menu` then `menu.<submodule>`).
- Test entry point: `uv run pytest tests/menu -q` (one pytest directory per invocation).

Keep changes scoped; match existing patterns in this directory.
