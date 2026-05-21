---
name: metainformant-src-metainformant-mcp
description: METAINFORMANT rules for directory src/metainformant/mcp. Use when editing, adding tests, or reviewing code under this path. Read the linked AGENTS.md first; use uv only, write outputs to output/, no mocks.
---

# METAINFORMANT — `src/metainformant/mcp`

Before editing files in this subtree:

- Read [`AGENTS.md`](../../../src/metainformant/mcp/AGENTS.md) for this folder (canonical technical context).
- Optional overview: [`README.md`](../../../src/metainformant/mcp/README.md).
- Global rules: [`CLAUDE.md`](../../../CLAUDE.md) at repo root (uv, `output/`, `.tmp/`, no mocks).
- Testing policy: [`docs/NO_MOCKING_POLICY.md`](../../../docs/NO_MOCKING_POLICY.md).
- Use `metainformant.core.io` for file I/O and `metainformant.core.utils.logging` for logs.

Keep changes scoped; match existing patterns in this directory.
