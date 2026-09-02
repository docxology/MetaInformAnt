---
name: metainformant-src-metainformant-information
description: METAINFORMANT rules for directory src/metainformant/information. Use when editing, adding tests, or reviewing code under this path. Read the linked AGENTS.md first; use uv only, write outputs to output/, real implementations.
---

# METAINFORMANT — `src/metainformant/information`

Before editing files in this subtree:

- Read [`AGENTS.md`](../../../src/metainformant/information/AGENTS.md) for this folder (canonical technical context).
- Optional overview: [`README.md`](../../../src/metainformant/information/README.md).
- Global rules: [`CLAUDE.md`](../../../CLAUDE.md) at repo root (uv, `output/`, `.tmp/`, real implementations).
- Testing policy: [`docs/REAL_IMPLEMENTATION_POLICY.md`](../../../docs/REAL_IMPLEMENTATION_POLICY.md).
- Use `metainformant.core.io` for file I/O and `metainformant.core.utils.logging` for logs.

## Module surface (generated, validated)
- Purpose: Information theory analysis module for METAINFORMANT.
- Public submodules: `integration`, `metrics`, `network_info`, `workflow`.
- Canonical import: `import metainformant.information` (submodules: `from metainformant import information` then `information.<submodule>`).
- Test entry point: `uv run pytest tests/information -q` (one pytest directory per invocation).

Keep changes scoped; match existing patterns in this directory.
