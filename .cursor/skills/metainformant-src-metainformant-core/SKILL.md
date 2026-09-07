---
name: metainformant-src-metainformant-core
description: METAINFORMANT rules for directory src/metainformant/core. Use when editing, adding tests, or reviewing code under this path. Read the linked AGENTS.md first; use uv only, write outputs to output/, real implementations.
---

# METAINFORMANT — `src/metainformant/core`

Before editing files in this subtree:

- Read [`AGENTS.md`](../../../src/metainformant/core/AGENTS.md) for this folder (canonical technical context).
- Optional overview: [`README.md`](../../../src/metainformant/core/README.md).
- Global rules: [`CLAUDE.md`](../../../CLAUDE.md) at repo root (uv, `output/`, `.tmp/`, real implementations).
- Testing policy: [`docs/REAL_IMPLEMENTATION_POLICY.md`](../../../docs/REAL_IMPLEMENTATION_POLICY.md).
- Use `metainformant.core.io` for file I/O and `metainformant.core.utils.logging` for logs.

## Module surface (generated, validated)
- Purpose: Core utilities for METAINFORMANT bioinformatics toolkit.
- Public submodules: `data`, `db`, `engine`, `execution`, `io`, `ncbi`, `ui`, `utils`.
- Canonical import: `import metainformant.core` (submodules: `from metainformant import core` then `core.<submodule>`).
- Test entry point: `uv run pytest tests/core -q` (one pytest directory per invocation).

Keep changes scoped; match existing patterns in this directory.
