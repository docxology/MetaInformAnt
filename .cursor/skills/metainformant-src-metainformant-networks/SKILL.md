---
name: metainformant-src-metainformant-networks
description: METAINFORMANT rules for directory src/metainformant/networks. Use when editing, adding tests, or reviewing code under this path. Read the linked AGENTS.md first; use uv only, write outputs to output/, real implementations.
---

# METAINFORMANT — `src/metainformant/networks`

Before editing files in this subtree:

- Read [`AGENTS.md`](../../../src/metainformant/networks/AGENTS.md) for this folder (canonical technical context).
- Optional overview: [`README.md`](../../../src/metainformant/networks/README.md).
- Global rules: [`CLAUDE.md`](../../../CLAUDE.md) at repo root (uv, `output/`, `.tmp/`, real implementations).
- Testing policy: [`docs/REAL_IMPLEMENTATION_POLICY.md`](../../../docs/REAL_IMPLEMENTATION_POLICY.md).
- Use `metainformant.core.io` for file I/O and `metainformant.core.utils.logging` for logs.

## Module surface (generated, validated)
- Purpose: Network analysis module for METAINFORMANT.
- Public submodules: `analysis`, `config`, `interaction`, `regulatory`, `workflow`.
- Canonical import: `import metainformant.networks` (submodules: `from metainformant import networks` then `networks.<submodule>`).
- Test entry point: `uv run pytest tests/networks -q` (one pytest directory per invocation).

Keep changes scoped; match existing patterns in this directory.
