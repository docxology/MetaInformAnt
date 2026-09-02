---
name: metainformant-src-metainformant-protein
description: METAINFORMANT rules for directory src/metainformant/protein. Use when editing, adding tests, or reviewing code under this path. Read the linked AGENTS.md first; use uv only, write outputs to output/, real implementations.
---

# METAINFORMANT — `src/metainformant/protein`

Before editing files in this subtree:

- Read [`AGENTS.md`](../../../src/metainformant/protein/AGENTS.md) for this folder (canonical technical context).
- Optional overview: [`README.md`](../../../src/metainformant/protein/README.md).
- Global rules: [`CLAUDE.md`](../../../CLAUDE.md) at repo root (uv, `output/`, `.tmp/`, real implementations).
- Testing policy: [`docs/REAL_IMPLEMENTATION_POLICY.md`](../../../docs/REAL_IMPLEMENTATION_POLICY.md).
- Use `metainformant.core.io` for file I/O and `metainformant.core.utils.logging` for logs.

## Module surface (generated, validated)
- Purpose: Protein sequence and structure analysis module for METAINFORMANT.
- Public submodules: `database`, `domains`, `function`, `sequence`, `structure`, `visualization`, `workflow`.
- Canonical import: `import metainformant.protein` (submodules: `from metainformant import protein` then `protein.<submodule>`).
- Test entry point: `uv run pytest tests/protein -q` (one pytest directory per invocation).

Keep changes scoped; match existing patterns in this directory.
