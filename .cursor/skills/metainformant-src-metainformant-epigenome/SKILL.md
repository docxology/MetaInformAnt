---
name: metainformant-src-metainformant-epigenome
description: METAINFORMANT rules for directory src/metainformant/epigenome. Use when editing, adding tests, or reviewing code under this path. Read the linked AGENTS.md first; use uv only, write outputs to output/, real implementations.
---

# METAINFORMANT — `src/metainformant/epigenome`

Before editing files in this subtree:

- Read [`AGENTS.md`](../../../src/metainformant/epigenome/AGENTS.md) for this folder (canonical technical context).
- Optional overview: [`README.md`](../../../src/metainformant/epigenome/README.md).
- Global rules: [`CLAUDE.md`](../../../CLAUDE.md) at repo root (uv, `output/`, `.tmp/`, real implementations).
- Testing policy: [`docs/REAL_IMPLEMENTATION_POLICY.md`](../../../docs/REAL_IMPLEMENTATION_POLICY.md).
- Use `metainformant.core.io` for file I/O and `metainformant.core.utils.logging` for logs.

## Module surface (generated, validated)
- Purpose: Epigenome analysis module for METAINFORMANT.
- Public submodules: `analysis`, `assays`, `chromatin_state`, `peak_calling`, `visualization`, `workflow`.
- Canonical import: `import metainformant.epigenome` (submodules: `from metainformant import epigenome` then `epigenome.<submodule>`).
- Test entry point: `uv run pytest tests/epigenome -q` (one pytest directory per invocation).

Keep changes scoped; match existing patterns in this directory.
