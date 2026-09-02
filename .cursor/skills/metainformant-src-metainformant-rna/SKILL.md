---
name: metainformant-src-metainformant-rna
description: METAINFORMANT rules for directory src/metainformant/rna. Use when editing, adding tests, or reviewing code under this path. Read the linked AGENTS.md first; use uv only, write outputs to output/, real implementations.
---

# METAINFORMANT — `src/metainformant/rna`

Before editing files in this subtree:

- Read [`AGENTS.md`](../../../src/metainformant/rna/AGENTS.md) for this folder (canonical technical context).
- Optional overview: [`README.md`](../../../src/metainformant/rna/README.md).
- Global rules: [`CLAUDE.md`](../../../CLAUDE.md) at repo root (uv, `output/`, `.tmp/`, real implementations).
- Testing policy: [`docs/REAL_IMPLEMENTATION_POLICY.md`](../../../docs/REAL_IMPLEMENTATION_POLICY.md).
- Use `metainformant.core.io` for file I/O and `metainformant.core.utils.logging` for logs.

## Module surface (generated, validated)
- Purpose: RNA transcriptomic analysis and workflow orchestration module for METAINFORMANT.
- Public submodules: `amalgkit`, `analysis`, `core`, `deconvolution`, `engine`, `retrieval`, `splicing`, `steps`.
- Canonical import: `import metainformant.rna` (submodules: `from metainformant import rna` then `rna.<submodule>`).
- Test entry point: `uv run pytest tests/rna -q` (one pytest directory per invocation).

Keep changes scoped; match existing patterns in this directory.
