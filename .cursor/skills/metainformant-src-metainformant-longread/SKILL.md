---
name: metainformant-src-metainformant-longread
description: METAINFORMANT rules for directory src/metainformant/longread. Use when editing, adding tests, or reviewing code under this path. Read the linked AGENTS.md first; use uv only, write outputs to output/, real implementations.
---

# METAINFORMANT — `src/metainformant/longread`

Before editing files in this subtree:

- Read [`AGENTS.md`](../../../src/metainformant/longread/AGENTS.md) for this folder (canonical technical context).
- Optional overview: [`README.md`](../../../src/metainformant/longread/README.md).
- Global rules: [`CLAUDE.md`](../../../CLAUDE.md) at repo root (uv, `output/`, `.tmp/`, real implementations).
- Testing policy: [`docs/REAL_IMPLEMENTATION_POLICY.md`](../../../docs/REAL_IMPLEMENTATION_POLICY.md).
- Use `metainformant.core.io` for file I/O and `metainformant.core.utils.logging` for logs.

## Module surface (generated, validated)
- Purpose: Long-read sequencing analysis module for METAINFORMANT.
- Public submodules: `analysis`, `assembly`, `io`, `methylation`, `phasing`, `quality`, `utils`, `visualization`, `workflow`.
- Canonical import: `import metainformant.longread` (submodules: `from metainformant import longread` then `longread.<submodule>`).
- Test entry point: `uv run pytest tests/longread -q` (one pytest directory per invocation).

Keep changes scoped; match existing patterns in this directory.
