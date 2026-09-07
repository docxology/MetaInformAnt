---
name: metainformant-tests-perf
description: METAINFORMANT rules for directory tests/perf. Use when editing, adding tests, or reviewing code under this path. Read the linked AGENTS.md first; use uv only, write outputs to output/, real implementations.
---

# METAINFORMANT — `tests/perf`

Before editing files in this subtree:

- Read [`AGENTS.md`](../../../tests/perf/AGENTS.md) for this folder (canonical technical context).
- Optional overview: [`README.md`](../../../tests/perf/README.md).
- Global rules: [`CLAUDE.md`](../../../CLAUDE.md) at repo root (uv, `output/`, `.tmp/`, real implementations).
- Testing policy: [`docs/REAL_IMPLEMENTATION_POLICY.md`](../../../docs/REAL_IMPLEMENTATION_POLICY.md).
- Use `metainformant.core.io` for file I/O and `metainformant.core.utils.logging` for logs.

Keep changes scoped; match existing patterns in this directory.
