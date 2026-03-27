---
name: metainformant-scripts-gwas-structure
description: METAINFORMANT rules for directory scripts/gwas/structure. Use when editing, adding tests, or reviewing code under this path. Read the linked AGENTS.md first; use uv only, write outputs to output/, no mocks.
---

# METAINFORMANT — `scripts/gwas/structure`

Before editing files in this subtree:

- Read [`AGENTS.md`](../../../scripts/gwas/structure/AGENTS.md) for this folder (canonical technical context).
- Optional overview: [`README.md`](../../../scripts/gwas/structure/README.md).
- Global rules: [`CLAUDE.md`](../../../CLAUDE.md) at repo root (uv, `output/`, `.tmp/`, no mocks).
- Testing policy: [`docs/NO_MOCKING_POLICY.md`](../../../docs/NO_MOCKING_POLICY.md).
- Use `metainformant.core.io` for file I/O and `metainformant.core.utils.logging` for logs.

Keep changes scoped; match existing patterns in this directory.
