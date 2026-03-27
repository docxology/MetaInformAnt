---
name: metainformant-scripts-gwas-pipelines
description: METAINFORMANT rules for directory scripts/gwas/pipelines. Use when editing, adding tests, or reviewing code under this path. Read the linked AGENTS.md first; use uv only, write outputs to output/, no mocks.
---

# METAINFORMANT — `scripts/gwas/pipelines`

Before editing files in this subtree:

- Read [`AGENTS.md`](../../../scripts/gwas/pipelines/AGENTS.md) for this folder (canonical technical context).
- Optional overview: [`README.md`](../../../scripts/gwas/pipelines/README.md).
- Global rules: [`CLAUDE.md`](../../../CLAUDE.md) at repo root (uv, `output/`, `.tmp/`, no mocks).
- Testing policy: [`docs/NO_MOCKING_POLICY.md`](../../../docs/NO_MOCKING_POLICY.md).
- Use `metainformant.core.io` for file I/O and `metainformant.core.utils.logging` for logs.

Keep changes scoped; match existing patterns in this directory.
