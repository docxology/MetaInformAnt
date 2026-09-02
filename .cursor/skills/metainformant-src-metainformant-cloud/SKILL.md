---
name: metainformant-src-metainformant-cloud
description: METAINFORMANT rules for directory src/metainformant/cloud. Use when editing, adding tests, or reviewing code under this path. Read the linked AGENTS.md first; use uv only, write outputs to output/, real implementations.
---

# METAINFORMANT — `src/metainformant/cloud`

Before editing files in this subtree:

- Read [`AGENTS.md`](../../../src/metainformant/cloud/AGENTS.md) for this folder (canonical technical context).
- Optional overview: [`README.md`](../../../src/metainformant/cloud/README.md).
- Global rules: [`CLAUDE.md`](../../../CLAUDE.md) at repo root (uv, `output/`, `.tmp/`, real implementations).
- Testing policy: [`docs/REAL_IMPLEMENTATION_POLICY.md`](../../../docs/REAL_IMPLEMENTATION_POLICY.md).
- Use `metainformant.core.io` for file I/O and `metainformant.core.utils.logging` for logs.

## Module surface (generated, validated)
- Purpose: Cloud deployment infrastructure for MetaInformAnt pipelines.
- Canonical import: `import metainformant.cloud` (submodules: `from metainformant import cloud` then `cloud.<submodule>`).
- Test entry point: `uv run pytest tests/cloud -q` (one pytest directory per invocation).

Keep changes scoped; match existing patterns in this directory.
