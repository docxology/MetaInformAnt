---
name: metainformant-src-metainformant-ontology
description: METAINFORMANT rules for directory src/metainformant/ontology. Use when editing, adding tests, or reviewing code under this path. Read the linked AGENTS.md first; use uv only, write outputs to output/, real implementations.
---

# METAINFORMANT — `src/metainformant/ontology`

Before editing files in this subtree:

- Read [`AGENTS.md`](../../../src/metainformant/ontology/AGENTS.md) for this folder (canonical technical context).
- Optional overview: [`README.md`](../../../src/metainformant/ontology/README.md).
- Global rules: [`CLAUDE.md`](../../../CLAUDE.md) at repo root (uv, `output/`, `.tmp/`, real implementations).
- Testing policy: [`docs/REAL_IMPLEMENTATION_POLICY.md`](../../../docs/REAL_IMPLEMENTATION_POLICY.md).
- Use `metainformant.core.io` for file I/O and `metainformant.core.utils.logging` for logs.

## Module surface (generated, validated)
- Purpose: Gene ontology and functional annotation module for METAINFORMANT.
- Public submodules: `annotation`, `core`, `pathway_enrichment`, `query`, `visualization`.
- Canonical import: `import metainformant.ontology` (submodules: `from metainformant import ontology` then `ontology.<submodule>`).
- Test entry point: `uv run pytest tests/ontology -q` (one pytest directory per invocation).

Keep changes scoped; match existing patterns in this directory.
