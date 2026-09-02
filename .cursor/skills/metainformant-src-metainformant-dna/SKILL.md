---
name: metainformant-src-metainformant-dna
description: METAINFORMANT rules for directory src/metainformant/dna. Use when editing, adding tests, or reviewing code under this path. Read the linked AGENTS.md first; use uv only, write outputs to output/, real implementations.
---

# METAINFORMANT — `src/metainformant/dna`

Before editing files in this subtree:

- Read [`AGENTS.md`](../../../src/metainformant/dna/AGENTS.md) for this folder (canonical technical context).
- Optional overview: [`README.md`](../../../src/metainformant/dna/README.md).
- Global rules: [`CLAUDE.md`](../../../CLAUDE.md) at repo root (uv, `output/`, `.tmp/`, real implementations).
- Testing policy: [`docs/REAL_IMPLEMENTATION_POLICY.md`](../../../docs/REAL_IMPLEMENTATION_POLICY.md).
- Use `metainformant.core.io` for file I/O and `metainformant.core.utils.logging` for logs.

## Module surface (generated, validated)
- Purpose: DNA sequence analysis and genomics module for METAINFORMANT.
- Public submodules: `alignment`, `annotation`, `expression`, `external`, `integration`, `io`, `phylogeny`, `population`, `sequence`, `transcription`, `translation`, `variation`.
- Canonical import: `import metainformant.dna` (submodules: `from metainformant import dna` then `dna.<submodule>`).
- Test entry point: `uv run pytest tests/dna -q` (one pytest directory per invocation).

Keep changes scoped; match existing patterns in this directory.
