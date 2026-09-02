---
name: metainformant-src-metainformant-metagenomics
description: METAINFORMANT rules for directory src/metainformant/metagenomics. Use when editing, adding tests, or reviewing code under this path. Read the linked AGENTS.md first; use uv only, write outputs to output/, real implementations.
---

# METAINFORMANT — `src/metainformant/metagenomics`

Before editing files in this subtree:

- Read [`AGENTS.md`](../../../src/metainformant/metagenomics/AGENTS.md) for this folder (canonical technical context).
- Optional overview: [`README.md`](../../../src/metainformant/metagenomics/README.md).
- Global rules: [`CLAUDE.md`](../../../CLAUDE.md) at repo root (uv, `output/`, `.tmp/`, real implementations).
- Testing policy: [`docs/REAL_IMPLEMENTATION_POLICY.md`](../../../docs/REAL_IMPLEMENTATION_POLICY.md).
- Use `metainformant.core.io` for file I/O and `metainformant.core.utils.logging` for logs.

## Module surface (generated, validated)
- Purpose: Metagenomics analysis module for METAINFORMANT.

This module provides comprehensive tools for metagenomic analysis, including
amplicon-based community profiling (16S/ITS), shotgun metagenomics (assembly,
binning, profiling), functional annotation (gene prediction, pathway reconstruction),
specialized visualization for microbial ecology data, community diversity metrics,
and comparative/differential abundance analysis.
- Public submodules: `amplicon`, `comparative`, `diversity`, `functional`, `shotgun`, `visualization`.
- Canonical import: `import metainformant.metagenomics` (submodules: `from metainformant import metagenomics` then `metagenomics.<submodule>`).
- Test entry point: `uv run pytest tests/metagenomics -q` (one pytest directory per invocation).

Keep changes scoped; match existing patterns in this directory.
