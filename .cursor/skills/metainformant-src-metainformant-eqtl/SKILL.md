---
name: metainformant-src-metainformant-eqtl
description: METAINFORMANT rules for directory src/metainformant/eqtl. Use when editing, adding tests, or reviewing code under this path. Read the linked AGENTS.md first; use uv only, write outputs to output/, real implementations.
---

# METAINFORMANT — `src/metainformant/eqtl`

Before editing files in this subtree:

- Read [`AGENTS.md`](../../../src/metainformant/eqtl/AGENTS.md) for this folder (canonical technical context).
- Optional overview: [`README.md`](../../../src/metainformant/eqtl/README.md).
- Global rules: [`CLAUDE.md`](../../../CLAUDE.md) at repo root (uv, `output/`, `.tmp/`, real implementations).
- Testing policy: [`docs/REAL_IMPLEMENTATION_POLICY.md`](../../../docs/REAL_IMPLEMENTATION_POLICY.md).
- Use `metainformant.core.io` for file I/O and `metainformant.core.utils.logging` for logs.

## Module surface (generated, validated)
- Purpose: eQTL and transcriptome-variant analysis module for METAINFORMANT.

This module provides reusable eQTL analysis components: synthetic and real
expression/genotype data preparation, transcriptome SNP-calling pipeline
methods (HISAT2 alignment + bcftools variant calling), and bcftools stats
parsing. Scripts under ``scripts/eqtl/`` are thin orchestrators that import
from this module.

Example:
    >>> from metainformant.eqtl.synthetic import create_synthetic_data
    >>> expr, geno, gene_pos, var_pos = create_synthetic_data(
    ...     n_genes=5, n_variants=25, n_samples=12
    ... )
    >>> expr.shape
    (5, 12)
- Public submodules: `pipeline`, `synthetic`, `variant_calling`, `variant_stats`.
- Canonical import: `import metainformant.eqtl` (submodules: `from metainformant import eqtl` then `eqtl.<submodule>`).
- Test entry point: `uv run pytest tests/eqtl -q` (one pytest directory per invocation).

Keep changes scoped; match existing patterns in this directory.
