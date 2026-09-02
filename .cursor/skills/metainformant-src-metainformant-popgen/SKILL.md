---
name: metainformant-src-metainformant-popgen
description: METAINFORMANT rules for directory src/metainformant/popgen. Use when editing, adding tests, or reviewing code under this path. Read the linked AGENTS.md first; use uv only, write outputs to output/, real implementations.
---

# METAINFORMANT — `src/metainformant/popgen`

Before editing files in this subtree:

- Read [`AGENTS.md`](../../../src/metainformant/popgen/AGENTS.md) for this folder (canonical technical context).
- Optional overview: [`README.md`](../../../src/metainformant/popgen/README.md).
- Global rules: [`CLAUDE.md`](../../../CLAUDE.md) at repo root (uv, `output/`, `.tmp/`, real implementations).
- Testing policy: [`docs/REAL_IMPLEMENTATION_POLICY.md`](../../../docs/REAL_IMPLEMENTATION_POLICY.md).
- Use `metainformant.core.io` for file I/O and `metainformant.core.utils.logging` for logs.

## Module surface (generated, validated)
- Purpose: Population genetics analysis module for METAINFORMANT.

Consolidates the population-genetics analysis surface: sequence-based
summary statistics, neutrality tests, two-population comparison, PCA /
kinship / HWE on genotype matrices, LD summaries, and demographic model
comparisons. The scenario analysis pipeline (as used by the synthetic
dataset workflow) is exposed via :func:`analyze_dataset`, with
``scripts/popgen/`` as its thin orchestrator.

Example:
    >>> from metainformant.popgen.analysis import summarize_scenario
    >>> result = summarize_scenario(
    ...     ["ATCG", "ATCG", "GCTA"], label="demo"
    ... )  # doctest: +SKIP
- Public submodules: `analysis`.
- Canonical import: `import metainformant.popgen` (submodules: `from metainformant import popgen` then `popgen.<submodule>`).
- Test entry point: `uv run pytest tests/popgen -q` (one pytest directory per invocation).

Keep changes scoped; match existing patterns in this directory.
