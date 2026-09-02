---
name: metainformant-src-metainformant-gwas
description: METAINFORMANT rules for directory src/metainformant/gwas. Use when editing, adding tests, or reviewing code under this path. Read the linked AGENTS.md first; use uv only, write outputs to output/, real implementations.
---

# METAINFORMANT — `src/metainformant/gwas`

Before editing files in this subtree:

- Read [`AGENTS.md`](../../../src/metainformant/gwas/AGENTS.md) for this folder (canonical technical context).
- Optional overview: [`README.md`](../../../src/metainformant/gwas/README.md).
- Global rules: [`CLAUDE.md`](../../../CLAUDE.md) at repo root (uv, `output/`, `.tmp/`, real implementations).
- Testing policy: [`docs/REAL_IMPLEMENTATION_POLICY.md`](../../../docs/REAL_IMPLEMENTATION_POLICY.md).
- Use `metainformant.core.io` for file I/O and `metainformant.core.utils.logging` for logs.

## Module surface (generated, validated)
- Purpose: Genome-Wide Association Studies (GWAS) module for METAINFORMANT.

The package exposes the historical convenience API lazily. Importing
``metainformant.gwas`` keeps subpackages and common functions discoverable
without eagerly importing the full workflow and visualization stack.
- Public submodules: `analysis`, `data`, `finemapping`, `heritability`, `reporting`, `simulation`, `validation`, `visualization`, `workflow`.
- Canonical import: `import metainformant.gwas` (submodules: `from metainformant import gwas` then `gwas.<submodule>`).
- Test entry point: `uv run pytest tests/gwas -q` (one pytest directory per invocation).

Keep changes scoped; match existing patterns in this directory.
