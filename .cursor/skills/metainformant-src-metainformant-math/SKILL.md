---
name: metainformant-src-metainformant-math
description: METAINFORMANT rules for directory src/metainformant/math. Use when editing, adding tests, or reviewing code under this path. Read the linked AGENTS.md first; use uv only, write outputs to output/, real implementations.
---

# METAINFORMANT — `src/metainformant/math`

Before editing files in this subtree:

- Read [`AGENTS.md`](../../../src/metainformant/math/AGENTS.md) for this folder (canonical technical context).
- Optional overview: [`README.md`](../../../src/metainformant/math/README.md).
- Global rules: [`CLAUDE.md`](../../../CLAUDE.md) at repo root (uv, `output/`, `.tmp/`, real implementations).
- Testing policy: [`docs/REAL_IMPLEMENTATION_POLICY.md`](../../../docs/REAL_IMPLEMENTATION_POLICY.md).
- Use `metainformant.core.io` for file I/O and `metainformant.core.utils.logging` for logs.

## Module surface (generated, validated)
- Purpose: Mathematical biology and theoretical modeling module for METAINFORMANT.
- Public submodules: `bayesian`, `core`, `decision_theory`, `epidemiology`, `evolutionary_dynamics`, `popgen`, `population_genetics`, `quantitative_genetics`.
- Canonical import: `import metainformant.math` (submodules: `from metainformant import math` then `math.<submodule>`).
- Test entry point: `uv run pytest tests/math -q` (one pytest directory per invocation).

Keep changes scoped; match existing patterns in this directory.
