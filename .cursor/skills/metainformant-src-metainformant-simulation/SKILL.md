---
name: metainformant-src-metainformant-simulation
description: METAINFORMANT rules for directory src/metainformant/simulation. Use when editing, adding tests, or reviewing code under this path. Read the linked AGENTS.md first; use uv only, write outputs to output/, real implementations.
---

# METAINFORMANT — `src/metainformant/simulation`

Before editing files in this subtree:

- Read [`AGENTS.md`](../../../src/metainformant/simulation/AGENTS.md) for this folder (canonical technical context).
- Optional overview: [`README.md`](../../../src/metainformant/simulation/README.md).
- Global rules: [`CLAUDE.md`](../../../CLAUDE.md) at repo root (uv, `output/`, `.tmp/`, real implementations).
- Testing policy: [`docs/REAL_IMPLEMENTATION_POLICY.md`](../../../docs/REAL_IMPLEMENTATION_POLICY.md).
- Use `metainformant.core.io` for file I/O and `metainformant.core.utils.logging` for logs.

## Module surface (generated, validated)
- Purpose: Evolutionary and Population Genetics Simulation module for METAINFORMANT.

This module provides comprehensive tools for simulating evolutionary processes,
including agent-based modeling, population genetics simulation, molecular evolution,
and sequence generation.
- Public submodules: `benchmark`, `methylation`, `models`, `workflow`.
- Canonical import: `import metainformant.simulation` (submodules: `from metainformant import simulation` then `simulation.<submodule>`).
- Test entry point: `uv run pytest tests/simulation -q` (one pytest directory per invocation).

Keep changes scoped; match existing patterns in this directory.
