---
name: metainformant-src-metainformant-life-events
description: METAINFORMANT rules for directory src/metainformant/life_events. Use when editing, adding tests, or reviewing code under this path. Read the linked AGENTS.md first; use uv only, write outputs to output/, real implementations.
---

# METAINFORMANT — `src/metainformant/life_events`

Before editing files in this subtree:

- Read [`AGENTS.md`](../../../src/metainformant/life_events/AGENTS.md) for this folder (canonical technical context).
- Optional overview: [`README.md`](../../../src/metainformant/life_events/README.md).
- Global rules: [`CLAUDE.md`](../../../CLAUDE.md) at repo root (uv, `output/`, `.tmp/`, real implementations).
- Testing policy: [`docs/REAL_IMPLEMENTATION_POLICY.md`](../../../docs/REAL_IMPLEMENTATION_POLICY.md).
- Use `metainformant.core.io` for file I/O and `metainformant.core.utils.logging` for logs.

## Module surface (generated, validated)
- Purpose: Life events and trajectory analysis module for METAINFORMANT.
- Public submodules: `analysis`, `core`, `models`, `survival`, `visualization`, `workflow`.
- Canonical import: `import metainformant.life_events` (submodules: `from metainformant import life_events` then `life_events.<submodule>`).
- Test entry point: `uv run pytest tests/life_events -q` (one pytest directory per invocation).

Keep changes scoped; match existing patterns in this directory.
