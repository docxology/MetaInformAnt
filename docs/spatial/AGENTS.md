# Agent Directives: docs/spatial

## Role

Maintain spatial transcriptomics documentation for tissue-coordinate expression
analysis.

## Scope

- Document current interfaces under `src/metainformant/spatial/`.
- Link to local spatial topic files such as [index.md](index.md), [io.md](io.md),
  [clustering.md](clustering.md), and [visualization.md](visualization.md).
- Keep examples runnable or clearly marked as pseudocode.

## Boundaries

- Use `metainformant.core.io` for domain data file I/O in examples.
- Direct stdlib parsing is acceptable only for protocol glue, CLI helpers, core
  utilities, or narrow parser internals covered by tests.
- Do not claim support for external spatial platforms unless source code or
  tests in this checkout exercise that path.
- Place generated figures and analysis output under `output/`.
