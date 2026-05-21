# Agent Directives: scripts/validate

## Role

Maintain repository validation and completeness-audit scripts.

## Rules

- Run with `uv`.
- Keep checks deterministic and repository-root aware.
- Prefer read-only validation; write reports only under `output/` unless a
  script explicitly documents another location.
- Use canonical package paths in examples and diagnostics.
