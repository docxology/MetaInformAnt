# Agent Directives: docs/tasks

## Role

Task-oriented quick-reference documentation for common METAINFORMANT workflows.

## Module Scope

End-to-end task execution guides mapping user goals to underlying module architectures and scripts.

## Key Source Files

| Path | Description |
|------|-------------|
| `scripts/` | Execution orchestrators referenced in tasks |
| `src/metainformant/` | Underlying library modules |

## Rules

- Name only commands that run in the current checkout, or mark future work as draft/planned.
- Prefer canonical imports and current CLI paths.
- Keep generated outputs under `output/`.
- Use `metainformant.core.io` for domain data file I/O in examples.
