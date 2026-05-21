# Agent Directives: scripts/reorganize

## Role

Support mechanical codebase reorganizations and import rewrites.

## Rules

- Use `uv` for all execution.
- Preserve public imports with compatibility facades when moving modules.
- Keep write sets narrow and review diffs after generated rewrites.
- Direct stdlib parsing is allowed for repository rewrite tooling when the
  behavior is deterministic and covered by checks.
