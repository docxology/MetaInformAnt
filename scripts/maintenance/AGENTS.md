# Agent Directives: scripts/maintenance

## Role

Maintain repository maintenance scripts.

## Rules

- Run commands with `uv`.
- Prefer structured parsers for Python, Markdown, YAML, TOML, and JSON when
  they are available.
- Direct stdlib parsing is allowed here for narrow maintenance glue when tested.
- Do not write generated runtime material into `src/`.
- Keep scripts idempotent where possible and document any broad rewrite.
