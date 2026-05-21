# Maintenance Scripts

Utilities in this directory support repository maintenance tasks such as docs
generation, import refactors, test metadata updates, and consistency checks.

Run scripts from the repository root with `uv run python` unless a script
documents a different invocation.

## Boundaries

- These scripts may parse repository metadata directly when they are narrow
  maintenance tools.
- Do not write runtime artifacts into `src/`; use `output/` or a documented
  temporary/cache directory.
- Keep generated docs and inventories truthful to the current checkout.
