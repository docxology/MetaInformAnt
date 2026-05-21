# Reorganization Scripts

Scripts in this directory support import-map construction and mechanical
repository reorganizations.

## Use

Run from the repository root and inspect generated reports before applying
rewrites:

```bash
uv run python scripts/reorganize/build_import_map.py
uv run python scripts/reorganize/rewrite_imports.py --help
```

## Boundaries

- Treat reports here as maintenance artifacts, not runtime package data.
- Keep public compatibility paths in place when moving implementations.
- Never rewrite unrelated user changes.
