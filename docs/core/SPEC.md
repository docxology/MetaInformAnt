# Specification: core

## Scope
Documentation for the core domain in MetaInformAnt.

## Architecture
- **Dependency Level**: Documentation
- **Component Type**: Guide
- **Location**: `docs/core/`

## Data Structures
- **Files**:
  - `README.md`: Overview
  - `AGENTS.md`: AI Attribution
  - `SPEC.md`: This file
  - `API_REFERENCE.md`: Generated public symbol and signature inventory
  - `CROSS_CHECK_REPORT.md`: Generated source/documentation integrity report
  - `*.md`: Topic-specific guides

## Integration
- **Source**: `src/metainformant/core/`
- **Tests**: `tests/core/test_core_*.py`

## Testing Policy
- **Real Implementation**: All tests must use real implementations. Mocks are strictly prohibited.

## Documentation Contract

`API_REFERENCE.md` is generated from public classes, functions, and methods in
the canonical implementation modules under `src/metainformant/core/`. The
compatibility facade files are intentionally excluded so re-exported names do
not create duplicate documentation. Run both commands below after changing
core source or signatures:

```bash
uv run python scripts/package/generate_core_api_reference.py
uv run python scripts/core_docs_cross_check.py --strict
```

The checker ignores fenced examples and ordinary prose headings, checks every
source-derived signature, and reports stale API-like headings in topic guides.
