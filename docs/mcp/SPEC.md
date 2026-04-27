# Specification: mcp

## Scope
Documentation for the mcp domain in MetaInformAnt.

## Architecture
- **Dependency Level**: Documentation
- **Component Type**: Infrastructure
- **Location**: `docs/mcp/`

## Data Structures
- **Files**:
  - `README.md`: Overview
  - `AGENTS.md`: AI Attribution
  - `SPEC.md`: This file
  - `*.md`: Topic-specific guides

## Integration
- **Source**: `src/metainformant/mcp/`
- **Tests**: `tests/test_mcp_*.py`

## Testing Policy
- **Zero Mock**: All tests must use real implementations. Mocks are strictly prohibited.
