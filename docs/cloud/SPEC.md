# Specification: cloud

## Scope
Documentation for the cloud domain in MetaInformAnt.

## Architecture
- **Dependency Level**: Documentation
- **Component Type**: Deployment
- **Location**: `docs/cloud/`

## Data Structures
- **Files**:
  - `README.md`: Overview
  - `AGENTS.md`: AI Attribution
  - `SPEC.md`: This file
  - `*.md`: Topic-specific guides

## Integration
- **Source**: `src/metainformant/cloud/`
- **Tests**: `tests/test_cloud_*.py`

## Testing Policy
- **Zero Mock**: All tests must use real implementations. Mocks are strictly prohibited.
