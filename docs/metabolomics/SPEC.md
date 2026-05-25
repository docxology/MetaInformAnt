# Specification: metabolomics

## Scope
Documentation for the metabolomics domain in MetaInformAnt.

## Architecture
- **Dependency Level**: Documentation
- **Component Type**: Analysis
- **Location**: `docs/metabolomics/`

## Data Structures
- **Files**:
  - `README.md`: Overview
  - `AGENTS.md`: AI Attribution
  - `SPEC.md`: This file
  - `*.md`: Topic-specific guides

## Integration
- **Source**: `src/metainformant/metabolomics/`
- **Tests**: `tests/metabolomics/test_metabolomics_*.py`

## Testing Policy
- **Zero Mock**: All tests must use real implementations. Mocks are strictly prohibited.
