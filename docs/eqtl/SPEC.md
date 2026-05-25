# Specification: eqtl

## Scope
Documentation for the eqtl domain in MetaInformAnt.

## Architecture
- **Dependency Level**: Documentation
- **Component Type**: Integration
- **Location**: `docs/eqtl/`

## Data Structures
- **Files**:
  - `README.md`: Overview
  - `AGENTS.md`: AI Attribution
  - `SPEC.md`: This file
  - `*.md`: Topic-specific guides

## Integration
- **Source**: `src/metainformant/gwas/finemapping/eqtl.py` (primary), `src/metainformant/gwas/finemapping/colocalization.py`, `src/metainformant/gwas/analysis/eqtl.py`
- **Scripts**: `scripts/eqtl/` (pipeline orchestrators)
- **Tests**: `tests/gwas/test_gwas_finemapping_eqtl.py`, `tests/gwas/test_gwas_analysis_eqtl.py`

## Testing Policy
- **Real Implementation**: All tests must use real implementations. Mocks are strictly prohibited.
