# Specification: longread

## 🎯 Scope
Documentation for the longread domain in MetaInformAnt.

## 🧱 Architecture
- **Dependency Level**: Documentation
- **Component Type**: Guide
- **Location**: `docs/longread/`

## 💾 Data Structures
- **Files**:
  - `README.md`: Overview
  - `AGENTS.md`: AI Attribution
  - `SPEC.md`: This file
  - `*.md`: Topic-specific guides

## 🔌 Integration
- **Source**: `src/metainformant/longread/`
- **Tests**: `tests/test_longread_*.py`

## 🧪 Testing Policy
- **Zero Mock**: All tests must use real implementations. Mocks are strictly prohibited.
