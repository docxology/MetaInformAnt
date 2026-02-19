# Specification: information

## 🎯 Scope
Documentation for the information domain in MetaInformAnt.

## 🧱 Architecture
- **Dependency Level**: Documentation
- **Component Type**: Guide
- **Location**: `docs/information/`

## 💾 Data Structures
- **Files**:
  - `README.md`: Overview
  - `AGENTS.md`: AI Attribution
  - `SPEC.md`: This file
  - `*.md`: Topic-specific guides

## 🔌 Integration
- **Source**: `src/metainformant/information/`
- **Tests**: `tests/test_information_*.py`

## 🧪 Testing Policy
- **Zero Mock**: All tests must use real implementations. Mocks are strictly prohibited.
