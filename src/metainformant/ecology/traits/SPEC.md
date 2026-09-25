# Specification: traits

## 🎯 Scope
Trait-based ecology sub-package.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 1 Python module
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`

> Consolidation note (2026-09): the former `functional.py` duplicate
> (numpy-based CWM/FRic/Rao with a non-Villeger FEve) was removed. The
> canonical functional diversity implementation is
> `metainformant.ecology.analysis.functional` (Villeger MST-based FEve/FDiv,
> hull-area FRic, list-based pure Python).
