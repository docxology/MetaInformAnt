# Specification: io

## 🎯 Scope
Core Input/Output utilities for METAINFORMANT.

## 🧱 Architecture
- **Dependency Level**: Core
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 13 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`
- `atomic.py`
- `cache.py`
- `checksums.py`
- `data_root.py`
- `disk.py`
- `download.py`
- `download_manager.py`
- `download_robust.py`
- `errors.py`
- ...

## Atomic JSON publication contract

`dump_json(..., atomic=True)` serializes plain or gzip JSON into a unique,
initially owner-only temporary file in the destination directory. Complete data
is closed before atomic replacement. Existing regular destination permission
bits are retained; new destinations keep the helper's owner-only mode.

Concurrent writers publish complete values with the last replacement winning.
Serialization failure preserves the previous destination and unrelated sibling
temporary files. This is publication atomicity, without read-modify-write
coordination or file/directory fsync durability.
