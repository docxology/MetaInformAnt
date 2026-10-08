# Specification: retrieval

## 🎯 Scope
RNA data retrieval modules.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 2 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`
- `ena_downloader.py`

## Transfer bounds

`ENADownloader(file_workers=1|2)` bounds concurrent file transfers per run.
Distinct output filenames are required before starting transfers. All submitted
transfers settle before fallback; gzip integrity, retry deadlines, partial resume
and retained invalid evidence apply independently to each file.
