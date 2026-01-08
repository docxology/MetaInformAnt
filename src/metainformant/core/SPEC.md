# Core Module Specification

## 1. Overview
The Core Module provides the foundational infrastructure for the METAINFORMANT toolkit. It adheres to the "No Mock" policy by providing real, thread-safe implementations for I/O, configuration, caching, and parallel processing.

## 2. Component Specifications

### 2.1 Configuration (`config.py`)
- **Requirement**: Must support loading from YAML/TOML and environment variable overrides.
- **Environment API**: `BIO_` prefix or module specific (e.g. `AK_` for RNA).
- **Type Safety**: strict typing for configuration dictionaries.

### 2.2 I/O Operations (`io.py`)
- **Supported Formats**: JSON, JSONL, CSV, TSV, Parquet.
- **Features**: 
  - Transparent compression (gzip) detection via file extensions.
  - Atomic writes to prevent partial file corruption.
  - Thread-safe file checking.

### 2.3 Caching (`cache.py`)
- **Structure**: Filesystem-based JSON cache.
- **Features**:
  - `JsonCache` class for object-oriented access.
  - Module-level helpers: `cache_json`, `get_cache_info`, `clear_cache_dir`.
  - Automatic TTL expiration and cleanup.

### 2.4 Parallel Processing (`parallel.py`)
- **Mechanism**: `concurrent.futures.ThreadPoolExecutor`.
- **Constraint**: Must preserve order for `thread_map`.
- **Safety**: Graceful handling of worker failures.

## 3. Testing Standards
- **Mocking**: STRICTLY PROHIBITED. All tests must use temporary directories (`tmp_path`) and real file operations.
- **Performance**: Benchmarks required for I/O heavy operations.

## 4. Dependencies
- **Required**: `pathlib`, `typing`, `json`, `yaml`.
- **Optional**: `pandas`, `pyarrow` (graceful degradation if missing).
