# IO

Unified file I/O (JSON, YAML, TOML, CSV, Parquet, JSONL), path management with security validation, disk space monitoring, caching, and robust downloads with progress tracking.

## Contents

| File | Purpose |
|------|---------|
| `io.py` | Load/dump for JSON, YAML, TOML, CSV, TSV, Parquet, and JSONL formats |
| `paths.py` | Path resolution, security checks, temp files, and output discovery |
| `download.py` | HTTP/FTP/file downloads with heartbeat, progress, and retry logic |
| `cache.py` | TTL-based JSON cache with directory management |
| `disk.py` | Disk usage monitoring, space checks, and cleanup utilities |
| `download_manager.py` | High-level download orchestration for batch operations |
| `download_robust.py` | Resilient download with automatic retry and fallback |
| `errors.py` | IO-specific exception types |

## Key Functions and Classes

| Symbol | Description |
|--------|-------------|
| `load_json()` / `dump_json()` | Read and write JSON with atomic writes and gzip support |
| `load_yaml()` / `dump_yaml()` | YAML I/O with auto-created parent directories |
| `read_csv()` / `write_csv()` | Pandas-based CSV read/write |
| `read_parquet()` / `write_parquet()` | Columnar Parquet format I/O |
| `expand_and_resolve()` | Expand ~ and resolve to absolute path |
| `is_within()` | Security check: verify a path stays under a parent |
| `download_with_progress()` | Download a file with speed, ETA, and heartbeat tracking |
| `JsonCache` | TTL-aware JSON object cache backed by filesystem |
| `get_disk_usage()` | Total, used, and free space for a given path |
| `ensure_disk_space()` | Abort early if insufficient disk space |
| `cleanup_old_files()` | Remove files older than a specified age |

## Usage

```python
from metainformant.core.io import io, paths

data = io.load_yaml("config/workflow.yaml")
io.dump_json(result, "output/result.json")
resolved = paths.expand_and_resolve("~/data/input.txt")
assert paths.is_within(resolved, base_path="/safe/dir")
```
