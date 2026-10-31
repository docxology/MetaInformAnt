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
| `atomic.py` | Atomic write/replace helpers and temp-directory contexts |
| verification.py | Full-hash snapshots with physical-identity checks limited to one operation |
| `checksums.py` | File and stream hashing (MD5/SHA-1/SHA-256) with verification |
| `data_root.py` | Canonical `AMALGKIT_DATA_ROOT` resolution and prefix remapping |
| `download_manager.py` | High-level download orchestration for batch operations |
| `download_robust.py` | Resilient download with automatic retry and fallback |
| `errors.py` | IO-specific exception types |
| `sra_environment.py` | Campaign-local SRA tool environments (NCBI_SETTINGS, TMPDIR, VDB_CONFIG) |

## Key Functions and Classes

| Symbol | Description |
|--------|-------------|
| `load_json()` / `dump_json()` | Read and write JSON with atomic writes and gzip support |
| `load_yaml()` / `dump_yaml()` | YAML I/O with auto-created parent directories |
| `read_csv()` / `write_csv()` | Pandas-based CSV read/write |
| `read_parquet()` / `write_parquet()` | Columnar Parquet format I/O |
| `expand_and_resolve()` | Expand ~ and resolve to absolute path |
| `is_within()` | Security check: verify a path stays under a parent |
| `read_csv()` / `write_csv()` | CSV read/write (pandas when available, native fallback otherwise) |
| `JsonCache` | TTL-aware JSON object cache backed by filesystem |
| `get_disk_usage()` | Total, used, and free space for a given path |
| `ensure_disk_space()` | Abort early if insufficient disk space |
| `cleanup_old_files()` | Remove files older than a specified age |
| `download_with_progress()` | HTTP(S)/FTP/file download with heartbeat file, retries and HTTP Range resume; a destination that already has the full size counts as downloaded, and a string ending in `/` names a directory (see [Downloads](../../../../docs/core/download.md)) |

## Atomic JSON publication

`dump_json(..., atomic=True)` uses the existing atomic helper's unique owner-only
temporary file in the target directory. Plain and gzip writers close a complete
value before replacing the destination. Existing regular destination permission
bits are retained at publication; new destinations remain owner-only. Concurrent
writes publish complete values with the last replacement winning; a serialization
failure preserves the existing destination and unrelated sibling temporary files.
This does not provide a read-modify-write transaction or a file/directory fsync durability guarantee.

## Usage

```python
from metainformant.core.io import io, paths

data = io.load_yaml("config/workflow.yaml")
io.dump_json(result, "output/result.json")
resolved = paths.expand_and_resolve("~/data/input.txt")
assert paths.is_within(resolved, parent="/safe/dir")
```

[Verification contract and usage](../../../../docs/core/verification.md) documents BoundReads and its limits.
