# Core I/O Module

File handling, serialization, caching, and downloads for MetaInformAnt.

## Purpose

This module provides:
- JSON, JSONL, and CSV I/O operations
- Robust file download with retries and progress tracking
- TTL-based caching
- Path manipulation and validation

## Key Components

| File | Description |
|------|-------------|
| [io.py](io.py) | `load_json`, `dump_json`, `read_csv`, `write_csv` |
| [download.py](download.py) | `download_file`, `download_json` with retries |
| [cache.py](cache.py) | `JsonCache` with TTL support |
| [paths.py](paths.py) | Path resolution, sanitization, and discovery |
| [disk.py](disk.py) | Disk space management and atomic writes |
| [download_manager.py](download_manager.py) | Concurrent download orchestration |

## Usage

```python
from metainformant.core.io import load_json, download_file

config = load_json("config.json")
download_file("https://example.com/data.tar.gz", "/tmp/data.tar.gz")
```

## Related Documentation

- **Parent**: [src/metainformant/core/README.md](../README.md)
- **SPEC**: [SPEC.md](SPEC.md)
- **AGENTS**: [AGENTS.md](AGENTS.md)
- **Data Module**: [../data/README.md](../data/README.md)
