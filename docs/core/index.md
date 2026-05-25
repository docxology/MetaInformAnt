# Core Infrastructure: Overview

Shared infrastructure and utilities used by all METAINFORMANT domain modules. The core package provides battle-tested components for I/O, configuration, logging, parallel execution, caching, database connectivity, and workflow orchestration.

## Quick Navigation

- **[Getting Started](./GETTING_STARTED.md)** — 5-minute tutorial with complete pipeline example
- **[Architecture](./ARCHITECTURE.md)** — System design, component interactions, and principles

## Core Components

| Component | Description | Documentation |
|-----------|-------------|---------------|
| [I/O Operations](./io.md) | File I/O, JSON/CSV/TSV/YAML, downloads, atomic writes | `core.io` |
| [Configuration](./config.md) | Config loading, environment overrides, merging | `core.utils.config` |
| [Path Handling](./paths.md) | Path resolution, security, sanitization | `core.io.paths` |
| [Logging](./logging.md) | Structured logging, metadata, environment config | `core.utils.logging` |
| [Caching](./cache.md) | JSON cache with TTL, thread-safe operations | `core.io.cache` |
| [Download](./download.md) | Robust HTTP/FTP downloads, retry, resume, heartbeat | `core.io.download` |
| [Parallel Execution](./parallel.md) | Thread/process pools, resource-aware workers | `core.execution.parallel` |
| [Database](./db.md) | PostgreSQL connectivity, connection pooling | `core.data.db` |
| [Hashing](./hash.md) | SHA256 file and content hashing | `core.utils.hash` |
| [Text Processing](./text.md) | Text cleaning, slugify, gene name standardization | `core.utils.text` |
| [Workflow](./workflow.md) | DAG orchestration, config-driven pipelines | `core.execution.workflow` |

## Module Structure

```
src/metainformant/core/
├── io/                    # Input/Output operations
│   ├── io.py             # Core file I/O (JSON, CSV, YAML, Parquet)
│   ├── paths.py          # Path utilities and security
│   ├── cache.py          # JSON caching with TTL
│   ├── download.py       # Download with retry/resume/heartbeat
│   ├── atomic.py         # Atomic file operations
│   ├── checksums.py      # Checksum verification
│   └── disk.py           # Disk space management
├── utils/                # Utility functions
│   ├── logging.py        # Structured logging
│   ├── config.py         # Configuration loader
│   ├── hash.py           # SHA256 hashing
│   ├── text.py           # Text processing
│   ├── errors.py         # Error hierarchy
│   └── timing.py         # Performance timing
├── execution/            # Execution engines
│   ├── parallel.py       # Parallel execution utilities
│   ├── workflow.py       # Workflow orchestration
│   └── discovery.py      # Symbol discovery
├── data/                 # Data layer
│   ├── db.py             # PostgreSQL integration
│   └── validation.py     # Validation utilities
├── engine/               # Pipeline engines
│   └── workflow_manager.py
└── ui/                   # User interfaces
    └── tui.py            # Terminal UI
```

## Core Design Principles

### 1. **Zero Mocking**
All tests use real implementations. No mock objects. This ensures production reliability.

### 2. **Atomicity**
All file writes use atomic replacement (temp file → rename) to prevent corruption.

### 3. **Observability**
- Consistent log format: `TIMESTAMP | LEVEL | MODULE | MESSAGE`
- Optional structured metadata via `log_with_metadata()`
- Download heartbeats for progress tracking

### 4. **Security**
- Path traversal prevention (`is_safe_path()`)
- Filename sanitization
- SQL injection protection (`sanitize_connection_params()`)

### 5. **Portability**
- Pure `pathlib.Path` (no `os.path`)
- UTF-8 everywhere
- Minimal external dependencies (optional)

## Usage Pattern

```python
# Standard import pattern
from metainformant.core import io
from metainformant.core.io import cache, paths
from metainformant.core.utils import logging, config

# Get logger
logger = logging.get_logger(__name__)

# Ensure directories
output = paths.ensure_directory(Path("output"))

# Load configuration
cfg = config.load_mapping_from_file("config.yaml")

# Download with caching
cached = cache.load_cached_json(cache_dir, "key", ttl_seconds=3600)
if cached is None:
    data = io.download_json(url)
    cache.cache_json(cache_dir, "key", data)

# Process files
for item in io.read_jsonl("data.jsonl"):
    process(item)

logger.info("Pipeline complete")
```

## Environment Variables

| Variable | Purpose | Default |
|----------|---------|---------|
| `CORE_LOG_LEVEL` | Logging level (DEBUG, INFO, WARNING, ERROR) | INFO |
| `AK_THREADS` | Override default thread count | CPU-dependent |
| `AK_WORK_DIR` | Working directory for outputs | `output/` |
| `AK_LOG_DIR` | Directory for log files | `logs/` |
| `PG_HOST` | PostgreSQL host | localhost |
| `PG_PORT` | PostgreSQL port | 5432 |
| `PG_DATABASE` | Database name | metainformant |
| `PG_USER` | Database user | postgres |
| `PG_PASSWORD` | Database password | (empty) |

## Sphinx Documentation

Build with:
```bash
uv run python scripts/package/uv_docs.sh
```

Or manually:
```bash
cd docs
sphinx-build -b html . _build
```

## Contributing

When modifying core components:

1. **Add tests** in `tests/core/test_core_*.py` (no mocks!)
2. **Update documentation** in `docs/core/*.md`
3. **Follow conventions**: `pathlib.Path`, type hints, `get_logger(__name__)`
4. **Check AGENTS.md**: `src/metainformant/core/AGENTS.md` has agent-specific rules

## Related Resources

- [Full API Reference](./SPEC.md) — Type signatures and data structures
- [Agent Directives](../AGENTS.md) — Documentation agent guidelines
- [Examples Directory](../../examples/core/) — Runnable code samples
- [Source Code](../../src/metainformant/core/) — Implementation details
