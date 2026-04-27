# Core Infrastructure Architecture

## Overview

The `metainformant.core` package provides foundational infrastructure shared across all METAINFORMANT domain modules. It follows a layered architecture with clear separation of concerns, enabling reusable bioinformatics pipeline components.

## Architectural Principles

### 1. **Layered Design**
- **I/O Layer** (`io/`): File operations, network downloads, path handling, caching
- **Utility Layer** (`utils/`): Logging, configuration, hashing, text processing, error handling
- **Execution Layer** (`execution/`): Parallel processing, workflow orchestration
- **Data Layer** (`data/`): Database connectivity, validation
- **Engine Layer** (`engine/`): Pipeline management and TUI

### 2. **Dependency Flow**
```
Domain Modules (rna, gwas, dna, ...)
    ↓ use
Core Infrastructure (metainformant.core)
    ↓ uses
Python Standard Library + Optional deps (requests, psycopg2, pandas, pyarrow, PyYAML)
```

**Rule**: Lower layers never import from higher layers. Domain modules can import any core component; core components should not import domain modules.

### 3. **Idiomatic Error Handling**
All core modules raise custom exceptions from `metainformant.core.utils.errors`:
- `ConfigError` for configuration issues
- `IOError` (aliased as `CoreIOError`) for file/network operations
- `ValidationError` for data validation failures
- `WorkflowError` for orchestration problems

### 4. **Structured Logging**
All modules use `metainformant.core.utils.logging.get_logger(__name__)` for consistent log formatting and metadata support.

## Component Interaction Diagram

```mermaid
graph TB
    subgraph "Domain Layer"
        RNA[RNA Module]
        GWAS[GWAS Module]
        DNA[DNA Module]
        PROTEIN[Protein Module]
    end

    subgraph "Core Infrastructure Layer"
        IO[io/<br/>I/O, Downloads,<br/>Paths, Cache]
        UTILS[utils/<br/>Config, Logging,<br/>Hash, Text, Errors]
        EXEC[execution/<br/>Parallel, Workflow,<br/>Discovery]
        DATA[data/<br/>Database,<br/>Validation]
        ENGINE[engine/<br/>Workflow Manager,<br/>TUI]
    end

    subgraph "External Dependencies"
        STD[Python Stdlib<br/>json, pathlib,<br/>logging, ...]
        OPT[Optional Packages<br/>requests, psycopg2,<br/>pandas, PyYAML]
    end

    RNA --> IO
    RNA --> UTILS
    RNA --> EXEC
    GWAS --> IO
    GWAS --> UTILS
    DNA --> DATA
    PROTEIN --> IO

    IO --> STD
    UTILS --> STD
    EXEC --> STD
    DATA --> STD

    IO --> OPT
    UTILS --> OPT
    EXEC --> OPT
    DATA --> OPT

    style IO fill:#e1f5ff
    style UTILS fill:#fff4e1
    style EXEC fill:#e8f5e9
    style DATA fill:#f3e5f5
    style ENGINE fill:#ffebee
```

## Module Relationships and Data Flow

### I/O Subsystem (`io/`)

The I/O subsystem provides a unified interface for all file and network operations:

```
┌─────────────────────────────────────────────────────────────────┐
│                         User Code                               │
│  from metainformant.core import io                              │
│  data = io.load_json("file.json")                               │
└────────────────────────────┬────────────────────────────────────┘
                             │
        ┌────────────────────┼────────────────────┐
        │                    │                    │
        ▼                    ▼                    ▼
┌───────────────┐   ┌───────────────┐   ┌───────────────┐
│ io.py         │   │ download.py   │   │ paths.py      │
│ - load_json   │   │ - download_   │   │ - expand_     │
│ - dump_json   │   │   with_progress│   │   and_resolve │
│ - read_csv    │   │ - retry logic │   │ - is_safe_path│
│ - write_jsonl │   │ - resume      │   │ - sanitize_   │
└───────┬───────┘   └───────┬───────┘   └───────┬───────┘
        │                   │                    │
        └───────────────────┼────────────────────┘
                            │
        ┌───────────────────┼────────────────────┐
        │                   │                    │
        ▼                   ▼                    ▼
┌───────────────┐   ┌───────────────┐   ┌───────────────┐
│ atomic.py     │   │ cache.py      │   │ disk.py       │
│ - atomic write│   │ - JsonCache   │   │ - space       │
│ - safe rename │   │ - TTL entries │   │   management  │
└───────────────┘   └───────────────┘   └───────────────┘
```

**Key Design Decisions**:
- All file writes use atomic operations (write to temp, then rename) to prevent corruption
- Gzip compression handled transparently via `.gz` suffix detection
- Progress tracking via heartbeat files for long-running downloads
- Cache keys are sanitized for filesystem safety

### Configuration Subsystem (`utils/config.py`)

```
┌──────────────┐     ┌──────────────┐     ┌──────────────┐
│ YAML/TOML/   │────▶│ load_mapping │────▶│ merge_configs│
│ JSON Files   │     │ from_file    │     │ (deep merge) │
└──────────────┘     └──────────────┘     └──────────────┘
                                                    │
┌──────────────┐     ┌──────────────┐              ▼
│ Environment  │────▶│ apply_env_  │──────────▶│ coerce_config │
│ Variables    │     │ overrides    │           │ types        │
└──────────────┘     └──────────────┘           └──────────────┘
```

**Environment Variable Conventions**:
- `AK_*`: Generic METAINFORMANT overrides (threads, work_dir, log_dir)
- `PG_*`: PostgreSQL database configuration (PG_HOST, PG_PORT, PG_DATABASE, PG_USER, PG_PASSWORD)
- `DB_*`: Alternate database variable names (DB_NAME, DB_USER, DB_PASSWORD)

### Logging Subsystem (`utils/logging.py`)

All loggers follow a consistent format:
```
YYYY-MM-DD HH:MM:SS | LEVEL | module.name | message
```

**Features**:
- Structured metadata logging via `log_with_metadata()`
- Environment-controlled log level via `CORE_LOG_LEVEL`
- Automatic console handler attachment on first `get_logger()` call
- Optional file handler support via `setup_logger()`

### Parallel Execution (`execution/parallel.py`)

Provides three mapping strategies:

```mermaid
graph LR
    A[Input Items] --> B{Task Type?}

    B -->|I/O-bound| C[thread_map<br/>2-4x cores]
    B -->|CPU-bound| D[process_map<br/>core count]
    B -->|Batch| E[parallel_batch<br/>chunked]

    C --> F[Results]
    D --> F
    E --> F
```

**Resource-Aware Worker Count**:
- I/O tasks: up to 4× CPU cores (capped at 32)
- CPU tasks: CPU count minus 1 (reserve for OS)
- Memory constrained: scales by available RAM (256 MB per worker default)

### Workflow Orchestration (`execution/workflow.py`)

DAG-based execution with topological sorting:

```mermaid
graph TD
    A[Config with<br/>steps] --> B[Validate DAG]
    B --> C{Valid?}
    C -->|No| D[Return errors]
    C -->|Yes| E[Topological<br/>Sort]
    E --> F[Execute Steps<br/>in Order]
    F --> G{Step<br/>Failed?}
    G -->|Yes| H[Skip Downstream<br/>Dependents]
    G -->|No| I[Store Result]
    H --> J[Continue]
    I --> J
    J --> K[All Done?]
    K -->|No| F
    K -->|Yes| L[Return<br/>Results]
```

**Features**:
- Dependency resolution via Kahn's algorithm
- Automatic function loading via `module.path:function_name` strings
- Failed step propagation (downstream steps skipped)
- Thread-safe result storage

### Database Integration (`data/db.py`)

Connection pooling with context manager support:

```python
with get_connection(...) as conn:
    with conn.connect() as db_conn:
        results = conn.execute_query(db_conn, "SELECT ...")
```

**Connection Pool**: `ThreadedConnectionPool` with configurable min/max connections.

### Caching (`io/cache.py`)

Two-tier caching strategy:
1. **Key-based** (`JsonCache` class): Thread-safe, per-key TTL
2. **File-based** (`get_json_cache`, `set_json_cache`): Simple file-based cache

Cache entries stored as JSON with metadata:
```json
{
  "value": { ... },
  "created_at": 1703023456.789,
  "ttl_seconds": 3600,
  "expires_at": 1703027056.789
}
```

## Cross-Cutting Concerns

### Atomicity

All write operations in `io/` use atomic file replacement:
1. Write to temporary file (`.tmp` suffix)
2. `fsync()` if available
3. Atomic rename (`os.replace()`)

This prevents partial writes on crash/power loss.

### Observability

- **Logging**: Every module uses `get_logger(__name__)`
- **Heartbeats**: Long downloads write progress to `.downloads/*.heartbeat.json`
- **Metrics**: Cache stats, directory sizes, workflow timing

### Security

- **Path safety**: `is_safe_path()` blocks path traversal (`..`) and absolute paths to sensitive locations
- **Filename sanitization**: Removes special characters, control codes, normalizes length
- **SQL injection prevention**: `sanitize_connection_params()` removes dangerous keywords and characters
- **Credential handling**: Database passwords not logged; URLs masked (`postgresql://user:***@host/db`)

### Portability

- **Path handling**: Pure `pathlib.Path` everywhere, never `os.path`
- **Text encoding**: UTF-8 default, auto-detection for gzip
- **Platform independence**: No shell commands; Python stdlib only where possible

## Performance Considerations

### I/O Optimization
- Streaming reads (`Iterator` for JSONL, delimited files)
- Chunked downloads (default 1 MB chunks)
- Gzip compression for large JSON datasets
- Lazy loading of optional dependencies (pandas, PyYAML, psycopg2)

### Memory Efficiency
- Generators for large file iteration (`read_jsonl`, `read_delimited`)
- Batch processing utilities (`parallel_batch`)
- Directory size scanning with error tolerance (skip unreadable files)

### Concurrency Model
- **ThreadPoolExecutor** for I/O-bound tasks (downloads, file operations)
- **ProcessPoolExecutor** for CPU-bound tasks (data transformation)
- **RLock** in `JsonCache` for thread-safe caching
- Thread-safe workflow result storage

## Extension Points

### Custom Download Protocols
Implement a function matching `_download_http`/`_download_file_url` signature and register in `download_with_progress`.

### Custom Configuration Formats
Add parser in `config.py`'s `load_mapping_from_file()` by extending suffix handling.

### Custom Workflow Steps
Pass `"module.path:function"` strings to `BaseWorkflowOrchestrator.add_step()` for dynamic loading.

## Testing Philosophy

**Zero Mocking Policy**: All tests use real implementations:
- Real file I/O with temporary directories
- Real network calls (or simulated via `httpbin`/local servers)
- Real database connections (PostgreSQL instance required)
- Real compression/decompression

This ensures production reliability and prevents false positives.

## Common Pitfalls

1. **Thread safety**: Don't share `JsonCache` across processes (only threads)
2. **Path traversal**: Always use `expand_and_resolve()` and `is_within()` for user-supplied paths
3. **Connection pooling**: Database connections must be released via context manager or `close()`
4. **Atomic writes**: Don't open the same file for writing concurrently from multiple processes
5. **Workflow dependency cycles**: Topological sort detects cycles; ensure DAG structure

## Future Directions

- **Async I/O**: Potential migration to `asyncio` for higher concurrency
- **Cache eviction policies**: LRU, LFU beyond TTL
- **Distributed execution**: Integration with Dask or Ray
- **Schema validation**: JSON Schema integration for configs
- **Metrics collection**: Prometheus exporters for production monitoring

## Related Documentation

- [I/O Operations](./io.md): File and network I/O deep dive
- [Configuration](./config.md): Config loading and environment overrides
- [Logging](./logging.md): Structured logging patterns
- [Parallel Execution](./parallel.md): Thread/process pool utilities
- [Workflow Orchestration](./workflow.md): DAG-based pipeline execution
