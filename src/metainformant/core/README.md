# Core Utilities Module

The `core` module provides essential shared utilities and infrastructure used across all METAINFORMANT domains. These utilities handle common operations like configuration management, I/O operations, logging, and parallel processing.

## Overview

This module contains foundational components that are used throughout the METAINFORMANT package. Each submodule provides specific functionality while maintaining consistency and reusability across the entire codebase.

## Submodules

### Configuration (`config.py`)
Centralized configuration management with environment variable overrides.

**Key Features:**
- YAML and TOML configuration file support (optional dependencies)
- Environment variable integration and validation
- Database connection configuration (PostgreSQL)
- Configuration inheritance and merging
- Runtime configuration validation

**Usage:**
```python
from metainformant.core import config

# Load configuration from file with environment overrides
cfg = config.load_config("config.yaml")

# Access configuration values
db_config = cfg.database
threads = cfg.compute.threads

# Validate configuration
config.validate_config(cfg)
```

### I/O Operations (`io.py`)
Robust file I/O utilities with support for multiple formats and compression.

**Key Features:**
- JSON/JSONL reading and writing
- CSV/TSV processing with type inference
- Gzip compression support
- Atomic write operations
- Directory creation and path validation

**Usage:**
```python
from metainformant.core import io

# JSON operations
data = io.read_json("data.json")
io.write_json(data, "output.json")

# CSV operations with type inference
df = io.read_csv("data.csv")
io.write_csv(df, "output.csv")

# Compressed file handling
io.write_json_gz(data, "compressed.json.gz")
```

### Logging (`logging.py`)
Consistent, structured logging across all modules.

**Key Features:**
- Structured logging with consistent formatting
- Multiple output targets (console, file, etc.)
- Log level management and filtering
- Context-aware logging for distributed operations

**Usage:**
```python
from metainformant.core import logging

# Setup logger
logger = logging.setup_logger("my_module")

# Log messages with context
logger.info("Processing started", extra={"items": 100})
logger.error("Processing failed", extra={"error": str(exc)})
```

### Parallel Processing (`parallel.py`)
Thread-based parallel execution with order preservation.

**Key Features:**
- Thread pool execution with result ordering
- Progress tracking and cancellation
- Memory-efficient processing of large datasets
- Exception handling and aggregation

**Usage:**
```python
from metainformant.core import parallel

# Parallel map with progress tracking
def process_item(item):
    return item * 2

results = parallel.thread_map(process_item, items, max_workers=4)
```

### Path Management (`paths.py`)
Path expansion, resolution, and containment validation.

**Key Features:**
- Path expansion with environment variables and user home
- Path resolution and normalization
- Containment checks for security
- Cross-platform path handling

**Usage:**
```python
from metainformant.core import paths

# Expand paths with environment variables
expanded = paths.expand_path("~/data/{species}")
resolved = paths.resolve_path("./relative/path")

# Security containment checks
if paths.is_contained(resolved, "/allowed/directory"):
    # Safe to use path
    pass
```

### Caching (`cache.py`)
JSON-based caching with TTL (Time To Live) support.

**Key Features:**
- Disk-based JSON caching
- TTL-based expiration
- Cache size management
- Thread-safe operations

**Usage:**
```python
from metainformant.core import cache

# Create cache with TTL
cache_obj = cache.JSONCache("cache_dir", ttl_seconds=3600)

# Cache operations
cache_obj.set("key", {"data": "value"})
value = cache_obj.get("key")
```

### Text Processing (`text.py`)
Text normalization and processing utilities.

**Key Features:**
- Unicode normalization and cleaning
- Case conversion with locale awareness
- Whitespace normalization
- Text encoding detection and conversion

**Usage:**
```python
from metainformant.core import text

# Text normalization
clean_text = text.normalize_text("  Mixed\tCase\nText  ")
title_case = text.title_case("species name")
```

### Hash Functions (`hash.py`)
Content and file hashing utilities.

**Key Features:**
- Multiple hash algorithms (MD5, SHA256, etc.)
- File content hashing
- String hashing with encoding support
- Hash comparison utilities

**Usage:**
```python
from metainformant.core import hash

# File hashing
file_hash = hash.hash_file("data.txt")
content_hash = hash.hash_string("content to hash")
```

### Database Integration (`db.py`) - Optional
Database client helpers for PostgreSQL integration.

**Key Features:**
- Connection pooling and management
- Query execution with parameter binding
- Transaction management
- Result set processing

**Usage:**
```python
from metainformant.core import db

# Database operations (requires database configuration)
with db.get_connection(db_config) as conn:
    results = db.execute_query(conn, "SELECT * FROM table WHERE id = %s", (123,))
```

## Architecture Principles

### Defensive Imports
Optional dependencies are imported defensively to avoid hard failures during unrelated operations. This ensures the package remains functional even when optional dependencies are missing.

### Consistency
All utilities follow consistent patterns for configuration, error handling, and return values.

### Performance
Operations are designed to be memory-efficient and suitable for large-scale biological data processing.

### Extensibility
Utilities are designed to be easily extended with new formats, algorithms, or integrations.

## Error Handling

All core utilities implement comprehensive error handling:
- Clear error messages with context
- Appropriate exception types
- Graceful degradation for optional features
- Input validation and sanitization

## Testing

Each submodule includes comprehensive tests covering:
- Normal operation scenarios
- Error conditions and edge cases
- Performance characteristics
- Integration with other modules

## Dependencies

- **Required**: Standard library only (pathlib, json, csv, etc.)
- **Optional**: PyYAML (for YAML config), tomllib (Python 3.11+ for TOML)
- **Database**: psycopg2 (optional, for PostgreSQL integration)

## Usage in Other Modules

Core utilities are designed to be imported and used throughout the METAINFORMANT package:

```python
from metainformant.core import config, io, logging, paths

# Typical usage pattern in other modules
logger = logging.get_logger(__name__)
cfg = config.load_config()
data = io.read_json(paths.expand_path(cfg.input_file))
```

This ensures consistency and reduces code duplication across the entire codebase.
