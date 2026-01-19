# Core Utilities Module

Configuration, logging, hashing, progress tracking, and text processing for MetaInformAnt.

## Purpose

This module provides foundational utilities used across all other modules:
- Configuration loading and merging
- Custom exception hierarchy
- Content and file hashing
- Structured logging
- Progress bars and tracking
- Symbol indexing for codebase navigation

## Key Components

| File | Description |
|------|-------------|
| [config.py](config.py) | `load_mapping_from_file`, `merge_configs`, `apply_env_overrides` |
| [logging.py](logging.py) | `get_logger`, `setup_logger` |
| [hash.py](hash.py) | `sha256_file`, `sha256_bytes`, `verify_file_integrity` |
| [progress.py](progress.py) | `progress_bar`, `task_context` |
| [text.py](text.py) | `slugify`, `safe_filename`, `standardize_gene_name` |
| [symbols.py](symbols.py) | `index_functions`, `find_symbol` |
| [errors.py](errors.py) | `retry_with_backoff`, `safe_execute` |
| [optional_deps.py](optional_deps.py) | `warn_optional_dependency` |

## Usage

All utility functions are re-exported from the package for easy access:

```python
from metainformant.core.utils import (
    get_logger,
    load_mapping_from_file,
    sha256_file,
    progress_bar,
    slugify,
)

logger = get_logger(__name__)
config = load_mapping_from_file("settings.yaml")
file_hash = sha256_file("/path/to/file.txt")
```

## Related Documentation

- **Parent**: [src/metainformant/core/README.md](../README.md)
- **SPEC**: [SPEC.md](SPEC.md)
- **AGENTS**: [AGENTS.md](AGENTS.md)
- **I/O Module**: [../io/README.md](../io/README.md)
