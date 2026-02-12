# Utils

Shared utility modules for logging, configuration, error handling, hashing, text processing, progress display, timing, and optional dependency management.

## Contents

| File | Purpose |
|------|---------|
| `logging.py` | Logger creation with file output and environment-based level config |
| `config.py` | YAML/TOML config loading, env overrides, merging, and schema discovery |
| `errors.py` | Exception hierarchy and retry/error-context helpers |
| `hash.py` | SHA-256 hashing for bytes, files, strings, and directory trees |
| `text.py` | Text normalization, slugification, gene name standardization |
| `progress.py` | Rich progress bars and step-logging for long-running tasks |
| `timing.py` | Execution timing decorator, Timer context manager, rate limiter |
| `optional_deps.py` | Suppress or warn about missing optional dependencies |
| `symbols.py` | AST-based symbol indexing, fuzzy search, and reference finding |

## Key Functions and Classes

| Symbol | Description |
|--------|-------------|
| `get_logger()` | Return a configured logger for a module name |
| `setup_logger()` | Create a logger with optional file handler and level |
| `load_mapping_from_file()` | Load YAML or TOML config into a dict |
| `apply_env_overrides()` | Override config keys from environment variables |
| `merge_configs()` | Deep-merge a base config with overrides |
| `METAINFORMANTError` | Base exception for all project errors |
| `retry_with_backoff()` | Retry a callable with exponential backoff |
| `error_context()` | Context manager that annotates exceptions with message |
| `sha256_file()` | Compute SHA-256 hash of a file |
| `slugify()` | Convert text to a URL/filename-safe slug |
| `standardize_gene_name()` | Normalize gene identifiers to a canonical form |
| `progress_bar()` | Display a rich progress bar for iterable processing |
| `timed()` | Decorator that logs function execution duration |
| `index_functions()` | Index all functions in a repo via AST parsing |
| `find_symbol()` | Locate a symbol definition by name |

## Usage

```python
from metainformant.core.utils.logging import get_logger
from metainformant.core.utils.config import load_mapping_from_file, apply_env_overrides
from metainformant.core.utils.errors import retry_with_backoff

logger = get_logger(__name__)
cfg = load_mapping_from_file("config/workflow.yaml")
cfg = apply_env_overrides(cfg, prefix="AK")
result = retry_with_backoff(fetch_data, max_retries=3)
```
