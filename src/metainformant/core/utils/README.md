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
| `progress.py` | tqdm progress bars (with stdlib fallback) and step-logging for long-running tasks |
| `timing.py` | Execution timing decorator, Timer context manager, rate limiter |
| `optional_deps.py` | Suppress or warn about missing optional dependencies |
| `symbols.py` | AST-based symbol indexing, fuzzy search, and reference finding |
| `batches.py` | Chunked iteration and DataFrame batching/memory reporting |
| `newick.py` | Newick tree parsing into label->children mappings |
| `seeds.py` | Deterministic SHA-256 replicate-seed derivation |
| `watchdog.py` | Process watchdog for stalled long-running pipelines |

## Key Functions and Classes

| Symbol | Description |
|--------|-------------|
| `get_logger()` | Return a configured logger for a module name |
| `setup_logger()` | Create a logger with optional file handler and level |
| `load_mapping_from_file()` | Load YAML or TOML config into a dict |
| `retry_with_backoff()` | Decorator factory that retries the decorated callable with exponential backoff |
| `merge_configs()` | Deep-merge a base config with overrides |
| `METAINFORMANTError` | Base exception for all project errors |
| `retry_with_backoff()` | Retry a callable with exponential backoff |
| `error_context()` | Context manager that annotates exceptions with message |
| `sha256_file()` | Compute SHA-256 hash of a file |
| `slugify()` | Convert text to a URL/filename-safe slug |
| `standardize_gene_name()` | Normalize gene identifiers to a canonical form |
| `progress_bar()` | Display a tqdm progress bar for iterable processing |
| `get_logger_with_level()` | Logger with explicit or env-based level |
| `timeout_after()` | Context manager that flags blocks exceeding a timeout |
| `chunked()` | Split an iterable into fixed-size batches |
| `deterministic_replicate_seeds()` | Stable per-replicate seeds from one base seed |
| `find_symbol_references()` | Find all references to a symbol across the repo |
| `fuzzy_find_symbol()` | Approximate symbol search with similarity scores |
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
cfg = apply_env_overrides(cfg, prefix="AMALGKIT")

@retry_with_backoff(max_attempts=3)
def fetch_data(url: str) -> dict:
    ...
```
