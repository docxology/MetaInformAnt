# Core Module Rules (`metainformant.core`)

## Purpose

Shared utilities for the domain modules. Most other modules build on it; `menu` and `cloud` do not import it.

## Dependencies

- **Required at import**: standard library and `psutil` (`utils/watchdog.py` imports it at module level and `core.utils` imports `watchdog`)
- **Optional**: guarded by `try/except ImportError` or imported inside the function that needs them: `pandas` (`read_csv` and `write_csv` fall back to native readers and writers; Parquet also needs `pyarrow`), `PyYAML` (YAML configs; `load_mapping_from_file` falls back to a minimal parser for simple mappings), `requests` (HTTP downloads), `psycopg2` (PostgreSQL helpers), `jsonschema` (`validate_json_schema`), `tqdm` (progress bars)

## Source Structure

```text
src/metainformant/core/
  db.py                 compatibility re-export of data/db.py
  ncbi.py               NCBI contact policy (resolve_ncbi_contact)
  sequence.py           shared sequence primitives (reverse_complement, validate_dna_sequence)
  data/
    db.py, validation.py
  engine/
    workflow_manager.py
  execution/
    discovery.py, parallel.py, workflow.py
  io/
    atomic.py, cache.py, checksums.py, data_root.py, disk.py
    download.py, download_manager.py, download_robust.py
    errors.py, io.py, paths.py, sra_environment.py
  ui/
    tui.py
  utils/
    batches.py, config.py, errors.py, hash.py, logging.py
    newick.py, optional_deps.py, progress.py, seeds.py, symbols.py
    text.py, timing.py, watchdog.py
```

## Package Management

- **ALWAYS use `uv`** for all Python package management and environment operations
- Use `uv venv` to create virtual environments
- Use `uv pip install` for installing packages
- Use `uv run` to execute commands in the project environment
- Use `uv sync` to sync dependencies from `pyproject.toml` and `uv.lock`
- Use `uv add <package>` to add new dependencies
- Use `uv remove <package>` to remove dependencies
- **NEVER use `pip`, `pip3`, `python -m pip`, or `python3 -m pip` directly**

## Key Submodules

### Configuration (`config`)

- `load_mapping_from_file` loads a mapping from YAML, TOML or JSON; it raises `FileNotFoundError` for a missing file and `ValueError` for an unsupported suffix or a non-mapping top level
- `apply_env_overrides(config, prefix=...)` reads only `{PREFIX}_THREADS`, `{PREFIX}_WORK_DIR` and `{PREFIX}_LOG_DIR`; `load_typed_env` reads other typed variables
- `merge_configs` (deep merge) and `coerce_config_types`
- `PostgresConfig` (frozen dataclass) and `load_postgres_config_from_env(prefix="PG")` for PostgreSQL settings

**Patterns**:

```python
from metainformant.core.utils.config import load_mapping_from_file, apply_env_overrides

raw = load_mapping_from_file("config/file.yaml")
config = apply_env_overrides(raw, prefix="MODULE")  # "MODULE" is a placeholder prefix
```

### I/O (`io`)

- JSON, JSONL, CSV, TSV, YAML and Parquet readers and writers, plus a TOML reader (`load_toml`)
- Gzip-aware file operations: a `.gz` suffix is handled transparently by `open_text_auto`, `load_json`, `dump_json`, `read_jsonl` and `write_jsonl`
- Atomic writes: `dump_json`, `write_jsonl` and `write_delimited` write to a temporary file and rename it by default (`atomic=True`); `atomic_write` in `core.io.atomic` does the same for other files
- Automatic directory creation: the writers create missing parent directories
- `read_csv` returns a pandas DataFrame when pandas is installed and a dict of lists otherwise
- `load_json` raises `FileNotFoundError` for a missing file and `metainformant.core.utils.errors.IOError` (a `METAINFORMANTError`, not the built-in) when parsing or reading fails

**CRITICAL**: Domain data file I/O goes through `metainformant.core.io`. Direct `json` or `csv` parsing is allowed only inside `core`, in protocol adapters, in subprocess/CLI glue, and in narrow parser internals covered by tests (the rule in `src/metainformant/core/AGENTS.md`). Outside those cases do not add `import json` or `import csv` to a domain module.

**Patterns**:

```python
from metainformant.core import io

# Reading
data = io.load_json("data/file.json")
records = list(io.read_jsonl("data/records.jsonl"))
df = io.read_csv("data/table.csv")

# Writing (creates parent dirs automatically)
io.dump_json(data, "output/result.json")
io.write_csv(df, "output/table.csv")
io.write_jsonl(records, "output/records.jsonl")

# Gzip-aware operations
with io.open_text_auto("data/large_file.txt.gz") as f:
    content = f.read()
```

**PROHIBITED** (outside the cases above):

```python
import json # WRONG - use metainformant.core.io instead
import csv # WRONG - use metainformant.core.io instead

data = json.load(f) # WRONG
json.dump(data, f) # WRONG
```

### Paths (`paths`)

- Path expansion and resolution: `expand_and_resolve` expands `~` and resolves to an absolute path without requiring the path to exist
- Containment checks (security): `is_within(path, parent)` resolves both paths before comparing
- Safe path validation: `is_safe_path` is a string check that rejects `..`, paths starting with `/etc/` or `/root/`, and the characters `;`, `|`, `&` and `$`; it is not a containment check, so use `is_within` for that

**Patterns**:

```python
from metainformant.core.io import paths
from pathlib import Path

# Expand and resolve
path = paths.expand_and_resolve("~/data/file.txt")

# Security check
if paths.is_within(user_path, base_dir):
    # Safe to use
    pass
```

### Logging (`logging`)

- One console format for every logger the helpers create: `time | level | name | message`
- Module-specific loggers from `get_logger(__name__)`; `setup_logger(name, log_file=...)` adds a file handler
- `log_with_metadata` appends a JSON metadata dict to the message (or logs it as a separate `METADATA:` line with `structured=True`)
- `CORE_LOG_LEVEL` sets the level through `configure_logging_from_env` and `get_logger_with_level`

**Patterns**:

```python
from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)
logger.info("Operation started")
logger.error("Error occurred", exc_info=True)
```

### Cache (`cache`)

- `JsonCache(cache_dir, ttl_seconds=3600)`: thread-safe JSON files with a time-to-live per entry
- The caller chooses the cache directory; keep it under `output/`

### Data root (`io.data_root`)

- `resolve_data_root(data_root=None, default="output/amalgkit")` returns an absolute path: the explicit argument first, then `AMALGKIT_DATA_ROOT`, then `default`
- `remap_amalgkit_prefix` maps config paths under `output/amalgkit` onto that root

### NCBI contact (`ncbi`)

- `resolve_ncbi_contact` uses an explicit `email`, then `NCBI_EMAIL`; with neither it raises `NCBIContactError` (a `ValueError`) unless `allow_anonymous=True`, which emits a `RuntimeWarning` and returns an anonymous contact

## Output Paths

By convention (`cursorrules/core.cursorrules`):

- Core outputs: `output/core/`
- Cache: `output/core/cache/`
- Logs: `output/core/logs/`

No core function creates these directories by itself: callers pass the directory (for example `JsonCache("output/core/cache")` or `setup_logger(name, log_file="output/core/logs/run.log")`). `paths.get_cache_dir()` and `paths.get_logs_dir()` return `<repository root>/cache` and `<repository root>/logs`, and the symbol-index cache is written to `output/.discovery_cache/`.

### Parallel (`parallel`)

- Thread-based `thread_map` (results in input order by default; `ordered=False` returns completion order) and `thread_map_unordered`
- Process-based `process_map` for CPU-bound work; thread-based `rate_limited_map` (maximum calls per second) and `parallel_batch`
- `resource_aware_workers` suggests a worker count from CPU cores and available memory
- When a task raises, `thread_map` finishes the remaining tasks and then re-raises the first error that completed

### Errors (`errors`)

- Exception hierarchy rooted at `METAINFORMANTError`: `ConfigError`, `IOError` (the core class, not the built-in), `ValidationError`, `NetworkError`, `CacheError`, `DownloadError` (a `NetworkError`), `TermNotFoundError`, `PipelineError`, `DependencyError`, `ResourceError`
- `ValidationError` is also a `ValueError`, so `except ValueError` catches it; the other classes (for example `ConfigError` and `NetworkError`) are not `ValueError`s
- Helpers: `retry_with_backoff`, `error_context`, `safe_execute`, `validate_not_none`, `validate_type`

### Progress (`progress`)

- `progress_bar` wraps `tqdm` when it is installed and falls back to a minimal logging counter otherwise
- `task_context` and `log_progress` for long-running operations

### Validation (`validation`)

- Validators raise `ValidationError`: `validate_type`, `validate_range`, `validate_not_none`, `validate_not_empty`, `validate_path_exists`, `validate_path_is_file`, `validate_path_is_dir`, `validate_path_within`
- `validate_schema` checks a dict against a field-to-type mapping; `validate_json_schema` validates against a JSON Schema file and needs `jsonschema`
- `validator` turns a predicate into a validator

### Workflow (`workflow`)

- `run_config_based_workflow(config_path)`, `validate_config_file`, `create_sample_config` and `download_and_process_data`
- `BaseWorkflowOrchestrator` and `WorkflowStep`: dependency-ordered execution with cycle detection; dependents of a failed step are skipped

### Database (`db`)

- Optional PostgreSQL helpers (`psycopg2` is imported defensively): `PostgresConnection` with a connection pool, `execute_query`, `execute_transaction`, `bulk_insert`, `create_table` and `table_exists`
- `get_db_client`, the `get_connection` context manager and `build_postgres_url`
- `bulk_insert` and `create_table` reject table or column names that are not simple SQL identifiers; like any other failure inside them, the rejection surfaces as `RuntimeError` (with the `ValueError` as its cause)
- `metainformant.core.db` re-exports `metainformant.core.data.db`

### Disk (`disk`)

- Disk space checks: `get_disk_usage`, `get_free_space`, `check_disk_space`, `ensure_disk_space`, `monitor_disk_space`
- Cleanup and sizing: `cleanup_temp_files`, `cleanup_old_files`, `safe_remove_directory`, `get_directory_size`, `get_largest_files`
- `get_recommended_temp_dir` prefers the location with more free space; `get_recommended_batch_size` derives a batch size from the drive's size category and the per-sample size (`detect_drive_size_category`)

### Hash (`hash`)

- `sha256_bytes`, `sha256_string`, `sha256_file` and `md5_file` return hex digests
- `hash_directory`, `file_hash_comparison` and `verify_file_integrity` compare content; `deterministic_seed` derives an integer seed from a string
- `core.io.checksums` adds `compute_md5`, `compute_sha256`, `verify_checksum` and checksum-file helpers

### Discovery (`discovery`)

- Repository introspection from parsed source: `discover_functions`, `discover_configs`, `discover_output_patterns`, `discover_workflows`
- `build_call_graph`, `find_symbol_usage` and `get_module_dependencies`

### Symbols (`symbols`)

- `index_functions` and `index_classes` build a symbol index of a source tree (cached under `output/.discovery_cache/`)
- `find_symbol`, `fuzzy_find_symbol`, `find_symbol_references`, `get_symbol_signature` and `get_symbol_metadata` query it

### Text (`text`)

- Whitespace and filename helpers: `normalize_whitespace`, `clean_whitespace`, `remove_control_chars`, `slugify`, `safe_filename`, `truncate_text`
- Biological-identifier helpers: `standardize_gene_name`, `format_species_name`, `clean_sequence_id`
- Extraction helpers: `extract_numbers`, `extract_email_addresses`, `count_words`

## Patterns

### Defensive Imports

```python-snippet
try:
    import optional_library
except ImportError:
    optional_library = None  # type: ignore

def function():
    if optional_library is None:
        raise ImportError("optional_library required")
    # Use library
```

### Pure Functions

- Prefer pure functions with minimal side effects
- Functions should be testable in isolation
- Side effects limited to I/O operations

### Error Handling

- Use `metainformant.core.utils.errors` types
- `ValidationError` is a `ValueError` subclass, so callers that catch `ValueError` also catch it
- Provide clear, actionable error messages
- Include biological context when relevant

### Type Hints

- Comprehensive type hints throughout
- Use `from __future__ import annotations` for forward references
- Union types: `str | Path | None`

## Testing

- Core utilities must have comprehensive tests
- **STRICTLY REAL IMPLEMENTATION**: Test real implementations only (see Global Project Rules (root `.cursorrules`) real-implementation policy)
- Use `tmp_path` fixture for test outputs
- Write test artifacts to `output/` only
- No test-double modules, test-double classes, or runtime function replacement
- Use real implementations with graceful skipping for unavailable dependencies

## Documentation

- Each submodule should have clear docstrings
- Examples in docstrings should use `output/` paths
- Document optional dependencies and their availability
- Reference module-specific rule files for domain patterns

## Integration with Other Modules

Most domain modules build on `core` utilities (`menu` and `cloud` do not import it):

- Import core utilities: `from metainformant.core import io`, `from metainformant.core.io import paths`, `from metainformant.core.utils import config`, and `from metainformant.core.utils import logging`
- Use core utilities for file I/O, logging, and path expansion and containment checks
- Do not reimplement core functionality in domain modules
- Use `metainformant.core.io` for domain data file I/O instead of direct `json` or `csv` imports (exceptions as listed under I/O)
- Build paths with `pathlib.Path`, and use `metainformant.core.io.paths` (`expand_and_resolve`, `is_within`) to expand user-supplied paths and check containment

## Reference

See Global Project Rules (root `.cursorrules`) for:

- Common directory structure and path handling
- Configuration patterns with env overrides
- Testing policy (real-implementation policy)
- Import patterns and code style
- Documentation guidelines
