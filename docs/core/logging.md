# Core: Structured Logging

The `logging` module provides consistent, structured logging across all METAINFORMANT modules. It builds on Python's standard `logging` module with sensible defaults, environment-based configuration, and support for metadata-rich log entries.

## Purpose

Logging is essential for:
- **Debugging**: Understanding what went wrong
- **Auditing**: Tracking pipeline execution and data lineage
- **Monitoring**: Detecting issues in production runs
- **Reproducibility**: Recording software versions, parameters, timestamps

The `core.utils.logging` module standardizes logging format, handlers, and log levels across all domain modules.

## Design Principles

### 1. **Get Logger by Module Name**
Always use `get_logger(__name__)` to inherit module hierarchy in logs:
```
metainformant.rna.amalgkit.quantifier  # Clear origin
```

### 2. **Structured Format**
All logs use consistent format:
```
YYYY-MM-DD HH:MM:SS | LEVEL | module.name | message
```
Example:
```
2026-04-26 20:15:01 | INFO | metainformant.rna.pipeline | Starting quantification for sample ERR123456
```

### 3. **Environment-Driven Configuration**
Log level controlled via `CORE_LOG_LEVEL` environment variable. No code changes required.

### 4. **Metadata Support**
Structured metadata can be appended to log messages or logged as separate JSON entries for machine parsing.

### 5. **Graceful Degradation**
If console output fails (unlikely), logging degrades to no-op rather than crashing.

## Component Overview

```
┌────────────────────────────────────────────┐
│         Application / Domain Module        │
│  from metainformant.core.utils.logging     │
│         import get_logger                  │
│                                            │
│  logger = get_logger(__name__)             │
│  logger.info("Processing ...")             │
└─────────────────┬──────────────────────────┘
                  │
                  ▼
┌────────────────────────────────────────────┐
│           get_logger() factory             │
│  • Creates logger if not exists            │
│  • Attaches StreamHandler (console)        │
│  • Sets default level (INFO)               │
│  • Applies CORE_LOG_LEVEL override         │
└─────────────────┬──────────────────────────┘
                  │
                  ▼
┌────────────────────────────────────────────┐
│        Configured Logger (singleton)       │
│  • Formatter: timestamp | level | name     │
│  • Handlers: Console (always), File (opt)  │
│  • Level: from env or explicit             │
│  • Propagation: False (no double-log)      │
└────────────────────────────────────────────┘
```

## API Reference

### Core Functions

#### `get_logger(name: str) -> logging.Logger`

Get or create a logger with default console handler.

**Parameters**:
- `name`: Logger name (typically `__name__`)

**Returns**: `logging.Logger` instance

**Behavior**:
- Checks `logging.getLogger(name)`. If no handlers attached, adds a `StreamHandler` with standard formatter.
- Sets level to `logging.INFO` by default (overridden by `CORE_LOG_LEVEL` or explicit `get_logger_with_level()`).
- Subsequent calls return same logger (singleton per name).

**Example**:
```python
from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)
logger.info("Module initialized")
logger.debug("Detailed trace (won't show unless CORE_LOG_LEVEL=DEBUG)")
```

#### `get_logger_with_level(name: str, level: str | int | None = None) -> logging.Logger`

Get logger with explicit or environment-based log level.

**Parameters**:
- `name`: Logger name
- `level`: Level as string (`"DEBUG"`, `"INFO"`, `"WARNING"`, `"ERROR"`, `"CRITICAL"`) or integer (e.g., `logging.DEBUG`). If `None`, uses `CORE_LOG_LEVEL` env var or `INFO`.

**Returns**: Configured logger

**Example**:
```python
# Explicit level
debug_logger = get_logger_with_level("my_module", "DEBUG")

# Environment-controlled
# $ export CORE_LOG_LEVEL=DEBUG
logger = get_logger_with_level("my_module")  # Will be DEBUG
```

#### `setup_logger(name: str, log_file: str | None = None, level: str = "INFO") -> logging.Logger`

Set up a logger with file output in addition to console.

**Parameters**:
- `name`: Logger name
- `log_file`: Optional file path for log file (directory created if missing)
- `level`: Logging level (default: `"INFO"`)

**Returns**: Configured logger with console + optional file handler

**Example**:
```python
logger = setup_logger(
    "pipeline",
    log_file="logs/pipeline_20260426.log",
    level="DEBUG"
)
logger.info("Pipeline started with file logging")
```

**Note**: This is typically called once at application startup. For per-module loggers, use `get_logger()`.

#### `configure_logging_from_env(default_level: str = "INFO") -> None`

Configure **root logger** from `CORE_LOG_LEVEL`. Affects all loggers that haven't been explicitly configured.

**Parameters**:
- `default_level`: Fallback if `CORE_LOG_LEVEL` unset

**Example**:
```python
from metainformant.core.utils.logging import configure_logging_from_env

# Call early in main()
configure_logging_from_env()

# Now all loggers use CORE_LOG_LEVEL or default
```

### Structured Metadata Logging

#### `log_with_metadata(logger: logging.Logger, message: str, metadata: dict, *, level: str = "INFO", structured: bool = False) -> None`

Attach structured data to log entries.

**Parameters**:
- `logger`: Logger instance
- `message`: Human-readable message
- `metadata`: Dictionary of key-value pairs (serializable)
- `level`: Log level (`"DEBUG"`, `"INFO"`, etc.)
- `structured`: If `True`, logs metadata as separate JSON line; if `False`, appends JSON to message.

**Examples**:

Unstructured (simple append):
```python
log_with_metadata(logger, "Sample processed", {
    "sample_id": "ERR123456",
    "genes": 20456,
    "time_sec": 45.2
})
# Log line: "2026-04-26 ... | INFO | ... | Sample processed | {'sample_id': 'ERR123456', ...}"
```

Structured (separate JSON entry):
```python
log_with_metadata(logger, "Batch completed", {
    "batch_id": "batch_001",
    "records": 1000,
    "success_rate": 0.99
}, structured=True)
# Two log lines:
#  "Batch completed"
#  "METADATA: {'batch_id': 'batch_001', ...}"
```

**Use structured** when feeding logs to systems like:
- ELK stack (Elasticsearch, Logstash, Kibana)
- Splunk
- Graylog
- JSON log parsers

## Supported Log Levels

Standard Python levels (increasing severity):

| Level | Numeric | When to Use |
|-------|---------|-------------|
| `DEBUG` | 10 | Fine-grained diagnostic info (variable values, SQL queries) |
| `INFO` | 20 | High-level progress (startup, completion, counts) |
| `WARNING` | 30 | Unexpected but non-critical (missing optional file, deprecated config) |
| `ERROR` | 40 | Operation failed but pipeline continues or fails gracefully |
| `CRITICAL` | 50 | Fatal error causing immediate shutdown |

**Default**: `INFO`

## Environment Variables

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `CORE_LOG_LEVEL` | `str` | `"INFO"` | Minimum log level to output. Case-insensitive. |
| `CORE_LOG_FORMAT` | `str` | `"default"` | Future: custom format strings (not yet implemented) |

**Examples**:
```bash
# Show all logs including DEBUG
export CORE_LOG_LEVEL=DEBUG

# Only warnings and errors
export CORE_LOG_LEVEL=WARNING

# In code (overrides env for specific logger)
logger = get_logger_with_level("my_module", "DEBUG")
```

## Format String

The default log format is:
```python
fmt="%(asctime)s | %(levelname)s | %(name)s | %(message)s"
datefmt="%Y-%m-%d %H:%M:%S"
```

Result:
```
2026-04-26 20:15:01 | INFO | metainformant.rna.pipeline | Starting sample ERR123456
```

**Customization**: For custom formats, use standard Python logging:
```python
import logging
from metainformant.core.utils.logging import get_logger

logger = get_logger("custom")
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(levelname)s:%(name)s:%(message)s"))
logger.addHandler(handler)
```

## Log Handler Architecture

### Console Handler (Always Attached)

Every logger from `get_logger()` gets a `StreamHandler` (stderr) with our formatter. This ensures logs appear without configuration.

### File Handler (Optional via `setup_logger`)

`setup_logger()` adds a `FileHandler` writing to specified path. Directory auto-created.

**Rotation**: Not built-in. For long-running pipelines, use external logrotate or implement `RotatingFileHandler`:
```python
from logging.handlers import RotatingFileHandler

handler = RotatingFileHandler("pipeline.log", maxBytes=10*1024*1024, backupCount=5)
logger.addHandler(handler)
```

### No Propagation

Loggers have `propagate = False` to prevent duplicate logs from ancestor loggers. Each module's logs appear only once.

## Usage Examples

### 1. Basic Module Logging

```python
from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

def process_sample(sample_id: str):
    logger.info(f"Processing sample {sample_id}")
    try:
        result = compute(sample_id)
        logger.debug(f"Sample {sample_id} result: {result}")
        logger.info(f"Completed sample {sample_id}")
        return result
    except Exception as e:
        logger.error(f"Sample {sample_id} failed: {e}", exc_info=True)
        raise
```

### 2. Conditional Debug Logging

Avoid expensive string formatting when debug disabled:
```python
logger = get_logger(__name__)

# BAD: Always formats even if DEBUG disabled
logger.debug(f"Large object: {expensive_to_string(huge_dict)}")

# GOOD: Guard with if
if logger.isEnabledFor(logging.DEBUG):
    logger.debug(f"Large object: {expensive_to_string(huge_dict)}")
```

### 3. Exception Logging

```python
try:
    run_analysis()
except Exception:
    logger.exception("Analysis failed")  # Includes traceback
    # Or:
    logger.error("Analysis failed", exc_info=True)
```

### 4. Metadata-Rich Logging

```python
from metainformant.core.utils.logging import log_with_metadata

log_with_metadata(
    logger,
    "Workflow step completed",
    {
        "step": "quantification",
        "sample": "ERR123456",
        "duration_sec": 45.3,
        "genes_detected": 20456,
        "mapping_rate": 0.87,
    },
    level="INFO"
)
```

JSON output (if structured=True):
```
METADATA: {"step": "quantification", "sample": "ERR123456", ...}
```

### 5. Pipeline Progress Tracking

```python
logger = get_logger("pipeline")

def process_batch(batch_id, samples):
    logger.info(f"Starting batch {batch_id} with {len(samples)} samples")
    for i, sample in enumerate(samples):
        process_sample(sample)
        if (i + 1) % 100 == 0:
            logger.info(f"Batch {batch_id}: {i+1}/{len(samples)} samples done")
    logger.info(f"Batch {batch_id} complete")
```

### 6. Component Initialization Logging

```python
class Quantifier:
    def __init__(self, index_path, threads=4):
        self.logger = get_logger(self.__class__.__name__)
        self.logger.info(f"Initializing Quantifier with index={index_path}, threads={threads}")
        # ...
```

### 7. Conditional File Logging

```python
import sys
from datetime import datetime

def main():
    # Setup root logger from env
    configure_logging_from_env()

    # Create timestamped log file
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = f"logs/pipeline_{timestamp}.log"

    # Get module-specific logger with file output
    logger = setup_logger("pipeline", log_file=log_file, level="DEBUG")

    logger.info("Pipeline started")
    # ...
```

## Best Practices

### Do's

✅ **Do** use `get_logger(__name__)` in every module.

✅ **Do** log at appropriate levels:
- `DEBUG` — internal state, variable values, SQL, API calls
- `INFO` — high-level progress, start/stop, counts
- `WARNING` — recoverable issues (missing optional file)
- `ERROR` — operation failures that are handled or cause graceful exit
- `CRITICAL` — fatal errors (out of disk space, database down)

✅ **Do** include relevant context: sample IDs, file paths, step names.

✅ **Do** use `exc_info=True` or `logger.exception()` for caught exceptions.

✅ **Do** use structured metadata for monitoring/aggregation.

### Don'ts

❌ **Don't** log sensitive data (passwords, PII) in plain text.

❌ **Don't** log huge objects (use summarized form):
```python
# BAD
logger.debug(f"Full response: {response.text}")  # Could be MBs

# GOOD
logger.debug(f"Response status: {response.status_code}, len={len(response.text)}")
```

❌ **Don't** create loggers inside tight loops (cache at module/class level).

❌ **Don't** change global logging config from library code (let application control).

❌ **Don't** use `print()` for production logs (use logger).

## Testing Logs

For test suites, you can capture log output:

```python
import logging
from metainformant.core.utils.logging import get_logger

def test_something(caplog):
    logger = get_logger("test_module")
    with caplog.at_level(logging.INFO):
        run_function()
        assert "Expected log message" in caplog.text
```

Or use `logging.handlers.MemoryHandler` to buffer logs in-memory.

## Performance Impact

Logging overhead is minimal when disabled:
- `logger.debug("...")` — string formatted **always**, but handler skips if level too high
- Use `if logger.isEnabledFor(logging.DEBUG):` guard for expensive string interpolation
- Metadata dict creation also has cost; guard with level check if expensive to build

**Benchmark** (simple call, no handlers):
```
logger.debug("msg")        # ~0.5 µs
logger.info("msg")         # ~0.5 µs
isEnabledFor(DEBUG) check  # ~0.1 µs
```

Negligible compared to I/O operations.

## Troubleshooting

### Issue: No Logs Appear

**Check**:
1. Is `CORE_LOG_LEVEL` set too high (e.g., `ERROR`)? Lower to `INFO` or `DEBUG`.
2. Did you call `get_logger(__name__)`? Check `logger.handlers` is non-empty.
3. Are logs going to stderr? May be captured by systemd or container. Check with `journalctl` or `docker logs`.

**Fix**: Explicitly attach a handler:
```python
logger = get_logger(__name__)
logger.addHandler(logging.StreamHandler())  # Force console
```

### Issue: Duplicate Log Lines

**Cause**: Logger has multiple handlers or propagation to ancestor logger also has handler.

**Fix**:
```python
logger.propagate = False  # Prevent ancestor logging
logger.handlers.clear()  # Remove duplicates, then add your own
```

### Issue: File Handler Not Writing

**Cause**: Directory doesn't exist or permissions.

**Fix**:
```python
Path("logs").mkdir(exist_ok=True)  # Ensure directory
```

### Issue: Log Level Not Changing

**Cause**: `CORE_LOG_LEVEL` set after logger creation, or logger level explicitly set to override.

**Fix**: Call `configure_logging_from_env()` before obtaining loggers, or set logger level directly:
```python
logger.setLevel(logging.DEBUG)
```

## Related Components

| Module | Relationship |
|--------|--------------|
| `core.io` | I/O errors logged via this module |
| `core.cache` | Cache hits/misses logged |
| `core.download` | Download progress logged |
| `core.parallel` | Worker task completion logged |
| `core.workflow` | Step status changes logged |

## Dependencies

- **Required**: Python stdlib `logging`, `json`
- **Optional**: None

## Further Reading

- Python logging HOWTO: https://docs.python.org/3/howto/logging.html
- Structlog (alternative): https://www.structlog.org/ (not used here, but similar philosophy)
