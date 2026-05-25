# Core: Input/Output Operations

The `io` module is the primary interface for all file and network operations in METAINFORMANT. It provides a unified, high-level API for reading and writing JSON, YAML, CSV, TSV, Parquet, and JSON Lines formats, with transparent gzip compression support, atomic file operations, and robust error handling.

## Purpose

Bioinformatics pipelines process diverse data formats from various sources (local files, HTTP, FTP, databases). The `io` module centralizes all I/O logic to ensure:
- Consistent error handling and logging
- Atomic writes to prevent corruption
- Transparent compression handling
- Progress tracking for long-running downloads
- Type-safe operations with clear contracts

## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                     User Code / Domain Modules                   │
│  from metainformant.core import io                               │
│  data = io.load_json("file.json")                                │
└────────────────────────────┬────────────────────────────────────┘
                             │
                ┌────────────┼────────────┐
                │            │            │
                ▼            ▼            ▼
        ┌─────────────┐ ┌─────────────┐ ┌─────────────┐
        │ io.py       │ │ download.py │ │ paths.py    │
        │ Core I/O    │ │ HTTP/FTP    │ │ Path utils  │
        └──────┬──────┘ └──────┬──────┘ └──────┬──────┘
               │               │               │
        ┌──────┴───────────────┴───────────────┴──────┐
        │                                             │
        ▼                                             ▼
┌──────────────┐                              ┌──────────────┐
│ atomic.py    │                              │ cache.py     │
│ Atomic write │                              │ JSON cache   │
│ + safe rename│                              │ with TTL     │
└──────────────┘                              └──────────────┘
```

## Core Philosophy

### 1. **Everything is a Path**
All functions accept `str | Path`. Use `pathlib.Path` throughout.

### 2. **Gzip Transparent**
Files ending in `.gz` are automatically compressed/decompressed. No special API.

### 3. **Atomic Writes by Default**
`dump_json(..., atomic=True)` writes to temp file then renames to prevent partial writes on crashes.

### 4. **Fail Early, Log Clearly**
Errors raise specific exceptions (`IOError`, `ConfigError`) with context. Never silent failures.

### 5. **Stream When Possible**
Functions like `read_jsonl()` and `read_delimited()` are generators, streaming data without loading entire file into memory.

## API Reference — io.py Core Functions

### File Discovery and Preparation

#### `ensure_directory(path: str | Path) -> Path`

Create directory (and parents) if missing.

**Parameters**:
- `path`: Directory path to create

**Returns**: `Path` object for convenience

**Raises**: `OSError` if directory creation fails (e.g., permission denied)

**Example**:
```python
output_dir = io.ensure_directory("output/analysis")
# Creates: output/analysis/ (if missing)
```

**Note**: Does nothing if directory already exists.

### Text File I/O

#### `open_text_auto(path: str | Path, mode: str = "rt", encoding: str = "utf-8") -> io.TextIOBase`

Open text file with automatic gzip detection based on `.gz` suffix.

**Parameters**:
- `path`: File path
- `mode`: Text mode (`"r"`, `"w"`, `"a"` with optional `"t"`). Binary modes (`"rb"`, `"wb"`) are **rejected**.
- `encoding`: Text encoding (default: `"utf-8"`)

**Returns**: Text file handle (like built-in `open()`)

**Raises**: `ValueError` if binary mode requested

**Examples**:
```python
# Read plain text
with io.open_text_auto("notes.txt") as fh:
    text = fh.read()

# Read gzipped
with io.open_text_auto("data.txt.gz") as fh:
    text = fh.read()  # Auto-decompressed

# Write gzipped
with io.open_text_auto("output.json.gz", mode="wt") as fh:
    json.dump(data, fh)
```

### JSON I/O

#### `load_json(path: str | Path) -> Any`

Load and parse JSON file.

**Parameters**:
- `path`: JSON file path (`.json` or `.json.gz`)

**Returns**: Parsed Python object (dict, list, or primitive)

**Raises**:
- `FileNotFoundError` if file missing
- `IOError` (aliased from `core.utils.errors.IOError`) on read error or JSON parse failure

**Example**:
```python
config = io.load_json("config/workflow.yaml")  # Wait, YAML? No! Use .yaml extension with load_yaml()
config = io.load_json("config/workflow.json")
```

**Note**: For YAML, use `load_yaml()`.

#### `dump_json(obj: Any, path: str | Path, *, indent: int | None = None, atomic: bool = True) -> None`

Write object to JSON file.

**Parameters**:
- `obj`: JSON-serializable object
- `path`: Output path (`.json` or `.json.gz` for gzip)
- `indent`: Indentation spaces for pretty-printing (default: compact)
- `atomic`: If `True` (default), write to temp file then rename

**Raises**: `IOError` on write failure

**Atomic write behavior**:
1. Write to `path + ".tmp"` (or `.json.gz.tmp` for gzipped)
2. `fsync()` if available
3. `os.replace(tmp, final)` — atomic on POSIX

This ensures readers never see partial writes.

**Example**:
```python
# Compact JSON (default)
io.dump_json(data, "results.json")

# Pretty-printed
io.dump_json(data, "results.json", indent=2)

# Gzipped
io.dump_json(data, "results.json.gz")

# Non-atomic (for temporary files, beware corruption risk)
io.dump_json(data, "temp.json", atomic=False)
```

#### `load_json_gz(path: str | Path) -> Any`
#### `dump_json_gz(obj: Any, path: str | Path, *, indent: int | None = None) -> None`

Convenience wrappers for gzipped JSON. Equivalent to using `.gz` suffix with `load_json`/`dump_json`.

## YAML and TOML

#### `load_yaml(path: str | Path) -> Any`

Load YAML file.

**Requires**: `PyYAML` (`pip install pyyaml`)

**Raises**: `ImportError` if unavailable

#### `load_toml(path: str | Path) -> Any`

Load TOML file.

**Requires**: Python 3.11+ (`tomllib` stdlib) or `tomli` backport on older Python.

**Raises**: `ImportError` if unavailable

#### `dump_yaml(obj: Any, path: str | Path) -> None`

Write YAML file.

**Requires**: `PyYAML`

**Fallback**: If `yaml` unavailable, writes JSON with `.yaml` extension (not ideal but preserves data).

**Note**: YAML dumping not heavily used in core (configs are read, not written). Consider using if generating config templates.

## Delimited Text (CSV/TSV)

#### `read_delimited(path: str | Path, *, delimiter: str = ",") -> Iterator[dict[str, str]]`

Read CSV/TSV file as a generator of row dictionaries.

**Parameters**:
- `path`: File path (`.csv`, `.tsv`, `.csv.gz`, `.tsv.gz`)
- `delimiter`: Field separator (`,` for CSV, `\t` for TSV)

**Yields**: `dict` mapping column names → values (all strings)

**Raises**:
- `FileNotFoundError` if missing
- `IOError` on read/parse failure

**Example**:
```python
# Read CSV
for row in io.read_delimited("samples.csv"):
    print(f"Sample: {row['sample_id']}, Condition: {row['condition']}")

# Read TSV
for row in io.read_delimited("genes.tsv", delimiter="\t"):
    print(row["gene_id"], row["gene_name"])
```

**Memory efficient**: Streams rows; suitable for large files.

#### `write_delimited(rows: Iterable[Mapping[str, Any]], path: str | Path, *, delimiter: str = ",", atomic: bool = True) -> None`

Write iterable of row dictionaries to CSV/TSV.

**Parameters**:
- `rows`: Iterable of `dict` (or Mapping) objects
- `path`: Output file path
- `delimiter`: Field separator
- `atomic`: Use atomic write (default: True)

**Raises**: `IOError` on write failure

**Example**:
```python
rows = [
    {"name": "Alice", "score": 95},
    {"name": "Bob", "score": 87},
]
io.write_delimited(rows, "results.csv")

# TSV
io.write_delimited(rows, "results.tsv", delimiter="\t")
```

**Note**: Fieldnames (column order) inferred from first row's keys.

#### `read_csv(path: str | Path, **kwargs) -> Any`
#### `write_csv(data: Any, path: str | Path, **kwargs) -> None`

Pandas-accelerated CSV I/O. If `pandas` is installed, delegated to `pd.read_csv`/`DataFrame.to_csv`. Otherwise falls back to `read_delimited`/`write_delimited`.

**Parameters**: Pass-through to pandas (e.g., `dtype`, `usecols`, `index=False`).

**Returns**: `pandas.DataFrame` if pandas available, else `dict` of columns (fallback).

#### `read_tsv(path: str | Path) -> list[list[str]]`
#### `write_tsv(data, path: str | Path) -> None`

Simple TSV readers returning list-of-lists (no header parsing). Less commonly used; prefer `read_delimited(..., delimiter="\t")`.

## Parquet Support

#### `read_parquet(path: str | Path, **kwargs) -> Any`

Read Parquet file using pandas.

**Requires**: `pandas` + `pyarrow` or `fastparquet`

**Raises**: `ImportError` if pandas unavailable; `ImportError` with install hint if parquet engine missing.

**Example**:
```python
df = io.read_parquet("data.parquet")
# Returns pandas DataFrame
```

#### `write_parquet(df: Any, path: str | Path, **kwargs) -> None`

Write DataFrame to Parquet.

**Parameters**:
- `df`: pandas DataFrame or DataFrame-like object
- `path`: Output `.parquet` file
- `**kwargs`: Passed to `DataFrame.to_parquet()` (e.g., `compression="snappy"`)

## JSON Lines (JSONL)

Newline-delimited JSON, one object per line. Ideal for streaming large datasets.

#### `read_jsonl(path: str | Path) -> Iterator[dict[str, Any]]`

Stream JSONL file.

**Yields**: Parsed `dict` per line

**Example**:
```python
total_records = 0
for record in io.read_jsonl("logs.jsonl"):
    process(record)
    total_records += 1
print(f"Processed {total_records} records")
```

#### `write_jsonl(rows: Iterable[Mapping[str, Any]], path: str | Path, *, atomic: bool = True) -> None`

Write iterable of dicts as JSONL.

**Parameters**:
- `rows`: Iterable of mappings (dict, namedtuple with `_asdict()`)
- `atomic`: Use atomic write (temp + rename)

**Example**:
```python
records = [{"id": 1, "value": "a"}, {"id": 2, "value": "b"}]
io.write_jsonl(records, "output.jsonl")
```

**Atomic**: Writes to `output.jsonl.tmp` then renames to prevent partial files.

## Download Functions

High-level download convenience wrappers. For robust, resumable downloads with progress/retry, see [Download](./download.md) module.

#### `download_file(url: str, dest_path: str | Path, *, chunk_size: int = 8192, timeout: int = 30) -> bool`

Download a file to local path.

**Parameters**:
- `url`: Source URL (http, https, ftp, file)
- `dest_path`: Local destination (file or directory; auto-names if directory)
- `chunk_size`: Streaming chunk size in bytes (default: 8192 = 8 KB)
- `timeout`: Request timeout seconds (default: 30)

**Returns**: `True` if successful, `False` otherwise

**Behavior**:
- Creates parent directories if needed
- Resumes partial downloads if server supports `Range` header
- Retries with exponential backoff (via `download_with_progress` under the hood)
- Logs errors (no exception raised)

**Example**:
```python
success = io.download_file(
    "https://example.com/genome.fa.gz",
    "data/genome.fa.gz",
    chunk_size=1024*1024,  # 1 MB chunks
    timeout=60
)
if success:
    print("Download complete")
else:
    print("Download failed")
```

#### `download_json(url: str, *, timeout: int = 30) -> Any | None`

Download and parse JSON in one call.

**Parameters**:
- `url`: JSON endpoint URL
- `timeout`: Request timeout

**Returns**: Parsed JSON object or `None` on failure

**Example**:
```python
data = io.download_json("https://api.example.com/v1/samples")
if data:
    for sample in data["samples"]:
        print(sample["id"])
```

#### `download_text(url: str, *, timeout: int = 30) -> str | None`

Download raw text content.

**Returns**: String or `None`

#### `download_csv(url: str, *, timeout: int = 30, **kwargs) -> Any | None`

Download CSV and parse via pandas (if available) or fallback. Returns DataFrame or `dict` of columns.

**Additional `**kwargs`**: Passed to `pandas.read_csv()` (e.g., `sep="\t"` for TSV).

#### `batch_download(urls: list[str], dest_dir: str | Path, *, timeout: int = 30) -> dict[str, bool]`

Download multiple files in batch.

**Parameters**:
- `urls`: List of URLs
- `dest_dir`: Directory to save files (filenames auto-derived from URLs)
- `timeout`: Per-request timeout

**Returns**: `{url: success_bool}` mapping

**Example**:
```python
urls = [
    "https://example.com/file1.csv",
    "https://example.com/file2.csv",
]
results = io.batch_download(urls, "data/downloads")
for url, ok in results.items():
    status = "OK" if ok else "FAILED"
    print(f"{url}: {status}")
```

## Path Utilities (re-exported from `io.paths`)

The `io` module re-exports path utilities for convenience:

| Re-export from `paths` | Description |
|------------------------|-------------|
| `expand_and_resolve` | Expand `~` and resolve symlinks to absolute path |
| `is_within` | Check if path is inside parent directory |
| `sanitize_filename` | Remove unsafe characters for filesystem |
| `create_temp_file` | Generate unique temporary file path |
| `find_files_by_extension` | Recursively find files by extension |
| `get_file_size` | Get file size in bytes (0 if missing) |
| `get_directory_size` | Total size of all files in directory |

See [Path Handling](./paths.md) for full documentation.

## Atomic Operations (`atomic.py`)

While `dump_json(atomic=True)` covers most use cases, `atomic.py` provides lower-level primitives:

#### `atomic_write(path: str | Path, content: bytes | str) -> None`

Write content to file atomically.

**Implementation**: Writes to `path + ".tmp"`, `fsync()`, then `os.replace()`.

#### `atomic_write_json(path: str | Path, data: Any, indent: int | None = None) -> None`

Convenience wrapper: serialize to JSON then atomic write.

## Caching Integration (`cache.py`)

The `io` module uses `cache` for its own internal download retry state (heartbeat files), but you can also use caching with I/O:

```python
from metainformant.core import io
from metainformant.core.io import cache

# Cache downloaded file contents
def load_or_download(url, cache_dir):
    key = f"download/{hash(url)}"
    data = cache.load_cached_json(cache_dir, key)
    if data:
        return data

    # Download and cache
    raw = io.download_text(url)
    cache.cache_json(cache_dir, key, raw)
    return raw
```

## Error Handling

All `io` functions raise exceptions from `metainformant.core.utils.errors`:

- `IOError` (distinct from built-in `IOError`, but similar interface): File/network errors
- `ConfigError`: Configuration-related failures (e.g., unsupported format)

**Pattern**:
```python
from metainformant.core.io import errors as io_errors

try:
    data = io.load_json("missing.json")
except io_errors.IOError as e:
    print(f"I/O error: {e}")
except FileNotFoundError as e:
    print(f"File not found: {e}")  # More specific
```

## Compression Handling

All functions that read/write files support transparent gzip via suffix check:

| Suffex | Behavior |
|--------|----------|
| `.gz` | Open with `gzip.open()` (compression level 6 by default) |
| No suffix | Regular `open()` |

**Autodetection**: No flags needed. The `open_text_auto()` function internally checks `Path.suffix`.

**Example**:
```python
# Read/write plain JSON
io.dump_json(data, "data.json")
d1 = io.load_json("data.json")

# Read/write compressed JSON (4–10× smaller for text data)
io.dump_json(data, "data.json.gz")
d2 = io.load_json("data.json.gz")  # Same API
assert d1 == d2
```

**Performance**:
- Compression adds ~50–200 ms for 10 MB JSON
- Decompression adds ~20–100 ms
- Trade-off: disk space vs CPU time

## Performance Tips

### 1. Stream Large Files

Don't load entire file if you can process row-by-row:

```python
# GOOD: Stream
for record in io.read_jsonl("big.jsonl"):
    process(record)  # One record at a time

# BAD: Load everything
records = list(io.read_jsonl("big.jsonl"))  # May OOM
```

### 2. Use Gzip for I/O-bound Pipelines

If network or disk is bottleneck, compress intermediate files:
```python
io.dump_json(large_data, "intermediate.json.gz")  # 90% size reduction
```

### 3. Batch Database Writes

If writing many rows to CSV for database bulk load:
```python
# Write in chunks to avoid holding all in memory
chunks = [records[i:i+10000] for i in range(0, len(records), 10000)]
for chunk_num, chunk in enumerate(chunks):
    io.write_delimited(chunk, f"batch_{chunk_num:04d}.csv")
```

### 4. Choose Right Format

| Format | Read Speed | Write Speed | Compression | Human-readable |
|--------|------------|-------------|-------------|----------------|
| JSON | Fast | Fast | Poor (text) | Yes |
| JSONL | Fast (stream) | Fast | Poor | Yes |
| CSV | Fast (pandas) | Fast | Poor | Yes |
| Parquet | Very fast | Very fast | Excellent (Snappy) | No |
| JSON.gz | Medium | Medium | Good | Yes (after decompress) |
| Parquet + Snappy | Fastest | Fastest | Excellent | No |

**Guideline**:
- Intermediate data between pipeline steps → **JSONL** (streamable, simple)
- Final tabular data for analysis → **Parquet** (fast random access, columnar)
- Config files → **YAML** (human-friendly)
- Archival storage → **Parquet + Snappy** or **JSON.gz**

### 5. Reuse Connections

`download_file()` opens a new connection per call. For many small downloads:
```python
# Use requests.Session() directly for connection pooling
import requests
session = requests.Session()
for url in urls:
    r = session.get(url)
    # ...
session.close()
```

## Dependencies

### Required
- Python stdlib: `json`, `csv`, `gzip`, `io`, `pathlib`, `urllib.parse`

### Optional
- `requests`: HTTP downloads. If missing, `download_json()` and `download_text()` fail.
- `pandas`: `read_csv`, `write_csv`, `read_parquet`, `write_parquet`. Falls back to native CSV if missing.
- `pyarrow` or `fastparquet`: Required for Parquet support (pandas uses one of these).
- `PyYAML`: Full YAML support. Without it: `.yaml` files may fail unless simple enough for built-in parser.
- `tomllib` (Python 3.11+) or `tomli`: TOML support.

## Common Pitfalls

### Pitfall 1: Forgetting to Create Parent Directories

**Symptom**: `FileNotFoundError: [Errno 2] No such file or directory: 'output/results.json'`

**Cause**: Parent directory `output/` doesn't exist.

**Fix**: Use `io.ensure_directory()` or `paths.ensure_directory()`:
```python
output_path = Path("output/results.json")
io.ensure_directory(output_path.parent)  # Creates "output/"
io.dump_json(data, output_path)
```

### Pitfall 2: Loading Large Files into Memory

**Symptom**: MemoryError or system slowdown.

**Fix**: Use streaming iterators:
```python
# Process 1M JSONL records without loading all into RAM
count = 0
for record in io.read_jsonl("huge.jsonl"):
    process(record)  # Process one at a time
    count += 1
```

### Pitfall 3: Corrupted JSON from Partial Writes

**Symptom**: `JSONDecodeError` reading a file you just wrote.

**Cause**: Non-atomic write interrupted (power loss, kill -9).

**Fix**: Always use `atomic=True` (default). If you disabled it for performance, accept risk.

### Pitfall 4: Wrong Delimiter for TSV

**Symptom**: CSV parser reads TSV as single column.

**Cause**: Using `read_delimited()` with default `delimiter=","` on TSV.

**Fix**: Specify `delimiter="\t"` or use `read_tsv()`:
```python
for row in io.read_delimited("data.tsv", delimiter="\t"):
    print(row)

# OR
rows = io.read_tsv("data.tsv")  # Returns list of lists (no header parsing)
```

### Pitfall 5: Silent Fallback to JSON When YAML Expected

**Symptom**: YAML config loaded but anchors/aliases missing.

**Cause**: `dump_yaml()` fell back to JSON because `yaml` module missing.

**Fix**: Install PyYAML. Check logs for warnings.

### Pitfall 6: Download Timeouts on Large Files

**Symptom**: Large file downloads fail with timeout.

**Fix**: Increase timeout:
```python
io.download_file(url, dest, timeout=300)  # 5 minutes
```

Consider using `download_robust.py`'s `Downloader` class for resumable multi-part downloads.

## Testing Considerations

- Use temporary paths (`tmp_path` in pytest)
- Test with both `.json` and `.json.gz`
- Simulate corrupted files to test error handling
- For download tests: use `httpbin.org` or local `http.server`

## Related Components

| Component | How It Relates |
|-----------|----------------|
| `core.io.cache` | Stores downloaded data for reuse |
| `core.io.paths` | Path manipulation before I/O |
| `core.download` | Robust download with retry, resume, heartbeat |
| `core.utils.logging` | I/O errors logged consistently |
| `core.utils.hash` | Verify downloaded file integrity |
| `core.utils.text` | Text cleaning for CSV/TSV parsing |

## Examples in Repository

See `examples/core/`:
- `example_io.py`: Demonstrates all I/O functions
- `example_config.py`: Loads YAML configs
- `example_logging.py`: Structured logging with I/O

## Next Steps

- Read [Download](./download.md) for robust network transfers
- Read [Paths](./paths.md) for path manipulation best practices
- Read [Caching](./cache.md) to avoid repeated downloads
- Explore [API Reference](./SPEC.md) for full type signatures
