# Core: Input/Output Operations

The `io` module provides comprehensive file I/O utilities with support for multiple formats, compression, and network operations.

## Core Functions

### Directory Management
- **`ensure_directory(path)`** → `Path`
  - Create directory and parents if missing
  - Returns Path object for the created directory

### File Opening
- **`open_text_auto(path, mode="rt", encoding="utf-8")`** → `TextIOBase`
  - Open text files with automatic gzip handling
  - Supports .gz files transparently
  - Text modes only ("rt", "wt", "at")

## JSON Operations

### Basic JSON
- **`load_json(path)`** → `Any`
  - Load JSON from file (handles gzip automatically)
- **`dump_json(obj, path, indent=None)`** → `None`
  - Write JSON to file with optional indentation
  - Creates parent directories automatically
  - Sorts keys for consistent output

### Compressed JSON
- **`load_json_gz(path)`** → `Any`
  - Load from gzipped JSON file
- **`dump_json_gz(obj, path, indent=None)`** → `None`
  - Write to gzipped JSON file

## JSON Lines (JSONL)

- **`read_jsonl(path)`** → `Iterator[Dict[str, Any]]`
  - Read JSON Lines file line by line
  - Skips empty lines automatically
- **`write_jsonl(rows, path)`** → `None`
  - Write iterable of dicts as JSON Lines
  - Each row becomes one JSON line

## CSV/TSV Operations

### Pandas-Based (Recommended)
- **`read_csv(path, **kwargs)`** → `DataFrame`
  - Read CSV with pandas (handles gzip)
  - Requires pandas installation
- **`write_csv(df, path, **kwargs)`** → `None`
  - Write DataFrame to CSV
- **`read_parquet(path, **kwargs)`** → `DataFrame`
  - Read Parquet files with pandas
- **`write_parquet(df, path, **kwargs)`** → `None`
  - Write DataFrame to Parquet

### Native Implementation
- **`read_delimited(path, delimiter=",")`** → `Iterator[Dict[str, str]]`
  - Read CSV/TSV without pandas dependency
  - Returns dict iterator with string values
- **`write_delimited(rows, path, delimiter=",")`** → `None`
  - Write rows to CSV/TSV
  - Auto-detects fieldnames from first row

### TSV Specific
- **`read_tsv(path)`** → `List[List[str]]`
  - Read tab-separated values as list of lists
- **`write_tsv(data, path)`** → `None`
  - Write list of lists as TSV

## Network Operations

### File Downloads
- **`download_file(url, dest_path, chunk_size=8192, timeout=30)`** → `bool`
  - Download file from URL to local path
  - Supports automatic filename detection
  - Returns success status

### Data Downloads
- **`download_json(url, timeout=30)`** → `Any`
  - Download and parse JSON from URL
- **`download_text(url, timeout=30)`** → `str | None`
  - Download text content from URL
- **`download_csv(url, timeout=30, **kwargs)`** → `DataFrame | None`
  - Download and parse CSV from URL

### Batch Operations
- **`batch_download(urls, dest_dir, timeout=30)`** → `Dict[str, bool]`
  - Download multiple files in batch
  - Returns success status for each URL

## Usage Examples

### Basic File Operations
```python
from metainformant.core import io

# Create directory structure
work_dir = io.ensure_directory("output/analysis")

# JSON operations
data = {"experiment": "test", "results": [1, 2, 3]}
io.dump_json(data, "output/analysis/results.json")
loaded = io.load_json("output/analysis/results.json")

# Compressed JSON
io.dump_json_gz(data, "output/analysis/results.json.gz")
compressed = io.load_json_gz("output/analysis/results.json.gz")
```

### JSON Lines Processing
```python
from metainformant.core import io

# Write JSONL
records = [
    {"id": 1, "name": "sample1", "value": 10.5},
    {"id": 2, "name": "sample2", "value": 20.3},
]
io.write_jsonl(records, "output/samples.jsonl")

# Read JSONL
for record in io.read_jsonl("output/samples.jsonl"):
    print(f"Sample {record['name']}: {record['value']}")
```

### CSV Processing
```python
from metainformant.core import io

# Using pandas (recommended)
try:
    df = io.read_csv("data/samples.csv")
    io.write_csv(df, "output/filtered_samples.csv")
except ImportError:
    # Fallback to native implementation
    rows = list(io.read_delimited("data/samples.csv"))
    # Process rows...
    io.write_delimited(rows, "output/filtered_samples.csv")
```

### File Downloads
```python
from metainformant.core import io

# Download single file
success = io.download_file(
    "https://example.com/data.csv",
    "output/downloaded_data.csv"
)

# Download JSON data
api_data = io.download_json("https://api.example.com/data")
if api_data:
    io.dump_json(api_data, "output/api_response.json")

# Batch download
urls = [
    "https://example.com/file1.txt",
    "https://example.com/file2.txt",
]
results = io.batch_download(urls, "output/downloads")
print(f"Downloaded {sum(results.values())}/{len(results)} files")
```

### TSV Operations
```python
from metainformant.core import io

# Read TSV
data = io.read_tsv("data/matrix.tsv")

# Write TSV
matrix = [
    ["gene", "sample1", "sample2"],
    ["GENE1", "10.5", "12.3"],
    ["GENE2", "8.7", "9.1"],
]
io.write_tsv(matrix, "output/expression_matrix.tsv")
```

## Compression Support

All text-based operations automatically handle gzip compression:
- Files ending in `.gz` are automatically compressed/decompressed
- Works with JSON, JSONL, CSV, TSV, and delimited files
- No code changes required - just use `.gz` in filenames

```python
from metainformant.core import io

# Automatic gzip handling
io.dump_json(data, "output/large_data.json.gz")  # Compressed
loaded = io.load_json("output/large_data.json.gz")  # Decompressed automatically
```

## Error Handling

The io module provides robust error handling:
- Automatic directory creation for output files
- Graceful handling of missing pandas dependency
- Network timeout and error handling for downloads
- Clear error messages for file operations

## Dependencies

- **Required**: Standard library only
- **Optional**: pandas (for DataFrame operations), requests (for downloads)
- **Network**: requests library for download functions
