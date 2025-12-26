# Core: Path Management

The `paths` module provides comprehensive path handling utilities with safety checks, file operations, and directory management for secure and predictable file system operations.

## Core Path Functions

### Path Resolution and Validation
- **`expand_and_resolve(path)`** → `Path`
  - Expand user home (~) and environment variables
  - Resolve to absolute path without requiring file existence

- **`is_within(path, parent)`** → `bool`
  - Check if path is contained within parent directory
  - Resolves both paths before comparison
  - Safe containment checking for security

### Directory Operations
- **`ensure_directory(path)`** → `None`
  - Create directory and all parent directories
  - Safe to call on existing directories

- **`prepare_file_path(file_path)`** → `None`
  - Ensure parent directories exist for file path
  - Prepares directory structure for file creation

### Safety and Validation
- **`is_safe_path(path)`** → `bool`
  - Check for path traversal attempts and dangerous patterns
  - Prevents directory traversal and injection attacks

## File Operations

### File Information
- **`get_file_extension(filename)`** → `str`
  - Extract file extension including the dot
  - Returns empty string if no extension

- **`get_file_size(path)`** → `int`
  - Get file size in bytes
  - Returns 0 for non-existent files

- **`get_directory_size(path)`** → `int`
  - Calculate total size of all files in directory recursively

### File Name Manipulation
- **`change_extension(path, new_extension)`** → `Path`
  - Change file extension (with or without dot)
  - Returns new Path object

- **`sanitize_filename(filename)`** → `str`
  - Remove dangerous characters for filesystem safety
  - Handles control characters and reserved names

### File Discovery
- **`find_files_by_extension(directory, extension)`** → `List[Path]`
  - Find all files with specific extension recursively
  - Extension can include or exclude the dot

## Temporary Files

- **`create_temp_file(suffix="", prefix="tmp", directory=None)`** → `Path`
  - Create unique temporary file path
  - Uses secure hashing for uniqueness
  - Optional custom directory

## Usage Examples

### Basic Path Operations
```python
from metainformant.core import paths

# Expand and resolve paths
home_config = paths.expand_and_resolve("~/config/settings.yaml")
absolute_path = paths.expand_and_resolve("./relative/path")

# Safety checks
if paths.is_within(user_path, "/safe/directory"):
    print("Path is safe to use")
else:
    raise ValueError("Unsafe path detected")
```

### Directory Management
```python
from metainformant.core import paths

# Ensure directories exist
output_dir = paths.ensure_directory("output/analysis")
nested_dir = paths.ensure_directory("output/deep/nested/structure")

# Prepare for file creation
paths.prepare_file_path("output/results/final_report.txt")
# Now safe to write to the file
```

### File Operations
```python
from metainformant.core import paths

# File information
size = paths.get_file_size("large_dataset.fasta")
ext = paths.get_file_extension("sequence.fasta")  # ".fasta"

# Change extensions
json_path = paths.change_extension("data.csv", ".json")
compressed = paths.change_extension("results.txt", ".txt.gz")

# Find files
fasta_files = paths.find_files_by_extension("data", ".fasta")
all_csv = paths.find_files_by_extension(".", "csv")
```

### Safe Filename Handling
```python
from metainformant.core import paths

# Sanitize potentially dangerous filenames
safe_name = paths.sanitize_filename("dangerous<file>name?.txt")
# Becomes: "dangerous_file_name_.txt"

# Check path safety
if paths.is_safe_path(user_input_path):
    process_path(user_input_path)
else:
    raise ValueError("Unsafe path detected")
```

### Directory Size Analysis
```python
from metainformant.core import paths

# Get sizes for cleanup decisions
data_size = paths.get_directory_size("data/")
output_size = paths.get_directory_size("output/")

print(f"Data directory: {data_size / 1024**3:.1f} GB")
print(f"Output directory: {output_size / 1024**3:.1f} GB")
```

### Temporary File Creation
```python
from metainformant.core import paths

# Create unique temporary files
temp_csv = paths.create_temp_file(".csv", "analysis")
temp_log = paths.create_temp_file(".log", "process", "/tmp/my_app")

# Safe to use without conflicts
with open(temp_csv, 'w') as f:
    f.write("temporary,data\n")
```

## Repository Conventions

The METAINFORMANT project follows strict path conventions:

### Directory Structure
- **`config/`**: Configuration files and run options
- **`data/`**: Input datasets and local databases (read-mostly)
- **`output/`**: All outputs from tests and real runs (ephemeral)
- **`src/`**: Source code
- **`tests/`**: Test files and data
- **`docs/`**: Documentation

### Path Safety
All path operations include security checks:
- No writing outside allowed directories without explicit user paths
- Containment validation for security
- Sanitization of user-provided filenames
- Prevention of path traversal attacks

### Integration with Core Modules
```python
from metainformant.core import paths, io, config

# Load config with path expansion
cfg_path = paths.expand_and_resolve("config/analysis.yaml")
cfg = config.load_mapping_from_file(cfg_path)

# Prepare output directory
output_dir = paths.ensure_directory(cfg.get("output_dir", "output"))

# Write results safely
result_path = output_dir / "analysis_results.json"
io.dump_json(results, result_path)
```

## Error Handling

Path operations provide robust error handling:
- Graceful handling of non-existent files/directories
- Clear error messages for invalid operations
- Safe fallbacks for edge cases

## Dependencies

- **Required**: Standard library `pathlib`, `os`, `tempfile`
- **Optional**: `re` module for filename sanitization
