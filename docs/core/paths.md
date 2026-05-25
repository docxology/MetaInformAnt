# Core: Path Management

The `paths` module provides comprehensive path handling utilities for safe, predictable, and secure filesystem operations across all METAINFORMANT modules. It centralizes path resolution, validation, directory management, and output structure discovery.

## Purpose

File system operations are ubiquitous in bioinformatics:
- **Input file discovery**: Locating FASTQ, BAM, VCF, annotation files across nested directories
- **Output organization**: Creating standardized directory hierarchies (`output/<module>/<analysis>/`)
- **Path safety**: Preventing path traversal attacks when handling user-provided paths
- **Cross-platform**: Handling Windows vs POSIX path differences (though project primarily targets Linux)
- **Cleanup**: Calculating directory sizes for disk space management

The `paths` module provides a vetted, consistent API across all domain modules to avoid ad-hoc path manipulation that leads to bugs and security issues.

## Design Principles

### 1. **Pure Functions**
All path functions are pure (no side effects) and return new `Path` objects or values. No module-level mutable state.

### 2. **Cross-Platform Compatibility**
Uses `pathlib.Path` exclusively (never `os.path` string manipulation). Handles Windows path separators if code ever runs on Windows (though Linux is primary target).

### 3. **Security-First**
`is_safe_path()` blocks path traversal patterns (`..`, absolute paths to sensitive locations). `sanitize_filename()` removes dangerous filesystem characters (`<>,:"/\\|?*`). Sanitization is applied at boundaries where user input enters path operations.

### 4. **Fail-Safe Defaults**
Functions return sensible defaults on error:
- `get_file_size()` returns 0 (not exception) for missing files
- `sanitize_filename()` returns `"untitled"` for empty/invalid input
- `get_directory_size()` returns 0 for non-existent directories

### 5. **Auto-Creation**
`ensure_directory()` and `prepare_file_path()` create missing parent directories automatically, preventing `FileNotFoundError` during writes.

### 6. **Output Pattern Discovery**
`discover_output_patterns()` and `find_output_locations()` provide introspection into the project's output directory structure, enabling automation and orchestration.

### 7. **No External Dependencies**
Uses only Python standard library (`pathlib`, `os`, `tempfile`, `re`). No third-party packages required.

### 8. **Project Root Inference**
`get_project_root()` computes the repository root relative to this file's location (`src/metainformant/core/io/paths.py` → 4 parents up). No environment variable or config needed.

## Module Organization

The `paths` module is in `src/metainformant/core/io/paths.py`. It's part of the `io` subpackage.

**Public API**:

**Path resolution & validation**:
- `expand_and_resolve()` — Expand `~` and env vars, resolve to absolute
- `is_within()` — Containment check (path-in-parent)
- `is_safe_path()` — Security validation

**Directory operations**:
- `ensure_directory()` — Create directory tree
- `prepare_file_path()` — Ensure parent dirs exist for a file

**File operations**:
- `get_file_extension()` — Get extension with dot
- `change_extension()` — Replace extension
- `get_file_size()` — Size in bytes
- `get_directory_size()` — Recursive directory total
- `sanitize_filename()` — Remove dangerous characters

**File discovery**:
- `find_files_by_extension()` — Recursive extension search
- `discover_output_patterns()` — Get standard output structure for a module
- `find_output_locations()` — Find existing output dirs
- `get_module_output_base()` — Get default output base path
- `list_output_structure()` — Map entire output tree

**Temporary files**:
- `create_temp_file()` — Create unique non-existent temp path

**Project structure**:
- `get_project_root()` — Repository root
- `get_data_dir()` — Data directory
- `get_cache_dir()` — Cache directory
- `get_logs_dir()` — Logs directory
- `get_temp_dir()` — Temp directory

## API Reference

### `expand_and_resolve(path: str | Path) -> Path`

Expand user home (`~`) and environment variables (`$VAR`), then resolve to absolute path.

**Parameters**:
- `path`: Path string or `Path` object

**Returns**: Absolute `Path` object (does not require file to exist)

**Example**:
```python
p1 = expand_and_resolve("~/data/sample.fastq")
# If HOME=/home/alice → PosixPath('/home/alice/data/sample.fastq')

p2 = expand_and_resolve("$PROJECT/output")
# If PROJECT=/mnt/raid → PosixPath('/mnt/raid/output')

p3 = expand_and_resolve("./relative/path")
# Converted to absolute based on current working directory
```

**Note**: This does not check existence. Use `Path.resolve(strict=True)` to require existence.

### `is_within(path: str | Path, parent: str | Path) -> bool`

Check whether `path` is contained within `parent` directory (after resolving both).

**Parameters**:
- `path`: Candidate path
- `parent`: Expected parent directory

**Returns**: `True` if `path` is within `parent`, `False` otherwise (including if either doesn't exist)

**Security use**: Prevent path traversal attacks where user input tries to access files outside allowed directory.

**Example**:
```python
safe_dir = "/safe/output"
user_path = "../../../etc/passwd"

if is_within(user_path, safe_dir):
    process(user_path)  # OK
else:
    raise ValueError("Path outside safe directory")
```

**Implementation**: Uses `Path.relative_to()` inside try/except catching `ValueError`.

### `ensure_directory(path: Path) -> None`

Create directory and all missing parents. Safe to call on existing directories.

**Parameters**:
- `path`: Directory to create

**Raises**: `OSError` if creation fails (permissions, invalid path)

**Example**:
```python
from metainformant.core.io import paths

# Create deep nested structure
paths.ensure_directory("/mnt/raid/output/analysis/2026-01-15/results")
# Equivalent to: mkdir -p /mnt/raid/output/analysis/2026-01-15/results

# Safe to call on existing
paths.ensure_directory("/mnt/raid/output")  # No error if already exists
```

### `prepare_file_path(file_path: Path) -> None`

Ensure parent directories exist for a file before writing.

**Parameters**:
- `file_path`: Target file path

**Example**:
```python
file_path = Path("output/analysis/results/sample_001.json")
paths.prepare_file_path(file_path)
# Creates output/analysis/results/ if not already present
with open(file_path, 'w') as fh:
    json.dump(data, fh)  # Now safe—parent dirs exist
```

**Combined with `expand_and_resolve`**:
```python
from metainformant.core.io import paths

user_provided = "~/analysis/../output/data.json"
safe_path = paths.expand_and_resolve(user_provided)
paths.prepare_file_path(safe_path)  # Creates parent safely
```

### `is_safe_path(path: str) -> bool`

Check for path traversal attempts and dangerous patterns. Blocks directory escapes outside intended root.

**Parameters**:
- `path`: Path string to validate

**Returns**: `True` if path appears safe, `False` otherwise

**Security checks**:
1. Disallows `..` components (path traversal)
2. Blocks absolute paths to sensitive system directories (`/etc/`, `/root/`, `/sys/`, `/proc/`)
3. Disallows shell metacharacters (`;`, `|`, `&`, `$`) that could enable command injection when used in shell contexts
4. Blocks backslash-escaped traversal on Windows (`..\`)

**Example**:
```python
safe_paths = [
    "output/data.csv",
    "output/subdir/file.txt",
    "./relative/path",
]

dangerous_paths = [
    "../../../etc/passwd",
    "/etc/shadow",
    "output/$(malicious)",
    "output; rm -rf /",
]

for p in safe_paths:
    assert paths.is_safe_path(p), f"Should be safe: {p}"

for p in dangerous_paths:
    assert not paths.is_safe_path(p), f"Should be blocked: {p}"
```

**Important**: `is_safe_path()` is a heuristic; it's not a replacement for proper sandboxing. Combine with `is_within()` for containment.

### `get_file_extension(filename: str) -> str`

Extract file extension including leading dot.

**Parameters**:
- `filename`: Filename or path

**Returns**: Extension with dot (e.g., `".fastq.gz"` returns `".gz"` — only last suffix); empty string if no extension

**Example**:
```python
ext1 = get_file_extension("reads.fastq")         # ".fastq"
ext2 = get_file_extension("archive.tar.gz")       # ".gz"
ext3 = get_file_extension("Makefile")             # ""
ext4 = get_file_extension("/path/to/file.txt")    # ".txt"
```

**Note**: Returns only the final suffix. For compound extensions (`.tar.gz`), you need extra logic:
```python
def get_all_extensions(path: str) -> list[str]:
    p = Path(path)
    exts = []
    while p.suffix:
        exts.append(p.suffix)
        p = p.with_suffix("")
    return list(reversed(exts))
# get_all_extensions("archive.tar.gz") → [".tar", ".gz"]
```

### `change_extension(path: str, new_extension: str) -> Path`

Replace file extension.

**Parameters**:
- `path`: Original file path
- `new_extension`: New extension (with or without leading dot)

**Returns**: New `Path` object with extension replaced

**Example**:
```python
p1 = change_extension("data/reads.fastq", ".fasta")
# PosixPath('data/reads.fasta')

p2 = change_extension("output.vcf.gz", "txt")  # No dot
# PosixPath('output.vcf.txt')

p3 = change_extension("report", ".pdf")  # No existing extension
# PosixPath('report.pdf')

# Chain
p4 = change_extension("data.csv", ".json")
io.dump_json(data, p4)
```

### `find_files_by_extension(directory: str | Path, extension: str) -> list[Path]`

Recursively find all files with given extension.

**Parameters**:
- `directory`: Root directory to search
- `extension`: Extension (with or without leading dot)

**Returns**: List of `Path` objects (sorted alphabetically; empty list if none)

**Example**:
```python
fasta_files = find_files_by_extension("data/refs", ".fasta")
# [PosixPath('data/refs/hg38.fasta'), PosixPath('data/refs/mm10.fasta'), ...]

all_json = find_files_by_extension(".", "json")
# Finds all JSON files in current directory tree
```

**Performance**: Uses `Path.rglob()` which recursively walks entire tree. For large directories (>100K files), consider `os.scandir()` recursive walk for better performance, or use `find` subprocess.

**Filter by pattern**:
```python
# Find FASTQ files (paired-end convention: _R1_, _R2_)
r1_files = [p for p in find_files_by_extension("data/", "fastq.gz")
            if "_R1" in p.name]
r2_files = [p for p in find_files_by_extension("data/", "fastq.gz")
            if "_R2" in p.name]
```

### `get_file_size(path: str | Path) -> int`

Get file size in bytes.

**Parameters**:
- `path`: File path

**Returns**: Size in bytes; 0 if file doesn't exist or inaccessible

**Example**:
```python
size = get_file_size("large.bam")
if size > 10_000_000_000:
    print(f"File is {size / 1e9:.2f} GB — consider splitting")

# Check before expensive operation
if get_file_size("input.fastq") == 0:
    raise ValueError("Input file is empty")
```

**Note**: Returns 0 for directories or missing files. To check existence, use `Path.exists()` separately.

**Compare two files**:
```python
def files_identical(a: str, b: str) -> bool:
    return get_file_size(a) == get_file_size(b) and \
           hashlib.md5(Path(a).read_bytes()).hexdigest() == \
           hashlib.md5(Path(b).read_bytes()).hexdigest()
```

### `get_directory_size(path: str | Path) -> int`

Calculate total size of all regular files recursively.

**Parameters**:
- `path`: Directory path

**Returns**: Total size in bytes (sum of all file sizes); 0 if directory doesn't exist

**Example**:
```python
size = get_directory_size("output/rna/")
print(f"Output directory: {size / 1024**3:.2f} GB")

# Check disk usage before cleanup
data_size = get_directory_size("data/")
if data_size > 1_000_000_000_000:  # 1 TB
    print("Data directory exceeding 1TB — review retention policy")
```

**Implementation**: Recursively walks tree with `rglob("*")`, sums `stat().st_size` for regular files, catches `OSError` on permission errors and continues.

### `sanitize_filename(filename: str) -> str`

Sanitize filename for safe filesystem use. Removes dangerous characters and normalizes to safe representation.

**Parameters**:
- `filename`: Original filename (may include extension)

**Returns**: Sanitized filename safe for POSIX filesystems

**Sanitization rules**:
1. Replace characters `<>:"/\\|?*` with underscores `_`
2. Remove Unicode control characters (C0/C1 control sets)
3. Strip leading/trailing whitespace and dots (Windows forbids trailing dots)
4. If result empty, return `"untitled"`
5. Truncate to 255 characters max (most filesystems' limit), preserving extension

**Example**:
```python
unsafe = 'report: 2023/analysis?<test>.pdf'
safe = sanitize_filename(unsafe)
# 'report-2023-analysis_test_.pdf'

# Query strings → underscores
sanitize_filename("file_name[version 2].txt")  # 'file_name[version 2].txt' → brackets kept (not dangerous)

# Control chars removed
sanitize_filename("file\x00with\x01control.txt")  # 'filewithcontrol.txt'
```

### `create_temp_file(suffix="", prefix="tmp", directory=None) -> Path`

Create a unique temporary file path that doesn't exist yet.

**Parameters**:
- `suffix`: File extension (default `""`)
- `prefix`: Filename prefix (default `"tmp"`)
- `directory`: Directory for temp file; if `None`, uses system temp dir

**Returns**: `Path` to non-existent file

**Note**: Unlike `tempfile.NamedTemporaryFile(delete=False)`, this only creates the pathname (not the file). Race condition window exists between path generation and file creation—caller should create file immediately or use `NamedTemporaryFile` directly if atomicity required.

**Example**:
```python
# Create temp CSV path in system temp directory
temp_csv = paths.create_temp_file(".csv", prefix="analysis")
# e.g., /tmp/tmp_1705331234_abcdef.csv

# Create in project's temp directory
project_tmp = get_temp_dir()
temp_log = paths.create_temp_file(".log", "job", directory=project_tmp)
# e.g., /path/to/project/tmp/tmp_job_xyz.log

# Use immediately
with open(temp_csv, 'w') as fh:
    writer = csv.writer(fh)
    writer.writerow(["sample", "value"])
```

**Collision handling**: The function loops with different hash until a non-existent path is found (very unlikely to collide).

### `discover_output_patterns(module_name: str) -> dict[str, Any]`

Get standard output directory pattern for a METAINFORMANT module.

**Parameters**:
- `module_name`: Module name (e.g., `'rna'`, `'gwas'`, `'dna'`)

**Returns**: Dictionary:
```python
{
    "base_pattern": "output/amalgkit/<species>/<step>/",
    "subdirs": ["quant/", "work/", "logs/"],
    "examples": ["output/amalgkit/Apis_mellifera/quant/", ...],
}
```

**Known modules**:
- `rna` → `output/amalgkit/<species>/<step>/` with subdirs `quant/`, `work/`, `logs/`
- `gwas` → `output/gwas/<analysis_type>/` with `association/`, `plots/`, `qc/`
- `dna` → `output/dna/<analysis_type>/` with `phylogeny/`, `population/`, `variants/`
- `protein` → `output/protein/<analysis_type>/` with `structures/`, `alignments/`
- Default → `output/<module>/`

**Example**:
```python
pattern = paths.discover_output_patterns("rna")
print(f"Base: {pattern['base_pattern']}")
# Base: output/amalgkit/<species>/<step>/

# Template variables like <species> indicate runtime parameterization
# Resolve for specific case:
species = "Homo_sapiens"
step = "quant"
output_dir = Path(pattern['base_pattern'].replace("<species>", species).replace("<step>", step))
# Path('output/amalgkit/Homo_sapiens/quant/')
```

### `find_output_locations(repo_root: str | Path, pattern: str | None = None) -> list[Path]`

Find existing output directories in `output/` subtree.

**Parameters**:
- `repo_root`: Repository root path
- `pattern`: Optional filter substring (case-insensitive); if `None`, returns all output subdirectories

**Returns**: Sorted list of unique directory `Path` objects

**Example**:
```python
# All output directories
all_outs = paths.find_output_locations(".")
# [Path('output/dna/phylogeny/'), Path('output/gwas/association/'), ...]

# Filtered
rna_outs = paths.find_output_locations(".", pattern="rna")
# [Path('output/amalgkit/Apis_mellifera/quant/'), ...]
```

### `get_module_output_base(module_name: str) -> str`

Get default output base path for a module (no template variables).

**Parameters**:
- `module_name`: Module name

**Returns**: Base path string like `"output/dna"` or `"output/amalgkit"` (amalgkit for RNA)

**Example**:
```python
base = get_module_output_base("gwas")
# "output/gwas"

full = Path(base) / "association"
# Path('output/gwas/association')
```

### `list_output_structure(repo_root: str | Path) -> dict[str, Any]`

Recursively map the entire `output/` directory structure with statistics.

**Parameters**:
- `repo_root`: Repository root

**Returns**: Nested dictionary:
```text
{
    "total_dirs": 150,
    "total_files": 12000,
    "total_size": 450_000_000_000,  # 450 GB
    "structure": {
        "type": "directory",
        "path": "output",
        "children": {
            "dna": { ... nested ... },
            "gwas": { ... },
        },
    },
}
```

**Example**:
```python
import json
stats = paths.list_output_structure(".")
print(f"Total output: {stats['total_files']} files, "
      f"{stats['total_size'] / 1e9:.1f} GB across {stats['total_dirs']} dirs")

# Pretty-print structure (first few levels)
print(json.dumps(stats["structure"], indent=2)[:1000])
```

**Use cases**:
- Cleanup scripts identify large output trees
- Audit trail of all analyses produced
- Workflow status dashboards

### Project Root Helpers

```python
get_project_root()  # → Path to repository root (via __file__ location)
get_data_dir()      # → Path(root / "data")
get_cache_dir()     # → Path(root / "cache")
get_logs_dir()      # → Path(root / "logs")
get_temp_dir()      # → Path(root / "tmp")
```

**Example**:
```python
root = paths.get_project_root()
print(f"Project root: {root}")
# /path/to/MetaInformAnt

cache = paths.get_cache_dir()
cache.mkdir(exist_ok=True)
# -> /path/to/MetaInformAnt/cache
```

## Usage Examples

### Standard Output Directory Setup

```python
from metainformant.core.io import paths
from datetime import datetime

def setup_run_output(module: str, run_id: str | None = None) -> Path:
    """Create timestamped output directory for a module run."""
    if run_id is None:
        run_id = datetime.now().strftime("%Y%m%d_%H%M%S")

    base = Path(paths.get_module_output_base(module))
    output_dir = base / run_id

    paths.ensure_directory(output_dir)
    print(f"Output directory: {output_dir}")
    return output_dir

# Usage:
rna_out = setup_run_output("rna")  # e.g., output/amalgkit/20260115_143022/
gwas_out = setup_run_output("gwas")  # e.g., output/gwas/20260115_143022/
```

### Safe User-Provided Path Handling

```python
def process_user_upload(user_path: str, allowed_root: str = "uploads/") -> Path:
    """Safely handle user-provided file path."""
    # Expand and resolve
    full_path = paths.expand_and_resolve(user_path)

    # Check containment
    root = paths.expand_and_resolve(allowed_root)
    if not paths.is_within(full_path, root):
        raise ValueError(f"Path {full_path} is outside allowed directory {root}")

    # Validate file exists
    if not full_path.is_file():
        raise FileNotFoundError(f"Not a file: {full_path}")

    # Sanitize for any downstream operations (logging, display)
    safe_name = paths.sanitize_filename(full_path.name)

    return full_path
```

### Finding Input Files by Pattern

```python
def find_fastq_pairs(data_dir: str) -> list[tuple[Path, Path]]:
    """Find paired-end FASTQ files."""
    all_fastq = paths.find_files_by_extension(data_dir, ".fastq.gz")

    pairs = []
    for r1 in all_fastq:
        if "_R1" not in r1.name:
            continue
        r2 = r1.parent / r1.name.replace("_R1", "_R2")
        if r2 in all_fastq:
            pairs.append((r1, r2))

    return pairs

pairs = find_fastq_pairs("data/raw/")
print(f"Found {len(pairs)} paired-end FASTQ files")
```

### Output Directory Cleanup by Age

```python
import shutil
from datetime import datetime, timedelta

def cleanup_old_outputs(max_age_days: int = 30) -> int:
    """Remove output directories older than max_age_days."""
    root = Path(get_project_root())
    output_root = root / "output"
    if not output_root.exists():
        return 0

    cutoff = datetime.now() - timedelta(days=max_age_days)
    removed = 0

    for subdir in output_root.iterdir():
        if not subdir.is_dir():
            continue
        mtime = datetime.fromtimestamp(subdir.stat().st_mtime)
        if mtime < cutoff:
            size = get_directory_size(subdir)
            print(f"Removing {subdir} ({size / 1e9:.2f} GB, older than {max_age_days} days)")
            shutil.rmtree(subdir)
            removed += 1

    return removed

removed = cleanup_old_outputs(30)
print(f"Cleaned up {removed} old output directories")
```

### Atomic File Move with Directory Creation

```python
def atomic_move(src: Path, dest_dir: str | Path, new_name: str | None = None) -> Path:
    """Move file to destination directory, creating parents as needed."""
    dest_parent = Path(dest_dir)
    paths.ensure_directory(dest_parent)

    dest_name = new_name or src.name
    dest = dest_parent / dest_name

    # Atomic rename within same filesystem
    src.rename(dest)
    return dest

src_temp = Path("/tmp/tmp_download.bam")
final = atomic_move(src_temp, "output/bam/", "sample001.bam")
```

### Output Structure Traversal for Reporting

```python
def summarize_outputs() -> dict[str, Any]:
    """Generate summary of all module outputs."""
    root = get_project_root()
    structure = list_output_structure(root)

    summary = {
        "total_size_gb": structure["total_size"] / 1e9,
        "total_files": structure["total_files"],
        "modules": {},
    }

    # Count per module
    for module in ["dna", "rna", "gwas", "protein", "singlecell"]:
        pattern = discover_output_patterns(module)
        base_pattern = pattern["base_pattern"]
        # Find all directories matching module pattern
        module_dirs = [
            p for p in find_output_locations(root)
            if base_pattern.split("<")[0] in str(p)
        ]
        summary["modules"][module] = {
            "directory_count": len(module_dirs),
            "size_gb": sum(get_directory_size(d) for d in module_dirs) / 1e9,
        }

    return summary

report = summarize_outputs()
print(json.dumps(report, indent=2))
```

### Cross-Platform Path Handling

```python
# Although METAINFORMANT targets Linux, pathlib handles Windows too
def make_platform_safe(path: str | Path) -> Path:
    """Normalize path separators for current OS."""
    p = Path(path)
    # On Windows, convert forward slashes to backslashes for APIs that require them
    # (Usually not needed—pathlib does it)
    return p

# Display paths consistently (forward slashes)
def display_path(p: Path) -> str:
    """Return path string with forward slashes for logs/docs."""
    return str(p).replace("\\", "/")
```

### Disk Space Check Before Large Operation

```python
def ensure_disk_space(needed_bytes: int, path: str | Path = "/") -> bool:
    """Check if enough free space exists at path."""
    stat = shutil.disk_usage(path)
    free = stat.free
    required = needed_bytes * 1.1  # 10% safety margin

    if free < required:
        raise OSError(
            f"Insufficient disk space at {path}: "
            f"need {required / 1e9:.1f} GB free, have {free / 1e9:.1f} GB"
        )
    return True

# Usage: Check before downloading 50GB file
ensure_disk_space(50_000_000_000, "output/")
result = download_with_progress(url, "output/large.bam")
```

### Temp Directory Management

```python
import tempfile

def process_with_temp_dir() -> Path:
    """Create isolated temporary working directory."""
    tmp_root = get_temp_dir()
    run_tmp = tmp_root / f"run_{int(time.time())}"
    ensure_directory(run_tmp)

    # Use for intermediate files
    intermediate = run_tmp / "step1.tmp"
    intermediate.write_text("data")

    # Clean up when done
    shutil.rmtree(run_tmp)
    return run_tmp
```

## Error Handling

Path operations are designed to be robust, but understanding failure modes helps debugging:

### Non-Existent Paths

**`expand_and_resolve()`**: Returns a Path even if file doesn't exist (uses `strict=False`). This is intentional—you often want the absolute path before checking existence.

Check explicitly:
```python
p = paths.expand_and_resolve("maybe_missing.txt")
if not p.exists():
    print(f"Path does not exist: {p}")
```

**`get_file_size()`**: Returns 0 for missing files instead of raising. Useful for:
```python
size = paths.get_file_size("input.fastq")
if size == 0:
    print("File missing or empty")
```

**`get_directory_size()`**: Returns 0 for non-existent directory.

### Permission Errors

**Symptom**: `PermissionError: [Errno 13] Permission denied` when calling `ensure_directory()` or accessing files.

**Causes**:
- Directory not writable by current user
- File owned by another user with restricted permissions
- Attempting to write to read-only filesystem

**Fix**:
```bash
# Check permissions
ls -ld /path/to/directory

# Change ownership (if you have sudo)
sudo chown -R $USER /path/to/directory

# Adjust permissions
chmod u+rwx /path/to/directory
```

Or in code:
```python
import os
os.chmod("/path", 0o755)  # Requires ownership or sudo
```

### Path Too Long

**Symptom**: `OSError: [Errno 36] File name too long`.

**Cause**: Path exceeds filesystem limit (typically 255 bytes for single component, 4096 for full path).

**Prevention**:
```python
def truncate_path_components(path: Path, max_len: int = 200) -> Path:
    """Truncate each path component to max_len characters."""
    parts = []
    for part in path.parts:
        if len(part) > max_len:
            parts.append(part[:max_len])
        else:
            parts.append(part)
    return Path(*parts)
```

`sanitize_filename()` already truncates to 255 chars.

### Invalid Characters

**Symptom**: `OSError: [Errno 22] Invalid argument` on Windows with `<>:"/\\|?*` in filename.

**Fix**: Use `sanitize_filename()` before creating files:
```python
unsafe = "con:file.txt"  # CON is reserved on Windows
safe = sanitize_filename(unsafe)  # "con_file.txt"
```

### Symlink Loops

**Symptom**: `RecursionError` or infinite loop when using `resolve()` on symlink cycle.

**Fix**: Use `resolve(strict=False)` which doesn't follow symlinks (already used in `expand_and_resolve`). If you need to detect cycles, use `Path.resolve()` in try/except catching `RuntimeError` (Python 3.6+).

### Disk Full During Operation

**Symptom**: `OSError: [Errno 28] No space left on device`.

**Pre-check**:
```python
def check_disk_space(path: Path, needed_bytes: int) -> bool:
    stat = shutil.disk_usage(path)
    return stat.free >= needed_bytes

# Before writing 10GB file:
if not check_disk_space(Path("output/"), 10_000_000_000):
    raise OSError("Insufficient space")
```

**Atomic writes** (used in `io.dump_json()`) write to temp file then rename. If disk fills during write, `OSError` raised and temp file cleaned up on next run.

### Unmounted Filesystem

**Symptom**: `FileNotFoundError` for path that should exist on mounted filesystem.

**Cause**: Network filesystem (NFS, SMB) not mounted or stale mount.

**Check**:
```bash
mount | grep /mnt/raid
# Or:
df -h /mnt/raid
```

## Integration with Other Core Modules

### With I/O Module

```python
from metainformant.core import io
from metainformant.core.io import paths

def safe_read_and_write(input_path: str, output_path: str):
    """Read input, transform, write to output with all safety checks."""
    # Resolve and validate input
    input_file = paths.expand_and_resolve(input_path)
    if not input_file.is_file():
        raise FileNotFoundError(f"Input not found: {input_file}")

    # Prepare output
    output_file = paths.expand_and_resolve(output_path)
    if not paths.is_within(output_file, paths.get_project_root() / "output"):
        raise ValueError("Output must be under output/ directory")
    paths.prepare_file_path(output_file)

    # Read and write using io module
    data = io.load_json(input_file)
    processed = transform(data)
    io.dump_json(processed, output_file, atomic=True)
```

### With Config Module

```python
from metainformant.core.io import paths
from metainformant.core.utils import config

def load_config_with_path_expansion(config_path: str) -> dict:
    """Load config file with ~ and env var expansion."""
    expanded = paths.expand_and_resolve(config_path)
    if not expanded.exists():
        raise FileNotFoundError(f"Config not found: {expanded}")

    return config.load_mapping_from_file(expanded)

cfg = load_config_with_path_expansion("~/config/pipeline.yaml")
```

### With Download Module

```python
from metainformant.core.io import download, paths

def safe_download_to_output(url: str, filename: str) -> Path:
    """Download file to output directory with sanitization."""
    # Sanitize filename to prevent path traversal in url
    safe_name = paths.sanitize_filename(filename)
    dest = Path(get_project_root()) / "output" / "downloads" / safe_name

    # Ensure output dir exists
    paths.prepare_file_path(dest)

    # Download
    result = download.download_with_progress(url, dest, show_progress=True)
    if not result.success:
        raise RuntimeError(result.error)

    return dest
```

### In a Workflow Context

```text
# In a Snakemake rule, use paths for consistent directory structure
rule align:
    output:
        bam="output/alignment/{sample}.bam",
        bai="output/alignment/{sample}.bam.bai"
    params:
        index="reference/genome",
        temp_dir=lambda wildcards: str(
            paths.get_temp_dir() / f"align_{wildcards.sample}"
        )
    shell:
        """
        bwa mem -t 8 {params.index} \
            input/{wildcards.sample}.fastq.gz | \
            samtools view -b - > {output.bam}

        # Ensure temp cleaned
        rm -rf {params.temp_dir}
        """
```

## Security Notes

### Path Traversal Prevention

Always validate user-provided paths before use:
```python
def validate_user_path(user_input: str, allowed_root: Path) -> Path:
    full = paths.expand_and_resolve(user_input)

    # Method 1: is_within() containment
    if not paths.is_within(full, allowed_root):
        raise ValueError(f"Path escapes allowed root: {full}")

    # Method 2: explicit relative_to()
    try:
        full.relative_to(allowed_root)
    except ValueError:
        raise ValueError(f"Path outside allowed root: {full}")

    return full
```

### Symlink Attacks

Beware of symlinks pointing outside allowed directories:
```python
def safe_access(path: Path, allowed_root: Path) -> bool:
    """Check that path and all its parents are within allowed_root."""
    try:
        resolved = path.resolve(strict=True)  # Follow symlinks
    except FileNotFoundError:
        return False

    # Check resolved real path
    try:
        resolved.relative_to(allowed_root.resolve())
        return True
    except ValueError:
        return False
```

### Filename Whitelisting

For strict validation, combine sanitization with allowlist:
```python
ALLOWED_EXTENSIONS = {".fastq.gz", ".fasta", ".fq", ".fa", ".bam", ".bai"}

def validate_download_filename(filename: str) -> str:
    safe = paths.sanitize_filename(filename)
    ext = "".join(Path(safe).suffixes)  # Get compound extension
    if ext not in ALLOWED_EXTENSIONS:
        raise ValueError(f"File type not allowed: {ext}")
    return safe
```

### Absolute Path Policy

Decide whether absolute paths are allowed:
```python
def normalize_to_relative(path: str | Path, base: Path) -> Path:
    """Convert absolute path to relative if under base; reject otherwise."""
    p = Path(path).resolve()
    try:
        rel = p.relative_to(base.resolve())
        return rel
    except ValueError:
        raise ValueError(f"Path {p} not under base {base}")
```

## Testing Guidelines

### Test with Temporary Directories

```python-snippet
from pathlib import Path

def test_expand_and_resolve(tmp_path, explicit restoration):
    explicit restoration.setenv("HOME", str(tmp_path / "home"))
    result = paths.expand_and_resolve("~/data")
    assert result.is_absolute()
    assert str(result).startswith(str(tmp_path / "home"))

def test_ensure_directory(tmp_path):
    new_dir = tmp_path / "a" / "b" / "c"
    paths.ensure_directory(new_dir)
    assert new_dir.is_dir()
    assert new_dir.exists()

def test_sanitize_filename_edge_cases():
    assert paths.sanitize_filename("") == "untitled"
    assert paths.sanitize_filename("   ") == "untitled"
    assert paths.sanitize_filename("a" * 300)  # Truncates to 255

def test_is_safe_path():
    assert paths.is_safe_path("output/data.csv")
    assert not paths.is_safe_path("../../../etc/passwd")
```

### Integration Tests with Real Filesystem

```python
def test_prepare_file_path_creates_nested_dirs(tmp_path):
    file_path = tmp_path / "deep" / "nested" / "file.txt"
    paths.prepare_file_path(file_path)
    assert file_path.parent.exists()
    assert file_path.parent.is_dir()

def test_find_files_by_extension_recursive(tmp_path):
    # Create nested structure
    (tmp_path / "a/b").mkdir(parents=True)
    (tmp_path / "a/file1.txt").touch()
    (tmp_path / "a/b/file2.txt").touch()
    (tmp_path / "a/b/file3.tsv").touch()

    txts = paths.find_files_by_extension(tmp_path, ".txt")
    assert len(txts) == 2
    assert all(f.suffix == ".txt" for f in txts)
```

### Test Output Pattern Discovery

```python
def test_discover_output_patterns_known_modules():
    rna_pattern = paths.discover_output_patterns("rna")
    assert "amalgkit" in rna_pattern["base_pattern"]
    assert "quant" in rna_pattern["subdirs"]

    gwas_pattern = paths.discover_output_patterns("gwas")
    assert "gwas" in gwas_pattern["base_pattern"]
    assert "association" in gwas_pattern["subdirs"]
```

## Dependencies

### Required

Standard library only:
- `pathlib` — Path objects
- `os` — Filesystem checks
- `tempfile` — Temporary file generation
- `re` — Regex for sanitization

### Optional

None.

## Troubleshooting

| Symptom | Likely Cause | Fix |
|---------|--------------|-----|
| `FileNotFoundError` after `ensure_directory()` | Permissions denied on parent | Check write permissions on parent directory; use `sudo` or change ownership if needed |
| `OSError: [Errno 36] File name too long` | Sanitized/auto-generated name exceeds 255 chars | Truncate name; use shorter prefix in `create_temp_file()` |
| `find_files_by_extension()` returns empty | Wrong starting directory or extension format | Verify `directory` exists; extension should be `".fastq"` or `"fastq"` (both accepted) |
| `is_within()` returns False for nested path | One path not resolved/resolvable | Call `expand_and_resolve()` on both inputs first |
| Disk full during operation | Insufficient space | Use `get_directory_size()` to audit; clean up old outputs |

## Repository Conventions

### Directory Layout

METAINFORMANT project root layout:
```
project_root/
├── config/          # Configuration files
├── data/            # Input datasets
├── output/          # All analysis outputs
│   ├── dna/
│   ├── rna/amalgkit/
│   ├── gwas/
│   └── ...
├── cache/           # Cached downloads, computed values
├── logs/            # Log files
├── tmp/             # Temporary files (ephemeral)
├── src/
│   └── metainformant/
└── docs/
```

All paths module functions understand and operate within these conventions.

### Path Safety Policy

- All user-facing file operations must use `paths` functions for validation
- Never trust raw string paths from environment variables, CLI args, or user input
- Always combine `expand_and_resolve()` + `is_within()` for containment checks
- Sanitize filenames before passing to filesystem

### Output Location Resolution

```python
from metainformant.core.io import paths

# Determine output location for module
module = "gwas"
base = paths.get_module_output_base(module)  # "output/gwas"
output_dir = Path(base) / "association" / run_id
paths.ensure_directory(output_dir)
```

## Performance Notes

- `expand_and_resolve()` uses `Path.resolve(strict=False)` — fast; doesn't hit disk unless symlink resolution needed
- `get_directory_size()` walks entire subtree — can be slow for deep trees with many files. Caching recommended.
- `find_files_by_extension()` uses `Path.rglob()` — efficient but still O(N) in directory entries. For repeated queries, cache results.

### Caching Directory Sizes

```python
from functools import lru_cache

@lru_cache(maxsize=32)
def get_cached_dir_size(path: str) -> int:
    return paths.get_directory_size(path)

# Repeated calls with same path hit cache
```

## Further Reading

- `pathlib` documentation: https://docs.python.org/3/library/pathlib.html
- OWASP Path Traversal: https://owasp.org/www-community/attacks/Path_Traversal
