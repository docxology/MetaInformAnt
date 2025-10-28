# Core: Hashing and Integrity Utilities

The `hash` module provides cryptographic hashing utilities for data integrity, deterministic seed generation, and file verification.

## Functions

### SHA-256 Hashing
- **`sha256_bytes(data)`** → `str`
  - Compute SHA-256 hex digest of byte data
  - Returns 64-character hexadecimal string

- **`sha256_file(path, chunk_size=1048576)`** → `str`
  - Compute SHA-256 hash of file without loading entirely into memory
  - Uses configurable chunk size for large files
  - Memory-efficient for files of any size

- **`sha256_string(s)`** → `str`
  - Compute SHA-256 hash of string content
  - UTF-8 encoding used internally

### Integrity Verification
- **`verify_file_integrity(file_path, expected_hash)`** → `bool`
  - Verify file hasn't been corrupted or tampered with
  - Compares computed hash against expected value
  - Returns False for unreadable files

- **`file_hash_comparison(file1, file2)`** → `bool`
  - Compare two files by computing their hashes
  - Efficient alternative to byte-by-byte comparison
  - Useful for duplicate detection

### Deterministic Operations
- **`deterministic_seed(data)`** → `int`
  - Generate deterministic integer seed from string data
  - Seed is in valid range for random number generators
  - Useful for reproducible simulations

### Batch Operations
- **`hash_directory(path, pattern="**/*")`** → `Dict[str, str]`
  - Compute hashes for all files in directory matching pattern
  - Returns relative paths as keys, hashes as values
  - Skips unreadable files gracefully

## Usage Examples

### Basic Hashing
```python
from metainformant.core import hash

# Hash byte data
data = b"Hello, World!"
digest = hash.sha256_bytes(data)
print(f"SHA256: {digest}")

# Hash string content
text_hash = hash.sha256_string("Hello, World!")
print(f"String hash: {text_hash}")

# Hash large files efficiently
file_hash = hash.sha256_file("large_dataset.fasta")
print(f"File hash: {file_hash}")
```

### File Integrity Verification
```python
from metainformant.core import hash

# Verify downloaded file integrity
expected_hash = "a665a45920422f9d417e4867efdc4fb8a04a1f3fff1fa07e998e86f7f7a27ae3"
is_valid = hash.verify_file_integrity("downloaded_file.zip", expected_hash)

if is_valid:
    print("File integrity verified ✓")
else:
    print("File corrupted or tampered with ✗")
```

### Deterministic Seeds
```python
from metainformant.core import hash
import random

# Generate reproducible random sequences
seed_value = hash.deterministic_seed("experiment_123")
random.seed(seed_value)

# Will always produce the same sequence
values = [random.random() for _ in range(5)]
print(f"Reproducible values: {values}")
```

### Directory Integrity Checking
```python
from metainformant.core import hash

# Hash all Python files in project
python_files = hash.hash_directory(".", "**/*.py")

# Check for changes
for file_path, file_hash in python_files.items():
    print(f"{file_path}: {file_hash[:8]}...")
```

### Duplicate File Detection
```python
from metainformant.core import hash

def find_duplicates(file_list):
    """Find duplicate files using hashing."""
    seen_hashes = {}
    duplicates = []

    for file_path in file_list:
        file_hash = hash.sha256_file(file_path)
        if file_hash in seen_hashes:
            duplicates.append((file_path, seen_hashes[file_hash]))
        else:
            seen_hashes[file_hash] = file_path

    return duplicates

# Usage
files = ["data1.fasta", "data2.fasta", "data3.fasta"]
dups = find_duplicates(files)
for dup, original in dups:
    print(f"Duplicate: {dup} matches {original}")
```

## Performance Considerations

### Memory Efficiency
- **`sha256_file`**: Uses chunked reading (default 1MB chunks)
- Large files processed without memory exhaustion
- Configurable chunk size for optimization

### Speed Optimization
- Cryptographic operations are CPU-intensive
- SHA-256 provides good balance of speed and collision resistance
- Consider caching hashes for frequently accessed files

## Security Notes

- SHA-256 is cryptographically secure for integrity checking
- Not suitable for password hashing (use dedicated functions)
- Hash collisions are theoretically possible but practically unlikely
- File integrity checking helps detect corruption and tampering

## Error Handling

Hash operations are designed to be robust:
- File I/O errors return appropriate default values
- Unreadable files are skipped in batch operations
- Deterministic seed generation handles edge cases
- All functions work with pathlib.Path and string paths

## Integration with Core Modules

```python
from metainformant.core import hash, cache, io

def cached_computation_with_integrity(cache_dir, key, func, *args):
    """Cache computation result with integrity verification."""
    cache_file = cache_dir / f"{key}.json"

    # Try to load cached result
    if cache_file.exists():
        cached_data = io.load_json(cache_file)
        # Verify cache file integrity
        if hash.verify_file_integrity(cache_file, cached_data.get('_integrity_hash')):
            return cached_data['result']

    # Compute and cache with integrity hash
    result = func(*args)
    data_to_cache = {
        'result': result,
        '_integrity_hash': hash.sha256_string(str(result))
    }
    io.dump_json(data_to_cache, cache_file)
    return result
```

## Dependencies

- **Required**: Standard library `hashlib`, `pathlib`
- **Optional**: None (all functionality uses stdlib)
