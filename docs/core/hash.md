# Core: Hashing Utilities

The `hash` module provides cryptographic and non-cryptographic hashing functions for file integrity verification, content fingerprinting, deterministic seeding, and deduplication.

## Purpose

Hashing serves multiple roles in bioinformatics pipelines:
- **Integrity verification**: Ensure downloaded files match expected checksums (MD5, SHA256)
- **Deduplication**: Identify identical sequences or files by hash
- **Deterministic randomness**: Generate reproducible random seeds from strings
- **Content addressing**: Unique identifiers for blob storage
- **Change detection**: Detect file modifications or corruption

The `hash` module focuses on SHA-256 (cryptographic strength) and simple content hashing, with utilities tailored for large files.

## Design Principles

### 1. **Streaming for Large Files**
Never load entire file into memory. Use chunked reading (default 1 MB chunks) for hash computation.

### 2. **SHA-256 as Default**
SHA-256 provides collision resistance suitable for integrity checks. For non-crypto dedup, MD5 could be faster but not provided (you can add).

### 3. **Deterministic Seeding**
`deterministic_seed()` converts arbitrary string (e.g., sample ID) into an integer suitable for `random.seed()` or `numpy.random.seed()`.

### 4. **Simple and Focused**
No complex hash table implementations; raw hash functions only. Use `core.cache` for memoization of hash results if needed.

## API Reference

### Core Hash Functions

#### `sha256_bytes(data: bytes) -> str`

Compute SHA256 hex digest of byte data.

**Parameters**:
- `data`: Input bytes

**Returns**: 64-character hexadecimal string

**Example**:
```python
from metainformant.core.utils import hash

digest = hash.sha256_bytes(b"hello world")
# "b94d27b9934d3e08a52e52d7da7dabfac484efe37a5380ee9088f7ace2efcde9"
```

#### `sha256_string(s: str) -> str`

Compute SHA256 of string content (UTF-8 encoded).

**Parameters**:
- `s`: Input string

**Returns**: Hex digest

**Example**:
```python
digest = hash.sha256_string("AGCTAGCTAGCT")
# Use as sequence fingerprint
```

#### `sha256_file(path: str | Path, *, chunk_size: int = 1024*1024) -> str`

Compute SHA256 of file without full memory load.

**Parameters**:
- `path`: File path
- `chunk_size`: Read chunk size in bytes (default: 1 MiB)

**Returns**: Hex digest

**Raises**: `OSError` if file unreadable

**Example**:
```python
# Verify downloaded file integrity
expected = "a1b2c3d4..."
actual = hash.sha256_file("download/genome.fa.gz")
if actual != expected:
    raise ValueError("Checksum mismatch!")
```

**Performance**:
- 1 GB file: ~2–3 seconds on SSD (limited by disk I/O)
- Chunk size tuning: 1 MB is a good default; larger (10 MB) may improve throughput on fast storage; smaller for memory-constrained systems.

### Comparison Utilities

#### `file_hash_comparison(file1: str | Path, file2: str | Path) -> bool`

Compare two files by SHA256 hash. Returns `True` if identical.

**Example**:
```python
if hash.file_hash_comparison("backup/genome.fa", "current/genome.fa"):
    print("Files identical — no backup needed")
else:
    print("Files differ — update backup")
```

#### `verify_file_integrity(file_path: str | Path, expected_hash: str) -> bool`

Check if file's hash matches expected hex digest.

**Parameters**:
- `expected_hash`: Expected SHA256 hex string (case-insensitive)

**Returns**: `True` if match, `False` otherwise (including file errors)

**Example**:
```python
if not hash.verify_file_integrity("data.tar.gz", expected_hash):
    raise RuntimeError("Download corrupted or tampered")
```

### Directory Hashing

#### `hash_directory(path: str | Path, pattern: str = "**/*") -> dict[str, str]`

Compute hashes for all files in directory, returned as mapping of relative path → hash.

**Parameters**:
- `path`: Root directory
- `pattern`: Glob pattern (default `"**/*"` = all files recursively)

**Returns**: `{relative_file_path: sha256_hex}`

**Example**:
```python
# Snapshot entire output directory
snapshot = hash.hash_directory("output/analysis_20260426")
for rel_path, digest in snapshot.items():
    print(f"{rel_path}: {digest}")

# Compare two directory snapshots
old_snapshot = hash.hash_directory("output/run1/")
new_snapshot = hash.hash_directory("output/run2/")
if old_snapshot != new_snapshot:
    print("Outputs differ between runs")
```

**Note**: Ignores directories, only regular files. Unreadable files skipped.

### Deterministic Seeding

#### `deterministic_seed(data: str) -> int`

Generate integer seed from string (range: 0 to 2³¹−2 inclusive, suitable for `random.seed()` or `numpy.random.seed()`).

**Algorithm**:
1. Compute SHA256 of UTF-8 bytes of `data`
2. Take first 8 hex characters (32 bits)
3. Modulo `2**31 - 1` (maxint for `random.seed()`)

**Parameters**:
- `data`: Arbitrary string (e.g., sample ID, file path)

**Returns**: Integer seed in `[0, 2**31 - 1]`

**Example**:
```python
seed = hash.deterministic_seed("sample_ERR123456")
import random
random.seed(seed)  # Reproducible across runs, Python versions

import numpy as np
np.random.seed(seed)  # Same seed works for numpy
```

**Common pattern**:
```python
def process_sample(sample_id: str, shuffle: bool = True):
    # Ensure reproducible per-sample randomness
    seed = hash.deterministic_seed(sample_id)
    rng = random.Random(seed)

    if shuffle:
        data = list(data)
        rng.shuffle(data)
    # ...
```

## Hash Formats and Standards

### SHA256
- **Digest size**: 256 bits (32 bytes)
- **Hex**: 64 characters (e.g., `"b94d27b9934d3e08a..."`)
- **Collision resistance**: Extremely high (birthday bound ~2¹²⁸)
- **Speed**: ~300–500 MB/s on modern CPU

### When to Use Which Hash

| Use Case | Recommended | Why |
|----------|-------------|-----|
| File integrity (downloads) | SHA256 | Industry standard (PGP, Git) |
| Large directory snapshot | SHA256 per file | Easy comparison |
| In-memory deduplication | SHA256 | Collision resistance important |
| Deterministic random seed | SHA256 first 32 bits | Uniform distribution |
| Fast non-crypto fingerprint | Not provided (MD5) | MD5 broken for security but fast; consider adding if needed |

## Integration Examples

### Example 1: Download with Integrity Check

```python
from metainformant.core import io, hash
import requests

def download_with_checksum(url: str, dest: Path, expected_sha256: str) -> None:
    """Download file and verify SHA256."""
    # Download
    io.download_file(url, dest)

    # Verify
    actual = hash.sha256_file(dest)
    if actual.lower() != expected_sha256.lower():
        dest.unlink(missing_ok=True)
        raise ValueError(f"Checksum mismatch: expected {expected_sha256}, got {actual}")

    print(f"✓ Verified {dest.name}")

# Use with known checksums
checksums = {
    "genome.fa.gz": "e3b0c44298fc1c149afbf4c8996fb924..."
}
download_with_checksum(url, Path("data/genome.fa.gz"), checksums["genome.fa.gz"])
```

### Example 2: Content Deduplication

```python
from metainformant.core.utils import hash as hash_module
from pathlib import Path

def find_duplicates(directory: Path) -> dict[str, list[Path]]:
    """Find duplicate files in directory by hash."""
    hash_to_files: dict[str, list[Path]] = {}

    for file_path in directory.rglob("*"):
        if file_path.is_file():
            file_hash = hash_module.sha256_file(file_path)
            hash_to_files.setdefault(file_hash, []).append(file_path)

    # Return only hashes with >1 file
    return {h: files for h, files in hash_to_files.items() if len(files) > 1}

dupes = find_duplicates(Path("output/"))
for file_hash, files in dupes.items():
    print(f"Duplicate set ({file_hash[:16]}...):")
    for f in files:
        print(f"  {f} ({f.stat().st_size} bytes)")
```

### Example 3: Pipeline Step Fingerprinting

```python
def fingerprint_config(config: dict) -> str:
    """Create fingerprint of configuration for cache key."""
    import json
    config_str = json.dumps(config, sort_keys=True)
    return hash.sha256_string(config_str)[:16]  # Short ID

config_id = fingerprint_config({"threads": 8, "reference": "hg38"})
cache_key = f"step_quant_{config_id}"
```

### Example 4: Change Detection for Incremental Pipeline

```python
import json
from pathlib import Path

def needs_reprocessing(input_files: list[Path], output_file: Path) -> bool:
    """Check if output is stale compared to inputs."""
    if not output_file.exists():
        return True

    output_mtime = output_file.stat().st_mtime
    for inp in input_files:
        if inp.stat().st_mtime > output_mtime:
            return True  # Input newer than output
    return False

def fingerprint_based_check(input_files: list[Path], output_file: Path, fingerprint_file: Path) -> bool:
    """More robust: compare content hashes of inputs."""
    if not fingerprint_file.exists():
        return True

    # Load previous fingerprints
    old_fingerprints = json.loads(fingerprint_file.read_text())
    new_fingerprints = {
        str(f): hash.sha256_file(f) for f in input_files
    }

    if old_fingerprints != new_fingerprints:
        # Inputs changed — reprocess
        fingerprint_file.write_text(json.dumps(new_fingerprints, indent=2))
        return True
    return False
```

## Common Pitfalls

### Pitfall 1: Using Hash as Canonical Identifier Without Context

**Symptom**: Two files have same hash but are actually different (collision).

**Cause**: SHA256 collisions are astronomically unlikely (1 in 2¹²⁸), but theoretically possible. For non-adversarial scenarios, safe.

**Mitigation**: For mission-critical systems (e.g., storing patient genomic data), use double-hash (SHA256 + SHA3-256) or include metadata (path, size) in comparison.

### Pitfall 2: Hash Values as Filenames

**Symptom**: Path length issues, filename collisions from truncation.

**Cause**: 64-char hex digest is long for filename; truncating loses uniqueness.

**Mitigation**:
```python
# Safe: Use full hash or check collisions
def hash_to_filename(file_hash: str, suffix: str = "") -> str:
    # Use first 16 chars + last 4 chars (20 total) with collision check
    short = file_hash[:16] + file_hash[-4:]
    return f"{short}{suffix}"  # e.g., "b94d27b9934d3e08aabcd1234.json"
```

### Pitfall 3: Modifying Files After Hashing

**Symptom**: Cached hash value stale.

**Cause**: File changed after `sha256_file()` computed hash; reuse of old hash.

**Fix**: Recompute hash each time or combine with mtime check:
```python
def cached_hash(file_path: Path, cache: dict) -> str:
    stat = file_path.stat()
    cache_key = (str(file_path), stat.st_mtime, stat.st_size)
    if cache_key in cache:
        return cache[cache_key]
    h = hash.sha256_file(file_path)
    cache[cache_key] = h
    return h
```

## Performance Considerations

### Chunk Size Tuning

Larger chunks reduce `read()` system calls but increase memory per call. Typical:

| Chunk Size | Memory per thread | Good for |
|------------|-------------------|-----------|
| 64 KB | Minimal | Many concurrent hashes (100+ files) |
| 1 MB | Moderate | Default, balanced |
| 10 MB | ~10 MB | Very large files on fast storage (NVMe) |
| 100 MB | ~100 MB | Sequential single-threaded hash of huge files |

**Recommendation**: Keep default (1 MB). Tune only if profiling shows I/O bottleneck.

### Parallel Hashing

Hash many files in parallel using `core.parallel`:

```python
from metainformant.core.execution import parallel
from metainformant.core.utils import hash

files = list(Path("data/").rglob("*.fastq"))
hashes = parallel.thread_map(hash.sha256_file, files, max_workers=8)
# Map file → hash
file_hashes = dict(zip(files, hashes))
```

### SSD vs HDD

- **SSD**: Hashing limited by CPU (SHA256 compute). Larger chunks help.
- **HDD**: Hashing limited by disk seek time; sequential read speed matters. Larger chunks reduce seek overhead.

## Security Notes

SHA-256 is **not** a password hash. Do not use for storing passwords. Use `bcrypt`, `scrypt`, or `argon2` for password hashing.

For download integrity, SHA256 is industry standard (Linux distributions, Git, etc.).

## Related Components

| Module | Relationship |
|--------|--------------|
| `core.io.download` | Download with optional post-hoc SHA256 verification |
| `core.cache` | Cache could be keyed by hash of inputs |
| `core.utils.text` | `clean_sequence_id()` may be combined with hashing for ID normalization |

## Testing

Test against known test vectors:
```python
def test_sha256():
    # NIST test vector
    assert sha256_bytes(b"abc") == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
```

## Future Enhancements

- **MD5 and SHA1** for legacy checksum compatibility (non-crypto)
- **xxHash** for ultra-fast non-crypto hashing (if `xxhash` package available)
- **Parallel directory hashing** with progress bar
- **Hash database**: SQLite backend for persistent hash cache (avoid recompute)
