# Core: Caching System

## Purpose

The `cache` module provides persistent, thread-safe JSON caching with Time-To-Live (TTL) support for expensive computations, API responses, and intermediate results. It's designed for bioinformatics pipelines where network requests and computations can take minutes or hours, making caching essential for productivity.

## Design Principles

### 1. **Persistence Over Ephemeral Memory**
Cache entries survive process restarts and system reboots. Unlike `functools.lru_cache` which is in-memory and per-process, `metainformant.core.io.cache` writes to disk.

### 2. **Automatic Expiration**
TTL-based eviction prevents stale data from being used indefinitely. Each cache entry records its creation time and expiry time.

### 3. **Corruption Resilience**
Corrupted JSON files are automatically detected and removed on read, falling back to cache miss.

### 4. **No External Dependencies**
Pure Python standard library + core modules (`core.io`, `core.utils.logging`). No Redis, Memcached, or database required.

### 5. **Hierarchical Key Namespace**
Keys like `"api/v1/users/sample"` map to sanitized filenames, allowing logical organization without manual directory management.

## Architecture

### Two Caching Strategies

The module provides two APIs:

| API Type | Primary Class/Function | Use Case |
|----------|----------------------|----------|
| **Key-based OO** | `JsonCache` class | Multiple related cache operations, fine-grained control |
| **Key-based functional** | `cache_json()`, `load_cached_json()` | Simple one-off caching |
| **File-based** | `get_json_cache()`, `set_json_cache()` | When cachefile name matters (e.g., pre-defined location) |
| **Bulk operations** | `clear_cache_dir()`, `get_cache_info()` | Maintenance and monitoring |

### Cache Entry Format

Each cache file (`<cache_dir>/<sanitized_key>.json`) stores:

```json
{
  "value": { ... },           // Arbitrary JSON-serializable data
  "created_at": 1703023456.78, // Unix timestamp (float)
  "ttl_seconds": 3600,         // TTL at time of storage (None → no expiry)
  "expires_at": 1703027056.78  // Computed: created_at + ttl_seconds
}
```

The `CacheEntry` class encapsulates this format and provides `is_expired()` checks.

### Thread Safety

`JsonCache` uses a `threading.RLock` (reentrant lock) so that:
- Multiple threads can safely read/write the same cache
- A thread already holding the lock can re-acquire it (useful for nested calls)
- No race conditions during cleanup operations

**Not process-safe**: For multi-process scenarios, use file locking or a shared cache backend (Redis, etc.).

## API Reference

### Classes

#### `CacheEntry`

Represents a single cache entry with TTL metadata.

```text
CacheEntry(value: Any, ttl_seconds: int | None = None)
```

**Methods**:

- `is_expired() -> bool`: Returns True if current time > expires_at or TTL has elapsed.
- `time_to_expiry() -> float | None`: Seconds until expiry; None if no TTL.
- `to_dict() -> Dict[str, Any]`: Serializes entry for JSON storage.
- `from_dict(data: Dict[str, Any]) -> CacheEntry`: Deserializes from JSON.

#### `JsonCache`

Thread-safe JSON cache manager.

```text
JsonCache(cache_dir: str | Path, ttl_seconds: int = 3600)
```

**Parameters**:
- `cache_dir`: Directory to store cache files (created if missing)
- `ttl_seconds`: Default TTL for entries stored without explicit TTL

**Methods**:

| Method | Signature | Description |
|--------|-----------|-------------|
| `get` | `get(key: str, default: Any = None) -> Any` | Retrieve value; returns default on miss/expiry/corruption |
| `set` | `set(key: str, value: Any, ttl_seconds: int | None = None) -> None` | Store value with optional TTL override |
| `delete` | `delete(key: str) -> bool` | Remove single entry; returns True if deleted |
| `clear` | `clear() -> int` | Remove all `.json` files; returns count removed |
| `cleanup_expired` | `cleanup_expired() -> int` | Scan and remove expired/corrupted entries |
| `size` | `size() -> int` | Number of `.json` files in cache directory |
| `stats` | `stats() -> Dict[str, Any]` | Detailed statistics (valid/expired/corrupted counts) |

### Functional API

#### `cache_json(cache_dir, key, data, ttl_seconds=3600) -> None`

Store a JSON-serializable object in cache.

**Parameters**:
- `cache_dir`: Path to cache directory
- `key`: Cache key (hierarchical, e.g., `"api/v1/users"`). Sanitized for filesystem.
- `data`: Any JSON-serializable Python object (dict, list, str, int, float, bool, None)
- `ttl_seconds`: Time-to-live in seconds (default 3600 = 1 hour). Set to `0` or `None` for no expiry.

**Raises**: None (errors logged, not raised)

**Returns**: None

**Notes**: Creates `cache_dir` if missing. Uses atomic write via `io.dump_json()`.

#### `load_cached_json(cache_dir, key, default=None, ttl_seconds=None) -> Any | None`

Load cached data if present and not expired.

**Parameters**:
- `cache_dir`: Path to cache directory
- `key`: Cache key (same as used for `cache_json()`)
- `default`: Value to return on cache miss/expiry/corruption (default: None)
- `ttl_seconds`: Optional override TTL check. If provided, entry is considered expired if age > `ttl_seconds`, regardless of stored TTL. Useful for forcing refresh.

**Returns**: Cached value or `default`

**Notes**:
- If `ttl_seconds=0`, forces expiry (always returns `default`)
- Checks both stored expiry and optional TTL override
- Removes corrupted files automatically

#### `get_json_cache(cache_file, default=None, max_age_seconds=None) -> Any`

File-based variant: operates on a specific cache file path rather than a key-to-file mapping.

**Parameters**:
- `cache_file`: Direct path to `.json` cache file
- `default`: Return value on miss/expiry/error
- `max_age_seconds`: If set, entry older than this is considered expired

**Returns**: Cached value or `default`

#### `set_json_cache(cache_file, data) -> None`

Write data to a specific cache file (no TTL metadata stored; you manage expiry externally via filename or separate metadata).

**Parameters**:
- `cache_file`: Path to cache file
- `data`: JSON-serializable data

**Note**: Simpler than `cache_json()`; doesn't embed TTL. Suitable for simple file-based caches.

#### `clear_cache_dir(cache_dir) -> int`

Delete all `.json` and `.json.gz` files in `cache_dir`. Doesn't recurse into subdirectories.

**Returns**: Number of files deleted

#### `get_cache_info(cache_dir) -> Dict[str, Any]`

Inspect cache directory without modifying it.

**Returns** dictionary:
```python
{
    "exists": bool,          # True if directory exists
    "total_files": int,       # Count of .json files
    "valid_entries": int,     # Non-expired, parseable
    "expired_entries": int,   # TTL elapsed
    "corrupted_files": int,   # JSON parse errors
    "cache_dir": str,         # Absolute path
    "default_ttl": int,       # Default TTL (from JsonCache default)
}
```

## Configuration Options

### Environment Variables

The cache module respects general core environment variables:

| Variable | Effect |
|----------|--------|
| `CORE_LOG_LEVEL` | Controls verbosity of cache-related logs |
| `AK_WORK_DIR` | If using `cache.get_cache_dir()` (affects default locations) |

No cache-specific environment variables exist.

### Default TTL Values

Common TTL patterns (recommended):

| Use Case | TTL (seconds) | Rationale |
|----------|---------------|-----------|
| API responses (public data) | 300–1800 (5–30 min) | Data changes infrequently |
| Expensive computation results | 3600–86400 (1–24 h) | Avoid re-compute within same day |
| Reference data (gene lists) | 604800 (7 days) or `None` | Static, rarely changes |
| Development/debugging | 0 (no caching) | Always fresh |

### Cache Directory Location

By default, you provide the path. Common patterns:

```python
# Project-wide cache
cache_dir = Path("cache")

# Module-specific cache
cache_dir = Path("cache") / "rna" / "quantification"

# Temporary cache (cleared on reboot)
cache_dir = Path("/tmp/metainformant_cache")

# User-level cache
cache_dir = Path.home() / ".cache" / "metainformant"
```

## Code Examples

### 1. Caching API Calls

```python
from metainformant.core import io
from metainformant.core.io import cache
import requests

def fetch_genome_info(species: str) -> dict:
    """Fetch genome metadata from Ensembl with caching."""
    cache_dir = Path("cache/ensembl")
    key = f"genome_info/{species}"

    # Check cache (1 week TTL)
    cached = cache.load_cached_json(cache_dir, key, ttl_seconds=604800)
    if cached:
        return cached

    # Fetch from API
    url = f"https://rest.ensembl.org/info/assembly/{species}?content-type=application/json"
    response = requests.get(url, timeout=30)
    response.raise_for_status()
    data = response.json()

    # Cache result
    cache.cache_json(cache_dir, key, data, ttl_seconds=604800)
    return data

# Usage
human_genome = fetch_genome_info("homo_sapiens")
```

### 2. Caching Expensive Computations

```python
from metainformant.core.io import cache
import numpy as np

def compute_pairwise_distances(sequences: list[str]) -> np.ndarray:
    """Compute pairwise distances with caching."""
    cache_dir = Path("cache/computations")
    key = f"pairwise_{hash(tuple(sequences))}"

    cached = cache.load_cached_json(cache_dir, key)
    if cached is not None:
        return np.array(cached)

    # Expensive computation
    n = len(sequences)
    distances = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            distances[i, j] = distances[j, i] = levenshtein(sequences[i], sequences[j])

    # Cache as list (JSON can't store numpy arrays directly)
    cache.cache_json(cache_dir, key, distances.tolist())
    return distances
```

### 3. Cache Invalidation Strategies

```python
from metainformant.core.io import cache
from datetime import datetime

def conditional_refresh(cache_dir, key, data_source, max_age_hours=24):
    """Refresh cache only if data is older than threshold."""
    cached = cache.load_cached_json(cache_dir, key)
    if cached:
        age_hours = (time.time() - cached.get("fetched_at", 0)) / 3600
        if age_hours < max_age_hours:
            return cached["data"]

    # Fetch fresh
    fresh_data = data_source()
    cache.cache_json(
        cache_dir, key,
        {"data": fresh_data, "fetched_at": time.time()},
        ttl_seconds=None  # No TTL; we manage age manually
    )
    return fresh_data
```

### 4. Cache Warming

```python
def warm_cache_from_manifest(cache_dir, manifest_path):
    """Pre-populate cache from a list of keys (e.g., on pipeline startup)."""
    from metainformant.core import io

    manifest = io.load_json(manifest_path)
    for key, url in manifest["endpoints"].items():
        # Pre-fetch API data
        try:
            data = requests.get(url, timeout=60).json()
            cache.cache_json(cache_dir, key, data, ttl_seconds=86400)
        except Exception as e:
            print(f"Failed to warm {key}: {e}")
```

### 5. Monitoring Cache Health

```python
def monitor_cache(cache_dir, max_size_mb=1000, max_age_days=7):
    """Monitor cache size and evict old entries."""
    import shutil

    info = cache.get_cache_info(cache_dir)

    # Size check
    size_mb = info["total_size"] / (1024 * 1024)
    if size_mb > max_size_mb:
        print(f"Cache size {size_mb:.1f} MB exceeds limit {max_size_mb} MB")
        # Evict: clear entirely or selectively by age
        cache.clear_cache_dir(cache_dir)

    # Age check (manual scan)
    cutoff = time.time() - (max_age_days * 86400)
    for cache_file in Path(cache_dir).glob("*.json"):
        try:
            mtime = cache_file.stat().st_mtime
            if mtime < cutoff:
                cache_file.unlink()
                print(f"Evicted old cache: {cache_file.name}")
        except OSError:
            pass
```

## Common Pitfalls and Troubleshooting

### Pitfall 1: Storing Non-JSON-Serializable Data

**Symptom**: `TypeError: Object of type ndarray is not JSON serializable`

**Cause**: Cache only stores JSON-serializable types (dict, list, str, int, float, bool, None). NumPy arrays, pandas DataFrames, and custom objects are not supported.

**Solutions**:
```python
# Convert to native types before caching
import numpy as np
data = {"matrix": my_array.tolist()}  # Convert ndarray → list
cache.cache_json(cache_dir, key, data)

# Or use pickle-based cache (external) for arbitrary objects
```

### Pitfall 2: Cache Key Collisions Due to Poor Sanitization

**Symptom**: `"api/v1"` and `"api/v1/"` both map to `api_v1.json` → collision

**Cause**: Cache sanitizes keys to filesystem-safe names:
- `/` → `_` (but first `/` preserved in some cases?) Actually implementation: `safe_key = "".join(c for c in key if c.isalnum() or c in "._-")` — so `/` removed entirely!
- `"api/v1"` → `"apiv1"`

**Recommendation**: Use only alphanumeric, `.`, `_`, `-` in keys. Avoid `/` if you rely on key distinction:
```python
# Good: distinct keys
key1 = "api_v1_users"
key2 = "api_v1_posts"

# Risky: both become "apiv1" after sanitization
key1 = "api/v1/users"
key2 = "api/v1/posts"
```

### Pitfall 3: Concurrent Writes from Multiple Processes

**Symptom**: Cache file corruption or JSON parse errors.

**Cause**: `JsonCache` uses thread locks (`threading.RLock`), which do **not** synchronize across processes. Two processes writing the same key simultaneously can corrupt the file.

**Solutions**:
- Use separate cache directories per process (include PID in path)
- Use file locking (`fcntl` on Unix, `msvcrt` on Windows) — not implemented
- Switch to Redis/Memcached for true multi-process safety

### Pitfall 4: Cache Never Expiring

**Symptom**: Stale data used indefinitely.

**Cause**: Setting `ttl_seconds=None` (no expiry) or extremely large TTL (years).

**Fix**: Set reasonable TTL based on data volatility:
- API responses (public): 1–4 hours
- Genome annotations: 1–7 days
- Computed protein structures: 30 days+

Use `cache.cleanup_expired()` in a periodic maintenance job.

### Pitfall 5: Large Cache Bloating Disk

**Symptom**: `cache/` directory grows to GBs.

**Causes**:
- Unbounded cache growth (no eviction policy beyond TTL)
- Storing large objects (whole genome sequences, dense matrices)

**Solutions**:
```python
# Monitor and prune
info = cache.get_cache_info(cache_dir)
if info["total_size"] > 10 * 1024**3:  # 10 GB
    cache.clear_cache_dir(cache_dir)

# Selective pruning: keep only recent entries
for f in Path(cache_dir).glob("*.json"):
    if time.time() - f.stat().st_mtime > 86400:  # >1 day old
        f.unlink()
```

### Pitfall 6: Ignoring Corrupted Cache Files

**Symptom**: Application works but warnings about "corrupted cache file" appear.

**Cause**: Partial writes (e.g., power loss during `cache_json()`) leave invalid JSON.

**Fix**: The module auto-removes corrupted files. But if corruption rate is high (>1%), investigate disk issues or concurrent writes.

## Performance Considerations

### TTL Duration Impact

| TTL | Cache Hit Rate | Data Freshness | Disk Usage |
|-----|---------------|----------------|------------|
| 0 (disabled) | 0% | Always fresh | 0 |
| 60 seconds | Low | Very fresh | Low |
| 3600 seconds (1h) | Medium-High | Reasonable | Medium |
| 86400 seconds (1d) | High | Acceptable for static data | Medium-High |
| ∞ (no expiry) | Very High | Stale | High |

Choose TTL based on acceptable staleness thresholds.

### Disk I/O vs Computation Trade-off

Caching is beneficial when:

```
Cost(computation) + Cost(load_from_source) >> Cost(load_from_cache)
```

For a 10 MB JSON file:
- Network download: 0.5–2 seconds
- Disk read: 0.01–0.05 seconds
- Speedup: 10×–100×

For 1 KB files, overhead may not be worth it.

### Compression

Cache files can be gzipped manually:
```python
import gzip
with gzip.open(cache_file.with_suffix(".json.gz"), "wt") as fh:
    json.dump(data, fh)
```

The cache module doesn't auto-compress, but `core.io.dump_json()` supports `.gz` suffix.

### Cache Directory Location

- **SSD**: Fast random reads; cache works well
- **Network mount (NFS)**: Latency may negate benefits; consider local cache
- **/tmp**: Fast, but cleared on reboot; use for ephemeral data

### Warm vs Cold Cache

First run (cold cache) will populate cache; subsequent runs (warm cache) should be faster. Benchmark by timing `load_cached_json()` calls:
```python
start = time.time()
cached = load_cached_json(cache_dir, "key")
elapsed = time.time() - start
print(f"Cache load: {elapsed*1000:.1f}ms")  # Typically <10ms for small objects
```

## Related Components and Dependencies

### Dependencies

- **Imports**: `metainformant.core.io` (for `load_json`, `dump_json` with atomic writes)
- **Imports**: `metainformant.core.utils.logging` (for `get_logger`)
- **Optional**: None

### Related Modules

| Module | Relationship |
|--------|--------------|
| `core.io` | Cache uses `io.dump_json()` and `io.load_json()` for atomic file operations |
| `core.utils.logging` | Cache operations log warnings for corruption/errors |
| `core.io.paths` | Cache directories often created via `paths.ensure_directory()` |
| `core.download` | Download heartbeat files use JSON cache format for state persistence |

### Typical Usage Patterns

```python
# Common combination: download → cache → process
# 1. Download
raw_data = io.download_json(url)

# 2. Cache
cache.cache_json(cache_dir, "data/latest", raw_data)

# 3. Process
result = process_data(raw_data)

# 4. Cache result
cache.cache_json(cache_dir, "results/latest", result)
```

## Security Considerations

### Key Sanitization

Cache keys are sanitized to prevent directory traversal:

```python
# Dangerous keys are sanitized
"../../etc/passwd"   # -> "etcpasswd" (safe but unusual)
"key; rm -rf /"      # -> "keyrmrf" (safe)
```

Sanitization: keeps only alphanumeric, `.`, `_`, `-`. Empty result → `"default"`.

### Symlink Attacks

If `cache_dir` is user-writable (e.g., `/tmp`), a malicious user could:
1. Create symlink `cache_dir/evil → /etc/cron.d/malicious`
2. Your `cache_json()` writes to that symlink → overwrite critical file

**Mitigation**: Use dedicated cache directory owned by the application user:
```python
cache_dir = Path("/var/cache/metainformant")  # Root-owned, not user-writable
```

### Sensitive Data in Caches

Cache files store plain JSON. Do **not** cache:
- Passwords, API keys, tokens
- Personally Identifiable Information (PII) unless encrypted
- Decrypted private keys

If needed, encrypt before caching:
```python-snippet
from cryptography.fernet import Fernet
encrypted = fernet.encrypt(json.dumps(data).encode())
cache.cache_json(cache_dir, "secure_key", encrypted.decode())
```

## Testing Guidelines

- Test with real filesystem (real implementation)
- Use temporary directories (`tmp_path` fixture in pytest)
- Verify TTL expiration with `time.sleep()` or explicit restoration `time.time()`
- Test concurrent access via `threading.Thread`

See `tests/core/test_core_cache.py` for examples.

## Future Enhancements

- **LRU eviction**: Beyond TTL, limit total number/size of entries
- **Compression**: Auto-gzip large entries (>1 MB)
- **Distributed cache**: Redis backend option
- **Cache statistics**: Histogram of hit rates, sizes, ages
- **Cache warming**: Load frequently-used keys at startup
- **Streaming cache**: For very large datasets, cache chunks separately
