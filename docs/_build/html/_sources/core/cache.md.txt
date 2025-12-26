# Core: Caching System

The `cache` module provides JSON-based caching with Time-To-Live (TTL) support for efficient data storage and retrieval, particularly useful for expensive operations and API responses.

## Functions

### Key-Based Caching
- **`cache_json(cache_dir, key, obj)`** → `Path`
  - Store JSON-serializable object with string key
  - Automatic directory creation and file management
  - Keys support hierarchical structure with "/"
  - Returns path to cached file

- **`load_cached_json(cache_dir, key, ttl_seconds)`** → `Any | None`
  - Load cached JSON object by key
  - TTL-based expiration (0 = no expiration)
  - Returns None if cache miss or expired

### File-Based Caching
- **`get_json_cache(cache_file, default=None, max_age_seconds=None)`** → `Any`
  - Load data from specific cache file
  - Optional TTL checking
  - Returns default value on miss/expiration/error

- **`set_json_cache(cache_file, data)`** → `None`
  - Store data in specific cache file
  - Automatic directory creation

### Cache Management
- **`clear_cache_dir(cache_dir)`** → `None`
  - Remove all JSON cache files from directory
  - Safe operation (only removes .json/.gz files)

- **`get_cache_info(cache_dir)`** → `Dict[str, Any]`
  - Get statistics about cache directory
  - File count, total size, existence status

## Usage Examples

### Basic Key-Based Caching
```python
from pathlib import Path
from metainformant.core import cache

# Set up cache directory
cache_dir = Path("output/cache")

# Cache expensive computation results
def expensive_analysis(data):
    # Simulate expensive operation
    import time
    time.sleep(2)
    return {"result": sum(data), "count": len(data)}

# Cache the result
data = [1, 2, 3, 4, 5]
cache.cache_json(cache_dir, "analysis/sum", expensive_analysis(data))

# Retrieve from cache (fast!)
cached_result = cache.load_cached_json(cache_dir, "analysis/sum", ttl_seconds=3600)
print(f"Cached result: {cached_result}")
```

### Hierarchical Cache Keys
```python
from metainformant.core import cache

# Use hierarchical keys for organization
cache.cache_json(cache_dir, "api/responses/users", user_data)
cache.cache_json(cache_dir, "api/responses/posts", post_data)
cache.cache_json(cache_dir, "computations/pca/results", pca_results)

# Load specific cached items
users = cache.load_cached_json(cache_dir, "api/responses/users", ttl_seconds=1800)
posts = cache.load_cached_json(cache_dir, "api/responses/posts", ttl_seconds=1800)
```

### File-Based Caching with TTL
```python
from metainformant.core import cache

# Cache API responses
api_cache = Path("output/api_cache.json")

# Cache with 30-minute TTL
api_data = cache.get_json_cache(api_cache, max_age_seconds=1800)
if api_data is None:
    # Fetch fresh data
    api_data = fetch_from_api()
    cache.set_json_cache(api_cache, api_data)

process_data(api_data)
```

### Cache Management and Monitoring
```python
from metainformant.core import cache

# Get cache statistics
info = cache.get_cache_info(cache_dir)
print(f"Cache contains {info['total_files']} files")
print(f"Total size: {info['total_size'] / 1024:.1f} KB")

# Clear old cache files
cache.clear_cache_dir(cache_dir)
print("Cache cleared")
```

## Cache Key Safety

Cache keys are automatically sanitized:
- Dangerous characters are replaced with safe alternatives
- Hierarchical keys use "/" separators
- File extensions are automatically managed (.json/.gz)

```python
# These keys become safe filenames:
"api/v1/users"     # -> api/v1/users.json
"analysis/PCA"     # -> analysis/PCA.json
"weird<key>name"   # -> weird_key_name.json
```

## Integration Patterns

### Expensive Computation Caching
```python
from metainformant.core import cache

def cached_computation(cache_dir, key, func, *args, ttl_seconds=3600):
    """Cache wrapper for expensive functions."""
    # Try to load from cache first
    result = cache.load_cached_json(cache_dir, key, ttl_seconds=ttl_seconds)
    if result is not None:
        return result

    # Compute and cache
    result = func(*args)
    cache.cache_json(cache_dir, key, result)
    return result

# Usage
result = cached_computation(
    cache_dir, "ml/model_training", train_model, X_train, y_train
)
```

### API Response Caching
```python
from metainformant.core import cache

class CachedAPIClient:
    def __init__(self, cache_dir, ttl_seconds=1800):
        self.cache_dir = Path(cache_dir)
        self.ttl = ttl_seconds

    def get(self, endpoint):
        # Try cache first
        cached = cache.load_cached_json(self.cache_dir, f"api/{endpoint}", ttl_seconds=self.ttl)
        if cached is not None:
            return cached

        # Fetch and cache
        response = requests.get(f"https://api.example.com/{endpoint}")
        data = response.json()
        cache.cache_json(self.cache_dir, f"api/{endpoint}", data)
        return data
```

## Performance Considerations

### TTL Strategy
- **Short TTL (seconds/minutes)**: API responses, temporary results
- **Medium TTL (hours)**: Computation results, analysis outputs
- **Long/No TTL (days/never)**: Reference data, static lookups

### Cache Size Management
```python
# Monitor cache size
info = cache.get_cache_info(cache_dir)
if info['total_size'] > 100 * 1024 * 1024:  # 100MB
    # Clear old cache
    cache.clear_cache_dir(cache_dir)
```

### Memory vs Disk Trade-offs
- JSON caching is disk-based (persistent across runs)
- Consider memory caching for frequently accessed data
- Use compression (.gz) for large cached objects

## Error Handling

Cache operations are designed to fail gracefully:
- Missing cache files return None/default values
- Corrupted cache files are treated as misses
- File system errors don't crash the application
- Automatic retry logic for transient failures

## Dependencies

- **Required**: Standard library `time`, `pathlib`
- **Optional**: None (uses core.io and core.text modules)
