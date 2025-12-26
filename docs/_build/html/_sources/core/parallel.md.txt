# Core: Parallel Processing

The `parallel` module provides thread-based parallel execution utilities with order preservation, batch processing, and timeout support.

## Functions

### Thread Mapping
- **`thread_map(func, items, max_workers=8, chunk_size=None, timeout=None, ordered=True)`** → `List[U]`
  - Map function across items using threads with order preservation
  - Configurable worker count and chunking
  - Optional timeout per task
  - Materializes input if not a sequence

- **`thread_map_unordered(func, items, max_workers=8, timeout=None)`** → `List[U]`
  - Map function across items without preserving order
  - Better performance for order-independent operations
  - Same timeout and worker configuration

### Batch Processing
- **`parallel_batch(func, items, batch_size=10, max_workers=4)`** → `List[U]`
  - Process items in batches for memory efficiency
  - Function receives entire batch, returns batch results
  - Useful for bulk operations and API calls

### System Information
- **`cpu_count()`** → `int`
  - Get available CPU core count
  - Fallback for systems without multiprocessing

## Usage Examples

### Basic Parallel Mapping
```python
from metainformant.core import parallel

def process_sequence(seq):
    # Simulate expensive computation
    return len(seq) * 2

sequences = ["ATCG", "GCTA", "TTAG", "CGAT"] * 100  # 400 items

# Process in parallel with 8 workers
results = parallel.thread_map(process_sequence, sequences, max_workers=8)
print(f"Processed {len(results)} sequences")
```

### Parallel Processing with Timeout
```python
from metainformant.core import parallel

def download_data(url):
    import time
    time.sleep(1)  # Simulate network delay
    return f"Downloaded {url}"

urls = [f"https://api.example.com/data/{i}" for i in range(50)]

# Process with 2-second timeout per request
results = parallel.thread_map(
    download_data,
    urls,
    max_workers=10,
    timeout=2.0
)
```

### Unordered Processing (Better Performance)
```python
from metainformant.core import parallel

def analyze_file(filepath):
    # File analysis that doesn't need ordered results
    return {"path": filepath, "size": len(filepath)}

files = [f"data/file_{i}.txt" for i in range(1000)]

# Faster unordered processing
results = parallel.thread_map_unordered(analyze_file, files, max_workers=8)
# Results order is not preserved but processing is faster
```

### Batch Processing for Efficiency
```python
from metainformant.core import parallel

def process_batch(batch):
    """Process multiple items as a batch."""
    results = []
    for item in batch:
        # Simulate batch processing (e.g., database inserts)
        results.append(f"processed_{item}")
    return results

items = list(range(1000))

# Process in batches of 50
results = parallel.parallel_batch(
    process_batch,
    items,
    batch_size=50,
    max_workers=4
)
print(f"Processed {len(results)} items in batches")
```

### CPU-Aware Worker Configuration
```python
from metainformant.core import parallel

# Use available CPU cores
cores = parallel.cpu_count()
optimal_workers = min(cores, 16)  # Cap at reasonable maximum

results = parallel.thread_map(
    expensive_function,
    data_items,
    max_workers=optimal_workers
)
```

## Error Handling

The parallel module provides robust error handling:

```python
from metainformant.core import parallel
import logging

logger = logging.getLogger(__name__)

def risky_operation(item):
    if item == "bad_data":
        raise ValueError(f"Invalid data: {item}")
    return f"processed_{item}"

items = ["good1", "bad_data", "good2"]

try:
    results = parallel.thread_map(risky_operation, items, max_workers=4)
    # Results will contain exceptions for failed items
    for i, result in enumerate(results):
        if isinstance(result, Exception):
            logger.error(f"Failed to process {items[i]}: {result}")
        else:
            print(f"Success: {result}")
except Exception as e:
    logger.error(f"Parallel processing failed: {e}")
```

## Performance Considerations

### Choosing the Right Function
- **`thread_map`**: Use when order matters and you have < 10,000 items
- **`thread_map_unordered`**: Use when order doesn't matter (better performance)
- **`parallel_batch`**: Use for memory-efficient batch processing

### Worker Count Optimization
```python
from metainformant.core import parallel

# Good: Scale with available cores
cores = parallel.cpu_count()
workers = min(cores * 2, 32)  # 2x cores, max 32

# Avoid: Hard-coded values
# workers = 8  # May be too many or too few
```

### Chunking for Large Datasets
```python
# For very large datasets, use chunking
large_data = list(range(100000))

results = parallel.thread_map(
    process_item,
    large_data,
    max_workers=8,
    chunk_size=1000  # Process 1000 items per task
)
```

## Timeouts and Cancellation

```python
from metainformant.core import parallel

def api_call(endpoint):
    # Simulate API call with potential hangs
    import time
    time.sleep(0.5)
    return f"result_from_{endpoint}"

endpoints = [f"/api/data/{i}" for i in range(100)]

# 5-second timeout per API call
results = parallel.thread_map(
    api_call,
    endpoints,
    max_workers=20,
    timeout=5.0
)
```

## Dependencies

- **Required**: Standard library `concurrent.futures`, `multiprocessing`
- **Optional**: None (signal module used for timeouts on Unix systems)
