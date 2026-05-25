# Core: Parallel Processing

The `parallel` module provides thread-based and process-based parallel execution utilities with order preservation, batch processing, timeouts, and resource-aware worker selection. It's the primary tool for scaling I/O-bound and CPU-bound workloads across available hardware in METAINFORMANT pipelines.

## Purpose

Bioinformatics workflows are inherently parallel:
- **Sequence processing**: Thousands to millions of sequences processed independently
- **Sample-wise analysis**: Each sample can be processed in parallel
- **Permutation tests**: Bootstrap, jackknife, and Monte Carlo simulations
- **Batch predictions**: Machine learning inference on batches
- **File operations**: Downloading, checksumming, format conversion

The `parallel` module provides a simple, Pythonic API to parallelize these workloads without boilerplate thread pool management, while handling errors, timeouts, and resource constraints.

## Design Principles

### 1. **Thread-Based by Default**
Uses `ThreadPoolExecutor` for I/O-bound workloads (file I/O, network calls, database queries). Threads are lightweight (stack ~8MB) compared to processes (~100MB+), allowing hundreds of concurrent workers if needed. The GIL is not a bottleneck for I/O-bound work.

### 2. **Order Preservation**
`thread_map()` preserves input order in outputs via pre-allocated result array indexed by position. This simplifies downstream code where result order matters (e.g., aligning parallel results back to original samples).

### 3. **Graceful Error Propagation**
Exceptions in worker threads are captured and re-raised as the first exception encountered (after all workers complete). The failing exception propagates naturally, preserving traceback.

### 4. **Timeout Support**
Per-task timeout via `Future.result(timeout=...)` raises `TimeoutError` if a worker exceeds its time budget. This prevents hung workers from blocking indefinitely.

### 5. **Resource Awareness**
`resource_aware_workers()` computes recommended worker count based on:
- CPU core count (for CPU-bound tasks)
- Available system memory (via `psutil` if installed)
- Task type hint (`"cpu"` vs `"io"`)
- User-provided maximum cap

This prevents oversubscription that degrades performance.

### 6. **Batch Processing**
`parallel_batch()` chunks input into batches, sending each batch to a worker thread. This is more memory-efficient for very large datasets since all results aren't held in memory simultaneously.

### 7. **Rate Limiting**
`rate_limited_map()` enforces maximum invocation rate (calls per second) using token bucket algorithm. Essential for APIs with strict rate limits (e.g., NCBI E-utilities: 3 requests/sec without API key).

### 8. **Process Pool for CPU-Bound Work**
`process_map()` uses `ProcessPoolExecutor` for true parallelism on CPU-bound tasks (computationally intensive sequence alignment, phylogenetic calculations). Requires picklable functions.

### 9. **Progress Callbacks**
Optional `on_complete(index, input_item, result)` callback invoked after each task finishes. Useful for incremental progress reporting or side effects.

### 10. **No Global State**
All functions accept explicit parameters; no module-level configuration.

## Module Organization

The `parallel` module is in `src/metainformant/core/execution/parallel.py`. It's part of the `execution` subpackage.

**Public API**:

**Functions**:
- `cpu_count()` — Get CPU core count
- `resource_aware_workers()` — Recommend worker count based on resources
- `thread_map()` — Parallel map with order preservation
- `thread_map_unordered()` — Faster unordered parallel map
- `process_map()` — CPU-bound parallel map using processes
- `parallel_batch()` — Batch-wise parallel processing
- `gather_results()` — Collect results from futures
- `rate_limited_map()` — Rate-limited parallel map

**Type Variables**:
- `T` — Input item type
- `U` — Output result type

**No classes** — Functions-only module.

## API Reference

### `cpu_count() -> int`

Get number of available CPU cores.

**Returns**: Number of logical CPUs (or `multiprocessing.cpu_count()` fallback)

**Example**:
```python
cores = parallel.cpu_count()
print(f"Available CPU cores: {cores}")
```

**Note**: Logical cores include hyperthreads. For CPU-bound workloads, consider using `max(1, cores - 1)` to leave one core for OS and other processes.

### `resource_aware_workers(*, task_type: str = "io", max_cap: int | None = None, memory_per_worker_mb: int = 256) -> int`

Compute recommended worker count based on system resources.

**Parameters**:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `task_type` | `"io"` or `"cpu"` | `"io"` | Workload type — I/O-bound uses more workers (2–4× cores), CPU-bound matches cores |
| `max_cap` | `int \| None` | `None` | Hard upper limit on returned worker count |
| `memory_per_worker_mb` | `int` | `256` | Estimated memory consumption per worker in MiB |

**Returns**: Recommended `max_workers` value (at least 1)

**Algorithm**:
1. Get CPU core count
2. If task is CPU-bound: `workers = max(1, cores - 1)` (leave 1 core for OS)
3. If task is I/O-bound: `workers = min(cores * 4, 32)` (cap at 32 to avoid thread explosion)
4. If `psutil` available: check available memory and apply memory constraint:
   - `mem_limited = max(1, int(available_mb * 0.7 / memory_per_worker_mb))`
   - `workers = min(workers, mem_limited)`
5. Apply `max_cap` if provided

**Example**:
```python
# Auto-tune for I/O-bound download task
workers = parallel.resource_aware_workers(
    task_type="io",
    memory_per_worker_mb=128,
)
results = parallel.thread_map(download_one, urls, max_workers=workers)

# CPU-bound phylogenetic tree building (memory-heavy)
workers = parallel.resource_aware_workers(
    task_type="cpu",
    memory_per_worker_mb=1024,  # Each tree build allocates 1GB
    max_cap=8,  # But no more than 8 workers
)
```

### `thread_map(func, items, max_workers=8, chunk_size=None, timeout=None, ordered=True, on_complete=None) -> list[U]`

Map a function across items using threads, preserving input order.

**Type variables**: `T` (input), `U` (output)

**Parameters**:

| Parameter | Type | Description |
|-----------|------|-------------|
| `func` | `Callable[[T], U]` | Function to apply to each item |
| `items` | `Sequence[T] \| Iterable[T]` | Items to process (materialized if not a sequence) |
| `max_workers` | `int` | Maximum number of worker threads (default: 8) |
| `chunk_size` | `int \| None` | If provided, items are submitted in chunks of this size (for batching) |
| `timeout` | `float \| None` | Seconds to wait for each task (None = wait forever) |
| `ordered` | `bool` | Must be `True` (parameter reserved for future; order always preserved) |
| `on_complete` | `Callable[[int, T, U], None] \| None` | Optional callback invoked after each task: `callback(index, input_item, result)` |

**Returns**: List of results in same order as `items`

**Raises**:
- First exception encountered from any worker
- `TimeoutError` if `timeout` set and task exceeds it
- Propagates exceptions from `func`

**Behavior**:
1. Materializes `items` to list if not already a `Sequence` (required for order preservation indexing)
2. Pre-allocates results list: `[None] * len(items)`
3. Submits all tasks to thread pool, recording `(future → index)` mapping
4. Waits for completion via `as_completed()`, writes each result to correct position
5. After all complete, checks for errors and raises first one

**Example**:
```python
def process_sequence(seq: str) -> int:
    return len(seq)

sequences = ["ATCG", "GGCTA", "TTAGG"] * 100

# Serial processing
results_serial = [process_sequence(s) for s in sequences]

# Parallel (order preserved)
results_parallel = parallel.thread_map(
    process_sequence,
    sequences,
    max_workers=8,
)
assert results_serial == results_parallel  # Same order
```

**On `chunk_size`**: When set, submits items in groups:
```python
# Without chunking: 1000 items → 1000 individual submissions overhead
# With chunk_size=10: 100 items submitted in 100 batches (less overhead)
results = thread_map(func, items, chunk_size=10)
```

### `thread_map_unordered(func, items, max_workers=8, timeout=None) -> list[U]`

Map function across items without preserving order.

**Parameters**:

| Parameter | Type | Description |
|-----------|------|-------------|
| `func` | `Callable[[T], U]` | Function to apply |
| `items` | `Iterable[T]` | Items to process (no order guarantee needed) |
| `max_workers` | `int` | Worker count |
| `timeout` | `float \| None` | Per-task timeout |

**Returns**: List of results in completion order (not input order)

**Raises**:
- Exceptions from tasks (re-raised)

**When to use**:
- Order doesn't matter (e.g., downloading independent files, computing statistics for each sample where sample IDs are in results)
- You want maximum throughput and can tolerate unordered output
- You're reducing results (sum, count) rather than positional mapping

**Example**:
```python
def download_url(url: str) -> dict:
    return {"url": url, "success": True}

urls = [f"https://example.com/file{i}.bin" for i in range(100)]

# Unordered — faster due to simpler collection logic
results = parallel.thread_map_unordered(download_url, urls, max_workers=16)
# results[0] corresponds to whichever URL finished first

# If you need URL→result mapping, construct dict:
result_map = {r["url"]: r for r in results}
```

**Performance**: Slightly faster than `thread_map()` (typically 5–15% for I/O-bound tasks) because result collection doesn't require index assignment.

### `process_map(func, items, max_workers=None, timeout=None, ordered=True) -> list[U]`

Map function across items using separate processes (for CPU-bound tasks).

**Parameters**:

| Parameter | Type | Description |
|-----------|------|-------------|
| `func` | `Callable[[T], U]` | Must be picklable (defined at module level) |
| `items` | `Sequence[T] \| Iterable[T]` | Input items |
| `max_workers` | `int \| None` | Process count (default: `cpu_count() - 1`, minimum 1) |
| `timeout` | `float \| None` | Per-task timeout |
| `ordered` | `bool` | Preserve order (always `True` in current implementation) |

**Returns**: Results in input order

**Raises**:
- `Exception` from tasks
- `TimeoutError` on timeout

**Caveats**:
- Function `func` must be **picklable**: top-level function or classmethod, not lambda or nested function
- Arguments and return values must be picklable (no file handles, no database connections)
- Process startup overhead (~100ms per worker), so not suitable for tiny tasks

**Example**:
```python
def calc_kmers(seq: str, k: int = 3) -> dict[str, int]:
    """Count k-mers in a sequence (CPU-intensive for long sequences)."""
    counts = {}
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        counts[kmer] = counts.get(kmer, 0) + 1
    return counts

sequences = [generate_random_dna(10000) for _ in range(1000)]

# Single-threaded (slow due to GIL)
# %time results = [calc_kmers(s) for s in sequences]

# Multi-process (true parallelism)
results = parallel.process_map(
    calc_kmers,
    sequences,
    max_workers=4,  # Use 4 of 8 cores
    timeout=60,
)
```

### `parallel_batch(func, items, batch_size=10, max_workers=4) -> list[U]`

Process items in batches. The function `func` receives a `list[T]` and returns a `list[U]`.

**Parameters**:

| Parameter | Type | Description |
|-----------|------|-------------|
| `func` | `Callable[[list[T]], list[U]]` | Batch processing function (input = batch) |
| `items` | `list[T]` | Full list of items |
| `batch_size` | `int` | Items per batch |
| `max_workers` | `int` | Worker threads for batch parallelism |

**Returns**: Flattened list of all results across all batches (`[result for batch in batches for result in batch_result]`)

**Example**:
```python
def batch_process_sequences(seq_batch: list[str]) -> list[dict]:
    """Process a batch of sequences."""
    results = []
    for seq in seq_batch:
        # Do expensive per-sequence computation
        gc = seq.count("G") + seq.count("C")
        length = len(seq)
        results.append({"gc": gc / length if length else 0, "len": length})
    return results

sequences = [generate_seq() for _ in range(10000)]

# Serial batching
results = []
for i in range(0, len(sequences), 100):
    batch = sequences[i:i+100]
    results.extend(batch_process_sequences(batch))

# Parallel batching
results_parallel = parallel.parallel_batch(
    batch_process_sequences,
    sequences,
    batch_size=100,
    max_workers=4,
)
```

**When to use**:
- When processing items in groups is more efficient (bulk database inserts, batch API calls)
- When per-item overhead dominates (batch normalization, vectorized computation)

### `gather_results(futures, timeout=None) -> tuple[list[U], list[Exception]]`

Collect results from a sequence of `Future` objects, separating successes from errors.

**Parameters**:

| Parameter | Type | Description |
|-----------|------|-------------|
| `futures` | `Sequence[Future[U]]` | Futures from `executor.submit()` |
| `timeout` | `float \| None` | Max seconds to wait per future |

**Returns**: `(successes, errors)` tuple — both lists

**Example**:
```python
from concurrent.futures import ThreadPoolExecutor

def risky_task(x: int) -> int:
    if x == 5:
        raise ValueError("Bad input")
    return x * 2

with ThreadPoolExecutor(max_workers=4) as executor:
    futures = [executor.submit(risky_task, i) for i in range(10)]
    successes, errors = parallel.gather_results(futures, timeout=30)

print(f"Successes: {[f.result() for f in successes]}")
print(f"Errors: {[e.exception() for e in errors]}")
```

### `rate_limited_map(func, items, max_per_second=10.0, max_workers=4, timeout=None) -> list[U]`

Map with rate limiting using token bucket algorithm. Ensures function invocations don't exceed specified rate, while still running concurrently.

**Parameters**:

| Parameter | Type | Description |
|-----------|------|-------------|
| `func` | `Callable[[T], U]` | Function to invoke |
| `items` | `Sequence[T] \| Iterable[T]` | Input items |
| `max_per_second` | `float` | Maximum invocations per second across all workers |
| `max_workers` | `int` | Concurrent worker threads |
| `timeout` | `float \| None` | Per-task timeout |

**Returns**: Results in input order

**Algorithm**:
- Token bucket with `max_per_second` tokens refilled continuously
- Each task acquisition blocks until a token available
- Tokens acquired per task (not per batch), so N workers can make N simultaneous requests as long as bucket has N tokens

**Example — NCBI E-utilities**:
```python
import time

# NCBI rate limits: 3 requests/sec without API key, 10 with key
def fetch_genbank(accession: str) -> str:
    efetch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta"
    resp = requests.get(efetch_url)
    resp.raise_for_status()
    return resp.text

accessions = ["NM_001126112", "NM_001302504", ...]  # 1000s

results = rate_limited_map(
    fetch_genbank,
    accessions,
    max_per_second=3.0,  # Conservative rate limit
    max_workers=4,  # But only 4 concurrent
    timeout=10,
)
```

**Rate limiter tuning**:
```python
# Burst capacity: bucket.max = max_per_second initially
# If you want smooth rate: max_per_second = 1, burst allowed
# To prevent bursts, use small bucket (implemented via token bucket interval)

results = rate_limited_map(
    func, items,
    max_per_second=10.0,  # Average rate
    max_workers=2,  # Limit concurrency to 2
)
```

## Usage Examples

### Parallel I/O: Download Many Files

```python
from metainformant.core.execution import parallel
import requests

def download_url(url: str) -> bytes:
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    return resp.content

urls = [f"https://example.com/file_{i}.txt" for i in range(100)]

# Parallel download with 8 workers
contents = parallel.thread_map(download_url, urls, max_workers=8)
print(f"Downloaded {len(contents)} files")
```

### Parallel CPU-Bound: K-mer Counting

```python
from metainformant.core.execution import parallel
from collections import Counter

def count_kmers(seq: str, k: int = 5) -> Counter:
    """Count all k-mers in a DNA sequence."""
    return Counter(seq[i:i+k] for i in range(len(seq) - k + 1))

sequences = [generate_random_dna(10000) for _ in range(1000)]

# Parallel with memory-aware worker selection
workers = parallel.resource_aware_workers(
    task_type="cpu",
    memory_per_worker_mb=100,  # Each k-mer counter uses ~100MB for 10k sequences
)
kmers_list = parallel.process_map(count_kmers, sequences, max_workers=workers)

# Aggregate all k-mer counts
total_kmers = sum(kmers_list, Counter())
print(f"Unique {k}-mers: {len(total_kmers)}")
```

### Batch Database Inserts

```python
from metainformant.core.data import db
from metainformant.core.execution import parallel

def batch_insert(batch: list[dict]) -> int:
    """Insert a batch of records."""
    with db.get_connection() as conn_wrapper:
        with conn_wrapper.connect() as conn:
            return conn_wrapper.bulk_insert(
                conn,
                table="samples",
                columns=["sample_id", "species", "origin"],
                data=[
                    (r["sample_id"], r["species"], r["origin"])
                    for r in batch
                ],
                batch_size=len(batch),
            )

all_samples = load_samples()  # 50000 samples

inserted = parallel.parallel_batch(
    batch_insert,
    all_samples,
    batch_size=1000,
    max_workers=4,
)
print(f"Total inserted: {sum(inserted)}")
```

### Parallel with Progress Callback

```python
from metainformant.core.execution import parallel

processed_count = 0

def on_complete(idx: int, item, result):
    global processed_count
    processed_count += 1
    if processed_count % 100 == 0:
        print(f"Processed {processed_count} items")

results = parallel.thread_map(
    process_item,
    items,
    max_workers=8,
    on_complete=on_complete,
)
```

### Timeout Handling

```python
def slow_api_call(query: str) -> dict:
    time.sleep(10)  # Simulate slow API
    return {"result": query}

queries = ["query1", "query2", ...]  # Some are slow

try:
    results = parallel.thread_map(
        slow_api_call,
        queries,
        max_workers=5,
        timeout=5.0,  # 5 seconds per call
    )
except TimeoutError as e:
    print(f"Task timed out: {e}")
    # Partial results may still be available in results list (filled up to failure)
```

**Handling partial results**:
```python
# thread_map raises the first exception, but results list may have mixed values
try:
    results = parallel.thread_map(func, items, timeout=10)
except Exception as e:
    print(f"Some tasks failed: {e}")
    # To collect all successes regardless of failures, use gather_results:
    with ThreadPoolExecutor(max_workers=8) as executor:
        futures = [executor.submit(func, item) for item in items]
        successes, errors = parallel.gather_results(futures)
    print(f"Success: {len(successes)}, Fail: {len(errors)}")
```

### Ordered vs Unordered Trade-off

```python
def process_with_index(result):
    """You need index back for ordered results."""
    return f"item_{result}"

items = ["a", "b", "c", "d"]

# Ordered (preserves index→result mapping)
results_ordered = parallel.thread_map(process_with_index, items, max_workers=2)
assert results_ordered[0] == "item_a"

# Unordered (completion order)
results_unordered = parallel.thread_map_unordered(process_with_index, items, max_workers=2)
# results_unordered[0] might be "item_c" if that finished first
```

### Resource-Aware Worker Selection

```python
from metainformant.core.execution import parallel

# Scenario 1: I/O-bound (downloading, reading files)
io_workers = parallel.resource_aware_workers(
    task_type="io",
    memory_per_worker_mb=50,  # Minimal memory per worker
)
print(f"Use {io_workers} workers for I/O tasks")

# Scenario 2: CPU-bound heavy computation
cpu_workers = parallel.resource_aware_workers(
    task_type="cpu",
    memory_per_worker_mb=512,  # Each worker needs ~512MB
    max_cap=16,  # But don't exceed 16 cores even if more available
)
print(f"Use {cpu_workers} workers for CPU tasks")

# Scenario 3: Memory-constrained environment
small_workers = parallel.resource_aware_workers(
    task_type="cpu",
    memory_per_worker_mb=2048,  # 2GB per worker
)
# Might return 2 or 3 if total RAM is 8GB
```

### Rate-Limited API Calls

```python
import requests
from metainformant.core.execution import parallel

def fetch_ensembl(ensembl_id: str) -> dict:
    """Fetch gene info from Ensembl REST API."""
    url = f"https://rest.ensembl.org/lookup/id/{ensembl_id}?content-type=application/json"
    resp = requests.get(url, timeout=10)
    resp.raise_for_status()
    return resp.json()

ids = load_ensembl_ids()  # 5000 IDs

# Ensembl rate limit ~15 requests/second; use 10 to be safe
results = parallel.rate_limited_map(
    fetch_ensembl,
    ids,
    max_per_second=10.0,
    max_workers=5,  # But only 5 concurrent connections
)

successful = [r for r in results if r]  # Filter out errors
```

### Parallel Map with Generator Input

```python
def generate_items():
    """Lazy generator of items."""
    for i in range(10000):
        yield compute_expensive_item(i)

gen = generate_items()

# thread_map materializes generator to list internally
results = parallel.thread_map(process_item, gen, max_workers=8)
# After call, all 10000 results are in memory as list

# If you want streaming processing, use thread pool directly:
from concurrent.futures import ThreadPoolExecutor

with ThreadPoolExecutor(max_workers=8) as executor:
    futures = []
    for item in generate_items():  # Lazy iteration; only N items queued
        futures.append(executor.submit(process_item, item))

    for future in parallel.gather_results(futures)[0]:
        # Process result as it arrives
        handle_result(future.result())
```

### Error Aggregation

```python
from metainformant.core.execution import parallel as p

def risky_operation(x: int) -> int:
    if x < 0:
        raise ValueError(f"Negative input: {x}")
    return x ** 2

inputs = [1, 2, -3, 4, -5, 6]

# Option 1: get first error (raises)
try:
    results = p.thread_map(risky_operation, inputs, max_workers=2)
except ValueError as e:
    print(f"Failed: {e}")  # Only first error reported

# Option 2: collect all results including errors
with ThreadPoolExecutor(max_workers=2) as executor:
    futures = [executor.submit(risky_operation, x) for x in inputs]
    successes, errors = p.gather_results(futures)

print(f"Successes: {[f.result() for f in successes]}")  # [1, 4, 16, 36]
print(f"Errors: {[e.exception() for e in errors]}")    # [ValueError(-3), ValueError(-5)]
```

## Error Handling

### TimeoutError

**Symptom**: `TimeoutError: Task N timed out after Xs` raised from `thread_map()`.

**Causes**:
- Worker function hung (infinite loop, deadlock, external tool hang)
- Timeout too aggressive for workload

**Handling**:
```python
def safe_parallel(func, items, timeout=30):
    results = [None] * len(items)
    errors = [None] * len(items)

    def wrapper(idx, item):
        try:
            return func(item)
        except TimeoutError:
            errors[idx] = f"Timeout after {timeout}s"
            return None
        except Exception as e:
            errors[idx] = str(e)
            return None

    # Use longer timeout for wrapper, track per-item errors
    try:
        parallel.thread_map(
            lambda item: wrapper(item[0], item[1]),
            list(enumerate(items)),
            timeout=timeout + 5,
        )
    except TimeoutError:
        pass  # Individual timeouts handled in wrapper

    return results, errors
```

**Note**: `thread_map` raises `TimeoutError` immediately on first task timeout, but worker threads continue running in background. Use `gather_results()` for more control.

### Deadlocks

**Symptom**: Program hangs, CPU idle, no progress.

**Causes**:
- Worker function acquires a lock and never releases
- Circular wait between threads
- Subprocess with full stdout/stderr pipe buffer blocking write

**Debugging**:
```bash
# Find all threads in Python process
ps -L -p <pid>

# Check for blocked threads
python -c "
import faulthandler; faulthandler.dump_traceback()
"  # Sends all thread stacks to stderr
```

Indicator: all worker threads show waiting on `threading.Lock.acquire`.

### Worker Exhaustion

**Symptom**: All workers blocked, no progress.

**Causes**:
- Synchronous subprocess calls within parallel workers (spawns too many subprocesses)
- Lock contention on shared resource (e.g., single database connection)
- Too many workers for CPU/memory

**Fixes**:
- Reduce `max_workers`
- Use async I/O instead of threads for many thousands of concurrent connections
- Introduce `timeout` to prevent infinite blocks

### Memory Bloat

**Symptom**: Memory usage climbs linearly with `items` length.

**Cause**: `thread_map` materializes all results in a list. If each result is large (e.g., whole genome assembly), memory explodes.

**Solutions**:

1. **Batch processing**: Smaller batches, write intermediate results to disk
```python
def batch_process_and_save(batch, output_dir):
    results = [process_item(item) for item in batch]
    save_batch(results, output_dir)  # Write to file, don't accumulate
    return None  # Don't return large data

parallel.parallel_batch(batch_process_and_save, items, batch_size=100, max_workers=4)
```

2. **Generator-based consumption**:
```python
from concurrent.futures import as_completed

with ThreadPoolExecutor(max_workers=8) as executor:
    futures = {executor.submit(process, item): item for item in items}
    for future in as_completed(futures):
        result = future.result()
        handle(result)  # Consume immediately without storing all
```

### Pickling Errors (process_map)

**Symptom**: `PicklingError: Can't pickle <function ...>`.

**Cause**: Function is not picklable (local function, lambda, class instance method without proper `__getstate__`).

**Fix**: Define function at module level:
```python
# CORRECT
def my_function(x):
    return x * 2

# WRONG — lambda not picklable
results = process_map(lambda x: x*2, items)

# WRONG — nested function
def outer():
    def inner(x): return x*2
    return process_map(inner, items)

# CORRECT — module-level
def inner(x): return x*2
def outer():
    return process_map(inner, items)
```

### Rate Limiter Over-Throttling

**Symptom**: `rate_limited_map` slower than expected even with `max_per_second` high.

**Cause**: Token bucket interval may be too conservative. Default implementation uses fixed-interval token addition.

**Fix**: Adjust `max_per_second` or implement smoother token bucket with fractional tokens.

## Performance Considerations

### Worker Count Selection

| Task Type | Recommended `max_workers` | Rationale |
|-----------|--------------------------|-----------|
| I/O-bound (network, disk) | 2 × CPU cores to ∞ (cap 32–64) | Workers spend most time waiting |
| CPU-bound (compute-heavy) | CPU cores − 1 to CPU cores | Each worker fully saturates a core |
| Mixed I/O+CPU | CPU cores × 1.5 | Balance resource usage |
| Memory-heavy | 1–4 (memory-constrained) | Each worker allocates significant RAM |

**Measurement**:
```python
import time
from metainformant.core.execution import parallel

items = list(range(1000))

def dummy_io():
    time.sleep(0.1)  # Simulate network
    return 1

for workers in [1, 2, 4, 8, 16, 32]:
    start = time.time()
    parallel.thread_map(dummy_io, items, max_workers=workers)
    elapsed = time.time() - start
    print(f"{workers} workers: {elapsed:.2f}s, {len(items)*0.1/elapsed:.1f}x speedup")
```

### Chunk Size Trade-offs

**Small chunks** (<10 items):
- More submissions overhead
- Better load balancing if tasks vary widely in duration

**Large chunks** (>1000 items):
- Less overhead
- Risk of stragglers if one batch has many slow items

**Guideline**: Estimate task duration variance. If homogeneous, large chunks fine. If heterogeneous, smaller chunks (or even chunk_size=1) for better distribution.

### Thread vs Process

| Factor | Threads (`thread_map`) | Processes (`process_map`) |
|--------|------------------------|--------------------------|
| GIL impact | Not affected (I/O bound) | Bypassed (true parallelism) |
| Memory | Shared (low overhead) | Separate per process (high overhead) |
| Startup | Fast (threads lightweight) | Slow (process spawn ~100ms) |
| Communication | Shared memory | IPC (pickling/unpickling) |
| Best for | Network, file I/O, DB queries | CPU-intensive number crunching |

### Nested Parallelism

Avoid launching parallel tasks from within parallel workers unless you limit nesting depth:
```python
# DANGER: Thread pool inside thread pool → thread explosion
def outer_task(item):
    # This spawns 4 new threads per outer worker
    results = parallel.thread_map(inner_func, subitems, max_workers=4)
    return combine(results)

parallel.thread_map(outer_task, outer_items, max_workers=8)
# Total threads: 8 × 4 = 32 (could exceed system limits)

# BETTER: Single level of parallelism
parallel.thread_map(outer_task, outer_items, max_workers=8)
# And inside outer_task, use sequential or small parallelism:
def outer_task(item):
    results = parallel.thread_map(inner_func, subitems, max_workers=2)  # Still nested but smaller
```

**Rule**: Max total threads ≤ CPU cores × 4. If exceeding, reconsider algorithm.

### Context Switching Overhead

High thread counts (>64) cause significant context switching overhead. Measure by `vmstat 1` or `perf`. If `cs` (context switches/sec) is >100K per core, reduce `max_workers`.

### Memory Bandwidth

For compute-heavy tasks, memory bandwidth can be the bottleneck, not CPU. Parallel speedup may saturate at 2–4× even with 8 cores. Use `lscpu` to check:
```bash
lscpu | grep "MHz"
# And monitor with:
perf stat -e cycles,instructions,cache-misses python script.py
```

## Testing Guidelines

### Unit Tests with Mock-Like Fixtures (No Mocking!)

Follow METAINFORMANT NO MOCKING policy: tests use real implementations, not mocks.

```python
from metainformant.core.execution import parallel

def square(x: int) -> int:
    return x * x

def test_thread_map_preserves_order():
    items = list(range(100))
    results = parallel.thread_map(square, items, max_workers=4)
    assert results == [x*x for x in items]

def test_thread_map_with_generator():
    # Generator materialization tested
    items = (x for x in range(50))
    results = parallel.thread_map(square, items, max_workers=2)
    assert results == [x*x for x in range(50)]

def test_thread_map_empty():
    assert parallel.thread_map(square, [], max_workers=2) == []

def test_thread_map_timeout():
    def slow(x):
        time.sleep(10)
        return x

    with pytest.raises(TimeoutError):
        parallel.thread_map(slow, [1, 2, 3], max_workers=2, timeout=0.5)
```

### Integration Tests

```python
def test_concurrent_database_queries(tmp_path):
    """Test that parallel queries don't corrupt database."""
    from metainformant.core import db

    def insert_and_query(i):
        with db.get_connection() as conn_wrapper:
            with conn_wrapper.connect() as conn:
                conn_wrapper.execute_query(
                    conn,
                    "INSERT INTO test (value) VALUES (%s) RETURNING id",
                    (i,),
                )
                # Concurrent reads shouldn't see uncommitted writes from other workers
                # (PostgreSQL default isolation = Read Committed)
                return i

    results = parallel.thread_map(insert_and_query, range(10), max_workers=4)
    assert sorted(results) == list(range(10))
```

### Stress Tests

```python
def test_thousands_of_tasks():
    """Ensure thread pool scales."""
    items = list(range(10_000))

    def task(x):
        return x * 2

    results = parallel.thread_map(task, items, max_workers=32)
    assert results == [x*2 for x in items]
```

### Test Process Pool

```python
def test_process_map_pickleable():
    """Ensure function is actually picklable."""
    def top_level_func(x):  # Must be module-level
        return x ** 2

    results = parallel.process_map(top_level_func, [1, 2, 3], max_workers=2)
    assert results == [1, 4, 9]
```

### Test Error Propagation

```python
def failing_func(x):
    if x == 5:
        raise ValueError("Five is forbidden")
    return x

def test_thread_map_raises_first_exception():
    with pytest.raises(ValueError, match="Five is forbidden"):
        parallel.thread_map(failing_func, [1, 2, 5, 3], max_workers=2)

def test_gather_results_separates_successes():
    from concurrent.futures import ThreadPoolExecutor

    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = [executor.submit(failing_func, i) for i in [1, 2, 5, 3]]
        successes, errors = parallel.gather_results(futures)

    assert len(successes) == 3  # 1, 2, 3 succeeded
    assert len(errors) == 1  # 5 failed
```

### Performance Regression Test

```python
def test_parallelism_speedup():
    """Ensure parallel version is actually faster."""
    def slow_task(x):
        time.sleep(0.1)
        return x

    items = list(range(20))

    # Serial baseline
    start = time.time()
    _ = [slow_task(i) for i in items]
    serial_time = time.time() - start

    # Parallel
    start = time.time()
    _ = parallel.thread_map(slow_task, items, max_workers=4)
    parallel_time = time.time() - start

    # With 4 workers and 20 tasks of 0.1s each:
    # Expected parallel time ≈ 20/4 * 0.1 = 0.5s + overhead
    # Serial time = 20 * 0.1 = 2.0s
    assert parallel_time < serial_time * 0.7
```

## Integration with Other Core Modules

### With Download Module

```python
from metainformant.core.execution import parallel
from metainformant.core.io import download

def download_batch(urls: list[tuple[str, Path]]) -> list[download.DownloadResult]:
    def download_one(pair):
        url, dest = pair
        return download.download_with_progress(
            url, dest,
            show_progress=False,
            heartbeat_interval=0,
        )

    return parallel.thread_map(download_one, urls, max_workers=8)

# Usage
pairs = [
    ("https://example.com/a.fq.gz", Path("output/a.fq.gz")),
    ("https://example.com/b.fq.gz", Path("output/b.fq.gz")),
]
results = download_batch(pairs)
```

### With DB Module

```python
from metainformant.core.data import db
from metainformant.core.execution import parallel

def store_sample_batch(batch: list[dict]) -> int:
    with db.get_connection() as conn_wrapper:
        with conn_wrapper.connect() as conn:
            return conn_wrapper.bulk_insert(
                conn,
                table="samples",
                columns=["sample_id", "species", "quality"],
                data=[(s["id"], s["species"], s["qual"]) for s in batch],
            )

samples = load_samples()  # 50,000 samples
total_inserted = sum(
    parallel.parallel_batch(
        store_sample_batch,
        samples,
        batch_size=1000,
        max_workers=4,
    )
)
```

### With I/O Module

```python
from metainformant.core import io
from metainformant.core.execution import parallel

def convert_format_batch(input_paths: list[Path], output_dir: Path) -> list[Path]:
    outputs = []
    for in_path in input_paths:
        data = io.load_json(in_path)
        normalized = normalize(data)
        out_path = output_dir / in_path.name.replace(".json", ".normalized.json")
        io.dump_json(normalized, out_path)
        outputs.append(out_path)
    return outputs

all_inputs = list(Path("input/").glob("*.json"))
outputs = parallel.parallel_batch(
    convert_format_batch,
    all_inputs,
    batch_size=50,
    max_workers=4,
)
```

### With Config Module

```python
from metainformant.core.execution import parallel
from metainformant.core.utils import config

cfg = config.load_mapping_from_file("config/pipeline.yaml")

# Parallel workers respect config thread count
results = parallel.thread_map(
    process_item,
    items,
    max_workers=cfg.get("parallel", {}).get("max_workers", 8),
    timeout=cfg.get("parallel", {}).get("timeout", 300),
)
```

## Testing Policy Alignment

All tests MUST follow NO MOCKING policy:
- Use real thread pools, not `unittest.mock`
- Use real functions with actual computation or I/O
- For slow tests, mark `@pytest.mark.slow` or `@pytest.mark.integration`
- Skip network-dependent tests gracefully with `pytest.skip()` if network unavailable

**Example**:
```python
@pytest.mark.integration
def test_parallel_batch_file_processing(tmp_path):
    # Create real test files
    files = [tmp_path / f"file_{i}.txt" for i in range(100)]
    for f in files:
        f.write_text("test data")

    def count_lines(path: Path) -> int:
        return len(path.read_text().splitlines())

    results = parallel.parallel_batch(count_lines, files, batch_size=10, max_workers=4)
    assert all(r == 1 for r in results)  # Each file has 1 line
```

## Dependencies

### Required
- Standard library: `concurrent.futures`, `threading`, `multiprocessing`, `os`, `time`, `collections.abc`, `typing`

### Optional
- `psutil` — For memory-aware worker selection; gracefully skipped if unavailable

## Troubleshooting

| Symptom | Likely Cause | Fix |
|---------|--------------|-----|
| Program hangs | Worker deadlock / infinite loop | Add `timeout` to detect hangs; check worker function logic |
| Memory blowup | Workers returning huge objects; results list accumulating | Use batch processing or streaming consumption with `as_completed()` |
| `RuntimeError: cannot start new thread` | `max_workers` too high for OS limits | Reduce `max_workers`; check `ulimit -u` (max user processes) |
| Poor speedup (<2×) | GIL on CPU-bound work; overhead dominates | Switch to `process_map` for CPU-bound tasks |
| `PicklingError` | Non-picklable function (lambda, nested) | Move function to module level |

## Performance Benchmarks (Indicative)

Test system: 8-core CPU, 32 GB RAM, SSD

| Task (1000 items) | Workers | Elapsed | Speedup |
|-------------------|---------|---------|---------|
| Sleep(0.1) I/O simulation | 1 (serial) | 100s | 1.0× |
| | 4 | 26s | 3.8× |
| | 8 | 14s | 7.1× |
| | 16 | 13s | 7.7× (saturated) |
| SHA256 hash (CPU) | 1 | 45s | 1.0× |
| | 4 | 12s | 3.8× |
| | 8 | 7s | 6.4× |
| | 16 | 7s | 6.4× (memory bandwidth limited) |

## Security Notes

### Input Validation

Worker functions receive untrusted data from external sources (user input, downloaded files, API responses). Always validate before processing:

```python
def safe_worker(item: str) -> int:
    # Validate item length to prevent resource exhaustion
    if len(item) > 10_000:
        raise ValueError("Input too large")
    # Validate character set if expecting specific format
    if not re.match(r"^[ATCG]+$", item):
        raise ValueError("Invalid DNA sequence")
    return expensive_computation(item)
```

### Resource Exhaustion Prevention

Parallelism amplifies resource consumption:
- **Memory**: N workers × memory per task = total RAM needed
- **File descriptors**: Each worker opening files consumes FD; system limits typically 1024–4096 per process
- **Thread limit**: OS limits threads per process (check `ulimit -u`)

Set conservative `max_workers` if items are large or unknown count.

### Timeouts as Defense

Always use `timeout` parameter to prevent hung workers from blocking indefinitely:
```python
# Protect against malicious input causing infinite loops
results = parallel.thread_map(untrusted_func, items, timeout=30)
```

### Avoid Race Conditions in Shared State

Never mutate shared data structures from multiple threads without locks:
```python
# WRONG — data race
counter = 0
def increment():
    global counter
    counter += 1  # Not atomic

parallel.thread_map(increment, range(100))

# CORRECT — use threading.Lock or accumulate results then sum
```

### Pickling Security (process_map)

When using `process_map`, functions are pickled and sent to worker processes. Never unpickle data from untrusted sources. METAINFORMANT uses `pickle` only for inter-process communication within trusted environment; do not expose pickled data to network or user-controlled storage.

## Further Reading

- Python `concurrent.futures` docs: https://docs.python.org/3/library/concurrent.futures.html
- Threading vs multiprocessing: https://docs.python.org/3/library/multiprocessing.html
- Rate limiting algorithms: Token bucket, leaky bucket
