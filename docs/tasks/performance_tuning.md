# Performance Tuning Quick Reference

Optimize METAINFORMANT for speed and memory efficiency across all modules.

## When to Use

Use `performance_tuning` when pipelines run too slowly or crash with `MemoryError` — not for basic usage (defaults optimized for 16GB RAM, 8 cores). Profile first with `python -m cProfile`.

## Table of Contents

- [Caching](#caching)
- [Parallel Processing](#parallel-processing)
- [Memory Optimization](#memory-optimization)
- [I/O Optimization](#io-optimization)
- [Advanced Examples](#advanced-examples)
- [Expected Output](#expected-output)
- [Common Pitfalls](#common-pitfalls)

---

```python
from metainformant.core.caching import memoize_disk

@memoize_disk(cache_dir=".cache/")
def expensive_computation(data):
    # Results cached to disk automatically
    return process(data)
```

## Parallel Processing

```python
from metainformant.core.parallel import parallel_map

# Process 1000 genomes in parallel (auto-chunks)
results = parallel_map(
    analyze_genome,
    genome_files,
    n_jobs=-1,  # Use all cores
    chunk_size=10
)
```

## Memory Optimization

- Enable streaming for large files (`stream=True` in read functions)
- Use binary formats (Parquet, HDF5) > CSV for tabular data
- Downcast NumPy arrays: `.astype(np.float32)` instead of float64
- Release references early: `del large_object; gc.collect()`

## I/O Optimization

### Batch reads/writes
```python
# Bad: one-by-one
for sample in samples:
    df = pd.read_csv(f"{sample}.tsv")  # Slow open/close per file

# Good: batch
dfs = [pd.read_csv(f) for f in all_files]  # Single op
# or use dask for out-of-core
import dask.dataframe as dd
ddf = dd.read_parquet("data/*.parquet")
```

### SSD vs HDD benchmarks
```bash
# Test disk speed
dd if=/dev/zero of=testfile bs=1G count=1 oflag=dsync
# SSD: ~500 MB/s, HDD: ~100 MB/s
# → If <200 MB/s, consider tmpfs (ramdisk) for scratch
```

### Compression tradeoffs
| Format | Compression | Speed | Use case |
|--------|-------------|-------|----------|
| gzip  | 6:1         | Slow  | Archived, infrequent access |
| bzip2 | 9:1         | Slower| Max compression, archival |
| xz    | 10:1        | Very slow | Long-term storage |
| zstd  | 4:1         | Fast  | Fast round-trip compression |
| none  | 1:1         | Fastest| Scratch space, SSD ample |

For intermediate files, use `zstd`; for final distribution, use `bgzf` (blocked gzip).

## Advanced Examples

### Numba JIT compilation for custom loops
```python
from numba import jit
import numpy as np

@jit(nopython=True, parallel=True, fastmath=True)
def compute_pairwise_distances(matrix):
    """Compute pairwise Euclidean distances in parallel."""
    n = matrix.shape[0]
    distances = np.empty((n, n), dtype=np.float32)
    for i in range(n):
        for j in range(i+1, n):
            dist = np.sqrt(((matrix[i] - matrix[j])**2).sum())
            distances[i, j] = dist
            distances[j, i] = dist
    return distances

# 100x speedup over pure Python
large_matrix = np.random.randn(5000, 100).astype(np.float32)
dist = compute_pairwise_distances(large_matrix)  # ~2s vs 200s
```

### Memory-mapped arrays for >RAM datasets
```python
import numpy as np

# Create memory-mapped file (doesn't load into RAM)
mmap = np.memmap("big_matrix.dat", dtype=np.float32, mode="w+", shape=(100000, 1000))

# Use like normal ndarray, but pages from disk
row_mean = mmap[0].mean()  # Only loads row 0

# Flush changes to disk
mmap.flush()
```

### Profiling bottleneck with cProfile + snakeviz
```bash
# Profile run
python -m cProfile -o profile.pstats scripts/rna/run_amalgkit_single.py --config config.yaml

# Visualize (install: pip install snakeviz)
snakeviz profile.pstats
# Opens browser: heatmap of time spent per function
# Find hot spots: e.g., alignment (45%), IO (30%), DESeq2 (15%)
```

### Parallel processing with joblib for CPU-bound tasks
```python
from joblib import Parallel, delayed
from metainformant.dna import alignment

# Align 1000 sequence pairs across all cores
seq_pairs = [(seq_a[i], seq_b[i]) for i in range(1000)]

results = Parallel(n_jobs=-1, backend="loky")(
    delayed(alignment.global_align)(a, b, match=2, mismatch=-1, gap=-2)
    for a, b in seq_pairs
)
# n_jobs=-1 uses all cores; backend="loky" avoids GIL
```

### Chunked processing with iterators
```python
from metainformant.core.io import read_csv_chunks

# Process 10M-row CSV on 8GB RAM laptop
for chunk in read_csv_chunks("huge.tsv", chunksize=100000):
    # chunk is DataFrame with 100k rows
    summary = chunk.groupby("sample").mean()
    write_results(summary)  # Flush per chunk
```

### Dask for distributed-like parallelism on single machine
```python
import dask.array as da

# Create lazy computation graph (1000×1000 matrix in chunks)
x = da.random.random((100000, 100000), chunks=(10000, 10000))

# Compute SVD lazily, using all cores and disk spill
u, s, vt = da.linalg.svd_compressed(x, k=10, compute=False)
result = u.compute()  # Evaluates graph
```

## Expected Output

### Cache hit rate report
```
Cache directory: .metainformant_cache/
  hits: 1247 (94.2%)
  misses: 77 (5.8%)
  disk usage: 2.3 GB
Top cached functions:
  1. alignment.global_align: 412 hits (avg 0.08s saved each)
  2. composition.kmer_frequencies: 298 hits (avg 0.12s saved each)
  3. population.tajima_d: 187 hits (avg 1.4s saved each)
Cache saved ≈ 9m 32s of compute time this session.
```

### Parallel speedup chart (console ASCII)
```
n_jobs | Time (s) | Speedup | Efficiency
-------|----------|---------|-----------
  1    |  124.5   |  1.00x  | 100%
  2    |   63.2   |  1.97x  | 98.5%
  4    |   32.1   |  3.88x  | 97.0%
  8    |   17.4   |  7.15x  | 89.4%
 16    |   10.2   |  12.2x  | 76.3%
Overhead: thread spawning 0.8s, chunking 0.3s
```

### Memory profile (before vs after)
```
Peak memory before optimization: 24.1 GB (OOM on 32GB machine)
Peak memory after streaming:   6.8 GB (comfortable)

Breakdown:
 - Sequence loading: 18.2 GB → 0.0 GB (streaming=True)
 - Alignment matrix: 4.1 GB → 4.1 GB (unchanged)
 - Intermediate objects: 1.8 GB → 0.9 GB (del + gc.collect)

Recommendation: Use streaming for all >1Gb FASTAs.
```

### I/O throughput comparison
```
Format      | Size   | Read Time | Write Time | Compression Ratio
------------|--------|-----------|------------|-----------------
CSV         | 8.2 GB | 42.1 s    | 38.7 s     | 1:1
Parquet     | 1.9 GB | 6.3 s     | 5.8 s      | 4.3:1
HDF5        | 2.1 GB | 8.2 s     | 7.9 s      | 3.9:1
zstd (lvl3) | 2.4 GB | 9.1 s     | 8.5 s      | 3.4:1
```

### Wall-clock time breakdown (RNA-seq 28 species)
```
Task                      | Local (48h) | Cloud (parallel) | Speedup
--------------------------|-------------|------------------|--------
Genome download           | 8h          | 0.5h (parallel)  | 16x
Index building (per sp)   | 12h         | 2h (80 cores)    | 6x
Alignment (8 samples)     | 24h         | 0.5h (96 cores)  | 48x
Quantification            | 4h          | 0.2h             | 20x
Total (wall)              | 48h         | ~3h              | 16x
Total (VM cost)           | $0          | $12              | N/A
```

## Common Pitfalls

| Problem | Likely Cause | Fix |
|---------|-------------|-----|
| `n_jobs=-1` slows down (not speeds up) | Too many processes, GIL contention, memory thrashing | Set `n_jobs=os.cpu_count()//2`; ensure data fits per-worker; use `backend="loky"` or `threading` for I/O-bound |
| Cache misses despite memoization | Different argument values (e.g., paths differ) or unhashable args | Normalize paths: `os.path.abspath()`; convert lists to tuples; use `functools.lru_cache` for in-memory (small objects) |
| `MemoryError` despite `stream=True` | Downstream code eagerly loads entire file | Check all function calls in pipeline; avoid `.read()` or `.load()` on results |
| Parallel speedup plateaus at 4x | Bottleneck elsewhere (disk I/O, not CPU) | Profile with `cProfile`; move data to RAM disk (/dev/shm); compress on-the-fly |
| Subprocesses hang after `parallel_map` | Child processes inherit locks, deadlock | Use `if __name__ == "__main__":` guard; avoid global state; prefer `multiprocess` (dill serialization) |
| Slow pandas on 1B rows | Pandas not designed for out-of-core | Use `dask`, `polars`, or `modin.pandas`; or chunk with `read_csv(chunksize=)` |
| Dask cluster overhead > gain | Task graph too fine-grained (10k tasks < 1ms each) | Coarsen tasks: increase chunk size; use `dask.delayed` for irregular tasks; persist intermediate results |

---

**Related:** [Core parallel processing](../core/parallel.md) | [Caching strategies](../core/cache.md) | [Profiling guide](../testing.md) | [I/O reference](../core/io.md)
