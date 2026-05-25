# Performance: structural_variants

Performance characteristics, optimization strategies, and scaling guidelines for structural variant analysis in METAINFORMANT.

## Benchmark Summary

**Test environment (unless noted):**

| Component | Specification |
|-----------|---------------|
| CPU | Intel Xeon Gold 6248R (2× 24 cores, 48 threads) @ 3.0 GHz |
| RAM | 256 GB DDR4 ECC |
| Storage | NVMe SSD (2 GB/s sequential) |
| OS | Ubuntu 22.04 LTS |
| Python | 3.11.8 |
| NumPy | 1.26.4 |
| pysam | 0.22.0 |

**Workload:** Human whole-genome sequencing (WGS), 30× mean depth, BAM size ~100 GB (≈700M reads).

| Stage | Input Size | Wall Time (1 thread) | Wall Time (8 threads) | Peak RAM | Notes |
|-------|------------|---------------------|-----------------------|----------|-------|
| **CNV detection** | Whole-genome depth | 8–12 min | 3–5 min | 4–6 GB | CBS segmentation is O(n) on number of bins (~250K for 1000 bp windows) |
| **SV calling** | 700M alignment records | 45–60 min | 12–18 min | 12–18 GB | Split + discordant evidence detection; single-pass streaming |
| **Breakpoint refinement** | 5K SV candidates | 3–5 min | 1–2 min | 3–4 GB | Refetches reads locally; memory scales with SV count |
| **Gene overlap annotation** | Up to 20K SVs vs 200K genes | 5–8 s | 2–4 s | 500 MB | `IntervalIndex` binary search; fast |
| **Functional impact** | 5K SVs | 1–2 min | 1 min | 1–2 GB | TAD lookup (linear search in TAD boundaries) |
| **Quality filtering** | 20K SVs | <1 s | <1 s | <100 MB | Simple thresholds |
| **Multi-caller merge** | 3 callers × 15K SVs = 45K | 30–90 s | 20–60 s | 2–4 GB | O(n²) worst-case but early termination on sorted chroms |
| **Whole pipeline (single sample)** | 100 GB BAM | 60–80 min | 16–25 min | 20–30 GB | End-to-end detection → filtering |

**Cohort scaling (N samples):**

| N samples | Single-threaded wall time | 8-thread (per-sample parallel) | Memory / sample |
|-----------|-------------------------|-------------------------------|-----------------|
| 10 | ~13 h | ~2 h (plus overhead) | ~25 GB max |
| 100 | ~5.5 days | ~20 h | ~25 GB (one at a time) |
| 1000 | ~55 days | ~8 days (batch-8) | — |

**Notes:**
- Per-sample steps are **embarrassingly parallel**. Use `parallel.map()` across samples.
- Memory is dominated by BAM loading (~100 GB for 700M reads in Python dicts). Consider streaming reads in chunks or downsampling for quick checks.
- Disk space: 100 GB BAM → SV JSON ~50–200 MB per sample (uncompressed); gzip saves ~60%.

---

## Optimization Checklist

### 1. Parallel Processing (Biggest Win)

```python
from metainformant import config

config.set("structural_variants.parallel", True)
config.set("structural_variants.max_workers", 8)  # ≤ physical core count
```

**Speedup:** 2.5–4× on 8 cores for SV calling (I/O-bound phases saturate earlier than CBS). CNV detection gains less (<2×) due to NumPy GIL release.

**Per-sample parallelism:**

```python
from metainformant import parallel

def process_sample(bam_path):
    # Full pipeline for one sample
    return results

results = parallel.map(
    process_sample,
    list_of_bams,
    workers=12,          # Max CPU cores
    chunk_size=10,       # Batch N tasks per worker
)
```

**Caveat:** Each worker loads its own BAM → RAM = N_workers × BAM_size. Limit workers to fit memory or use `workers=1` + `chunk_size` streaming.

### 2. Chunk Size Tuning (SV Calling)

SV calling processes reads in memory. Increase `chunk_size` to reduce Python loop overhead:

```python
config.set("structural_variants.chunk_size", 20000)  # default 10000
```

**Trade-off:** Larger chunks = faster but more RAM per worker. Benchmark 1% sample to find sweet spot.

### 3. Threshold Tuning for Speed vs. Sensitivity

| Setting | Fast | Balanced | Accurate |
|---------|------|----------|----------|
| `SV_MIN_MAPQ` | 30 | 20 | 10 |
| `SV_MIN_SUPPORT` | 5 | 3 | 2 |
| `SV_CNV_SIGNIFICANCE` | 0.05 | 0.01 | 0.001 |

Lower thresholds increase runtime (more candidates) but improve recall.

### 4. Compressed I/O

```python
from metainformant.core import io

# Result files
io.dump_json(results, "results.json.gz")   # gzip compression (default if .gz extension)
# or
io.dump_json(results, "results.json", compress="gzip")
```

**Impact:** ~70% size reduction; adds ~10–15% CPU overhead but cuts disk I/O (esp. network FS).

### 5. Caching Intermediate Results

```python
from metainformant.core.io import cache

cache.enable(ttl=3600)  # Cache for 1 hour
cache.disk_path = "/ssd/cache"  # Fast SSD

# Subsequent runs with same BAM + config skip CNV recompute
```

**Best use:** Intermediate depth arrays (expensive to recompute from BAM). Not for final SV lists (already small).

### 6. I/O Optimization

**Disk:** Use local SSD, not network-mounted home directory. If using NFS, set `TMPDIR=/local/ssd` for temp files.

**Streaming large BAMs:** Use `pysam.fetch()` with coordinate intervals rather than `fetch(until_eof=True)`. Process chromosome-by-chromosome:

```python
for chrom in bam.references:
    chrom_reads = [r for r in bam.fetch(chrom)]  # Stream one chrom at a time
    process(chrom_reads)
    # Clear memory
    del chrom_reads
    import gc; gc.collect()
```

**Memory-mapped arrays:** If depth arrays become big (>1 GB), use `np.memmap`:

```python
depth = np.memmap("depth.dat", dtype=np.float64, mode="w+", shape=(n_bins,))
# Write in chunks, then flush to disk
```

### 7. Algorithm Selection

- **CNV CBS**: Only algorithm currently; O(n). Accepts GC correction (adds O(n) passes).
- **SV genotyping (cohort)**: Use `"depth"` method for speed; `"split"` for accuracy (requires read-level evidence).
- **Multi-caller merging**: `merge_callsets()` (reciprocal overlap) more accurate; `survivor_merge()` faster for many variants.

### 8. Memory Reduction Strategies

- **Release large objects:** `del depth_data; import gc; gc.collect()`
- **Use float32 instead of float64:** Set `dtype=np.float32` in internal arrays if precision acceptable (controlled by config flag `SV_FLOAT_PRECISION`).
- **Process chromosomes sequentially:** Do not accumulate `all_variants` across chroms; write per-chrom intermediate files then combine.
- **Filter early:** Apply size/quality filters during SV calling to reduce candidate list before refinement/annotation.

---

## Profiling & Bottleneck Diagnosis

### Profile with `cProfile`

```python
import cProfile
import pstats
from metainformant.structural_variants import pipeline

profiler = cProfile.Profile()
profiler.enable()

results = pipeline.run("sample.bam", ...)

profiler.disable()
stats = pstats.Stats(profiler)
stats.sort_stats("cumtime")  # or "tottime"
stats.print_stats(20)  # Top 20 functions
```

**Typical hotspots:**

| Function | % Time | Optimization |
|----------|--------|--------------|
| `call_structural_variants()` | 40–60% | Reduce `min_support`; increase `min_mapq` |
| `detect_split_reads()` | 20–30% | CIGAR parsing is CPU-bound; consider Cython (future) |
| `detect_discordant_pairs()` | 10–15% | Fewer paired reads → faster |
| `segment_coverage()` (CBS) | 5–10% | Already O(n); not a bottleneck |
| Annotation (`IntervalIndex`) | <5% | Already optimized |

### Enable Debug Timing

```bash
export CORE_LOG_LEVEL=DEBUG
```

Logs include:
```
[DEBUG] Loading BAM: 45.2 s
[DEBUG] Computing depth: 124.6 s
[DEBUG] CBS segmentation: 8.3 s
[DEBUG] SV calling: 2145.1 s  ← biggest cost
[DEBUG] Annotation: 14.2 s
[DEBUG] Filtering: 1.2 s
[DEBUG] Total: 2329.6 s
```

### Memory Profiling

```bash
pip install memory-profiler
```

```python
from memory_profiler import profile

@profile
def run_pipeline():
    ...

run_pipeline()
```

**Expected peak:** ~25 GB for one 30× WGS sample full pipeline. Exceeds 32 GB limit → swap → crash. Solution: reduce workers to 1–2, or sample down depth for testing.

---

## Scaling Beyond a Single Machine

### Distributed Computing (Cluster)

Use `metainformant.cloud` or custom Dask:

```python
from dask.distributed import Client
from metainformant import parallel

client = Client(n_workers=20)  # SLURM cluster

# Each sample is a task
futures = []
for bam in bam_list:
    future = client.submit(process_sample, bam)
    futures.append(future)

results = client.gather(futures)
```

### Cloud (GCP/AWS)

See `docs/cloud/` for preemptible VM instructions.

- **Storage:** Use regional SSD persistent disks (PD-SSD) for BAMs
- **Compute:** Use `n1-standard-16` or `n2-highmem-32` depending on dataset size
- **Parallelize:** One sample per VM; merge results in Cloud Storage

---

## Large Cohort Best Practices (N ≥ 1000)

1. **Pre-create reference panels** (not per-sample):
   - Compute depth once per sample, store as binary `.npz` (compressed) → reuse for CNV, genotyping.
   - Example: `depth_{sample}.npz`

2. **Downsample for quick SV discovery:**
   - Extract subset of reads (`samtools view -s 0.1`) to calibrate parameters quickly.

3. **Run stages separately:**
   - **Stage 1:** Detection only (CNV + SV call) → store intermediates.
   - **Stage 2:** Annotation (deterministic, can parallelize across samples).
   - **Stage 3:** Multi-caller merge (if needed) → requires all callers for all samples; batch-process.
   - **Stage 4:** Population analysis (single process; memory-bound by genotype matrix).

4. **Genotype matrix sparsity:**
   - Genotype matrix (100K SVs × 10K samples) is sparse (<5% alt allele). Store as `scipy.sparse.csr_matrix` or `uint8` dense if dense enough.
   - PCA can use randomized SVD on sparse matrices.

5. **Intermediate storage format:**
   - Use Parquet (columnar) for per-sample variant tables: `pd.DataFrame.to_parquet()`.
   - More efficient than JSON for large N.

6. **Sample sharding for parallel map:**
   - Divide samples into shards of ~100, process each shard, then merge.

---

## Algorithmic Complexity Summary

| Component | Time Complexity | Space Complexity | Scaling |
|-----------|----------------|------------------|---------|
| CBS (`segment_coverage`) | O(n) per chromosome | O(n) | Linear in bin count |
| SV calling | O(m log m) from sorting | O(m) for alignments | Linear in read count |
| Evidence clustering | O(e log e) from sort | O(e) | e = evidence items; ≤ m |
| Gene overlap (`IntervalIndex`) | O(log g + k) per query | O(g) | g = # genes |
| Multi-caller merge | O(v²) worst-case; ~O(v log v) with early break | O(v) | v = variants per chrom |
| SV population PCA | O(n × m²) full SVD; O(n × m × k) randomized | O(n × m) | n=samples, m=SVs |

**Linear or near-linear in input data size. Biggest cost is SV calling (read-level inference); consider downsampling or distributed processing for large cohorts.**

---

## Hardware Recommendations

| Scenario | RAM | CPU Cores | Storage |
|----------|-----|-----------|---------|
| Development / testing (few samples) | 32 GB | 4–8 | Local SSD |
| Production single-sample (30× WGS) | 64–128 GB | 16–32 | NVMe SSD |
| Cohort processing (100 samples) | 256 GB | 32–64 | RAID 0 NVMe or parallel network FS |
| Population-scale (1000+ samples) | 512 GB – 1 TB | 64+ | Lustre/GPFS parallel FS |

**Note:** Peak memory = loading full BAM + intermediate arrays. Estimate ~2× BAM file size for full pipeline.

---

## Common Pitfalls & Fixes

| Symptom | Likely Cause | Fix |
|---------|--------------|-----|
| Runtime too slow (single sample > 4 h) | Low `min_mapq` → too many low-quality reads | Raise `SV_MIN_MAPQ` to 30 |
| MemoryError: cannot allocate array | Workers × BAM memory | Reduce `max_workers` to 1–2, or stream by chromosome |
| Disk fills with intermediate JSON | Compressed output not enabled | Use `.json.gz` extension or set `CORE_COMPRESS=true` |
| Very few SVs detected | Over-filtering (thresholds too strict) | Lower `min_support`, `min_qual`; inspect intermediate counts |
| Too many SVs (noisy) | Under-filtering or low coverage | Raise thresholds; check data quality |
| Merge step hangs | Graph building O(v²) too large | Filter more before merge; reduce variant count per chromosome with `pre_merge_filter` |
| PCA crashes | Genotype matrix too large for memory | Use randomized SVD (`svd_solver='randomized'`), or subset SVs |

---

## Configuration for Speed

**For quick exploration (speed priority):**

```bash
export SV_MIN_MAPQ=30
export SV_MIN_SUPPORT=5
export SV_CNV_SIGNIFICANCE=0.05
export SV_FILTER_MIN_QUAL=30
```

**For publication-quality results (accuracy priority):**

```bash
export SV_MIN_MAPQ=10
export SV_MIN_SUPPORT=2
export SV_CNV_SIGNIFICANCE=0.001
export SV_FILTER_MIN_QUAL=10
```

---

## Monitoring Progress

Use structured logging to a file:

```python
import logging
from metainformant.core.utils.logging import setup_logging

setup_logging(
    level="INFO",
    log_file="pipeline.log",
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
```

Log file shows per-stage timings on every `logger.info()` call containing timing.

For real-time monitoring, `tail -f pipeline.log | grep "TIMING"`.

---

## Future Optimizations (Roadmap)

- **Cython rewrite** of CIGAR parsing + evidence clustering (2–3× speedup)
- **GPU-accelerated CBS** for CNV detection (experimental)
- **K-d tree indexing** for gene overlap (faster for massive annotation databases)
- **Batch PCA** for incremental cohort addition
- **Compressed sparse genotype matrix** as default for N > 100

---

## Related

- **[GETTING_STARTED.md](GETTING_STARTED.md)** — Basic pipeline setup
- **[CONFIGURATION.md](CONFIGURATION.md)** — All performance-related settings
- **[TROUBLESHOOTING.md](TROUBLESHOOTING.md)** — Diagnose slow runs, OOM
