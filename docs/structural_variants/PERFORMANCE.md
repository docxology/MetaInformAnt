# Performance: structural_variants

Performance characteristics, optimization strategies, and scaling guidelines for structural variant analysis in METAINFORMANT.

## Benchmark Summary

This module ships no benchmark suite, and the repository holds no whole-genome timing, memory or cohort-scaling measurements. The functions work on in-memory Python lists and NumPy arrays, which sets what drives the cost:

| Stage | What drives the cost |
|-------|----------------------|
| **CNV segmentation** (`segment_coverage`, called by `detect_cnv_from_depth`) | Every arc of at most half the segment is scored, vectorised over the arcs of one length with NumPy and repeated at each recursion level: still grows with the square of the number of bins per chromosome |
| **SV calling** (`call_structural_variants`) | Linear in the reads: one pass for split/discordant evidence, one pass and a sort to build the `ReadIndex`, then genotyping asks the index for the few reads near each call, so the cost per called SV is small and does not depend on the number of reads |
| **Gene overlap annotation** (`annotate_gene_overlap`) | Converts a gene list of dicts to `GenomicInterval` objects and builds an `IntervalIndex` on every call (a sort plus a linear pass); each query then costs time proportional to the number of genes it returns, wherever on the chromosome it falls |
| **Multi-caller merge** (`merge_callsets`, `survivor_merge`) | Sort, then a pairwise sweep within each chromosome that only compares variants whose starts fall inside the current variant's span plus 10 kb (`merge_callsets`) or `max_distance` (`survivor_merge`) |
| **Population functions** | NumPy on dense matrices; PCA is an eigendecomposition of the smaller of the sample × sample and SV × SV covariance matrices |

**Indicative timings** (synthetic data, single thread, Python 3.12.15 on an arm64 Mac; CPU time, best of three runs unless noted, so read the scaling rather than the absolute values):

| Operation | Input | Time |
|-----------|-------|------|
| `segment_coverage` | 1,000 / 2,000 / 4,000 / 8,000 / 16,000 bins; from a flat noise profile to sixteen segments | 0.003–0.011 s / 0.006–0.020 s / 0.02–0.06 s / 0.07–0.23 s / 0.27–0.64 s |
| `segment_coverage` | 32,000 / 64,000 / 128,000 / 256,000 bins; three-segment profile, one run each | 1.2 s / 4.5 s / 19 s / 68 s |
| `detect_cnv_from_depth` | one chromosome, 2,490 bins with no copy-number change (chr1 at 100 kb windows) | 0.01 s |
| `call_structural_variants` | 50k reads, 50 SVs | 0.14 s |
| | 100k reads, 50 SVs | 0.30 s |
| | 100k reads, 100 SVs | 0.33 s |
| | 200k reads, 100 SVs | 0.62 s |
| `IntervalIndex` build | 20k / 80k genes on one chromosome | 16 ms / 154 ms |
| `IntervalIndex.query_overlap` | a 10 kb window near the end of the chromosome, 20k / 80k genes | 2 µs / 7 µs |
| `annotate_gene_overlap` | 60,000 genes on 20 chromosomes; 100 / 2,000 / 20,000 variants | 0.10 s / 0.12 s / 0.30 s |

The synthetic inputs were: CBS profiles of Gaussian noise (standard deviation 0.15) with a +0.8 step over the middle third (three segments) or +0.7 steps over every second sixteenth (sixteen segments); reads of 150 bp at random positions on a 50 Mb chromosome, properly paired with a 300 bp insert (`insert_size_stats` passed), plus six split reads with an SA tag for each planted deletion; genes of 1–60 kb at random positions on chromosomes of 100–150 Mb (one chromosome for the index timings, twenty for `annotate_gene_overlap`). Each doubling of the bins multiplies the CBS time by about 3 to 4 (the first row gives the range over the three profile shapes: a flat noise profile at the low end, sixteen segments at the high end). A human chromosome at 1 kb windows (about 50k–250k bins) therefore takes seconds to a minute or more per chromosome: pick `window_size` so that each chromosome has at most a few tens of thousands of bins. `call_structural_variants` grows linearly with the reads (about 3 µs of CPU per read here) and the number of called SVs adds little. Memory is dominated by the alignment records you hold: an alignment dict with 11 keys and no sequence took about 0.6 KB when its string values were shared objects and about 0.7 KB when each record owned its strings (measured with `tracemalloc`), so 100M reads are about 60-70 GB before any sequence strings.

**Notes:**
- Per-sample steps are **embarrassingly parallel**. Use `parallel.process_map()` across samples; the cost of a cohort is then the per-sample cost times N divided by the number of workers.
- Work per chromosome or region rather than loading a whole genome: the detection and annotation functions have no streaming mode.
- JSON results shrink with gzip when written to a `.json.gz` path (see Compressed I/O below).

---

## Optimization Checklist

### 1. Parallel Processing (Biggest Win)

The functions themselves are single-threaded and have no thread or worker settings (there is no `config.set("structural_variants.parallel", ...)`); the practical form of parallelism is one worker per sample.

**Per-sample parallelism:**

```python-snippet
from metainformant.core.execution import parallel

def process_sample(bam_path):
    # Full pipeline for one sample (a top-level function, so the process pool can pickle it)
    return results

results = parallel.process_map(
    process_sample,
    list_of_bams,
    max_workers=12,      # worker processes; None = cpu_count() - 1
)
```

`process_map` runs a process pool, submits one task per item and has no chunk-size option. On macOS and Windows, where workers are spawned, call it under `if __name__ == "__main__":`.

**Caveat:** each worker builds its own in-memory alignment list, so peak RAM is roughly the number of workers times the memory of one sample's records. Limit `max_workers` to fit memory, or process each sample by chromosome (next section).

### 2. Processing in Slices (SV Calling)

`call_structural_variants` works on the list of alignment dicts you pass it and has no chunking parameter. To bound memory, build the list and call it per chromosome or region (for example one `bam.fetch(chrom)` at a time), then concatenate the results. Two things change with slicing: the insert-size distribution is estimated from each slice unless you pass the same `insert_size_stats` to every call, and an inter-chromosomal event is reported once per slice (chr1→chr8 and chr8→chr1), so drop or merge the mirrored calls afterwards.

**Trade-off:** smaller slices cost less memory (the run time is linear in the reads either way, because genotyping uses a `ReadIndex`), at the price of the two caveats above. Try a 1% region first to see how time and memory behave.

### 3. Threshold Tuning for Speed vs. Sensitivity

| Setting | Quick exploration | Balanced | Maximum sensitivity |
|---------|-------------------|----------|---------------------|
| `min_mapq` (argument of `call_structural_variants`) | 30 | 20 | 10 |
| `min_support` argument or `SV_MIN_SUPPORT` | 5 | 3 | 2 |
| `significance` argument or `SV_CNV_SIGNIFICANCE` | 0.001 | 0.01 | 0.05 |

Lower `min_mapq` and `min_support` admit more reads and more candidate calls (more work, more false positives, better recall in repeats and at low coverage). `significance` works the other way round: it is the p-value threshold for accepting a CBS change point, so a larger value accepts more change points and a smaller value gives fewer, more confident segments. It hardly changes the run time (CPU time 0.012 s at 0.001, 0.011 s at 0.01 and 0.012 s at 0.05 on a 2,000-bin profile with three segments), because the scans at the top levels of the recursion dominate the cost. There is no environment variable for `min_mapq`.

### 4. Compressed I/O

```python
from metainformant.core import io

# Result files: a ".gz" extension makes dump_json write gzip-compressed JSON
io.dump_json(results, "results.json.gz")
results = io.load_json("results.json.gz")   # load_json reads .gz transparently
```

`dump_json` has no `compress=` argument; the extension alone selects gzip. Compression trades CPU time for smaller files, which helps most on network filesystems. The data must already be plain JSON types (convert dataclasses with `dataclasses.asdict`).

### 5. Caching Intermediate Results

Nothing in the module is cached automatically (`core.io.cache` has no `enable()` or `disk_path`). To avoid recomputing an expensive intermediate such as the per-bin depth arrays, store it yourself with `JsonCache`:

```python
import numpy as np
from metainformant.core.io.cache import JsonCache

cache = JsonCache("/ssd/cache", ttl_seconds=3600)   # cache directory on a fast SSD, 1 hour default lifetime

key = f"depth:{bam_path}:1000"                       # include everything the result depends on
cached = cache.get(key)
if cached is None:
    depth_data = compute_depth(bam_path, window_size=1000)   # your BAM-reading helper
    cache.set(key, {chrom: depth.tolist() for chrom, depth in depth_data.items()})
else:
    depth_data = {chrom: np.asarray(values) for chrom, values in cached.items()}
```

**Best use:** Intermediate depth arrays (expensive to recompute from BAM). Not for final SV lists (already small). `JsonCache` stores JSON, which is bulky for large arrays; for those, `np.save` / `np.load` on a path derived from the key is the better fit. `JsonCache` names the file after the letters, digits and `._-` of the key only, so keys that differ in other characters alone (for example `depth:/data/a/b.bam:1000` and `depth:/data/ab.bam:1000`) share one entry; make the key unique in those characters, for example by including a hash of the full path.

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

- **CNV CBS**: `method="segmentation"` is the only method, and its cost grows with the square of the bins per chromosome (see above), so control cost with `window_size`. The input is read depth, not log2 ratios: negative or non-finite values raise `ValueError`. GC correction is a separate step (`calculate_log2_ratio(..., gc_content=...)`) and is linear in the number of bins.
- **SV genotyping (one sample)**: `call_structural_variants` builds one `ReadIndex` and passes it to `genotype_sv`. When you call `genotype_sv` yourself for many variants against the same reads, build the index once (`ReadIndex(reads)`) and pass it in place of the list; the genotype is the same and each call stops being a pass over every read.
- **SV genotyping (cohort)**: `genotype_sv_population` turns precomputed per-sample evidence into genotypes with a constant-time rule per cell, so `"depth"` (uses `depth_ratio`) and `"split"` (uses `split_reads`, `discordant_pairs`, `total_reads`) cost the same; choose by which evidence you have. The expensive part is computing that evidence from the BAMs.
- **Multi-caller merging**: `merge_callsets()` matches by reciprocal overlap and `survivor_merge()` by breakpoint distance; both sort and sweep each chromosome in the same way, so choose between them by matching criterion, not by speed.

### 8. Memory Reduction Strategies

- **Release large objects:** `del depth_data; import gc; gc.collect()`
- **Keep the genotype matrix compact:** store it as `np.uint8` (values 0/1/2). The matrix-based population functions convert it internally (`sv_population_structure` to `float64`, making several copies) and do not accept `scipy.sparse` matrices.
- **Process chromosomes sequentially:** Do not accumulate `all_variants` across chroms; write per-chrom intermediate files then combine.
- **Keep records small:** drop keys the functions do not read (for example the read sequence `seq`, which `call_structural_variants` ignores) before holding millions of alignment dicts.
- **Filter early:** `call_structural_variants` already applies `min_mapq` before clustering and `min_support` / `min_sv_size` afterwards; run `filter_by_quality` / `filter_by_size` / `apply_blacklist` before annotation so the slower per-query annotation steps see fewer variants.

---

## Profiling & Bottleneck Diagnosis

### Profile with `cProfile`

```python-snippet
import cProfile
import pstats
from metainformant.structural_variants.detection.sv_calling import call_structural_variants

profiler = cProfile.Profile()
profiler.enable()

variants = call_structural_variants(alignments, min_support=3)

profiler.disable()
stats = pstats.Stats(profiler)
stats.sort_stats("cumtime")  # or "tottime"
stats.print_stats(20)  # Top 20 functions
```

**Where the time goes in `call_structural_variants`** (example profile on synthetic data: 100k reads, 100 SVs; shares of the call's cumulative time under `cProfile`, as the range over three profiled runs):

| Function | Share | Notes |
|----------|-------|-------|
| `detect_split_reads()` | ≈ 23–34% | CIGAR parsing in pure Python |
| `ReadIndex(...)` construction | ≈ 26–28% | One pass over the reads plus a sort by start position; built once per call |
| `detect_discordant_pairs()` | ≈ 22–26% | One pass over reads with mate information |
| `_estimate_insert_size()` | ≈ 12–17% | Skipped when you pass `insert_size_stats` |
| `genotype_sv()` | < 1% | Called once per SV; a binary search plus the few reads near the breakpoint |

For CNV detection the cost is almost entirely `segment_coverage()` at fine bin sizes. For annotation the split depends on the sizes: with a large gene list and few variants most of the call goes into converting the gene dicts and building the index (60,000 genes on 20 chromosomes: 0.10 s for 100 variants, 0.12 s for 2,000), and the queries take a growing share as the number of variants grows (0.30 s for 20,000). Profile your own data before optimizing; the split changes with the number of SVs, reads and bins.

### Timing the Stages

The modules log stage summaries but no timings. With logging configured (`configure_logging_from_env()`, see [GETTING_STARTED.md](GETTING_STARTED.md)), INFO shows lines such as

```
CNV detection complete: 3 segments across 1 chromosomes
Called 3 structural variants from 28 evidence items
Annotated 3 variants: 2 overlap with genes
```

and `CORE_LOG_LEVEL=DEBUG` adds details such as `Detected N split-read evidence items from M reads`. To time stages, wrap them yourself:

```python
import time
from metainformant.structural_variants.detection.sv_calling import call_structural_variants

t0 = time.perf_counter()
variants = call_structural_variants(alignments)
print(f"SV calling: {time.perf_counter() - t0:.1f} s")
```

### Memory Profiling

```bash
uv pip install memory-profiler
```

```python-snippet
from memory_profiler import profile

@profile
def run_pipeline():
    ...

run_pipeline()
```

**Rule of thumb:** peak memory is dominated by the alignment records you hold, about 0.6-0.7 KB per dict without a sequence (plus the sequence length if you keep `seq`). If a run swaps or hits `MemoryError`, process by chromosome, drop unused keys, or reduce the number of `process_map` workers.

---

## Scaling Beyond a Single Machine

### Distributed Computing (Cluster)

The structural variant functions are plain Python functions over in-memory data, so any scheduler that can run one function per sample works, for example Dask:

```python-snippet
from dask.distributed import Client

client = Client(n_workers=20)  # a local cluster; use dask_jobqueue.SLURMCluster for SLURM

# Each sample is a task
futures = []
for bam in bam_list:
    future = client.submit(process_sample, bam)
    futures.append(future)

results = client.gather(futures)
```

`metainformant.cloud` (`CloudConfig`, `GCPDeployer`) is GCP deployment infrastructure for pipelines; it has no structural-variant job runner.

### Cloud (GCP)

See `docs/cloud/` for the GCP deployment and preemptible VM notes.

- **Parallelize:** One sample per VM; collect the per-sample result files in object storage and merge them afterwards
- **Size VMs by memory:** use the per-record estimate above and measure one region before choosing a machine type

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

4. **Genotype matrix size:**
   - Store the genotype matrix (e.g. 100K SVs × 10K samples) as dense `np.uint8`. `sv_allele_frequency`, `sv_population_structure` and the other population functions call `np.asarray` on it, so `scipy.sparse` matrices are rejected, and there is no randomized-SVD option.
   - `sv_population_structure` converts to `float64` (8 bytes per genotype, so about 8 GB per copy at this size) and makes several copies; subset the SVs (common, LD-pruned) before PCA.

5. **Intermediate storage format:**
   - Use Parquet (columnar) for per-sample variant tables: `pd.DataFrame.to_parquet()`.
   - More efficient than JSON for large N.

6. **Sample sharding for parallel map:**
   - Divide samples into shards of ~100, process each shard, then merge.

---

## Algorithmic Complexity Summary

| Component | Time Complexity | Space Complexity | Scaling |
|-----------|----------------|------------------|---------|
| CBS (`segment_coverage`) | O(n²) per chromosome (every arc of at most half the segment is scored, vectorised over the arcs of one length, repeated per recursion level) | O(n) | n = bins per chromosome; each doubling cost 3–4× in the measurements above |
| Evidence detection (`detect_split_reads`, `detect_discordant_pairs`) | O(m) | O(e) | m = reads, e = evidence items |
| Evidence clustering | O(e log e) sort, plus a comparison of each item against the existing clusters, newest first, until one matches (worst case quadratic in the number of clusters) | O(e) | e = evidence items, roughly one per split or discordant read |
| SV genotyping (`genotype_sv` with a `ReadIndex`, as `call_structural_variants` calls it) | O(m log m) to build the index once, then O(log m + r) per called SV | O(m) for the index | r = reads that start within one longest-read length before the genotyping window; under 1% of the example profile. With a plain read list each call costs O(m) |
| Gene overlap (`IntervalIndex`) | O(g log g) build, O(k log g) per query, where k is the number of intervals returned; every `annotate_*` / `find_nearest_gene` call additionally converts every gene dict and builds the index | O(G) | g = genes on the chromosome, G = genes in the list |
| Multi-caller merge | O(v log v) sort plus pairwise comparisons inside each variant's span (+10 kb or `max_distance`); O(v²) worst case | O(v) | v = variants |
| SV population PCA | eigendecomposition of an s × s (s < v) or v × v covariance matrix: O(min(s, v)² · max(s, v) + min(s, v)³) | O(s × v) | s = samples, v = SVs; dense `float64`, no randomized solver |

**Time grows faster than linearly in CBS (bins per chromosome), in the worst case of evidence clustering and in multi-caller merging (variants within each other's span); everything else on the SV-calling and annotation path is linear in its input or logarithmic per query. Keep bins per chromosome in the tens of thousands at most, call per chromosome or region to bound memory, and filter before merging.**

---

## Hardware Recommendations

There are no tested hardware tiers for this module, so size the machine from the data you actually hold in memory:

- **RAM per worker** ≈ (alignment records held at once) × ≈0.6-0.7 KB, plus the sequence length per record if you keep `seq`, plus the depth arrays (8 bytes per bin) and the genotype matrix (1 byte per cell as `uint8`, 8 bytes in `float64` during PCA).
- **CPU:** the functions are single-threaded; extra cores only help through per-sample parallelism (`process_map`), so cores beyond (RAM ÷ RAM per worker) are wasted.
- **Storage:** fast local storage for the BAMs and the cache/intermediate files (`TMPDIR` on a local SSD) matters more than capacity.

---

## Common Pitfalls & Fixes

| Symptom | Likely Cause | Fix |
|---------|--------------|-----|
| SV calling is slow or runs out of memory | One list holding a whole genome of alignment dicts: the time is linear in the reads (about 3 µs each in the timings above) and the memory is 0.6–0.7 KB per read | Call per chromosome or region; raise `min_mapq` (argument) to drop reads before clustering and `min_support` to reduce the number of calls |
| CNV detection takes minutes per chromosome | Fine `window_size`: CBS grows with the square of the bins per chromosome (about a minute at 256,000 bins) | Use a coarser `window_size` so a chromosome has at most a few tens of thousands of bins |
| `ValueError: depth_data['chr1'] has ... negative ... depth value(s)` | log2 ratios (or other negative values) were passed to `detect_cnv_from_depth`, which takes read depths | Pass depths, or `2**ratio` for a tumor/normal ratio; segment log2 ratios with `segment_coverage` and `call_cnv_states` |
| MemoryError: cannot allocate array | Alignment dicts × `process_map` workers | Reduce `max_workers`, process by chromosome, drop unused keys such as `seq` |
| Disk fills with intermediate JSON | Uncompressed output | Write to a `.json.gz` path (`io.dump_json` gzips on that extension) |
| Very few SVs detected | Over-filtering (thresholds too strict), or `SV_MIN_SUPPORT` / `SV_MIN_QUAL` set in the environment (they override the arguments) | Lower `min_support`, `min_qual`; check `env \| grep SV_`; inspect intermediate counts |
| Too many SVs (noisy) | Under-filtering or low coverage | Raise thresholds; check data quality |
| Merge step slow | Large variants: each is compared with every variant that starts inside its span (+10 kb or `max_distance`) | Filter (size, quality, blacklist) before merging; merge per chromosome |
| PCA runs out of memory | Several `float64` copies of the genotype matrix | Subset the SVs (common, LD-pruned) before `sv_population_structure`; there is no randomized or sparse option |

---

## Configuration for Speed

Only `SV_MIN_SUPPORT`, `SV_MIN_QUAL` (read by `filter_by_quality`), `SV_CNV_SIGNIFICANCE`, `SV_CNV_DEL_THRESHOLD` and `SV_CNV_DUP_THRESHOLD` exist; each overrides the matching function argument. `min_mapq` has no environment variable (pass it to `call_structural_variants`). Of the settings below, `SV_MIN_SUPPORT` and `min_mapq` change how much work SV calling does; `SV_CNV_SIGNIFICANCE` and `SV_MIN_QUAL` change which segments and calls are kept, not the run time.

**For quick exploration (speed priority):**

```bash
export SV_MIN_SUPPORT=5
export SV_CNV_SIGNIFICANCE=0.001
export SV_MIN_QUAL=30
```
and `call_structural_variants(..., min_mapq=30)`.

**For maximum sensitivity (slower, more false positives):**

```bash
export SV_MIN_SUPPORT=2
export SV_CNV_SIGNIFICANCE=0.05
export SV_MIN_QUAL=10
```
and `call_structural_variants(..., min_mapq=10)`. Lower thresholds raise recall; they do not by themselves make results publication-quality, so inspect the extra calls.

---

## Monitoring Progress

Send the module's log records to a file:

```python-snippet
from metainformant.core.utils.logging import setup_logger

setup_logger(
    "metainformant.structural_variants",
    log_file="pipeline.log",
    level="INFO",
)
```

The file receives each stage summary once (for example `CNV detection complete: 3 segments across 1 chromosomes`); the console shows every line twice, because each module logger has its own console handler and the package logger adds another. The lines carry timestamps but no stage timings, so measure those yourself (see "Timing the Stages").

For real-time monitoring, `tail -f pipeline.log`.

---

## Future Optimizations (Roadmap)

These are ideas, not implemented features:

- **Cython rewrite** of CIGAR parsing and evidence clustering
- **GPU-accelerated CBS** for CNV detection
- **A sub-quadratic CBS scan** (the vectorised scan still scores every arc), so that 1 kb windows become cheap
- **Reusing the gene index across calls** (every `annotate_*` call currently converts the gene list and builds the index again)
- **Batch PCA** for incremental cohort addition
- **Sparse genotype matrices** in the population functions (they are currently rejected)

The measurements above point to the remaining larger costs: the quadratic CBS scan at fine bin sizes, and in `call_structural_variants` the per-read CIGAR parsing and the index build.

---

## Related

- **[GETTING_STARTED.md](GETTING_STARTED.md)** — Basic pipeline setup
- **[CONFIGURATION.md](CONFIGURATION.md)** — All performance-related settings
- **[TROUBLESHOOTING.md](TROUBLESHOOTING.md)** — Diagnose slow runs, OOM
