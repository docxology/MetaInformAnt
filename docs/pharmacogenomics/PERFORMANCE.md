# Performance: Pharmacogenomics

## Per-call cost

| Operation | Mean time per call | Notes |
|-----------|--------------------|-------|
| `call_star_alleles(variants, 'CYP2D6', allele_definitions=...)` | ~4 µs | 1 to 3 observed variants against the 12 CYP2D6 definitions |
| `determine_diplotype()` | ~1.4 µs | two table lookups and a sum |
| `classify_phenotype('*1/*4', 'CYP2D6')` | ~2.6 µs | string diplotype in, band or allele-function rule applied |
| `predict_metabolizer()` | ~5 µs | dict front end over the same classification |
| variants → phenotype (`call_star_alleles` → `determine_diplotype` → `classify_phenotype`) | ~11 µs | the whole per-sample pipeline |
| `get_dosing_recommendation('codeine', 'Poor Metabolizer')` | ~1 µs | list scan of the 9 built-in guidelines |
| `generate_clinical_report(...)` | ~0.26 ms | two genes, all CPIC drugs for them |
| `export_report(report, format='text')` | ~11 µs | |

Measured with the benchmark script below on an Apple-silicon laptop (Python 3.12, one core,
default logging): 10,000 calls per row, 1,000 for the report functions. Expect different
absolute numbers on other hardware; the per-sample cost does not depend on cohort size, so
cohort time grows linearly with the number of samples.

## Cohort processing

The per-sample functions are plain Python, so a cohort is processed by mapping a function over
its samples with `metainformant.core.execution.parallel`:

- `parallel.thread_map(func, items, max_workers=8)`: thread pool, results in input order
- `parallel.process_map(func, items, max_workers=None)`: process pool for CPU-bound work;
  `func` must be a top-level (picklable) function
- `parallel.parallel_batch(func, items, batch_size=10, max_workers=4)`: splits `items` into
  batches; `func` takes a list and returns a list

There is no `chunk_size` setting in the module. Because one sample costs microseconds, process
start-up and pickling dominate small jobs; give each worker a batch of many samples (for
example with `parallel_batch`, or `np.array_split` plus `process_map`) rather than one task per
sample.

```python-snippet
from metainformant.core.execution import parallel
workers = parallel.resource_aware_workers(task_type="cpu", max_cap=8)
results = parallel.process_map(func, items, max_workers=workers)
```

The matching code is pure Python, so threads share the GIL; use `process_map` for CPU-bound
cohorts.

## What a call costs

| Stage | Cost model |
|-------|------------|
| Allele tables | Python literals created at import: no parsing on the first call, no warm-up per process |
| `call_star_alleles` | O(n × m): n definitions of the gene, m defining variants each; set operations on the observed variants |
| `determine_diplotype` | two dictionary lookups |
| `classify_phenotype` | a walk over at most five score bands, or one lookup in the allele-function rules |
| CPIC lookup | a scan of the 9 built-in guidelines |
| Report | classification plus lookup per gene and drug, then string formatting |

The largest built-in table is CYP2D6 with 12 calling definitions; the others have 3 to 7
(CYP2C19 5, CYP2C9 7, CYP3A5 4, DPYD 5, TPMT 5, NUDT15 3, SLCO1B1 4). A custom definition file
with many alleles raises the matching cost linearly with its size, so load it once and pass the
list to every call:

```python-snippet
definitions = load_allele_definitions('CYP2E1', filepath='cyp2e1.json')
alleles = call_star_alleles(variants, 'CYP2E1', allele_definitions=definitions)
```

The cost also grows linearly with the number of observed variants per sample, since each is
tested for membership in each definition's variant set. Pre-filtering a sample's variants to
those in the gene's definitions (`{v for a in definitions for v in a.defining_variants}`) shortens
the sets when a sample carries many unrelated identifiers.

## Profiling recipe

```bash
python -m cProfile -o pgx.prof myscript.py
snakeviz pgx.prof
```

Look for accumulated time in `alleles/star_allele.py:call_star_alleles` (the matching loop).

## Caching strategies

Repeated calls on identical variant sets can be memoized:

```python-snippet
from metainformant.core.io.cache import JsonCache
from metainformant.pharmacogenomics.alleles.star_allele import call_star_alleles

cache = JsonCache('output/cache/pgx', ttl_seconds=86400)

def cached_allele_calls(variants, gene='CYP2D6'):
    key = f"{gene}_" + "_".join(sorted(variants))
    names = cache.get(key)
    if names is None:
        names = [a.name for a in call_star_alleles(set(variants), gene)]
        cache.set(key, names)
    return names
```

At a few microseconds per call, a file-backed cache only pays off for the slow steps around the
call (a large report batch, or I/O); for the matching itself, an in-memory
`functools.lru_cache` keyed on `frozenset(variants)` is cheaper.

## Scaling to >1M samples

For national-scale biobank processing:

- Use Dask distributed (`dask.distributed.Client`): each worker is a Python process that
  imports the module once.
- Store results in Parquet with `df.to_parquet('out/', partitioning_on=['phenotype'])` to reduce downstream I/O.
- Disable `parallel` in child workers to avoid nested pools (guard with `if __name__ == '__main__'`).

## Limitations & edge cases

- CYP2D6 copy-number variants are not inferred from variants: pass the copy number to `handle_cyp2d6_cnv`; every other function treats the gene as diploid.
- A variant outside every allele definition is ignored; it contributes nothing to matching (`detect_novel_alleles` lists such variants).

## Runtime checklist

- [ ] Are you CPU-bound? → `parallel.process_map` with `parallel.resource_aware_workers(task_type="cpu")` workers and a top-level function
- [ ] Is each task a single sample? → give each worker a batch of samples
- [ ] Loading a custom definition file per sample? → load it once with `load_allele_definitions(gene, filepath=...)` and pass `allele_definitions=definitions`
- [ ] Memory spiking? → do not accumulate results in a list inside the mapped function; let `parallel.thread_map` collect

## Benchmark script (reproducible)

```python
import time

import numpy as np

from metainformant.pharmacogenomics.alleles.diplotype import determine_diplotype
from metainformant.pharmacogenomics.alleles.phenotype import classify_phenotype
from metainformant.pharmacogenomics.alleles.star_allele import call_star_alleles, load_allele_definitions

N = 10_000
rng = np.random.default_rng(42)
definitions = load_allele_definitions('CYP2D6')
# Random variant sets drawn from the identifiers the CYP2D6 table knows
pool = sorted({v for a in definitions for v in a.defining_variants})
variants_list = [
    set(rng.choice(pool, size=int(rng.integers(1, 4)), replace=False).tolist()) for _ in range(N)
]

def phenotype_of(variants):
    names = [a.name for a in call_star_alleles(variants, 'CYP2D6', allele_definitions=definitions)]
    names += ['*1', '*1']  # pad: a single called allele pairs with *1
    return classify_phenotype(determine_diplotype(names[0], names[1], 'CYP2D6'), 'CYP2D6')

start = time.perf_counter()
for variants in variants_list:
    phenotype_of(variants)
elapsed = time.perf_counter() - start
print(f'{N} samples → {elapsed:.2f}s  ({N/elapsed:.0f} samples/s)')
```

## I/O considerations

If your cohort is stored in a CSV/Parquet file, I/O often dominates.

```python-snippet
import polars as pl  # faster than pandas for this

df = pl.read_csv('cohort.csv', columns=['rsids'])
# use `filter` to limit to genes of interest early
```

## Tips for HPC clusters

- Size the pool from the scheduler allocation:
  `parallel.process_map(func, items, max_workers=int(os.environ['SLURM_CPUS_PER_TASK']))`.
- Call `load_allele_definitions(gene)` once per gene in the master and pass the result to
  `call_star_alleles(..., allele_definitions=definitions)` so child processes do not rebuild it.
- If using MPI (`mpi4py`), wrap the map in `if rank == 0:` to avoid duplicated work.
