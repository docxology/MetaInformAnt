# Performance baselines — tests/perf/

Recorded, reproducible baselines for the campaign-scale benchmark scripts.
All runs use **synthetic data only** (deterministic seeds, see each script)
on the local checkout; never benchmark against the live campaign data root.

Environment for the baselines below: macOS (arm64), Python 3.12 via `uv`,
pandas 2.x / numpy 1.x, external-drive checkout under live campaign I/O
(2026-09-01). Absolute times are machine- and load-dependent — treat them as
order-of-magnitude references; the **CONSISTENCY: OK** line is the portable
invariant.

## How to run

```bash
uv run python tests/perf/bench_tau_expression.py
uv run python tests/perf/bench_ortholog_classification.py
uv run python tests/perf/bench_batches_memory.py
```

Each script prints a `RESULT-LINE:` suitable for pasting into the tables
below after a re-run.

## 1. Tau tissue-specificity (`bench_tau_expression.py`)

Frame: 8,000 genes x 120 tissue columns (~7.7 MiB deep) — the campaign's
per-species upper envelope. Default: lowest-10% mean-expression filter
(Xu & Colgan 2025) + log2-free tau.

| Metric | Baseline (2026-09-01) |
| --- | --- |
| `compute_tau` (whole frame) | 0.010 s (~96–98 M cells/s) |
| `tau_summary` | 0.0009 s |
| Global-filter + chunked pass (2k-row chunks) | 0.012 s |
| Chunk/full consistency | OK (exact: 3125.709325 both paths) |

**Measured pitfall (load-bearing):** `compute_tau` applies its
`lowest_fraction` filter **per call**, over whatever rows it receives.
Chunking raw rows therefore selects different genes per chunk and shifts
tau results (~0.14% drift observed on the default frame before the fix).
The correct chunked pattern is: filter **globally** first, then chunk the
filtered frame (pure tau arithmetic is row-local), passing
`lowest_fraction=0.0` per chunk. The benchmark asserts this consistency
and fails if the pattern regresses.

## 2. Ortholog classification (`bench_ortholog_classification.py`)

Bridge table: 150,000 orthogroups x 27 species (campaign roster scale),
~20% multicopy cells. Plus a 200k-transcript expression join for one
species.

| Metric | Baseline (2026-09-01) |
| --- | --- |
| bridge generation (150k x 27) | 4.228 s |
| `classify_orthogroups` (whole table) | 4.916 s (~30,500 OG/s) |
| `join_expression_with_orthology` (200k transcripts) | 8.529 s |
| chunked classification (25k-OG chunks) | 4.872 s (counts identical) |
| class-count consistency | OK (multicopy 141,315 / single_copy 8,685) |
| bridge deep memory | 68.5 MiB |

**Profiled headroom (for the rna/analysis lane — outside perf-lane scope):**
`classify_orthogroups` is bound by `iterrows()` + per-cell `_parse_ids`
(`tissue_specificity.py:190`). A vectorized equivalent
(`bridge_table.applymap`-free: split all cells at once via
`df.stack().str.split(",")`, count lengths per row) projects well above the
2x threshold given ~30.5k OG/s current throughput. Left unimplemented here
by lane-scope discipline; measured evidence recorded in
PROJECT_STATE_REPORT_2026-09-01_R5_PERF.md.

Note: classification is row-local (per-orthogroup), so chunked
classification is exactly consistent; class counts must match the
whole-table pass.

## 3. Memory hygiene (`bench_batches_memory.py`)

Frame: 500,000 rows x 3 cols (~16 MiB deep). Whole-frame vs chunked
sum over `metainformant.core.utils.batches.iter_frame_chunks`.

| Metric | Baseline (2026-09-01) |
| --- | --- |
| whole-frame sum (500k rows) | 0.0004 s |
| chunked sum (50k-row chunks) | 0.0005 s (0.92x — same speed) |
| consistency | OK (exact) |
| deep memory | 16.2 MiB (34 B/row) |

The chunked path costs ~nothing at this scale because vectorized pandas
sum is per-chunk; its value is bounded peak memory for operations that
materialize intermediates (joins, sorts, dtype conversions), not for
single-pass reductions.

## Scope and boundary

- Benchmarks are time-real only: real pandas/numpy on synthetic frames.
  No mock frameworks (zero-mocks policy).
- Synthetic distributions (log-normal expression, exponential counts)
  are shape-plausible for TPM/count tables but are not biological data.
- The live campaign (`run_all_species.py`, data root) is never touched.
