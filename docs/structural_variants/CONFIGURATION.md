# Configuration: structural_variants

Reference for configuring structural variant analysis: function arguments, the five `SV_` environment variables the code reads, and how to drive the functions from your own settings file.

## Configuration Hierarchy

Settings are resolved in this order (later overrides earlier):

1. **Code defaults** — the default values in each function signature
2. **Function arguments** — what you pass to the call
3. **Environment variables** — five `SV_*` variables replace the matching argument when the function runs (see below)

There is no `config.yaml`, no `config.get()`/`config.set()` registry and no `structural_variants.*` key namespace: nothing in this module reads a configuration file. To keep settings in a file, load it yourself and pass the values as arguments ([Driving Settings From a File](#driving-settings-from-a-file)).

Example:
```bash
# .env or shell
export SV_MIN_SUPPORT=5          # Replaces min_support (default 3), even when passed explicitly
export SV_CNV_SIGNIFICANCE=0.001 # More stringent CBS
```

```python
# Python code: the environment variable above still wins over this argument
from metainformant.structural_variants.detection.sv_calling import call_structural_variants

variants = call_structural_variants(alignments, min_support=10)
```

---

## Environment Variables

Exactly five variables are read. Each is read when the function is called (not at import), parsed as the type shown, and replaces the matching argument. A value that does not parse as that type (for example `5.5` for an int) raises `ValueError`.

| Variable | Type | Read by | Replaces argument | Default |
|----------|------|---------|-------------------|---------|
| `SV_MIN_SUPPORT` | int | `sv_calling.call_structural_variants()`, `quality_filter.filter_by_quality()` | `min_support` | `3` |
| `SV_MIN_QUAL` | float | `quality_filter.filter_by_quality()` | `min_qual` | `20.0` |
| `SV_CNV_SIGNIFICANCE` | float | `cnv.detect_cnv_from_depth()` | `significance` | `0.01` |
| `SV_CNV_DEL_THRESHOLD` | float | `cnv.call_cnv_states()` | `del_threshold` | `-0.3` |
| `SV_CNV_DUP_THRESHOLD` | float | `cnv.call_cnv_states()` | `dup_threshold` | `0.3` |

- The quality threshold is `SV_MIN_QUAL`. An earlier spelling in these docs, `SV_FILTER_MIN_QUAL`, is not read and has no effect.
- `SV_CNV_DEL_THRESHOLD` and `SV_CNV_DUP_THRESHOLD` affect `call_cnv_states()` only. `detect_cnv_from_depth()` assigns states with its built-in thresholds (deletion `-0.3`, duplication `0.3`, amplification `1.0`, homozygous deletion `-1.5`) and ignores both variables.
- `SV_MIN_SUPPORT` acts on both stages, so one export changes the detector and the quality filter together.
- Any other `SV_` variable is ignored.

---

## Function Arguments

Everything else is an ordinary keyword argument. The defaults below are the ones in the signatures. Modules are under `metainformant.structural_variants`.

### Detection

| Argument | Default | Function | Description |
|----------|---------|----------|-------------|
| `min_support` | `3` | `detection.sv_calling.call_structural_variants` | Minimum supporting reads (split reads plus discordant pairs); replaced by `SV_MIN_SUPPORT` |
| `min_mapq` | `20` | `call_structural_variants`, `genotype_sv` | Minimum mapping quality of a read |
| `min_sv_size` | `50` | `call_structural_variants` | Minimum SV size in bp (translocations are exempt) |
| `insert_size_stats` | `None` | `call_structural_variants` | `InsertSizeStats`; estimated from the reads when omitted (mean 400, std 100 if it cannot be estimated) |
| `min_clip` | `20` | `detect_split_reads` | Minimum soft-clip length for split-read evidence. `call_structural_variants` always uses 20 |
| `n_std` | `4.0` | `detect_discordant_pairs` | Standard deviations from the mean insert size that count as discordant. `call_structural_variants` always uses 4.0 |
| `window_size` | `1000` | `detection.cnv.detect_cnv_from_depth` | Bin size in bp; it only scales segment coordinates, the depth arrays must already be binned |
| `method` | `"segmentation"` | `detect_cnv_from_depth` | The only supported method |
| `significance` | `0.01` | `detect_cnv_from_depth`, `segment_coverage` | CBS p-value threshold; replaced by `SV_CNV_SIGNIFICANCE` in `detect_cnv_from_depth` |
| `ploidy` | `2` | `detect_cnv_from_depth`, `call_cnv_states` | Expected ploidy |
| `min_segment_bins` | `3` | `detect_cnv_from_depth` | Minimum bins in a segment |
| `merge_distance` | `1000` | `detect_cnv_from_depth` | Maximum gap (bp) when merging same-state segments |
| `max_gap` | `1000` | `cnv.merge_adjacent_segments` | Same, when called directly |
| `max_iterations`, `min_width` | `100`, `2` | `cnv.segment_coverage` | Recursion depth limit and minimum width (bins) of the pieces tested for a change point; `detect_cnv_from_depth` then drops segments narrower than `min_segment_bins` |
| `del_threshold`, `dup_threshold` | `-0.3`, `0.3` | `cnv.call_cnv_states` | Log2 ratio thresholds; replaced by `SV_CNV_DEL_THRESHOLD` / `SV_CNV_DUP_THRESHOLD` |
| `amp_threshold`, `homodel_threshold` | `1.0`, `-1.5` | `call_cnv_states` | Log2 ratio for amplification and homozygous deletion (arguments only) |
| `pseudocount`, `gc_content` | `1.0`, `None` | `cnv.calculate_log2_ratio` | Pseudocount; per-bin GC fractions switch on GC correction |
| `window`, `min_clip` | `200`, `5` | `detection.breakpoints.refine_breakpoints` | Search window around each breakpoint and minimum soft-clip length |
| `max_distance` | `100` | `breakpoints.cluster_breakpoints` | Maximum distance (bp) between clustered breakpoints |
| `max_homology` | `50` | `breakpoints.detect_microhomology` | Longest microhomology reported |

The evidence clustering distance inside `call_structural_variants` (500 bp) is fixed.

### Annotation

| Argument | Default | Function | Description |
|----------|---------|----------|-------------|
| `max_distance` | `100000` | `annotation.overlap.find_nearest_gene` | Maximum distance (bp) to assign a non-overlapping gene |
| `max_distance` | `1000000` | `annotation.overlap.IntervalIndex.query_nearest` | Maximum search distance (bp) |
| `boundary_window` | `10000` | `annotation.functional_impact.predict_tad_disruption` | Window (bp) around a TAD boundary midpoint. `predict_functional_impact` has no such argument and always uses 10 kb |

The dosage-sensitivity cut-offs inside `assess_dosage_sensitivity` are fixed in the code: a gene is haploinsufficient when its HI score is at least 0.5, its pLI at least 0.9, or its LOEUF below 0.35, and triplosensitive when its TS score is at least 0.5.

### Filtering and Merging

| Argument | Default | Function | Description |
|----------|---------|----------|-------------|
| `min_qual` | `20.0` | `filtering.quality_filter.filter_by_quality` | Minimum quality (Phred); replaced by `SV_MIN_QUAL` |
| `min_support` | `3` | `filter_by_quality` | Minimum supporting reads; replaced by `SV_MIN_SUPPORT` |
| `min_mapq`, `min_gq` | `0.0`, `0.0` | `filter_by_quality` | Minimum mean mapping quality and genotype quality (0 disables) |
| `min_size`, `max_size` | `50`, `None` | `quality_filter.filter_by_size` | Size window in bp |
| `max_af` | `0.01` | `quality_filter.filter_by_frequency` | Remove variants whose population frequency is above this |
| `match_window` | `200` | `filter_by_frequency` | Position window (bp) used as a pre-filter on database entries; the reciprocal-overlap test below decides matches, so the window only matters when `min_reciprocal_overlap` is 0 |
| `min_reciprocal_overlap` | `0.5` | `filter_by_frequency` | Reciprocal overlap needed to match a database entry |
| `min_overlap_fraction` | `0.0` | `quality_filter.apply_blacklist` | A variant is removed when its overlap with a single blacklist region is greater than this fraction of the variant (0.0 removes on any overlap) |
| `min_overlap`, `type_match` | `0.5`, `True` | `filtering.merge.merge_callsets` | Reciprocal overlap threshold; require matching SV types |
| `max_distance`, `min_callers`, `type_match`, `strand_match` | `1000`, `2`, `True`, `False` | `merge.survivor_merge` | Breakpoint distance (bp), minimum callers, type and strand matching |
| `max_distance`, `min_reciprocal_overlap`, `type_match` | `100`, `0.8`, `True` | `merge.deduplicate_variants` | Duplicate collapsing within one callset |

### Population

| Argument | Default | Function | Description |
|----------|---------|----------|-------------|
| `method` | `"depth"` | `population.sv_population.genotype_sv_population` | `"depth"` or `"split"` |
| `n_components` | `10` | `sv_population.sv_population_structure` | Number of principal components |
| `covariates` | `None` | `sv_population.sv_association_test` | Covariate matrix (samples x covariates) |
| `min_overlap`, `max_breakpoint_distance` | `0.5`, `500` | `sv_population.merge_sv_callsets` | Reciprocal overlap and breakpoint distance for cross-sample merging |

`sv_ld_analysis` has no window argument: it computes r² against every SNP you pass.

### Visualization

| Argument | Default | Function |
|----------|---------|----------|
| `figsize` | `(12, 12)` | `visualization.plots.plot_circos` |
| `figsize`, `bin_size` | `(14, 5)`, `1000` | `plots.plot_coverage_track` |
| `figsize`, `log_scale` | `(10, 6)`, `True` | `plots.plot_sv_size_distribution` |
| `figsize` | `(12, 5)` | `plots.plot_sv_type_summary` |
| `figsize`, `flanking` | `(14, 8)`, `500` | `plots.plot_breakpoint_detail` |
| `figsize`, `ploidy` | `(16, 5)`, `2` | `plots.plot_cnv_profile` |

Figures are saved at 150 dpi (fixed). Colours come from the module constants `SV_COLORS` and `CNV_COLORS` in `visualization.plots`; they are not arguments.

### Not Implemented

These settings appeared in earlier versions of this page. Nothing in `src/` reads them, so setting them has no effect:

- Environment variables: `SV_MIN_MAPQ`, `SV_MIN_SV_SIZE`, `SV_CNV_WINDOW_SIZE`, `SV_CNV_AMP_THRESHOLD`, `SV_CNV_HOMODEL_THRESHOLD`, `SV_CNV_MIN_SEGMENT_BINS`, `SV_CNV_MERGE_DISTANCE`, `SV_MIN_CLIP`, `SV_DISCORDANT_N_STD`, `SV_EVIDENCE_CLUSTER_DISTANCE`, `SV_MAX_GENE_DISTANCE`, `SV_TAD_WINDOW`, `SV_IMPACT_*`, `SV_FILTER_*`, `SV_MERGE_*`, `SV_PCA_COMPONENTS`, `SV_ASSOCIATION_METHOD`, `SV_LD_WINDOW`, `SV_VIS_*`, `SV_LOG_LEVEL`, `SV_OUTPUT_DIR`, `SV_CACHE_TTL`
- A `structural_variants:` YAML tree, `config.yaml` / `~/.hermes/config.yaml` lookup, `config.get()`, `config.set()`, `config.validate_schema()`
- `METAINFORMANT_DATA` data-directory resolution and `structural_variants.features.*` flags (`experimental_merge`, `gc_correction`, `microhomology_detection`)

For logging, use the core helpers (`configure_logging_from_env()` reads `CORE_LOG_LEVEL`); there is no module-specific log level.

---

## Driving Settings From a File

The module does not read configuration files, but a settings file of your own is easy to feed in with the core loader. `load_mapping_from_file` reads YAML, TOML or JSON into a dict:

```yaml
# sv_settings.yaml (your own file and schema; nothing in structural_variants reads it)
detection:
  min_support: 3
  min_mapq: 20
  min_sv_size: 50
quality:
  min_qual: 20.0
  min_support: 3
```

```python
from metainformant.core.utils.config import load_mapping_from_file
from metainformant.structural_variants.detection.sv_calling import call_structural_variants
from metainformant.structural_variants.filtering.quality_filter import filter_by_quality

cfg = load_mapping_from_file("sv_settings.yaml")
variants = call_structural_variants(alignments, **cfg["detection"])
filtered, stats = filter_by_quality(sv_dicts, **cfg["quality"])
```

The five environment variables above still override these values when they are set.

---

## Python API

### Reading the Defaults

```python
import inspect
from metainformant.structural_variants.detection.cnv import detect_cnv_from_depth

sig = inspect.signature(detect_cnv_from_depth)
for name, param in sig.parameters.items():
    if param.default is not inspect.Parameter.empty:
        print(f"{name} = {param.default!r}")
```

### Setting the Environment Variables From Python

```python
import os

# Takes effect for every later call: each function reads the variable when it runs
os.environ["SV_MIN_SUPPORT"] = "5"
os.environ["SV_CNV_SIGNIFICANCE"] = "0.001"

# Remove the override again
os.environ.pop("SV_MIN_SUPPORT", None)
```

### Seeing Which Values a Call Used

Results record the values that were applied after any environment override:

```python
import numpy as np
from metainformant.structural_variants.detection.cnv import detect_cnv_from_depth
from metainformant.structural_variants.filtering.quality_filter import filter_by_quality

depth = {"chr1": np.concatenate([np.full(40, 30.0), np.full(20, 15.0), np.full(40, 30.0)])}
results = detect_cnv_from_depth(depth)
print(results["chr1"].parameters)    # {'significance': 0.01, 'ploidy': 2, 'min_segment_bins': 3, 'merge_distance': 1000}

kept, stats = filter_by_quality([{"quality": 45.0, "support": 6}])
print(stats.parameters)              # {'min_qual': 20.0, 'min_support': 3, 'min_mapq': 0.0, 'min_gq': 0.0}
```

---

## Real-World Settings by Use Case

### Use Case 1: Germline SV Discovery (High Sensitivity)

```bash
# Maximize recall, sacrifice some precision
export SV_MIN_SUPPORT=2           # detector and quality filter
export SV_CNV_SIGNIFICANCE=0.05   # Less stringent CBS
export SV_MIN_QUAL=10             # Low quality threshold in filter_by_quality
```

Expect: More candidates → more manual review; higher false positive rate. The soft-clip threshold is fixed at 20 bp inside `call_structural_variants`.

### Use Case 2: Somatic CNV Calling (Tumor vs Normal)

```bash
# Paired tumor-normal with matched normal depth
# More stringent because somatic expected to be clear
export SV_CNV_SIGNIFICANCE=0.001  # Stricter CBS change-point test
export SV_CNV_DEL_THRESHOLD=-0.4  # call_cnv_states only
export SV_CNV_DUP_THRESHOLD=0.4
```

`detect_cnv_from_depth` keeps its built-in state thresholds. To classify with the stricter thresholds, pass the segments to `call_cnv_states` (it does not recompute `confidence`, which comes back as 0.0):

```python
from metainformant.structural_variants.detection.cnv import call_cnv_states

# results: the dict returned by detect_cnv_from_depth()
segments = [seg for result in results.values() for seg in result.segments]
recalled = call_cnv_states(segments, del_threshold=-0.4, dup_threshold=0.4)
```

### Use Case 3: Clinical Diagnostic Pipeline

```python
from metainformant.structural_variants.detection.sv_calling import call_structural_variants
from metainformant.structural_variants.filtering.merge import merge_callsets
from metainformant.structural_variants.filtering.quality_filter import filter_by_frequency, filter_by_quality

# Settings that favour specificity
variants = call_structural_variants(alignments, min_support=10, min_sv_size=100)
filtered, _ = filter_by_quality(sv_dicts, min_qual=50.0, min_support=10)
rare, _ = filter_by_frequency(filtered, population_db=gnomad_sv, max_af=0.001)  # Rare variants only
merged, _ = merge_callsets(callsets, min_overlap=0.8)  # Stricter merging
```

### Use Case 4: Large Cohort (10,000+ samples)

```python
from metainformant.structural_variants.detection.sv_calling import call_structural_variants
from metainformant.structural_variants.filtering.merge import merge_callsets
from metainformant.structural_variants.filtering.quality_filter import filter_by_quality
from metainformant.structural_variants.population.sv_population import sv_population_structure

# Fewer low-quality reads and candidates, fewer components kept
variants = call_structural_variants(alignments, min_mapq=30)
filtered, _ = filter_by_quality(sv_dicts, min_qual=30)
merged, _ = merge_callsets(callsets, min_overlap=0.7)
pca = sv_population_structure(genotype_matrix, n_components=5)
```

---

## Database & Reference Data

The annotation and filtering functions take their reference data as in-memory Python objects. Nothing in the module locates or loads reference files for you (the only file reading is the VCF input of `survivor_merge`), and no data-directory environment variable is read. Load the files yourself, for example with `core.io.load_json`:

```python
from metainformant.core import io

gene_db = io.load_json("genes.json")              # [{"chrom": "chr1", "start": 1000, "end": 5000, "name": "GENE1"}, ...]
hi_db = io.load_json("clingen_hi.json")           # {gene: haploinsufficiency score (0-1)}
tad_bounds = io.load_json("tad_boundaries.json")  # [{"chrom": "chr1", "start": 990000, "end": 1010000}, ...]
pop_db = io.load_json("gnomad_sv.json")           # [{"chrom": "chr1", "start": 1000, "end": 5000, "sv_type": "DEL", "af": 0.003}, ...]
```

| Argument | Accepted layouts |
|----------|------------------|
| `gene_db` (`annotate_gene_overlap`, `find_nearest_gene`), `regulatory_db` | List of dicts with `chrom`, `start`, `end`, `name`, or a list of `GenomicInterval` (any `is_coding` key is ignored) |
| `gene_annotations` (`predict_functional_impact`) | List of dicts with `chrom`, `start`, `end`, `name` and optionally `is_coding` (default coding); dicts only, not `GenomicInterval` |
| `haploinsufficiency_db` | `{gene: score}` |
| `tad_boundaries` | List of dicts with `chrom`, `start`, `end` (and optionally `genes`) |
| `population_db` | List of dicts with `chrom`, `start`, `end`, `af` (and optionally `sv_type`), or `{"chrom:start-end": af}` |
| `blacklist_regions` | List of dicts with `chrom`, `start`, `end`, or a list of `(chrom, start, end)` tuples |

---

## Optional Behaviour

There are no feature flags. The behaviours the old flag table described are plain function-level choices:

- **GC-bias correction** — pass per-bin GC fractions as `gc_content=` to `cnv.calculate_log2_ratio()`.
- **Microhomology scan** — `breakpoints.detect_microhomology()` can be called on its own; `refine_breakpoints()` runs it once, on the first read that carries a `seq` on the variant's chromosome, and reports the result on the first breakpoint only.
- **Experimental merging** — no such option exists; use `merge_callsets()` or `survivor_merge()`.

---

## Troubleshooting Configuration

### Problem: SV_ env var not taking effect

1. Confirm the name is one of the five in [Environment Variables](#environment-variables); any other `SV_` variable is ignored.
2. Confirm variable is exported: `echo $SV_MIN_SUPPORT`
3. Verify the process inherits the environment (check `ps eww <pid>`).
4. The variable is read each time the function runs, so no restart is needed after setting it from Python. A value that does not parse as the expected type raises `ValueError` instead of falling back to the default.

### Problem: an argument seems to be ignored

The environment variable wins over the argument. List what is set and unset it if it should not apply:

```bash
env | grep '^SV_'
unset SV_MIN_SUPPORT
```

### Problem: a settings file is not being used

The module reads no settings file. Load it yourself and pass the values as arguments ([Driving Settings From a File](#driving-settings-from-a-file)).

---

## Complete Example: Configuration-Driven Run

```python-snippet
#!/usr/bin/env python3
"""Run an SV workflow with environment-driven configuration."""

import os
from dataclasses import asdict
from metainformant.core import io
from metainformant.core.utils import logging as ml_logging
from metainformant.structural_variants.annotation.overlap import annotate_gene_overlap
from metainformant.structural_variants.detection.cnv import detect_cnv_from_depth
from metainformant.structural_variants.detection.sv_calling import call_structural_variants

# ── Configure via environment (read when each function is called; the value
#    overrides the corresponding function argument) ──────────────────────────
os.environ["SV_MIN_SUPPORT"] = "5"
os.environ["SV_CNV_SIGNIFICANCE"] = "0.001"

# Setup logging
ml_logging.setup_logger("metainformant.structural_variants", log_file="results/sv.log", level="INFO")

# Run (alignments: read-record dicts and depth_data: {chrom: per-bin depth}, see GETTING_STARTED.md)
variants = call_structural_variants(alignments)                                   # SV_MIN_SUPPORT applies
cnvs = detect_cnv_from_depth(depth_data)                                          # SV_CNV_SIGNIFICANCE applies
annotated = annotate_gene_overlap([asdict(v) for v in variants], io.load_json("genes.json"))

print(f"Completed: {len(annotated)} SVs, {sum(len(r.segments) for r in cnvs.values())} CNV segments")
```

---

## Validation & Debug

Check which of the five overrides are active:

```python
import os

for key in ("SV_MIN_SUPPORT", "SV_MIN_QUAL", "SV_CNV_SIGNIFICANCE", "SV_CNV_DEL_THRESHOLD", "SV_CNV_DUP_THRESHOLD"):
    if key in os.environ:
        print(f"{key} = {os.environ[key]}")
```

To see what a call actually used, read the `parameters` recorded on its result (`CNVResult.parameters`, `FilterStats.parameters`; see [Seeing Which Values a Call Used](#seeing-which-values-a-call-used)). `call_structural_variants` returns a plain list of `StructuralVariant` objects and records nothing.

---

## What's Next?

- **[GETTING_STARTED.md](GETTING_STARTED.md)** — See these settings in action
- **[PERFORMANCE.md](PERFORMANCE.md)** — Tune these for large datasets
- **[TROUBLESHOOTING.md](TROUBLESHOOTING.md)** — Configuration-related errors
