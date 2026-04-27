# Configuration: structural_variants

Complete reference for structural variant analysis configuration through environment variables, YAML files, and the Python API. All settings use the `SV_` prefix for clarity.

## Configuration Hierarchy

Settings are resolved in this order (later overrides earlier):

1. **Code defaults** — Built-in constants in each function signature
2. **YAML config file** — `config.yaml` or `~/.hermes/config.yaml`
3. **Environment variables** — `SV_*` variables in shell
4. **Runtime Python API** — `config.set()` during execution

Example:
```bash
# .env or shell
export SV_MIN_SUPPORT=5          # Overrides default 3
export SV_CNV_SIGNIFICANCE=0.001 # More stringent CBS
```

```python
# Python code (highest priority)
from metainformant import config
config.set("structural_variants.min_support", 10)
```

---

## Environment Variables

All structural variant settings use the `SV_` prefix to avoid collisions.

### Detection Settings

| Variable | Default | Description | Used By |
|----------|---------|-------------|---------|
| `SV_MIN_SUPPORT` | `3` | Minimum supporting reads for an SV call | `sv_calling.call_structural_variants()` |
| `SV_MIN_MAPQ` | `20` | Minimum mapping quality to consider a read | `sv_calling` |
| `SV_MIN_SV_SIZE` | `50` | Minimum SV size in base pairs | `sv_calling`, `quality_filter.filter_by_size()` |
| `SV_CNV_SIGNIFICANCE` | `0.01` | CBS p-value threshold for change-point detection | `cnv.segment_coverage()` |
| `SV_CNV_WINDOW_SIZE` | `1000` | CNV bin size in bp | `cnv.detect_cnv_from_depth()` |
| `SV_CNV_DEL_THRESHOLD` | `-0.3` | Log2 ratio threshold for deletion | `cnv.call_cnv_states()` |
| `SV_CNV_DUP_THRESHOLD` | `0.3` | Log2 ratio threshold for duplication | `cnv.call_cnv_states()` |
| `SV_CNV_AMP_THRESHOLD` | `1.0` | Log2 ratio for amplification | `cnv.call_cnv_states()` |
| `SV_CNV_HOMODEL_THRESHOLD` | `-1.5` | Log2 ratio for homozygous deletion | `cnv.call_cnv_states()` |
| `SV_CNV_MIN_SEGMENT_BINS` | `3` | Minimum bins for a valid CNV segment | `cnv.detect_cnv_from_depth()` |
| `SV_CNV_MERGE_DISTANCE` | `1000` | Max gap to merge same-state segments | `cnv.merge_adjacent_segments()` |
| `SV_MIN_CLIP` | `20` | Minimum soft-clip length for split-read evidence | `sv_calling.detect_split_reads()` |
| `SV_DISCORDANT_N_STD` | `4.0` | σ threshold for insert size discordance | `sv_calling.detect_discordant_pairs()` |
| `SV_EVIDENCE_CLUSTER_DISTANCE` | `500` | Max bp between evidence items to cluster | `sv_calling._cluster_evidence()` |

### Annotation Settings

| Variable | Default | Description |
|----------|---------|-------------|
| `SV_MAX_GENE_DISTANCE` | `100_000` | Max distance to assign nearest gene (not overlapping) |
| `SV_TAD_WINDOW` | `10_000` | Window around TAD boundary midpoint for disruption | |
| `SV_IMPACT_HI_THRESHOLD` | `0.5` | HI score threshold for haploinsufficiency |
| `SV_IMPACT_PLI_THRESHOLD` | `0.9` | pLI threshold for haploinsufficiency |
| `SV_IMPACT_LOEUF_THRESHOLD` | `0.35` | LOEUF upper fraction threshold (lower = more constrained) |

### Filtering Settings

| Variable | Default | Description |
|----------|---------|-------------|
| `SV_FILTER_MIN_QUAL` | `20.0` | Minimum variant quality (Phred) |
| `SV_FILTER_MIN_SUPPORT` | `3` | Minimum supporting reads (overrides detector) |
| `SV_FILTER_MIN_MAPQ` | `0.0` | Minimum mean mapping quality |
| `SV_FILTER_MIN_SIZE` | `50` | Minimum variant size |
| `SV_FILTER_MAX_SIZE` | *(none)* | Maximum variant size |
| `SV_FILTER_MAX_AF` | `0.01` | Maximum allele frequency to keep (remove common) |
| `SV_FILTER_MATCH_WINDOW` | `200` | Position matching window for frequency DB lookup |
| `SV_FILTER_MIN_RECIPROCAL_OVERLAP` | `0.5` | Overlap fraction for frequency matching |

### Multi-Caller Merge Settings

| Variable | Default | Description |
|----------|---------|-------------|
| `SV_MERGE_MIN_OVERLAP` | `0.5` | Reciprocal overlap threshold for merging |
| `SV_MERGE_TYPE_MATCH` | `true` | Require matching SV types |
| `SV_MERGE_MAX_BP_DISTANCE` | `1000` | Breakpoint distance for SURVIVOR merge |

### Population Analysis Settings

| Variable | Default | Description |
|----------|---------|-------------|
| `SV_PCA_COMPONENTS` | `10` | Number of principal components |
| `SV_ASSOCIATION_METHOD` | `"linear"` | Default regression method |
| `SV_LD_WINDOW` | `1_000_000` | Max distance for LD calculation (bp) |

### Visualization Settings

| Variable | Default | Description |
|----------|---------|-------------|
| `SV_VIS_FIGSIZE` | `(12, 12)` | Default figure size inches |
| `SV_VIS_DPI` | `150` | Figure DPI |
| `SV_VIS_COLOR_SCHEME` | `default` | Color palette name |

### General & Logging

| Variable | Default | Description |
|----------|---------|-------------|
| `SV_LOG_LEVEL` | *(inherits `CORE_LOG_LEVEL`)* | Module-specific logging |
| `SV_OUTPUT_DIR` | `./output` | Default output directory |
| `SV_CACHE_TTL` | `3600` | Result cache TTL in seconds |

---

## YAML Configuration

Create a project-local `config.yaml` in the repository root:

```yaml
# Structural Variants Configuration
structural_variants:
  # Detection
  detection:
    min_support: 3
    min_mapq: 20
    min_sv_size: 50
    cnv:
      significance: 0.01
      window_size: 1000
      del_threshold: -0.3
      dup_threshold: 0.3

  # Annotation
  annotation:
    max_gene_distance: 100000
    tad_window: 10000
    impact:
      hi_threshold: 0.5
      pli_threshold: 0.9
      loeuf_threshold: 0.35

  # Filtering
  filtering:
    quality:
      min_qual: 20.0
      min_support: 3
      min_mapq: 30.0
    size:
      min_size: 50
      max_size: 50000000  # 50 Mb
    frequency:
      max_af: 0.01
      match_window: 200
      min_reciprocal_overlap: 0.5
    blacklist:
      enabled: true

  # Multi-caller merging
  merging:
    min_overlap: 0.5
    type_match: true
    max_breakpoint_distance: 1000

  # Population
  population:
    pca_components: 10
    association_method: linear
    ld_window: 1000000

  # Visualization
  visualization:
    figsize: [12, 12]
    dpi: 150
    color_scheme: default

  # I/O
  output_dir: ./output
  cache_ttl: 3600
```

Alternatively, place at `~/.hermes/config.yaml` for user-wide defaults.

---

## Python API

### Accessing Configuration

```python
from metainformant import config

# Get entire module config as dict
sv_config = config.get("structural_variants")
print(sv_config)
# {'detection': {'min_support': 3, ...}, 'filtering': {...}, ...}

# Get specific nested key
min_support = config.get("structural_variants.detection.min_support")
print(min_support)  # 3 (or env-overridden value)

# Get with default fallback
chunk_size = config.get("structural_variants.chunk_size", default=5000)
```

### Setting Configuration

```python
# Set top-level key
config.set("structural_variants.min_support", 5)

# Set nested keys
config.set("structural_variants.detection.cnv.significance", 0.001)
config.set("structural_variants.filtering.quality.min_qual", 30.0)

# Override entire sub-dictionary
new_detection = {"min_support": 10, "min_mapq": 50}
config.set("structural_variants.detection", new_detection)

# Current values accessible anywhere
current_min_support = config.get("structural_variants.min_support")
```

### Validation

```python
from metainformant import config

# Validate schema (raises if invalid)
errors = config.validate_schema("structural_variants")
if errors:
    for err in errors:
        print(f"  - {err.path}: {err.message}")
    raise ValueError("Invalid configuration")
else:
    print("Configuration valid ✓")
```

### Module-Specific Configuration Function (Optional)

Some modules expose `configure()` function:

```python
from metainformant.structural_variants import configure

# High-level configuration
configure(
    detection=dict(method="fast", min_support=3),
    filtering=dict(min_quality=25, min_size=100),
    output_dir="results/",
)
```

(Implementation-dependent; check `__init__.py` for availability.)

---

## Real-World Settings by Use Case

### Use Case 1: Germline SV Discovery (High Sensitivity)

```bash
# Maximize recall, sacrifice some precision
export SV_MIN_SUPPORT=2
export SV_CNV_SIGNIFICANCE=0.05   # Less stringent CBS
export SV_MIN_CLIP=15             # Shorter clips acceptable
export SV_FILTER_MIN_QUAL=10      # Low quality threshold
```

Expect: More candidates → more manual review; higher false positive rate.

### Use Case 2: Somatic CNV Calling (Tumor vs Normal)

```bash
# Paired tumor-normal with matched normal depth
# More stringent because somatic expected to be clear
export SV_CNV_SIGNIFICANCE=0.001  # Bonferroni for multiple testing
export SV_CNV_DEL_THRESHOLD=-0.4  # More conservative DEL threshold
export SV_CNV_DUP_THRESHOLD=0.4
```

### Use Case 3: Clinical Diagnostic Pipeline

```python
# Certified clinical settings: maximum specificity
config.set("structural_variants.min_support", 10)
config.set("structural_variants.min_quality", 50.0)
config.set("structural_variants.min_sv_size", 100)
config.set("structural_variants.filtering.frequency.max_af", 0.001)  # Rare variants only
config.set("structural_variants.merging.min_overlap", 0.8)  # Stricter merging
```

### Use Case 4: Large Cohort (10,000+ samples)

```python
# Focus on runtime efficiency
config.set("structural_variants.detection.min_mapq", 30)  # Fewer reads to process
config.set("structural_variants.filtering.quality.min_qual", 30)  # Reduce candidates
config.set("structural_variants.merging.min_overlap", 0.7)  # Faster graph (fewer edges)
config.set("structural_variants.population.pca_components", 5)  # Less PCA overhead
```

---

## Database & Reference File Paths

Many annotation steps require external reference files. Paths can be absolute or relative to `METAINFORMANT_DATA` environment variable.

```bash
export METAINFORMANT_DATA=/path/to/metainformant_data

# Directory layout expected:
# $METAINFORMANT_DATA/
#   annotations/
#     genes.bed           # BED6: chrom,start,end,name,score,strand
#     enhancers.bed       # Regulatory elements
#   dosage/
#     clingen_hi.json     # {gene: HI_score}
#     gnomad_pli.json     # {gene: pLI}
#   tad/
#     hg38_tad_boundaries.json  # [{"chrom":..., "start":..., "end":...}, ...]
#   population/
#     gnomad_sv.json      # {chr: [{pos, sv_type, af}]}
```

References loaded via:

```python
from metainformant.core import paths

# Resolve with environment override
gene_path = paths.resolve_data_path("annotations/genes.bed")
hi_db = paths.resolve_data_path("dosage/clingen_hi.json")
```

**Register custom data location:**

```python
config.set("core.data_dir", "/my/data")
# or
export METAINFORMANT_DATA="/my/data"
```

---

## Advanced: Feature Flags

Some module behavior controlled by less-documented flags:

| Setting | Values | Purpose |
|---------|--------|---------|
| `structural_variants.features.experimental_merge` | `true`/`false` | Enable experimental edge-weight merging (future) |
| `structural_variants.features.gc_correction` | `true`/`false` | Enable GC-bias correction in CNV |
| `structural_variants.features.microhomology_detection` | `true`/`false` | Enable microhomology scan (extra cost) |

These may change between releases; use caution.

---

## Troubleshooting Configuration

### Problem: SV_ env var not taking effect

1. Confirm variable is exported: `echo $SV_MIN_SUPPORT`
2. Verify process inherits env (check `ps eww <pid>`)
3. Environmental read happens at module import. Restart Python after setting vars.

### Problem: YAML file not loading

```bash
# Check file location
ls -la config.yaml           # Repository root
ls -la ~/.hermes/config.yaml # User home

# Validate YAML syntax
python3 -c "import yaml; yaml.safe_load(open('config.yaml'))"
```

### Problem: Config value not being used

Some functions read environment variables only once at import. To force reread:

```python
import importlib
import metainformant.structural_variants.detection.cnv as cnv_module
importlib.reload(cnv_module)  # Re-reads env vars on next call
```

Better: Use Python API `config.set()` for dynamic overrides.

### Problem: Nested key not found

Use dot-notation for nested:

```python
# Correct
config.get("structural_variants.detection.min_support")

# Wrong (no nested lookup)
config.get("structural_variants")["detection"]["min_support"]  # works but fragile
```

---

## Complete Example: Configuration-Driven Run

```python
#!/usr/bin/env python3
"""Run SV pipeline with full configuration control."""

import os
from pathlib import Path
from metainformant import config, logging as ml_logging
from metainformant.structural_variants import pipeline

# ── Configure via environment (highest priority after code defaults) ────────
os.environ["SV_MIN_SUPPORT"] = "5"
os.environ["SV_CNV_SIGNIFICANCE"] = "0.001"
os.environ["SV_FILTER_MAX_AF"] = "0.005"

# Or via Python API (overrides env)
config.set("structural_variants.filtering.quality.min_qual", 30.0)
config.set("structural_variants.output_dir", "results/")

# Setup logging
ml_logging.setup_logging(level="INFO")

# Run
results = pipeline.run(
    bam="sample.bam",
    gene_annotations="genes.bed",
    dosage_db="dosage_scores.json",
    tad_boundaries="tads.json",
)

print(f"Completed: {len(results['filtered_svs'])} variants")
```

---

## Default Values Reference

Resolved effective configuration can be printed:

```python
from metainformant import config
from metainformant.structural_variants import get_default_config

defaults = get_default_config()  # Module-provided helper
print("\n".join(f"{k}: {v}" for k, v in defaults.items()))
```

Or inspect module source:

```python
import inspect
from metainformant.structural_variants.detection.cnv import detect_cnv_from_depth
sig = inspect.signature(detect_cnv_from_depth)
for name, param in sig.parameters.items():
    if param.default is not inspect.Parameter.empty:
        print(f"{name} = {param.default!r}")
```

---

## Validation & Debug

Check configuration validity at startup:

```python
from metainformant import config

cfg = config.get("structural_variants")
print(f"Detection min_support: {cfg['detection']['min_support']}")
print(f"CNV significance: {cfg['detection']['cnv']['significance']}")
```

To debug what's being used inside a function:

```python
import os
from metainformant.structural_variants.detection.cnv import _ENV_PREFIX

# Show all SV_ env vars that would be read
for key, val in os.environ.items():
    if key.startswith("SV_"):
        print(f"{key} = {val}")
```

---

## What's Next?

- **[GETTING_STARTED.md](GETTING_STARTED.md)** — See these settings in action
- **[PERFORMANCE.md](PERFORMANCE.md)** — Tune these for large datasets
- **[TROUBLESHOOTING.md](TROUBLESHOOTING.md)** — Configuration-related errors
