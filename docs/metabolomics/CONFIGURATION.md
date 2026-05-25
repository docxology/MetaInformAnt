# Metabolomics Configuration

Module-specific configuration options, environment variables, and resource requirements for metabolomics analysis.

## Configuration Overview

The metabolomics module does not currently have a dedicated configuration file system like the RNA amalgkit module. Configuration is primarily **function-parameter-driven** with limited environment variable support through the generic METAINFORMANT config system.

---

## Environment Variables

Metabolomics respects the **generic METAINFORMANT environment variables** defined in `metainformant.core.utils.config`:

| Variable | Purpose | Default | Applies To |
|----------|---------|--------|------------|
| `METAINFORMANT_OUTPUT` | Base output directory for generated files | `output/` | All I/O functions when given relative paths |
| `META_LOG_LEVEL` | Logging level (`DEBUG`, `INFO`, `WARNING`, `ERROR`) | `INFO` | All module loggers |
| `META_CACHE_DIR` | Cache directory for downloaded databases | `~/.cache/metainformant` | Future: database caching |

**Example**:
```bash
export METAINFORMANT_OUTPUT="/data/metabolomics_results"
export META_LOG_LEVEL=DEBUG
python3 analyze.py
```

**Access in code**:
```python
import os
from metainformant.core.utils import config

output_base = config.get_output_base()  # reads METAINFORMANT_OUTPUT
log_level = config.get_log_level()      # reads META_LOG_LEVEL
```

---

### Future Module-Specific Variables (Planned)

These are not yet implemented but are planned for future releases:

| Variable | Purpose | Planned |
|----------|---------|---------|
| `METABOLOMICS_PPM_TOLERANCE` | Default ppm tolerance for `identify_metabolites()` | Q3 2026 |
| `METABOLOMICS_DEFAULT_DB` | Path to default metabolite database (HMDB, KEGG) | Q3 2026 |
| `METABOLOMICS_NORM_METHOD` | Default normalization method | Q4 2026 |
| `METABOLOMICS_MAX_MATCHES` | Maximum matches to return per m/z | Q4 2026 |

---

## Function Parameters vs. Configuration

All current settings are passed directly to functions rather than read from config files:

### Identification

```python
# All controlled by function arguments:
matches = analysis.identification.identify_metabolites(
    observed_mz=mz_array,
    database=my_db,
    ppm_tolerance=10.0,           # ← configurable
)
```

**Recommendation**: For reproducible pipelines, create a simple YAML/JSON config and load it:

```yaml
# config/metabolomics.yaml
metabolomics:
  identification:
    ppm_tolerance: 10.0
    max_matches: 10
  normalization:
    method: "total_ion_count"
  differential:
    alpha: 0.05
    fdr_method: "bh"
```

```python
import yaml
cfg = yaml.safe_load(open("config/metabolomics.yaml"))
matches = identify_metabolites(..., ppm_tolerance=cfg["metabolomics"]["identification"]["ppm_tolerance"])
```

---

## Resource Requirements

### Memory (RAM)

| Data Size | Estimate | Notes |
|-----------|----------|-------|
| 100 metabolites × 20 samples | < 10 MB | Negligible |
| 5000 metabolites × 200 samples | ~8 MB (float64: 8 bytes per value) | `5000 × 200 × 8 ≈ 8 MB` |
| 10000 metabolites × 1000 samples | ~80 MB | Still modest; NumPy arrays are memory-mapped friendly |

**Note**: Numbers are for intensity matrices only. Peak-picked m/z tables from upstream tools (MZmine, XCMS) are similarly sized.

### CPU

| Operation | Complexity | Typical Time (n=5000, s=100) | Parallelizable |
|-----------|------------|-------------------------------|---------------|
| `identify_metabolites` | O(n × m) | 50 ms (n=1000 peaks, m=5000 DB) | Yes — per-peak |
| `identify_with_adducts` | O(n × m × a) | 350 ms (a=7 adducts) | Yes — per-peak |
| `normalize_intensities` | O(n × s) | 5 ms | Yes — per-row |
| `differential_abundance` | O(n × s) | 10 ms | Yes — per-metabolite |
| `metabolite_set_enrichment` | O(q × p) | 20 ms (q=100, p=500) | Yes — per-pathway |

**Parallelization**: All per-metabolite or per-peak operations can be parallelized across CPU cores using `concurrent.futures` or `joblib`. Not yet built-in but easily scriptable.

---

## Database Configuration

### Metabolite Reference Database

**Format**: Plain dictionary mapping metabolite name → monoisotopic neutral mass (Da).

**Example** (`my_database.py`):
```python
HMDB_DB = {
    "Glucose": 180.0634,
    "Lactate": 89.0240,
    "Pyruvate": 87.0083,
    "Citrate": 191.0197,
    # ... thousands more
}
```

**Sources**:

| Database | Format | Access | Size |
|----------|--------|--------|------|
| **HMDB** | CSV/JSON | Free download | ~114,000 metabolites |
| **KEGG** | KEGG COMPOUND database | Free (academic) | ~20,000 compounds |
| **METLIN** | Proprietary | Commercial | ~1,000,000 metabolites |
| **LipidMaps** | TSV | Free | ~50,000 lipids |
| **Reactome** | SBML | Free | Pathway-to-metabolite mappings |

**Loading from file** (future feature):
```python
from metainformant.metabolomics.io import load_database

db = load_database("hmdb.csv")  # auto-detects format, parses m/z column
```

**Currently**: You must parse database files yourself and pass dict to functions.

---

### Pathway Database

**Format**: `dict[str, list[str]]` mapping pathway name → list of metabolite names (exact matches to database keys).

**Example**:
```python
KEGG_GLYCOLYSIS = [
    "Glucose", "Glucose-6-phosphate", "Fructose-6-phosphate",
    "Fructose-1,6-bisphosphate", "Glyceraldehyde-3-phosphate",
    "Pyruvate", "Lactate",
]
PATHWAY_DB = {
    "Glycolysis": KEGG_GLYCOLYSIS,
    "TCA Cycle": TCA_METABOLITES,
}
```

**Building pathway_db from KEGG** (custom script, not built-in):
```python
import requests
from metainformant.core import io

# Fetch KEGG pathway maps (requires account)
# Or download pre-built JSON from MetaboAnalyst
import json
pathway_db = json.load(open("kegg_pathways.json"))
```

**Recommendation**: Start with a small curated pathway set (20–50 pathways) for exploratory analysis; expand to full KEGG (~300 pathways) for comprehensive enrichment.

---

## Performance Configuration

### Tuning `ppm_tolerance`

**Trade-off**:
- **Smaller** (e.g., 5 ppm): Fewer false positives, but may miss true matches if mass error is high.
- **Larger** (e.g., 20–50 ppm): Higher recall, but more false matches.

**Guidelines**:

| Instrument | Resolving Power | Recommended ppm |
|------------|----------------|-----------------|
| Quadrupole | ~1000 | 20–50 ppm (low res) |
| TOF | 20,000–40,000 | 10–20 ppm |
| Orbitrap | 60,000–140,000 | 5–10 ppm |
| FT-ICR | >200,000 | <5 ppm |

**Calibration**: Run quality control (QC) pooled samples first to assess mass accuracy. Compute median ppm error of known internal standards; set tolerance to 2–3× that error.

---

### Normalization Method Selection

| Data Type | Recommended | Reason |
|-----------|-------------|--------|
| Untargeted LC-MS (positive) | `"total_ion_count"` → `"log2"` | TIC corrects for injection volume; log2 stabilizes variance |
| Targeted quantification (MRM) | `"none"` (raw) or `"median"` | Targeted assays are already quantitative |
| NMR | `"none"` or `"pqn"` (not implemented) | NMR inherently quantitative |
| Large cohort with batches | `"median"` + ComBat | Median less sensitive to batch outliers |

**Multi-step normalization** (common pattern):
```python
# Step 1: Sample-wise normalization (TIC)
X_norm = normalize_intensities(X, method="total_ion_count")

# Step 2: Feature-wise transformation (log2)
X_log = normalize_intensities(X_norm, method="log2")

# Step 3: Batch correction (external library)
from pycombat import Combat
X_corrected = Combat().fit_transform(X_log.T).T
```

---

## I/O Configuration

### CSV Format Recommendations

**Column delimiter**: Use comma (`,`). For European formats with semicolon, pass `delimiter=";"`.

**Header row**: First row must contain sample names. First column must contain metabolite identifiers.

**Missing values**: Represent as blank cells or `0`. Zeros are treated as missing during imputation; use `missing_value_imputation()` before statistics.

**Example well-formed CSV**:
```csv
metabolite,QC_1,QC_2,Sample_A,Sample_B,Sample_C
Glucose,1200,1150,980,1050,1100
Lactate,450,430,380,420,395
```

### MGF Format Recommendations

**Required fields**: At least `PEPMASS` (precursor m/z) and peak list (m/z intensity pairs). Optional: `RTINSECONDS`, `SCANS`, `CHARGE`.

**Encoding**: UTF-8 plain text.

**Expected structure**:
```
BEGIN IONS
PEPMASS=180.0634 100.0
RTINSECONDS=245.6
SCANS=1
104.0527 45.2
118.0651 12.8
180.0634 100.0
END IONS
```

---

## Advanced: Custom Adduct Definitions

While `identify_with_adducts()` has built-in adducts, you can provide custom ones:

```python
CUSTOM_ADDUCTS = {
    "[M+2H]2+": 1.007276 / 2,   # For doubly charged
    "[M+3H]3+": 1.007276 / 3,   # Triply charged
    "[M+Cl]+": 34.969402,
}
matches = identify_with_adducts(
    mz_array, db,
    adducts=CUSTOM_ADDUCTS,
    ppm_tolerance=5.0,
    ion_mode="positive",
)
```

**Caveat**: For multiply charged ions, the mass calculation is `m/z = (M + n×H) / n`, so the mass delta per charge is `(1.007276 / n)`, not 1.007.

---

## Troubleshooting Configuration Issues

### "Config file not found"

**Problem**: Want to load external database but file missing.

**Fix**:
```bash
# Check path
ls -lh /path/to/database.csv

# Use absolute path
db = load_database(Path("/absolute/path/to/hmdb.csv"))
```

### "Invalid ppm tolerance"

**Problem**:
```python
identify_metabolites(..., ppm_tolerance=0)
ValueError: ppm_tolerance must be > 0
```

**Fix**: Use a reasonable tolerance (5–50 ppm depending on instrument).

---

## Performance Tuning Parameters

No dedicated "performance config" — optimize via algorithmic choices:

| Lever | Effect | How to adjust |
|-------|--------|---------------|
| **Database size** | O(n × m) time | Filter DB to relevant m/z range or metabolite classes |
| **ppm tolerance** | Slight filter | Use smallest tolerance that still captures true matches (calibrate with standards) |
| **Adduct set** | O(a) multiplier | Limit to 2–3 most common adducts per mode if certain of ionization |
| **Parallelization** | Speedup ∝ cores | Wrap per-metabolite/per-peak loops in `concurrent.futures.ProcessPoolExecutor` |

Example parallel identification:
```python
from concurrent.futures import ProcessPoolExecutor

def identify_single_peak(mz):
    matches = identify_metabolites(np.array([mz]), db, ppm_tolerance=10)
    return matches[0]

with ProcessPoolExecutor(max_workers=8) as ex:
    all_matches = list(ex.map(identify_single_peak, observed_mz))
```

---

## Module Versions and Compatibility

| Component | Version | Python |
|-----------|---------|--------|
| metainformant core | 0.1.0+ | 3.11+ |
| NumPy | ≥1.24 | — |
| SciPy | ≥1.10 (optional) | — |

**Check versions**:
```bash
python3 -c "import metainformant, numpy, scipy; print(metainformant.__version__)"
```

---

## Summary

The metabolomics module currently uses **parameter-driven configuration** rather than config files. Key knobs:

- `ppm_tolerance`: Identification sensitivity (default 10)
- `method`: Normalization method (TIC / median / log2 / pareto)
- `adducts`: Custom adduct dictionary or default built-in

Environment variables only affect output paths and logging through the core system. Future versions may add dedicated config files for metabolite databases and pathway collections.

For reproducible research:
1. Pin all package versions in `requirements.txt` (`uv pip freeze > requirements.txt`)
2. Store analysis parameters in a YAML/JSON alongside your script
3. Log session info:
   ```python
   import metainformant, sys
   print(f"METAINFORMANT {metainformant.__version__}, Python {sys.version}")
   ```
