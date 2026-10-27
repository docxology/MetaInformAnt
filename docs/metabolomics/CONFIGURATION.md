# Metabolomics Configuration

Module-specific configuration options, environment variables, and resource requirements for metabolomics analysis.

## Configuration Overview

The metabolomics module has no configuration system: no configuration file, no environment variables and no module-level settings, unlike the RNA amalgkit module. All configuration is **function-parameter-driven**: every setting is an argument of the function you call. The package imports only NumPy (and SciPy, optionally) and never calls `metainformant.core`, so the generic core configuration helpers are not involved either.

---

## Environment Variables

**The metabolomics module reads no environment variables.** Nothing under `src/metainformant/metabolomics` calls `os.environ`, `os.getenv`, `apply_env_overrides`, `load_typed_env`, `get_env_or_default` or a config-file loader. The following names are **not implemented** anywhere in the source (neither in the metabolomics package nor in `metainformant.core`), so setting them has no effect:

| Variable | Status |
|----------|--------|
| `METAINFORMANT_OUTPUT` | Not implemented. Output locations are whatever paths you pass to `write_csv`, `write_mgf`, `fig.savefig`, ... |
| `META_LOG_LEVEL` | Not implemented. Core logging reads `CORE_LOG_LEVEL` (via `configure_logging_from_env()`), but metabolomics functions emit no log records |
| `META_CACHE_DIR` | Not implemented. Nothing downloads or caches databases |
| `METABOLOMICS_PPM_TOLERANCE`, `METABOLOMICS_DEFAULT_DB`, `METABOLOMICS_NORM_METHOD`, `METABOLOMICS_MAX_MATCHES` | Not implemented |

**Environment-driven defaults in your own scripts**: the core helpers can read variables for you, but you must pass the values to the metabolomics functions yourself:

```python
from metainformant.core.utils.config import load_typed_env
from metainformant.metabolomics.analysis.identification import identify_metabolites

# Your own convention: `export METABOLOMICS_PPM_TOLERANCE=5` before running the script
env = load_typed_env(prefix="METABOLOMICS", keys={"PPM_TOLERANCE": float})
ppm = env.get("PPM_TOLERANCE", 10.0)  # variables that are not set are simply absent

matches = identify_metabolites(mz_array, my_db, ppm_tolerance=ppm)
```

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

**Recommendation**: For reproducible pipelines, keep the settings in a YAML/JSON file of your own and pass them to the functions. The module never reads this file; the keys below are only a convention that maps one-to-one onto real function arguments:

```yaml
# config/metabolomics.yaml (read by your script, not by the module)
metabolomics:
  identification:
    ppm_tolerance: 10.0
  normalization:
    method: "total_ion_count"
  enrichment:
    fdr_threshold: 0.05
```

```python
from metainformant.core.utils.config import load_mapping_from_file
from metainformant.metabolomics.analysis.identification import identify_metabolites, normalize_intensities
from metainformant.metabolomics.pathways.enrichment import enrichment_with_fdr

# YAML, TOML or JSON, chosen by the file suffix
cfg = load_mapping_from_file("config/metabolomics.yaml")["metabolomics"]

matches = identify_metabolites(mz_array, my_db, ppm_tolerance=cfg["identification"]["ppm_tolerance"])
normalized = normalize_intensities(intensities, method=cfg["normalization"]["method"])
enriched = enrichment_with_fdr(query_metabolites, pathway_db, fdr_threshold=cfg["enrichment"]["fdr_threshold"])
```

---

## Resource Requirements

### Memory (RAM)

| Data Size | Estimate | Notes |
|-----------|----------|-------|
| 100 metabolites × 20 samples | < 10 MB | Negligible |
| 5000 metabolites × 200 samples | ~8 MB (float64: 8 bytes per value) | `5000 × 200 × 8 ≈ 8 MB` |
| 10000 metabolites × 1000 samples | ~80 MB | Still modest once loaded |

**Note**: Numbers are for intensity matrices only. Peak-picked m/z tables from upstream tools (MZmine, XCMS) are similarly sized. `io.formats.read_csv` reads the whole file as text and builds Python lists before creating the array, so peak memory while loading is several times the matrix size (about 7× in a measurement: 59 MB peak for an 8 MB, 2000 × 500 matrix).

### CPU

| Operation | Complexity | Measured time | Parallelizable |
|-----------|------------|---------------|---------------|
| `identify_metabolites` | O(n × m) | 15–60 ms (n=1000 peaks, m=5000 DB) | Yes — per-peak |
| `identify_with_adducts` | O(n × m × a) | ~150 ms (a=7 adducts, same sizes) | Yes — per-peak |
| `normalize_intensities` | O(n × s) | 1–10 ms (5000 metabolites × 100 samples, by method) | Yes — per-sample (TIC, median) or per-row (pareto) |
| `differential_abundance` | O(n × s) | 2–10 ms (5000 metabolites × 100 samples) | Yes — per-metabolite |
| `metabolite_set_enrichment` | O(q × p) | 60–95 ms (q=100 query metabolites, p=500 pathways of 20 metabolites) | Yes — per-pathway |

Times are best-to-median wall-clock values from one arm64 machine (Python 3.12, NumPy 2.5); expect different absolute numbers elsewhere.

**Parallelization**: All per-metabolite or per-peak operations can be parallelized across CPU cores using `concurrent.futures` or `joblib`. Not yet built-in but easily scriptable.

---

## Database Configuration

### Metabolite Reference Database

**Format**: Plain `dict[str, float]` mapping metabolite name → mass. What the value means depends on the function: `identify_with_adducts` subtracts an adduct offset from each observed m/z first, so store the **monoisotopic neutral mass** (Da); `identify_metabolites` compares the observed m/z with the stored value directly (no adduct handling), so store the **exact m/z of the ion you expect to observe**.

**Example** (`my_database.py`, neutral monoisotopic masses for `identify_with_adducts`):
```python
HMDB_DB = {
    "Glucose": 180.0634,
    "Lactate": 90.0317,
    "Pyruvate": 88.0160,
    "Citrate": 192.0270,
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

**Loading from file**: There is no built-in database loader. Parse database files yourself (for example with `csv` or `pandas`) and pass the resulting dict to the functions.

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
TCA_METABOLITES = ["Citrate", "Isocitrate", "Alpha-ketoglutarate", "Succinate", "Fumarate", "Malate"]
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

`normalize_intensities` accepts exactly `"total_ion_count"`, `"median"`, `"log2"` and `"pareto"`; any other value (including `"none"`) raises `ValueError`. To leave the data unnormalized, skip the call.

| Data Type | Recommended | Reason |
|-----------|-------------|--------|
| Untargeted LC-MS (positive) | `"total_ion_count"` → `"log2"` | TIC corrects for injection volume; log2 stabilizes variance |
| Targeted quantification (MRM) | No normalization (skip the call) or `"median"` | Targeted assays are already quantitative |
| NMR | No normalization (skip the call); PQN is not implemented | NMR inherently quantitative |
| Large cohort with batches | `"median"` + batch correction (`correct_batch_combat`) | Median less sensitive to batch outliers |

**Multi-step normalization** (common pattern):
```python-snippet
from metainformant.quality.batch.detection import correct_batch_combat

# Step 1: Sample-wise normalization (TIC)
X_norm = normalize_intensities(X, method="total_ion_count")

# Step 2: Element-wise transformation (log2)
X_log = normalize_intensities(X_norm, method="log2")

# Step 3: Batch correction (simplified ComBat-style, in metainformant.quality;
# takes samples × features, so transpose in and out). batch_labels: one label per sample.
X_corrected = correct_batch_combat(X_log.T, batch_labels).T
```

---

## I/O Configuration

### CSV Format Recommendations

**Column delimiter**: Use comma (`,`). For European formats with semicolon, pass `delimiter=";"`; decimal commas are not understood, so a cell such as `1,5` is read as `0.0` (see Missing values below).

**Header row**: First row must contain sample names. First column must contain metabolite identifiers. (`read_csv` takes the sample names from the header cells starting at `data_start_col` and the identifiers from column `metabolite_col`; both default to the layout above.)

**Missing values**: Represent as blank cells or `0`. `read_csv` silently reads a cell that is not a number (blank, `NA`, text) as `0.0`, while cells that Python's `float()` accepts, including `nan` and `inf`, are kept as they are. Zeros and NaNs are treated as missing during imputation; use `missing_value_imputation()` before statistics.

**Example well-formed CSV**:
```csv
metabolite,QC_1,QC_2,Sample_A,Sample_B,Sample_C
Glucose,1200,1150,980,1050,1100
Lactate,450,430,380,420,395
```

### MGF Format Recommendations

**Fields read**: `PEPMASS` (precursor m/z; the first number is used, and the spectrum's `precursor_mz` is `None` if the line is missing), `RTINSECONDS`, `SCANS` and the peak list (one `m/z intensity` pair per line). Every other header line (`CHARGE`, `TITLE`, ...) is ignored, and each block must be closed with `END IONS`. All spectra are returned with `ms_level=2`.

**Encoding**: Plain text, read with Python's default text encoding (UTF-8 on most systems).

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

`identify_with_adducts()` ships with default adduct sets (positive: `[M+H]+`, `[M+Na]+`, `[M+K]+`, `[M+NH4]+`; negative: `[M-H]-`, `[M+Cl]-`, `[M+FA-H]-`, selected by `ion_mode`). To use your own, pass a dict of adduct name → m/z offset:

```python
CUSTOM_ADDUCTS = {
    "[M+H]+": 1.007276,
    "[M+Li]+": 7.015454,         # lithium adduct
    "[M-H2O+H]+": -17.003289,    # in-source water loss
}
matches = identify_with_adducts(
    mz_array, db,
    adducts=CUSTOM_ADDUCTS,
    ppm_tolerance=5.0,
)
```

**Caveats**:
- A custom `adducts` dict replaces the defaults and is used exactly as given; `ion_mode` only selects the default set when `adducts` is `None`.
- Each value is a single-charge m/z offset: the neutral mass is computed as `observed m/z - offset`. Charge states are not modelled, so multiply charged ions such as `[M+2H]2+` (`m/z = (M + 2×1.007276) / 2`) are not supported and would match the wrong masses.

---

## Troubleshooting Configuration Issues

### "File not found"

**Problem**: `read_csv` / `read_mgf` (or your own database loader) raises `FileNotFoundError` because the path is wrong. The module reads no config files and has no database loader, so the only files it opens are the ones you pass in.

**Fix**:
```bash
# Check path
ls -lh /path/to/database.csv

# Use an absolute path in the code that reads the file
```

### "No matches returned" / invalid ppm tolerance

**Problem**: `ppm_tolerance` is not validated, so a bad value does not raise: a negative tolerance can never match, and `ppm_tolerance=0` only matches bit-for-bit identical masses (with score 0.0 and a NumPy `RuntimeWarning`). A tolerance that is too tight for your mass accuracy behaves the same way and returns empty match lists.

```python-snippet
identify_metabolites(mz_array, db, ppm_tolerance=-1)
# -> [[], [], ...]  (one empty list per observed m/z, no exception)
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
| **Parallelization** | Process start-up and per-task overhead dominate small jobs | Give each `concurrent.futures.ProcessPoolExecutor` task a chunk of many peaks (one task per peak is far slower than a single call); worth trying only for very large peak lists |

Example parallel identification (saved as a script):
```python
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from metainformant.metabolomics.analysis.identification import identify_metabolites

db = {"Glucose": 180.0634, "Lactate": 89.0240}  # your {name: m/z} dict; module level, so worker processes can import it

def identify_chunk(mz_chunk):
    return identify_metabolites(mz_chunk, db, ppm_tolerance=10)

if __name__ == "__main__":  # required where processes are spawned (macOS, Windows)
    observed_mz = np.array([180.0634, 89.0240, 150.0])  # in practice: a very long peak list
    chunks = np.array_split(observed_mz, 8)  # one task per chunk of peaks, not per peak
    with ProcessPoolExecutor(max_workers=8) as ex:
        all_matches = [matches for part in ex.map(identify_chunk, chunks) for matches in part]
```

---

## Module Versions and Compatibility

| Component | Version | Python |
|-----------|---------|--------|
| metainformant | 0.4.0 (`metainformant.__version__`) | 3.11+ |
| NumPy | ≥1.24 declared; ≥2.0 for `detect_chromatographic_peaks` (`np.trapezoid`) | — |
| SciPy | ≥1.11 declared; optional for the metabolomics module itself (normal-approximation p-values without it) | — |

**Check versions**:
```bash
python3 -c "import metainformant, numpy, scipy; print(metainformant.__version__, numpy.__version__, scipy.__version__)"
```

---

## Summary

The metabolomics module currently uses **parameter-driven configuration** rather than config files. Key knobs:

- `ppm_tolerance`: Identification sensitivity (default 10)
- `method`: Normalization method (TIC / median / log2 / pareto)
- `adducts`: Custom adduct dictionary or default built-in

The module reads no environment variables and no configuration files; there are no dedicated config files for metabolite databases or pathway collections either (both are plain dicts you load yourself).

For reproducible research:
1. Pin all package versions in `requirements.txt` (`uv pip freeze > requirements.txt`)
2. Store analysis parameters in a YAML/JSON alongside your script
3. Log session info:
   ```python
   import metainformant, sys
   print(f"METAINFORMANT {metainformant.__version__}, Python {sys.version}")
   ```
