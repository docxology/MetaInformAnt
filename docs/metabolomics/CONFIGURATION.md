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

**Missing values.** `read_csv` returns `NaN` for a missing cell. `"total_ion_count"` and `"median"` rescale each sample by a statistic of that sample, so one `NaN` would distort the scale of its sample and, through the shared median, of every other sample; they raise `ValueError` when the matrix contains `NaN`. Impute first with `missing_value_imputation` (it treats `NaN` and zero as missing). `"log2"` keeps a `NaN` where it is. `"pareto"` does not check: a `NaN` makes the whole row of its metabolite `NaN`, so impute before it too.

| Data Type | Recommended | Reason |
|-----------|-------------|--------|
| Untargeted LC-MS (positive) | `"total_ion_count"` → `"log2"` | TIC corrects for injection volume; log2 stabilizes variance |
| Targeted quantification (MRM) | No normalization (skip the call) or `"median"` | Targeted assays are already quantitative |
| NMR | No normalization (skip the call); PQN is not implemented | NMR inherently quantitative |
| Large cohort with batches | `"median"` + batch correction (`correct_batch_combat`, parametric empirical-Bayes ComBat) | Median less sensitive to batch outliers |

**Multi-step normalization** (common pattern):
```python-snippet
from metainformant.quality.batch.detection import correct_batch_combat

# Step 0: Impute missing cells (read_csv returns NaN; TIC/median normalization and ComBat refuse NaN)
X = missing_value_imputation(X, method="min_half")

# Step 1: Sample-wise normalization (TIC)
X_norm = normalize_intensities(X, method="total_ion_count")

# Step 2: Element-wise transformation (log2)
X_log = normalize_intensities(X_norm, method="log2")

# Step 3: Batch correction with parametric empirical-Bayes ComBat (Johnson, Li & Rabinovic 2007,
# no covariates), in metainformant.quality. It takes samples × features, so transpose in and out.
# batch_labels: one label per sample.
X_corrected = correct_batch_combat(X_log.T, batch_labels).T
```

`correct_batch_combat` standardizes the data once, estimates each batch's location and scale effect per feature, shrinks those estimates toward each other with empirical-Bayes priors fitted across features, and adjusts the data. The result does not depend on the order or the labels of the batches or on the order of the samples. Requirements and behaviour:

- the data must be finite (no `NaN` or `inf`) and every batch needs at least 2 samples, otherwise `ValueError`;
- a feature that is constant within some batch cannot have its variance estimated and is returned unchanged; if exactly one feature varies within every batch the priors cannot be fitted and `ValueError` is raised;
- a single batch comes back as an unchanged copy;
- there are no covariates, so biological groups that coincide with batches are removed along with them.

---

## I/O Configuration

### CSV Format Recommendations

**Column delimiter**: Use comma (`,`). For files with a semicolon delimiter, pass `delimiter=";"`. Decimal commas are not understood and are reported rather than repaired: in a comma-delimited file `1,5` splits one number into two cells, so the row has more cells than the header and `read_csv` raises `ValueError`; in a semicolon-delimited file the cell `1,5` raises `ValueError` ("a decimal comma is not supported; use a decimal point"). Convert the file to decimal points before reading.

**Header row**: First row must contain sample names. First column must contain metabolite identifiers. (`read_csv` takes the sample names from the header cells starting at `data_start_col` and the identifiers from column `metabolite_col`; both default to the layout above.)

**Missing values**: Represent them as blank cells or `NA`. A blank cell or one of `NA`, `N/A`, `#N/A`, `NaN`, `null` or `None` (any case) is returned as `NaN`, never as `0.0`; a true zero stays `0.0`. Any other cell that is not a number (text, a decimal comma) raises `ValueError` naming the line, the metabolite and the sample, and a row with a different number of cells than the header raises `ValueError` naming the line. Zeros and `NaN`s are both treated as missing by `missing_value_imputation()`; run it before `normalize_intensities` with `"total_ion_count"` or `"median"`, which refuse `NaN`, and before statistics that use variance.

**Annotation columns**: columns before `data_start_col` other than `metabolite_col` (for example `mz` and `rt` of an XCMS or MZmine table) are per-metabolite annotations. `read_csv` returns them in `dataset.metadata` as `{metabolite: {column: value}}` (numbers as floats, missing-value tokens as `NaN`, anything else as text). A repeated metabolite name or annotation column name raises `ValueError` when annotation columns are present, because the metadata would be ambiguous. `write_csv` writes `NaN` as `nan` (which `read_csv` reads back as `NaN`) and does not write `metadata`.

**Example well-formed CSV**:
```csv
metabolite,QC_1,QC_2,Sample_A,Sample_B,Sample_C
Glucose,1200,1150,980,1050,1100
Lactate,450,430,380,420,395
```

**Example with annotation columns and a missing value** (`read_csv(path, metabolite_col=0, data_start_col=3)`):
```csv
name,mz,rt,QC_1,QC_2,Sample_A
Glucose,180.0634,245.6,1200,1150,980
Lactate,89.0240,150.2,450,NA,380
```
This gives `dataset.samples == ["QC_1", "QC_2", "Sample_A"]`, `dataset.metadata["Glucose"] == {"mz": 180.0634, "rt": 245.6}` and `dataset.intensities[1, 1]` equal to `NaN`.

### MGF Format Recommendations

**Fields read**: `PEPMASS` (precursor m/z; the first number is used, and the spectrum's `precursor_mz` is `None` if the line is missing), `RTINSECONDS`, `SCANS` (an integer) and the peak list. Every other header line (`CHARGE`, `TITLE`, ...) is ignored. All spectra are returned with `ms_level=2`.

**Line rules inside a `BEGIN IONS` / `END IONS` block**: blank lines and comment lines (starting with `#`, `;`, `!` or `/`) are skipped; a line containing `=` is a header; every other line is a peak and must start with two finite numbers, `m/z intensity` (further columns are ignored). A peak line is never skipped silently, so `.5 100` and `+200.0 20` are read and a malformed one is an error. An empty `PEPMASS=`, a `PEPMASS`, `RTINSECONDS` or `SCANS` value that is not a number, and a malformed peak line raise `ValueError` naming the file and the spectrum (counted from 1 in file order). Each block must be closed with `END IONS`: a `BEGIN IONS` block that is never closed (for example the last spectrum of a truncated file) is dropped without an error, so check the number of spectra returned.

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

`identify_with_adducts()` ships with default adduct sets (positive: `[M+H]+`, `[M+Na]+`, `[M+K]+`, `[M+NH4]+`; negative: `[M-H]-`, `[M+Cl]-`, `[M+FA-H]-`, selected by `ion_mode`, which is `"positive"` or `"negative"` in any case). To use your own, pass a dict of adduct name → m/z offset:

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
- A custom `adducts` dict replaces the defaults and is used exactly as given; `ion_mode` only selects the default set when `adducts` is `None`, but it is validated either way: any value other than `"positive"` or `"negative"` (in any case) raises `ValueError`, so a typo cannot go unnoticed.
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

**Problem**: `ppm_tolerance` must be a positive, finite number of ppm. Zero, a negative value, `NaN`, infinity or a non-number raises `ValueError` in both `identify_metabolites` and `identify_with_adducts`:

```python-snippet
identify_metabolites(mz_array, db, ppm_tolerance=-1)
# ValueError: ppm_tolerance must be a positive finite number, got -1
```

A tolerance that is valid but too tight for your mass accuracy returns empty match lists (one empty list per observed m/z, no exception).

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
| NumPy | ≥1.24 declared; every function runs on any supported version (`detect_chromatographic_peaks` integrates peaks with its own trapezoidal rule, because neither `np.trapz` nor `np.trapezoid` exists on all supported versions) | — |
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
