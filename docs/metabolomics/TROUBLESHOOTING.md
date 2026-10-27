# Metabolomics Troubleshooting

Common issues, error messages, and solutions for metabolomics analysis with METAINFORMANT.

## Installation Issues

### "No module named 'metainformant.metabolomics'"

**Symptom**:
```pycon
>>> from metainformant.metabolomics import analysis
ModuleNotFoundError: No module named 'metainformant.metabolomics'
```

**Cause**: Metabolomics module not installed or virtual environment not activated.

**Fix**:
```bash
# Ensure you're in the MetaInformAnt repo root
cd /path/to/MetaInformAnt

# Activate virtual environment
source .venv/bin/activate

# Reinstall in editable mode
uv pip install -e .

# Verify
python3 -c "from metainformant.metabolomics import analysis; print('OK')"
```

---

### "ImportError: numpy is required"

**Symptom**:
```
ImportError: numpy is required for metabolomics analysis
```

**Cause**: NumPy not installed in the active environment.

**Fix**:
```bash
# Activate venv first
source .venv/bin/activate

# Install numpy and scipy
uv pip install numpy scipy

# Verify
python3 -c "import numpy; print(numpy.__version__)"
```

---

### "AttributeError: module 'metainformant.metabolomics.analysis' has no attribute 'identification'"

**Cause**: Circular import or incomplete installation.

**Fix**:
```bash
# Reinstall cleanly
uv pip uninstall metainformant -y
uv pip install -e .
```

---

## Data Loading Errors

### "ValueError: ... cannot read 'abc' as a number" from `read_csv()`

**Symptom**:
```python
dataset = io.formats.read_csv("data.csv")
# ValueError: data.csv: line 2, metabolite 'Glucose', sample 'S1': cannot read 'abc' as a number;
#             write a missing value as a blank cell or NA
```

**Cause**: a cell is neither a number nor a missing-value token. `read_csv` reports it instead of replacing it with a value. Accepted missing-value tokens (any case) are a blank cell, `NA`, `N/A`, `#N/A`, `NaN`, `null` and `None`; they are returned as `NaN`.

**Fix**:
- Write missing measurements as a blank cell or `NA`
- A decimal comma (`1,5`) is not understood: the message says so. Convert to a decimal point
- Check the delimiter (use `delimiter="\t"` for TSV, `delimiter=";"` for semicolon files)

### "ValueError: ... line 2 has 4 cells but the header has 3" from `read_csv()`

**Cause**: a row has a different number of cells than the header. A decimal comma in a comma-delimited file is a common reason, because it splits one number into two cells; so are a comma inside a metabolite name and a missing trailing column. `read_csv` splits each line on the delimiter and does not interpret quotes.

**Fix**: make every row as wide as the header; remove the stray delimiter from the cell, or write the file with a delimiter (such as a tab) that the cells do not contain and pass it as `delimiter=`.

### "ValueError: ... do not fit the 3 header columns", "duplicate annotation column names" or "occurs more than once" from `read_csv()`

**Cause**: `metabolite_col` or `data_start_col` points outside the header, or annotation columns (those before `data_start_col`) are present and two of them have the same name, or a metabolite name occurs twice, so `dataset.metadata` would be ambiguous.

**Fix**: pass indices that match the header, and make the annotation column names and metabolite names unique.

### An empty dataset from `read_csv()`

An empty file (or one with no non-blank line) returns a dataset with no metabolites and no samples; it does not raise. Check `len(dataset.metabolites)` after loading.

### "ValueError: ... PEPMASS has no value" or "malformed peak line" from `read_mgf()`

**Cause**: inside a `BEGIN IONS` block an empty `PEPMASS=`, a `PEPMASS`, `RTINSECONDS` or `SCANS` value that is not a number (`SCANS` must be an integer), or a line that is not a header and does not start with two finite numbers `m/z intensity`. The message names the file and the spectrum (counted from 1).

**Fix**: correct the named spectrum. Comment lines (starting with `#`, `;`, `!` or `/`) and blank lines are skipped, and peak lines such as `.5 99` and `+200.0 20` are read.

### `read_mgf()` returns fewer spectra than the file contains

**Cause**: a `BEGIN IONS` block that is never closed by `END IONS` (for example the last spectrum of a truncated file) is dropped without an error.

**Fix**: check that every block ends with `END IONS`.

**Example correct CSV**:
```csv
metabolite,Sample1,Sample2,Sample3
Glucose,1000,1200,950
Lactate,450,380,520
```

---

### "FileNotFoundError" when reading CSV/MGF

**Symptom**:
```
FileNotFoundError: [Errno 2] No such file or directory: 'data.csv'
```

**Cause**: Wrong working directory or relative path.

**Fix**:
```python
from pathlib import Path

# Use absolute paths
data_path = Path("/full/path/to/data.csv")
dataset = io.formats.read_csv(data_path)

# Or change to data directory first
import os
os.chdir("/path/to/data")
dataset = io.formats.read_csv("data.csv")
```

---

## Identification Errors

### "ValueError: ppm_tolerance must be a positive finite number" or "ion_mode must be 'positive' or 'negative'"

**Cause**: `identify_metabolites` and `identify_with_adducts` reject a `ppm_tolerance` that is zero, negative, `NaN`, infinite or not a number, and `identify_with_adducts` rejects an `ion_mode` other than `"positive"` or `"negative"` (any case), even when you pass your own `adducts`.

**Fix**: pass a tolerance such as 5–50 ppm and `ion_mode="positive"` or `"negative"`.

---

### "No metabolites identified — all match lists empty"

**Symptom**:
```python
matches = identify_metabolites(observed_mz, db, ppm_tolerance=10)
# All inner lists are empty []
```

**Causes & Diagnosis**:

1. **ppm_tolerance too tight**:
   ```python
   # Check actual ppm errors
   for mz in observed_mz[:5]:
       for name, db_mz in db.items():
           ppm = abs(db_mz - mz) / mz * 1e6
           if ppm < 50:  # see what error range you're getting
               print(f"{name}: {ppm:.1f} ppm")
   # If all >10, increase tolerance: ppm_tolerance=20 or 50
   ```

2. **Database units mismatch**: Database might be in Da but data is m/z (same for z=1) or database is neutral mass while observed includes adduct.
   ```python
   # If positive mode, observed m/z = neutral + 1.007 ([M+H]+)
   observed_neutral = observed_mz - 1.007276
   matches = identify_metabolites(observed_neutral, db)
   ```

3. **Wrong mass type**: Database contains exact monoisotopic masses but observed are average masses (or vice versa). Ensure both are monoisotopic.

4. **Mass range outside database**: Observed m/z might be outside database range (e.g., database only has 100–500 Da, data has 50–1000 Da). Filter or expand database.

**Solution**: Use `identify_with_adducts()` for ESI data to automatically account for common adducts.

---

### "Too many matches per peak — low specificity"

**Symptom**: Every observed m/z matches dozens of database entries, even with 10 ppm tolerance.

**Cause**: Database is too dense (many metabolites with very similar masses) or ppm tolerance too loose.

**Diagnosis**:
```python
# Count matches per peak with current tolerance
for i, mz in enumerate(observed_mz[:10]):
    matches = [name for name, db_mz in db.items()
               if abs(db_mz - mz) / mz * 1e6 <= 10]
    print(f"m/z {mz:.4f}: {len(matches)} matches")
```

**Solutions**:
1. **Increase resolution**: Use a smaller ppm tolerance (5 ppm for high-res data like Orbitrap).
2. **Restrict database**: Filter to metabolites plausible for your sample type (e.g., exclude lipids for polar extracts).
3. **Add adduct context**: Use `identify_with_adducts()` which computes neutral mass first, reducing candidate pool.
4. **Add retention time**: If RT data available, filter by RT window (future: RT prediction integration).

---

### "All identified metabolites are false positives"

**Symptom**: Identified compounds make no biological sense (e.g., plant metabolites in bacterial sample).

**Cause**: Database contamination or incorrect database choice.

**Fix**:
- Use a curated, relevant database (HMDB for human, KEGG for general, custom for specific taxa)
- Verify database format: ensure m/z values are monoisotopic neutral masses (for `identify_metabolites`) or neutral masses (for `identify_with_adducts`)
- Add RT filters manually if RT data available

---

## Normalization Errors

### "ValueError: normalize_intensities(method='total_ion_count') needs complete data"

**Symptom**:
```python
normalize_intensities(dataset.intensities, method="total_ion_count")
# ValueError: normalize_intensities(method='total_ion_count') needs complete data: a NaN would silently
#             distort the scale factor of its sample. Impute missing values first with
#             missing_value_imputation, or use method='log2', which keeps a NaN in place
```

**Cause**: the matrix contains `NaN`, for example a cell that `read_csv` returned from a blank cell or `NA`. `"total_ion_count"` and `"median"` rescale each sample by a statistic of that sample, so one `NaN` would distort the scale of its sample and, through the shared median, of every other sample.

**Fix**:
```python
imputed = missing_value_imputation(dataset.intensities, method="min_half")  # NaN and zeros are filled
normalized = normalize_intensities(imputed, method="total_ion_count")
```
`method="log2"` keeps a `NaN` where it is; `fold_change()` returns `NaN` for a metabolite with a missing value in either group.

---

### "Zero division warning in normalize_intensities"

**Symptom**:
```
RuntimeWarning: divide by zero encountered in true_divide
```

**Cause**: Sample has zero total ion count (all zeros), creating division by zero.

**Diagnosis**:
```python
# Check for empty samples
zero_samples = np.where(dataset.intensities.sum(axis=0) == 0)[0]
print(f"Empty samples: {[dataset.samples[i] for i in zero_samples]}")
```

**Fix**:
1. **Remove empty samples**: Filter out samples with no signal before normalization
   ```python
   nonzero = dataset.intensities.sum(axis=0) > 0
   intensities = dataset.intensities[:, nonzero]
   samples = [s for i, s in enumerate(dataset.samples) if nonzero[i]]
   ```
2. **Handle in normalization**: The function already handles zero totals internally (`totals = np.where(totals > 0, totals, 1.0)`), so warnings are harmless but indicate data quality issues.

---

### "Normalized values are negative or NaN"

**Symptom**: After `normalize_intensities(method="pareto")`, some values are NaN.

**Cause**: Row has zero standard deviation (all values identical).

**Fix**: The code handles `std=0` by replacing with 1.0, but NaN might propagate if row also has NaN. Use `missing_value_imputation()` first:
```python
imputed = missing_value_imputation(intensities, method="median")
normalized = normalize_intensities(imputed, method="pareto")
```

---

## Statistical Testing Errors

### "t-test returns NaN or Inf p-values"

**Symptom**:
```python
t, p = differential_abundance(X, group_a, group_b)
# p contains nan or inf
```

**Causes**:

1. **Zero variance within groups**: All values in a group are identical → zero variance → infinite SE → NaN t.
   ```python
   # Check variances
   var_a = X[:, group_a].var(axis=1)
   var_b = X[:, group_b].var(axis=1)
   zero_var = (var_a == 0) | (var_b == 0)
   print(f"Metabolites with zero variance in one group: {zero_var.sum()}")
   ```
   **Fix**: Filter out features with zero variance before testing.

2. **Degrees of freedom collapse**: If n=1 in a group, variance undefined.
   **Fix**: Ensure at least 2 samples per group (n ≥ 2).

3. **Extreme outliers**: One huge outlier inflates variance → huge SE → t ≈ 0 → p ≈ 1 (not NaN). NaN indicates denominator was exactly zero.

**Fix**:
```python
# Filter features with sufficient variance
var_threshold = 1e-6
keep = (X.var(axis=1) > var_threshold)
X_filtered = X[keep, :]
# Then run differential_abundance on X_filtered
```

---

### "All p-values are 1.0 (nothing significant)"

**Symptom**: Differential analysis finds no significant features, even with large effect sizes.

**Possible causes**:

1. **Insufficient sample size**: n=2 per group gives low statistical power.
   **Solution**: Increase n or use effect size filtering (|logFC| > 1) without p-value.

2. **High variance within groups**: Biological variation too large relative to group differences.
   **Solution**: Check variance; consider `limma` (empirical Bayes shrinkage) via R or `rpy2`.

3. **Data not normalized**: Raw intensities have batch effects.
   **Solution**: Normalize first (TIC, then log2).

4. **Multiple testing too strict**: 5000 metabolites with BH-FDR is harsh.
   **Solution**: Use less stringent threshold (q < 0.1) or pre-filter to likely metabolites only.

---

### "Welch's approximation gave negative variance"

**Symptom**:
```
RuntimeWarning: invalid value encountered in sqrt
```

**Cause**: Numerical instability when group variances are tiny and n is small; denominator in df formula rounds to negative.

**Fix**: Already handled in code with `np.where(pooled_se > 0, pooled_se, 1e-10)`. If still seeing warnings, your data may have extreme precision issues. Ensure data is in appropriate units (not extremely small floats). Scale data (multiply by 1000) if needed.

---

## Pathway Enrichment Errors

### "Enrichment p-values are all 1.0 (nothing significant)"

**Symptom**: `metabolite_set_enrichment()` returns no significant pathways even with BH-FDR.

**Causes**:

1. **Query set too small**: With <5 metabolites, enrichment is underpowered.
   **Solution**: Use `pathway_activity_scoring()` instead, which uses continuous scores.

2. **Pathway database too large**: Pathways with hundreds of members are hard to enrich.
   **Solution**: Filter pathway_db to size 10–200 metabolites.

3. **Background size wrong**: If `background_size` vastly overestimates universe, p-values become too conservative.
   **Solution**: Let function auto-compute: `background_size=None` (default). A `background_size` smaller than the query or than a pathway raises `ValueError` (see below).

4. **Query metabolites not in pathway database**: Check overlap.
   ```python
   query_set = set(query_metabolites)
   for name, members in pathway_db.items():
       overlap = query_set & set(members)
       print(f"{name}: {len(overlap)}/{len(members)}")
   # If all overlaps are 0 or 1, enrichment underpowered
   ```

---

### "ValueError: inconsistent hypergeometric parameters"

**Symptom**:
```
ValueError: inconsistent hypergeometric parameters: population N=2, successes K=2, draws n=3
(need 0 <= K <= N and 0 <= n <= N). In metabolite_set_enrichment N is background_size, K the
pathway size and n the number of distinct query metabolites, so background_size must be at least
as large as the query and as every pathway
```

**Cause**: `background_size` is smaller than the number of distinct query metabolites or than a pathway, so it cannot be a universe that contains them.

**Fix**: pass a `background_size` that covers every metabolite you measured (at least as large as the query and the largest pathway), or leave it at `None` to use the union of the pathway members and the query.

**Also check** that `pathway_db` values are lists of metabolite names (strings), not counts; the function counts overlaps itself:
```python
pathway_db = {
    "Glycolysis": ["Glucose", "Pyruvate", "Lactate"],  # list of names, not a count such as 3
}
```

---

## Visualization Errors

### "Plot not showing in Jupyter notebook"

**Symptom**: `fig.show()` or `fig.write_html()` produces no output.

**Fix**: For Plotly figures in Jupyter:
```python-snippet
# In Jupyter, use:
fig.show(renderer="notebook")

# Or if using plotly.express:
import plotly.io as pio
pio.renderers.default = "notebook"
```

For static plots (matplotlib):
```python
import matplotlib.pyplot as plt
plt.show()  # Only works in IPython with %matplotlib inline
```

---

## Performance Issues

### "Identification is slow with large database (10k+ metabolites)"

**Symptom**: `identify_metabolites()` takes >5 seconds for 1000 peaks against 50k DB.

**Solutions**:

1. **Pre-filter database by m/z range**:
   ```python
   # Only keep DB entries near observed range (with buffer)
   obs_min, obs_max = observed_mz.min(), observed_mz.max()
   buffer = 50  # Da
   db_filtered = {k: v for k, v in db.items()
                  if obs_min - buffer <= v <= obs_max + buffer}
   ```

2. **Use integer m/z binning** (approx nearest neighbor):
   ```python
   from collections import defaultdict
   bins = defaultdict(list)
   bin_size = 1  # 1 Da bin
   for name, mz in db.items():
       bins[int(mz // bin_size)].append((name, mz))
   # Then only search same bin and adjacent bins
   ```

3. **Vectorize with NumPy** (future): Vectorized ppm calculation across all DB entries at once uses SIMD.

4. **Parallelize** (future): Split observed m/z array across processes.

---

### "Enrichment is slow with 1000 pathways"

**Symptom**: `metabolite_set_enrichment()` takes >1 second.

**Cause**: Loop over all pathways is O(n_query × n_pathways). For 500+ pathways, Python loops are slow.

**Solutions**:

1. **Pre-filter pathways**: Only test pathways that contain at least one query metabolite.
   ```python
   query_set = set(query_metabolites)
   pathway_db_filtered = {
       name: members for name, members in pathway_db.items()
       if query_set & set(members)  # at least one overlap
   }
   results = metabolite_set_enrichment(query, pathway_db_filtered)
   ```

2. **Batch similar queries**: If running many enrichment analyses (e.g., per time point), pre-compute pathway member sets as Python `set` objects for fast intersection.

3. **Parallelize by pathway** (advanced): Use `concurrent.futures` to distribute pathway tests across CPU cores. Each pathway test is independent:
   ```python-snippet
   from concurrent.futures import ThreadPoolExecutor
   with ThreadPoolExecutor() as ex:
       futures = [ex.submit(_test_one_pathway, ...) for ...]
   ```

---

## Configuration Errors

### "METAINFORMANT_OUTPUT not respected"

**Symptom**: Files saved to unexpected location despite setting `METAINFORMANT_OUTPUT`.

**Cause**: The metabolomics module reads no environment variables, so `METAINFORMANT_OUTPUT` has no effect. `write_csv` and `write_mgf` write exactly to the path you pass and do not create missing directories.

**Fix**: Within your code, construct the path explicitly:
```python
import os
from pathlib import Path

from metainformant.metabolomics.io.formats import write_csv

output_base = Path(os.getenv("METAINFORMANT_OUTPUT", "output"))  # your own convention
output_path = output_base / "metabolomics" / "results.csv"
output_path.parent.mkdir(parents=True, exist_ok=True)
write_csv(dataset, output_path)
```

---

## Common Pitfalls

### Pitfall 1: Using raw intensities for statistics

**Wrong**:
```python
# Don't use raw intensities for t-tests
t, p = differential_abundance(raw_intensities, ...)
```

**Why**: Raw MS intensities are not normally distributed, have heteroscedastic variance, and are prone to batch effects.

**Correct**:
```python
norm = normalize_intensities(raw, method="total_ion_count")
log_norm = normalize_intensities(raw, method="log2")  # often best for stats
t, p = differential_abundance(log_norm, ...)
```

### Pitfall 2: Ignoring multiple testing

**Wrong**:
```python
sig = pvals < 0.05  # 5000 tests → expect 250 false positives!
```

**Correct**:
```python
from metainformant.metabolomics.pathways.enrichment import benjamini_hochberg
qvals = benjamini_hochberg(pvals)
sig = np.array(qvals) < 0.05  # FDR-controlled
```

Or use `enrichment_with_fdr()` which handles this internally.

### Pitfall 3: Forgetting to impute missing values

**Wrong**:
```python
# read_csv returns NaN for blank/NA cells; zeros may also stand for "not detected"
normalize_intensities(intensities, method="total_ion_count")  # ValueError: needs complete data
var = intensities.var(axis=1)  # NaN, or zeros reducing the variance artificially
```

**Correct**:
```python
imputed = missing_value_imputation(intensities, method="min_half")  # fills NaN and zeros
normalized = normalize_intensities(imputed, method="total_ion_count")
var = imputed.var(axis=1)
```

Note: Imputation is needed before TIC/median normalization (which refuse `NaN`) and for statistics that use variance (PCA, clustering, differential analysis). For presence/absence tests (e.g., detection rate comparison), zeros have meaning, so count them before imputing.

---

## Diagnostic Checklist

When analysis produces unexpected results, run through this checklist:

- [ ] **Data loaded correctly?** `print(dataset.intensities.shape, dataset.metabolites[:5])`
- [ ] **No empty samples?** `assert (intensities.sum(axis=0) > 0).all()`
- [ ] **No NaN/Inf in matrix?** `assert np.isfinite(intensities).all()`
- [ ] **Reasonable m/z range?** `print(observed_mz.min(), observed_mz.max())`
- [ ] **Database units match?** Check if database m/z are neutral or [M+H]+
- [ ] **Using adduct-aware ID?** For ESI, prefer `identify_with_adducts()`
- [ ] **ppm tolerance appropriate?** 5–10 ppm for high-res; 10–20 for low-res
- [ ] **Normalized before stats?** Always normalize before differential tests
- [ ] **Multiple testing corrected?** Use BH-FDR
- [ ] **Sufficient replicates?** n ≥ 3 per group recommended; n=2 gives no df for variance

---

## Getting Help

If you encounter an issue not covered here:

1. **Check logs**: Enable debug logging:
   ```python
   import logging
   logging.basicConfig(level=logging.DEBUG)
   ```

2. **Minimal reproducible example**: Reduce your data to a tiny subset that triggers the bug.

3. **File an issue**: Include:
   - OS and Python version
   - METAINFORMANT commit hash (`git rev-parse HEAD`)
   - Full error traceback
   - Minimal code snippet that reproduces the bug
   - Sample data (or synthetic data that triggers issue)

4. **Alternative resources**:
   - [MetaboAnalyst documentation](https://www.metaboanalyst.ca/) for general metabolomics guidance
   - [HMDB 5.0 reference](https://pubmed.ncbi.nlm.nih.gov/34986597/) for metabolite information
   - Mass spectrometry community forums (ResearchGate, SEQanswers)
