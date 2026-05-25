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

### "ValueError: Empty intensity matrix" from `read_csv()`

**Symptom**:
```python
dataset = io.formats.read_csv("data.csv")
# ValueError: Intensity matrix is empty
```

**Cause**: CSV file has no data rows or all intensity values are non-numeric.

**Check**:
```bash
# Inspect file
head -10 data.csv
cat data.csv | wc -l
```

**Fix**:
- Ensure CSV has header row and at least one data row
- Check delimiter (use `delimiter="\t"` for TSV)
- Verify intensity columns contain numbers, not strings

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
   **Solution**: Let function auto-compute: `background_size=None` (default).

4. **Query metabolites not in pathway database**: Check overlap.
   ```python
   query_set = set(query_metabolites)
   for name, members in pathway_db.items():
       overlap = query_set & set(members)
       print(f"{name}: {len(overlap)}/{len(members)}")
   # If all overlaps are 0 or 1, enrichment underpowered
   ```

---

### "Hypergeometric calculation raised ValueError"

**Symptom**:
```
ValueError: k must be integer or exact integer
```

**Cause**: Input contains non-integer counts (bug in pathway_db where sizes not integers).

**Fix**: Ensure pathway_db values are lists of metabolite names (strings), not counts. The function counts overlaps internally.

**Correct**:
```python
pathway_db = {
    "Glycolysis": ["Glucose", "Pyruvate", "Lactate"],  # list of names
}
```

**Incorrect**:
```python
pathway_db = {
    "Glycolysis": 3,  # WRONG: should be list, not count
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

**Cause**: Metabolomics I/O functions don't automatically respect `METAINFORMANT_OUTPUT`; this is handled at higher level (scripts/CLI).

**Fix**: Within your code, construct paths explicitly:
```python
import os
from pathlib import Path

output_base = Path(os.getenv("METAINFORMANT_OUTPUT", "output"))
output_path = output_base / "metabolomics" / "results.csv"
dataset.to_csv(output_path)
```

Future: Direct support for `META_` prefixed environment variables in module config (planned).

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
# If intensity matrix has zeros (missing) and you compute variance:
var = intensities.var(axis=1)  # zeros reduce variance artificially
```

**Correct**:
```python
imputed = missing_value_imputation(intensities, method="min_half")
var = imputed.var(axis=1)
```

Note: Imputation is only needed for statistics that use variance (PCA, clustering, differential analysis). For presence/absence tests (e.g., detection rate comparison), zeros have meaning.

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
   - [HMDB](https://hmdb.ca/) for metabolite information
   - Mass spectrometry community forums (ResearchGate, SEQanswers)
