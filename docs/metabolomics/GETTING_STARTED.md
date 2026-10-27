# Getting Started with Metabolomics Analysis

Complete guide for setting up and running metabolomics workflows with METAINFORMANT, from data preparation to pathway analysis.

## Quick Start (10-Minute Tutorial)

### 1. Installation

```bash
cd /path/to/MetaInformAnt

# Create virtual environment
uv venv
source .venv/bin/activate

# Install metainformant
uv pip install -e .

# Install dependencies
uv pip install numpy scipy  # For numerical operations
```

### 2. Prepare Your Data

Metabolomics data should be in CSV format with metabolites as rows and samples as columns:

**`metabolites.csv`**:
```csv
metabolite,Sample1,Sample2,Sample3,Sample4,Sample5
Glucose,1250,1340,980,1120,1450
Lactate,450,380,520,410,390
Pyruvate,120,135,95,110,140
Citrate,85,92,78,88,95
```

Write a missing measurement as a blank cell or `NA`: `read_csv` returns it as `NaN` (never as `0.0`), and impute it before normalizing (see below). Any other non-numeric cell, a decimal comma (`1,5`), or a row with a different number of cells than the header raises `ValueError` naming the line instead of being repaired silently.

**Or MGF format** for MS/MS spectra:
```
BEGIN IONS
PEPMASS=180.0634
RTINSECONDS=245.6
SCANS=1
104.0527  45.2
118.0651  12.8
180.0634  100.0
END IONS
```

### 3. Run Your First Analysis

Create a simple Python script `analyze_metabolites.py`:

```python
#!/usr/bin/env python3
"""Quick metabolomics analysis demo."""

from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

from metainformant.metabolomics import (
    io,
    analysis,
    pathways,
    visualization,
)

# ── 1. Load data ──────────────────────────────────────────────────────
print("Loading metabolomics data...")
dataset = io.formats.read_csv("metabolites.csv")
print(f"  Loaded {len(dataset.metabolites)} metabolites across {len(dataset.samples)} samples")

# ── 2. Normalize intensities ──────────────────────────────────────────
print("Normalizing intensities (total ion count)...")
normalized = analysis.identification.normalize_intensities(
    dataset.intensities, method="total_ion_count"
)
print(f"  Normalized matrix shape: {normalized.shape}")

# ── 3. Differential abundance (group A: samples 0-1, group B: samples 2-4) ──
print("Computing differential abundance...")
group_a = [0, 1]   # e.g., control
group_b = [2, 3, 4]  # e.g., treatment
t_stats, pvals = analysis.identification.differential_abundance(
    normalized, group_a, group_b
)

# Find significant metabolites (unadjusted p < 0.05)
significant = np.where(pvals < 0.05)[0]
print(f"  Found {len(significant)} significant metabolites (p < 0.05)")
for idx in significant:
    print(f"    {dataset.metabolites[idx]}: t={t_stats[idx]:.3f}, p={pvals[idx]:.4g}")

# ── 4. Pathway enrichment ────────────────────────────────────────────
print("Running pathway enrichment...")
# Build a simple pathway database (in practice, load from KEGG/Reactome)
pathway_db = {
    "Glycolysis": ["Glucose", "Lactate", "Pyruvate"],
    "TCA Cycle": ["Citrate", "Alpha-ketoglutarate", "Succinate"],
}
query_metabolites = [dataset.metabolites[i] for i in significant]
enriched = pathways.enrichment.metabolite_set_enrichment(
    query_metabolites, pathway_db, background_size=50
)
print(f"  Enriched pathways (p < 0.05):")
for result in enriched:
    if result.p_value < 0.05:
        print(f"    {result.pathway_name}: {result.overlap}/{result.pathway_size} metabolites, "
              f"fold={result.fold_enrichment:.2f}, p={result.p_value:.4g}")

# ── 5. Visualize ────────────────────────────────────────────────────
print("Generating volcano plot...")
# fold_change() is the log2 ratio of group means, so pass linear-scale intensities
log2_fc = analysis.identification.fold_change(normalized, group_a, group_b)

# volcano_plot_data() only prepares the points (it draws nothing); plot them with matplotlib
points = visualization.plots.volcano_plot_data(
    dataset.metabolites, log2_fc, pvals, fc_threshold=1.0, p_threshold=0.05
)
fig, ax = plt.subplots(figsize=(5, 4))
ax.scatter(
    [p.log2_fc for p in points],
    [p.neg_log10_p for p in points],
    c=["tab:red" if p.significant else "tab:gray" for p in points],
)
for p in points:
    ax.annotate(p.name, (p.log2_fc, p.neg_log10_p), fontsize=8)
ax.axhline(-np.log10(0.05), linestyle="--", color="black", linewidth=0.8)
ax.set_xlabel("log2 fold change (group A / group B)")
ax.set_ylabel("-log10 p-value")
fig.savefig("volcano.png", dpi=150, bbox_inches="tight")
print("  Saved volcano plot to volcano.png")

print("\nDone!")
```

Run it:

```bash
python3 analyze_metabolites.py
```

Expected output (for the five-sample example table above):
```
Loading metabolomics data...
  Loaded 4 metabolites across 5 samples
Normalizing intensities (total ion count)...
  Normalized matrix shape: (4, 5)
Computing differential abundance...
  Found 0 significant metabolites (p < 0.05)
Running pathway enrichment...
  Enriched pathways (p < 0.05):
Generating volcano plot...
  Saved volcano plot to volcano.png

Done!
```

With only 2 vs. 3 samples, no metabolite reaches p < 0.05, so the significant-metabolite and enrichment lists stay empty; with a real cohort the loops print one line per hit.

### 4. Explore the Output

- **Volcano plot** (`volcano.png`): Static matplotlib scatter of log2 fold change vs. -log10 p-value, with metabolites that pass both thresholds highlighted
- **Console output**: List of significant metabolites and enriched pathways

## System Requirements

| Component | Version | Purpose |
|-----------|---------|---------|
| Python | 3.11+ | Core runtime |
| NumPy | ≥1.24 (declared); every function runs on any supported version | Numerical arrays |
| SciPy | ≥1.11 (declared) | Welch t-test p-values from Student's t distribution; without SciPy `differential_abundance` falls back to a normal approximation |
| matplotlib | ≥3.8 (declared) | Drawing the plots in the examples |
| uv | latest | Package manager |

The metabolomics package itself imports only NumPy (and SciPy when available). It has no plotting backend: `visualization.plots` prepares plot-ready data, and the examples draw it with matplotlib. SciPy, matplotlib and pandas (used for result tables in the examples) are installed with METAINFORMANT.

## Detailed Setup

### Environment Configuration

#### Required Environment Variables

None.

#### Optional Configuration

None. The metabolomics module reads no environment variables and no configuration files: input and output paths, tolerances and methods are all function arguments (see [CONFIGURATION.md](CONFIGURATION.md)).

### Installation from Source

```bash
# Clone repository
git clone https://github.com/docxology/MetaInformAnt.git
cd MetaInformAnt

# Install with uv (editable mode)
uv pip install -e .

# Verify installation
python3 -c "from metainformant.metabolomics import analysis; print('OK')"
```

## Understanding the Core Functions

### Metabolite Identification

**`identify_metabolites(observed_mz, database, ppm_tolerance=10.0)`**

Matches observed m/z values to a reference database within a specified mass tolerance (ppm). Returns a list of matches per query m/z, sorted by confidence score.

**Parameters**:
- `observed_mz`: 1D numpy array of m/z values from your MS run
- `database`: Dictionary `{metabolite_name: exact_mz}`
- `ppm_tolerance`: Maximum allowed mass error (default 10 ppm)

**Returns**: `list[list[MetaboliteMatch]]` — matches for each observed m/z

**Example**:
```python
db = {"Glucose": 180.0634, "Lactate": 89.0240, "Pyruvate": 87.0083}
observed = np.array([180.062, 89.0245, 150.000])
matches = analysis.identification.identify_metabolites(observed, db, ppm_tolerance=10)
for i, mz_matches in enumerate(matches):
    print(f"m/z {observed[i]}:")
    for m in mz_matches[:3]:  # top 3
        print(f"  {m.matched_name} (Δ{m.delta_ppm:.2f} ppm, score={m.score:.2f})")
# m/z 180.062:
#   Glucose (Δ7.78 ppm, score=0.22)
# m/z 89.0245:
#   Lactate (Δ5.62 ppm, score=0.44)
# m/z 150.0:          <- nothing within 10 ppm
```

Each `MetaboliteMatch` has `query_mz`, `matched_name`, `matched_mz`, `delta_ppm` and `score` (`1 - delta_ppm / ppm_tolerance`, so 1.0 is an exact mass match). The database values are compared with the observed m/z directly, so store the exact m/z of the ion you expect to observe (`identify_metabolites` does no adduct handling). `ppm_tolerance` must be a positive, finite number: zero, a negative value, `NaN`, infinity or a non-number raises `ValueError`.

### Adduct-Aware Identification

**`identify_with_adducts(observed_mz, database, adducts=None, ppm_tolerance=10.0, ion_mode="positive")`**

Accounts for common ionization adducts. For positive mode, tries [M+H]+, [M+Na]+, [M+K]+ and [M+NH4]+ (negative mode: [M-H]-, [M+Cl]-, [M+FA-H]-); subtracts the adduct mass to infer the neutral mass before database lookup, so `database` values are neutral monoisotopic masses. `ion_mode` is `"positive"` or `"negative"` in any case and anything else raises `ValueError`, as does a `ppm_tolerance` that is not a positive finite number.

**Example**:
```python
# Observed m/z 148.0604 could be glutamate [M+H]+ (148.0604 - 1.0073 = 147.0531)
matches = analysis.identification.identify_with_adducts(
    observed_mz=np.array([148.0604]),
    database={"Glutamate": 147.0532},
    ion_mode="positive",
    ppm_tolerance=10,
)
# Returns list[list[AdductMatch]] (one inner list per observed m/z), here:
# [[AdductMatch(query_mz=148.0604, neutral_mass=147.053124, adduct_type='[M+H]+',
#                matched_name='Glutamate', delta_ppm=0.5168...)]]
```

### Normalization Methods

| Method | Formula | When to Use |
|--------|---------|-------------|
| `total_ion_count` | `x / sum(x) * median(total)` | Standard, corrects for total signal |
| `median` | `x / median(x) * median(all medians)` | Robust to outliers |
| `log2` | `log2(x + 1)` | Variance stabilization, before PCA |
| `pareto` | `(x - mean) / sqrt(std)` | Mean-center + scale by sqrt(std) |

**Recommendation**: Start with `"total_ion_count"`; use `"log2"` for downstream statistical modeling (linear models assume normality).

**Missing values**: `"total_ion_count"` and `"median"` raise `ValueError` if the matrix contains `NaN` (for example a cell that `read_csv` returned from `NA`), because one missing value would distort the scale of its whole sample. Impute first:

```python-snippet
imputed = analysis.identification.missing_value_imputation(dataset.intensities, method="min_half")
normalized = analysis.identification.normalize_intensities(imputed, method="total_ion_count")
```

`"log2"` keeps a `NaN` in place, `fold_change()` returns `NaN` (not a floored value) for a metabolite with a missing value in either group, and `"pareto"` turns the whole row of a metabolite that has a `NaN` into `NaN`.

### Differential Abundance

**`differential_abundance(intensities, group_a, group_b)`**

Performs a two-sample t-test (Welch's method, unequal variance) per metabolite.

**Returns**: `(t_statistics, p_values)` — both 1D arrays length n_metabolites; a positive t means a higher mean in `group_a`. There is no built-in multiple-testing step: pass the p-values to `benjamini_hochberg` as below.

**Example**:
```python
t, p = analysis.identification.differential_abundance(
    normalized,
    group_a=[0, 1, 2],    # 3 control samples
    group_b=[3, 4, 5, 6], # 4 treatment samples
)

# Multiple testing correction
from metainformant.metabolomics.pathways.enrichment import benjamini_hochberg
qvals = benjamini_hochberg(list(p))
sig = np.array(qvals) < 0.05
```

### Pathway Enrichment

**`metabolite_set_enrichment(query_metabolites, pathway_db, background_size=None)`**

Over-representation analysis using a one-sided hypergeometric (Fisher's exact) test.

**Returns**: list of `EnrichmentResult` (`pathway_name`, `pathway_size`, `overlap`, `fold_enrichment` = observed / expected overlap, `p_value`, `metabolites_in_pathway`), sorted by p-value.

`background_size` is the size of the universe the query was drawn from. It must be at least as large as the number of distinct query metabolites and as every pathway; a smaller value raises `ValueError`, because it cannot describe a universe that contains them.

**Example**:
```python
from metainformant.metabolomics import pathways

pathway_db = {
    "Glycolysis": ["Glucose", "Fructose-6-P", "Pyruvate", "Lactate"],
    "TCA": ["Citrate", "Isocitrate", "Alpha-KG", "Succinate"],
}
query = ["Glucose", "Pyruvate", "Citrate"]
results = pathways.enrichment.metabolite_set_enrichment(
    query_metabolites=query,
    pathway_db=pathway_db,
    background_size=100,  # total metabolites in reference universe
)
for r in results:
    if r.p_value < 0.05:
        print(f"{r.pathway_name}: fold={r.fold_enrichment:.2f}, p={r.p_value:.4g}")
# Glycolysis: fold=16.67, p=0.003587
```

**With FDR correction**:
```python
significant = pathways.enrichment.enrichment_with_fdr(
    query, pathway_db, background_size=100, fdr_threshold=0.05
)
# Returns the EnrichmentResult entries whose Benjamini-Hochberg q-value is <= 0.05,
# sorted by p-value. The q-values themselves are not returned; compute them with
# pathways.enrichment.benjamini_hochberg() if you need them.
```

## Common Workflows

### Workflow 1: Case-Control Comparison

1. Load intensity matrix from CSV
2. Normalize (TIC)
3. Log2 transform for normality
4. t-tests per metabolite
5. Multiple testing correction (BH-FDR)
6. Volcano plot
7. Pathway enrichment on significant hits

### Workflow 2: Time Series or Multiple Groups

For >2 groups, use one-way ANOVA per metabolite (not implemented in the metabolomics module; use external `scipy.stats.f_oneway` or `statsmodels`):

```python
from scipy.stats import f_oneway

# groups: one list of sample indices per group, e.g. [[0, 1, 2], [3, 4, 5], [6, 7, 8]]
# Split intensities by group (metabolites x samples-in-group)
group_data = [normalized[:, group_indices] for group_indices in groups]

# axis=1 runs one ANOVA per metabolite (the default axis=0 would run it across metabolites)
f_stats, pvals = f_oneway(*group_data, axis=1)
```

### Workflow 3: Spectral Library Matching

`cosine_spectral_similarity()` compares two intensity vectors on the same m/z grid, so bin the spectra onto a common grid first (raw peak lists differ in length and m/z positions):

```python
import numpy as np
from metainformant.metabolomics import io, analysis

def bin_spectrum(spec, min_mz=50.0, max_mz=1000.0, bin_width=0.1):
    """Sum peak intensities into fixed m/z bins."""
    edges = np.arange(min_mz, max_mz + bin_width, bin_width)
    binned, _ = np.histogram(spec.mz_array, bins=edges, weights=spec.intensity_array)
    return binned

# Load query spectra (MGF)
query_spectra = io.formats.read_mgf("sample.mgf")

# Load reference library and bin it once
library = io.formats.read_mgf("reference_library.mgf")
library_bins = [bin_spectrum(ref) for ref in library]

# For each query spectrum, find best library match
for query in query_spectra:
    query_bins = bin_spectrum(query)
    best_score = 0.0
    best_match = None
    for ref, ref_bins in zip(library, library_bins):
        score = analysis.identification.cosine_spectral_similarity(query_bins, ref_bins)
        if score > best_score:
            best_score = score
            best_match = ref
    if best_match is not None:
        print(f"{query.scan_number} -> {best_match.scan_number} (cosine={best_score:.3f})")
```

## Next Steps

- Read [ARCHITECTURE.md](ARCHITECTURE.md) for system design details
- See [EXAMPLES.md](EXAMPLES.md) for real-world analysis scenarios
- Consult [CONFIGURATION.md](CONFIGURATION.md) for module-specific settings
- Review [CAPABILITIES.md](CAPABILITIES.md) for complete function reference
- Check [TROUBLESHOOTING.md](TROUBLESHOOTING.md) for common issues

## External References

- **MS Data Formats**: mzML (HUPO-PSI), mzXML (ISB), MGF (spectral library)
- **Preprocessing Tools**: XCMS (R), MZmine (Java), MS-DIAL
- **Databases**: HMDB, KEGG, METLIN, LipidMaps, GNPS
- **Normalization**: NormalyzerDE,statistical methods for LC-MS
- **Enrichment**: MSEA (MetaboAnalyst), GSEA for metabolites
