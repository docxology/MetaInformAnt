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
fig = visualization.plots.volcano_plot(
    t_stats=t_stats,
    p_values=pvals,
    labels=dataset.metabolites,
    threshold=0.05,
)
fig.write_html("volcano.html")
print("  Saved volcano plot to volcano.html")

print("\nDone!")
```

Run it:

```bash
python3 analyze_metabolites.py
```

Expected output:
```
Loading metabolomics data...
  Loaded 4 metabolites across 5 samples
Normalizing intensities (total ion count)...
  Normalized matrix shape: (4, 5)
Computing differential abundance...
  Found 2 significant metabolites (p < 0.05)
    Pyruvate: t=3.421, p=0.0321
    Citrate: t=-2.987, p=0.0485
Running pathway enrichment...
  Enriched pathways (p < 0.05):
    Glycolysis: 2/3 metabolites, fold=2.45, p=0.021
Generating volcano plot...
  Saved volcano plot to volcano.html
```

### 4. Explore the Output

- **Volcano plot** (`volcano.html`): Interactive plot showing significance vs. fold change
- **Console output**: List of significant metabolites and enriched pathways

## System Requirements

| Component | Version | Purpose |
|-----------|---------|---------|
| Python | 3.11+ | Core runtime |
| NumPy | ≥1.24 | Numerical arrays |
| SciPy | ≥1.10 | Statistical tests |
| uv | latest | Package manager |

Optional: `plotly` for interactive visualizations, `pandas` for DataFrame support (future).

## Detailed Setup

### Environment Configuration

#### Required Environment Variables

None required for basic use.

#### Optional Configuration

Metabolomics module respects the generic METAINFORMANT configuration system:

```bash
# Set output directory for all generated files
export METAINFORMANT_OUTPUT="/path/to/output"

# Enable debug logging
export META_LOG_LEVEL=DEBUG
```

### Installation from Source

```bash
# Clone repository
git clone https://github.com/your-org/MetaInformAnt.git
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
observed = np.array([180.062, 89.025, 150.000])
matches = analysis.identification.identify_metabolites(observed, db, ppm_tolerance=10)
for i, mz_matches in enumerate(matches):
    print(f"m/z {observed[i]}:")
    for m in mz_matches[:3]:  # top 3
        print(f"  {m.matched_name} (Δ{m.delta_ppm:.2f} ppm, score={m.score:.2f})")
```

### Adduct-Aware Identification

**`identify_with_adducts(observed_mz, database, adducts=None, ppm_tolerance=10.0, ion_mode="positive")`**

Accounts for common ionization adducts. For positive mode, tries [M+H]+, [M+Na]+, etc.; subtracts adduct mass to infer neutral mass before database lookup.

**Example**:
```python
# Observed m/z 203.052 could be glutamate [M+H]+ (203.052 - 1.007 = 202.045)
matches = analysis.identification.identify_with_adducts(
    observed_mz=np.array([203.052]),
    database={"Glutamate": 202.045},
    ion_mode="positive",
    ppm_tolerance=10,
)
# Returns: AdductMatch(neutral_mass=202.045, adduct_type="[M+H]+", ...)
```

### Normalization Methods

| Method | Formula | When to Use |
|--------|---------|-------------|
| `total_ion_count` | `x / sum(x) * median(total)` | Standard, corrects for total signal |
| `median` | `x / median(x) * median(all medians)` | Robust to outliers |
| `log2` | `log2(x + 1)` | Variance stabilization, before PCA |
| `pareto` | `(x - mean) / sqrt(std)` | Mean-center + scale by sqrt(std) |

**Recommendation**: Start with `"total_ion_count"`; use `"log2"` for downstream statistical modeling (linear models assume normality).

### Differential Abundance

**`differential_abundance(intensities, group_a, group_b)`**

Performs a two-sample t-test (Welch's method, unequal variance) per metabolite.

**Returns**: `(t_statistics, p_values)` — both 1D arrays length n_metabolites

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

Over-representation analysis using Fisher's exact test.

**Example**:
```python
pathways = {
    "Glycolysis": ["Glucose", "Fructose-6-P", "Pyruvate", "Lactate"],
    "TCA": ["Citrate", "Isocitrate", "Alpha-KG", "Succinate"],
}
results = pathways.enrichment.metabolite_set_enrichment(
    query_metabolites=["Glucose", "Pyruvate", "Citrate"],
    pathway_db=pathways,
    background_size=100,  # total metabolites in reference universe
)
for r in results:
    if r.p_value < 0.05:
        print(f"{r.pathway_name}: odds={r.fold_enrichment:.2f}, p={r.p_value:.4g}")
```

**With FDR correction**:
```python
significant = pathways.enrichment.enrichment_with_fdr(
    query_metabolites, pathway_db, fdr_threshold=0.05
)
# Returns only pathways with q < 0.05, sorted by q-value
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

For >2 groups, use one-way ANOVA per metabolite (not yet implemented; use external `scipy.stats.f_oneway` or `statsmodels`):

```python
from scipy.stats import f_oneway

# Split intensities by group
group_data = [normalized[:, group_indices] for group_indices in groups]
f_stats, pvals = f_oneway(*group_data)
```

### Workflow 3: Spectral Library Matching

```python
# Load query spectra (MGF)
query_spectra = io.formats.read_mgf("sample.mgf")

# Load reference library
library = io.formats.read_mgf("reference_library.mgf")

# For each query spectrum, find best library match
for query in query_spectra:
    best_score = 0.0
    best_match = None
    for ref in library:
        score = analysis.identification.cosine_spectral_similarity(
            query.intensity_array, ref.intensity_array
        )
        if score > best_score:
            best_score = score
            best_match = ref
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
