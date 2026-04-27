# Metabolomics Integration Guide

How the metabolomics module integrates with other METAINFORMANT modules and external tools for multi-omic analysis and extended functionality.

## Table of Contents

1. [Module Integration](#module-integration)
2. [Multi-Omics Integration](#multi-omics-integration)
3. [External Tool Integration](#external-tool-integration)
4. [Data Exchange Formats](#data-exchange-formats)
5. [Workflow Orchestration](#workflow-orchestration)

---

## Module Integration

### With `core` Module

The metabolomics module relies on core infrastructure for I/O, logging, and configuration:

```python
from metainformant.core import io as core_io
from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Use core_io for JSON/config files (not MS data)
config = core_io.load_json("config.json")
```

**Key touchpoints**:
- **`core.io`** — Reading/writing JSON, YAML config files; MS data uses own `io.formats`
- **`core.config`** — Environment variable handling (`METAINFORMANT_OUTPUT`, etc.)
- **`core.utils.logging`** — Structured logging across all functions

---

### With `visualization` Module

The `metainformant.visualization` module provides advanced plotting utilities that complement metabolomics-specific plots.

**Integration pattern**:

```python
from metainformant.metabolomics import io, analysis
from metainformant.visualization import plots as viz  # shared plotting module

# Load data
dataset = io.formats.read_csv("data.csv")
norm = analysis.identification.normalize_intensities(dataset.intensities, method="log2")

# Use shared visualization
pca_fig = viz.pca(norm.T, sample_labels=dataset.samples)
heatmap_fig = viz.heatmap(norm, row_labels=dataset.metabolites)

# Metabolomics-specific plot (volcano) from its own visualization submodule
from metainformant.metabolomics.visualization import plots as metab_viz
volcano = metab_viz.volcano_plot(t_stats, pvals, labels=dataset.metabolites)
```

**Recommended**: Use `metainformant.visualization` for PCA, heatmaps, clustering dendrograms; use `metabolomics.visualization` for domain-specific plots (volcano, enrichment bar charts).

---

### With `multiomics` Module

**Primary integration point**: Metabolite-gene correlation and integrated network analysis.

```python
from metainformant.metabolomics import analysis as meta
from metainformant.multiomics import integration as multi

# Metabolite matrix (metabolites × samples)
M = meta.normalize_intensities(intensities, method="log2")

# Gene expression matrix (genes × samples)
E = rna.merged_abundance  # from RNA module

# Correlate across samples
correlations = multi.correlate_layers(M, E, method="spearman")

# Identify metabolite-gene pairs with |ρ| > 0.7
strong_pairs = [(m, g, r) for (m, g), r in correlations.items() if abs(r) > 0.7]
```

**See also**: [multiomics documentation](../multiomics/README.md) for canonical correlation analysis (CCA), MOFA-style factor analysis.

---

### With `protein` Module

Similar to RNA integration but with proteomics data:

```python
from metainformant.protein import quant as prot
from metainformant.metabolomics import analysis as metab

# Protein abundances (samples × proteins, transposed to metabolites×samples for consistency)
P = prot.normalize intensities).T  # shape: proteins × samples

# Correlation
from metainformant.multiomics import correlate_layers
corr = correlate_layers(metab_matrix, P)
```

**Use case**: Identify metabolites correlated with enzyme abundance (e.g., substrate-product relationships).

---

### With `quality` Module

QC metrics from MS runs (e.g., total ion current, number of features, retention time drift) can be analyzed using `metainformant.quality`:

```python
from metainformant.quality import qc_metrics
from metainformant.metabolomics import io

samples = io.formats.read_mgf("batch_run.mgf")
qc_report = qc_metrics.assess_ms_qc(samples)
# Checks: TIC consistency, peak count, mass accuracy drift
```

**Future**: Dedicated metabolomics QC module planned (2027).

---

## Multi-Omics Integration

### Scenario A: Metabolomics + Transcriptomics

**Goal**: Identify genes whose expression correlates with metabolite levels (suggesting transcriptional regulation).

**Method**: Sample-wise correlation across matched metabolomics + RNA-seq samples.

```python
import numpy as np
from scipy.stats import spearmanr
from metainformant.metabolomics import analysis as meta
from metainformant.rna import amalgkit

# Load both datasets (assay matrices must be aligned by sample)
M = meta.io.formats.read_csv("metabolites.csv").intensities  # (n_metabs, n_samples)
R = amalgkit.load_merged_matrix("rna_merged.tsv").T  # (n_genes, n_samples)

# Ensure same sample order
# (load metadata to match, or use sample barcodes)
assert M.shape[1] == R.shape[1], "Sample count mismatch"

# Correlation for each metabolite-gene pair
# Vectorized approach: compute correlation matrix
def corr_matrix(A, B):
    """Pearson correlation between rows of A and rows of B."""
    n = A.shape[1]
    A_centered = A - A.mean(axis=1, keepdims=True)
    B_centered = B - B.mean(axis=1, keepdims=True)
    cov = A_centered @ B_centered.T / n
    varA = (A_centered**2).sum(axis=1, keepdims=True) / n
    varB = (B_centered**2).sum(axis=1, keepdims=True) / n
    return cov / np.sqrt(varA @ varB.T)

corr = corr_matrix(M, R)  # shape: (n_metabs, n_genes)

# Find strong correlations
metab_idx, gene_idx = np.where(np.abs(corr) > 0.8)
for m_i, g_i in zip(metab_idx, gene_idx):
    print(f"{metabolites[m_i]} ↔ {genes[g_i]}: r={corr[m_i, g_i]:.3f}")
```

**Visualization**: Heatmap of top correlated pairs, or network graph.

---

### Scenario B: Multi-Omic Factor Analysis

Use `multiomics` module's latent factor models to find joint signals:

```python
from metainformant.multiomics import factor_analysis
from metainformant.metabolomics import analysis as meta
from metainformant.rna import amalgkit as rna
from metainformant.protein import quant as prot

# Load all omics (standardize to features × samples)
M = meta.normalize_intensities(...)
R = rna.normalize_tpm(...)
P = prot.normalize_intensity(...)

# Joint factor analysis (e.g., MOFA+ style)
factors = factor_analysis.joint_nmf([M, R, P], n_factors=5)

# Factor 1 might represent "cell proliferation" affecting all omics layers
```

**See**: `docs/multiomics/METHODS.md` for mathematical details.

---

### Scenario C: Metabolite Set Enrichment with Gene Sets

Combine metabolite enrichment with gene set enrichment for pathway-level cross-omics:

```python
# Metabolite pathway enrichment (KEGG)
meta_enriched = metabolite_set_enrichment(meta_sig, kegg_metab_pathways)

# Gene set enrichment (GO BP, KEGG)
gene_enriched = gse(genes_sig, gene_pathways)

# Concordance analysis: are the same pathways hit in both omics?
shared_pathways = set(meta_enriched.pathway) & set(gene_enriched.pathway)
print(f"Concordant pathways: {shared_pathways}")
```

**Interpretation**: Concordance strengthens biological conclusion (e.g., "Glycolysis pathway upregulated at both transcript and metabolite levels").

---

## External Tool Integration

### Export to MetaboAnalyst

MetaboAnalyst is a widely-used web-based platform for metabolomics. Export METAINFORMANT results for upload:

```python
from metainformant.metabolomics import analysis, io

# After normalization and differential analysis:
norm = analysis.identification.normalize_intensities(raw, method="log2")
t, p = analysis.identification.differential_abundance(norm, group_a, group_b)

# Create MetaboAnalyst-format CSV
df = pd.DataFrame({
    "metabolite": dataset.metabolites,
    "log2FC": analysis.identification.fold_change(norm, group_a, group_b),
    "p-value": p,
    "q-value": benjamini_hochberg(list(p)),
})
df.to_csv("for_metaboanalyst.csv", index=False)

# Also: normalized intensity matrix
pd.DataFrame(norm, index=dataset.metabolites, columns=dataset.samples)\
  .to_csv("normalized_matrix.csv")
```

**Upload**: Go to [MetaboAnalyst.ca](https://www.metaboanalyst.ca/) → "Pathway Analysis" → upload CSV.

---

### Import from XCMS / MZmine

XCMS (R) and MZmine (Java) are popular peak picking and alignment tools. Their outputs can be read by METAINFORMANT:

**XCMS** outputs `**xcmsexport**` CSV with columns: `name`, `rt`, `mz`, `<sample1>`, `<sample2>`, ...

```python
# XCMS CSV: first non-metadata rows contain metabolite peaks
# Column layout: name | mz | rt | sample1 | sample2 | ...
# We need to extract the intensity matrix (skip mz/rt cols)

dataset = io.formats.read_csv(
    "xcms_output.csv",
    metabolite_col=0,       # name column
    data_start_col=3,       # first sample column (after name, mz, rt)
)
# After loading, you may want to filter by mz if needed
```

**MZmine** exports CSV similarly. Use `pandas` for more control:

```python
import pandas as pd
df = pd.read_csv("mzmine_output.csv")
intensities = df.iloc[:, 3:].values  # skip first 3 metadata columns
metabolite_names = df["feature_id"].tolist()
```

---

### GNPS Molecular Networking

Export feature table for GNPS networking:

```python
# GNPS requires: feature ID, m/z, RT, intensity per sample
feature_table = []
for i, name in enumerate(dataset.metabolites):
    # Parse name to extract mz and rt if encoded in feature ID
    # Or have separate metadata dict from loading
    feature_table.append({
        "feature_id": name,
        "mz": metadata[name]["mz"],
        "rt": metadata[name]["rt"],
        **{sample: dataset.intensities[i, j] for j, sample in enumerate(dataset.samples)},
    })

import json
with open("gnps_feature_table.json", "w") as f:
    json.dump(feature_table, f, indent=2)
print("Upload to GNPS: https://gnps.ucsd.edu")
```

---

### Spectral Library Search (MS2)

For MS/MS data, export to MGF for external search:

```python
from metainformant.metabolomics.io import formats

spectra = formats.read_mgf("query.mgf")
# Already done — now search against NIST, MassBank, or GNPS libraries
# Use cosine similarity as in EXAMPLES.md Example 4
```

**Integration with Spectral Libraries**:
- **NIST 20**: Commercial, gold standard
- **MassBank**: Open, community-contributed
- **GNPS**: Public repository + networking platform

---

## Data Exchange Formats

### Standardized Metabolomics Tabular (mzTab-M)

**mzTab-M** is a community standard for metabolomics data exchange (HUPO-PSI).

**Status**: Not yet implemented in METAINFORMANT. Future: `io.formats.read_mztab()`.

**Workaround**: Convert CSV to mzTab using `pymztab` (external):

```bash
pip install pymztab
python -m pymztab.convert metabolites.csv output.mzTab
```

---

### ISA-Tab (Investigation / Study / Assay)

For complete metadata + data packages (multi-omics), use ISA-Tab structure.

**Structure**:
```
investigation.txt
study.txt
assay.txt    # metabolomics assay
.../metabolites.txt  # data file
```

**Conversion**: Not built-in; use `isatools` Python package to generate ISA from METAINFORMANT outputs.

---

## Workflow Orchestration

### Chaining with RNA Pipeline

Real-world study: same samples underwent both metabolomics and RNA-seq. Want integrated analysis.

```python
#!/usr/bin/env python3
"""Integrated metabolomics + RNA pipeline."""

import subprocess
from pathlib import Path
from metainformant.metabolomics import io as meta_io, analysis as meta_ana
from metainformant.rna import amalgkit as rna

# Step 1: Run RNA pipeline (assuming samples already downloaded/quantified)
print("=== RNA Analysis ===")
rna_matrix = rna.load_merged("output/rna/merged_abundance.tsv")

# Step 2: Load metabolomics
print("=== Metabolomics Analysis ===")
meta_dataset = meta_io.read_csv("metabolomics/intensities.csv")
meta_norm = meta_ana.identification.normalize_intensities(
    meta_dataset.intensities, method="log2"
)

# Step 3: Align samples (must be same biological replicates)
# Assuming both use same sample IDs in filenames
samples_meta = set(meta_dataset.samples)
samples_rna = set(rna_matrix.columns)
common = sorted(samples_meta & samples_rna)

# Subset to common samples
meta_idx = [meta_dataset.samples.index(s) for s in common]
rna_idx = [rna_matrix.columns.get_loc(s) for s in common]
M_common = meta_norm[:, meta_idx]
R_common = rna_matrix.iloc[:, rna_idx].values

# Step 4: Integrated analysis
from metainformant.multiomics import integrate
corr = integrate.sample_wise_correlation(M_common, R_common)

# Step 5: Save
pd.DataFrame(corr).to_csv("integrated_correlation_matrix.csv")
```

**Benefits**: Joint analysis reveals biology missed by single-omics (e.g., post-transcriptional regulation, allosteric control).

---

### Using `snakemake` or `nextflow` with METAINFORMANT

METAINFORMANT functions are **library calls**, not CLI-first. To integrate into workflow engines:

**Snakemake** rule:
```python
rule metabolomics_analysis:
    input:
        csv="data/metabolites.csv"
    output:
        results="results/differential.csv",
        volcano="figures/volcano.html"
    params:
        group_a=[0,1,2],
        group_b=[3,4,5]
    script:
        "scripts/analyze_metabolomics.py"  # Your script using METAINFORMANT
```

Your `analyze_metabolomics.py` reads `input.csv`, writes `output.results`, `output.volcano`.

**Nextflow** similar: call Python script.

---

## Cross-Module Data Flow

### Example: Spatial Metabolomics (Future)

**Potential**: Spatial transcriptomics + spatial metabolomics (imaging MS) integration via `spatial` module.

```python
from metainformant.spatial import io as spatial_io
from metainformant.metabolomics import io as meta_io

# Load spatial transcriptomics (Visium)
st = spatial_io.visium.load("visium_folder/")

# Load spatial metabolomics (MALDI imaging MS)
msi = meta_io.formats.read_imzml("sample.imzML")  # future extension

# Align by spatial coordinates; correlate gene expression with metabolite intensity at each spot
results = integrate.spatial_correlation(st, msi)
```

**Not yet implemented**; placeholder for future spatial metabolomics submodule.

---

## Interoperability Best Practices

### 1. Consistent Sample Naming

When combining datasets, sample IDs must match exactly.

**Best practice**: Use a master `sample_metadata.csv`:

```csv
sample_id,condition,batch,omics
S001,control,A,metabolomics+rna
S002,control,A,metabolomics+rna
S003,treatment,B,metabolomics
```

Programmatically align:
```python
meta = pd.read_csv("metadata.csv")
meta_omics = meta[meta["omics"].str.contains("metabolomics")]
rna_omics = meta[meta["omics"].str.contains("rna")]

common = sorted(set(meta_omics["sample_id"]) & set(rna_omics["sample_id"]))
```

### 2. Feature Matching Across Omics

Genes and metabolites don't have 1:1 mapping. Use pathway/network databases:

- **KEGG**: Links genes → reactions → metabolites
- **Reactome**: Curated pathways with gene-metabolite pairs
- **HMDB**: Metabolite → protein interactions

**Implementation**: Build a mapping dict:

```python
# gene_id -> list[metabolite]
GENE_TO_METAB = {
    "PFKL": ["Fructose-6-phosphate", "Fructose-1,6-bisphosphate"],
    "LDHA": ["Pyruvate", "Lactate"],
}
```

Then gene → metabolite correlation can be tested at pathway level rather than individual pairs.

### 3. Cross-Platform Normalization

Different omics have different scale/distributions:
- **Metabolomics**: Log-normal, heavy right tail
- **RNA-seq**: Count data, negative binomial
- **Proteomics**: Log-normal after normalization

**Solution**: Transform each to comparable scale before integration:

```python
# Rank-based inverse normal transform (common in QTL)
from scipy.stats import rankdata, norm

def int_normalize(x):
    ranks = rankdata(x)
    return norm.ppf((ranks - 0.5) / len(ranks))

M_z = int_normalize(M, axis=1)  # per-metabolite
R_z = int_normalize(R, axis=1)  # per-gene
# Now M_z and R_z are approximately N(0,1) per feature
```

---

## Limitations and Future Directions

### Current Gaps

1. **Limited direct integration**: No built-in `metabolomics.integration` submodule (unlike `spatial.integration`); integration handled via `multiomics`.
2. **No pathway mapping API**: Must construct `pathway_db` manually or load external JSON.
3. **No direct export to multi-omics formats** (e.g., MultiAssayExperiment R class).
4. **Batch correction not built-in**: External library (`pycombat`, `scanpy`) required.

### Planned Enhancements

| Feature | Description | ETA |
|---------|-------------|-----|
| `metabolomics.integration` submodule | Metabolite-gene correlation network builder | Q4 2026 |
| Built-in KEGG/Reactome pathway loader | `load_kegg_pathways()` function | Q3 2026 |
| COM-BAT batch correction wrapper | Native integration | 2027 |
| Multi-omics factor analysis (MOFA) | Joint decomposition of metabolomics+omics | 2027 |
| Spatial metabolomics I/O | imzML, Analyze format readers | 2028 |

---

## FAQ

**Q: Can I use METAINFORMANT metabolomics for targeted (MRM) data?**

A: Yes. Targeted data are just intensity matrices. Use `normalize_intensities(method="median")` or internal standard scaling (custom). The identification functions are less relevant (targets are pre-specified).

**Q: How do I correlate metabolites with clinical phenotypes (not omics)?**

A: Same as gene correlation — phenotype is just another vector:
```python
from scipy.stats import spearmanr
rho, p = spearmanr(metabolite_vector, phenotype_vector)
```

Use `differential_abundance()` for group comparisons (case/control).

**Q: Is there integration with R packages (e.g., `mixOmics`)?**

A: Not directly. Export METAINFORMANT-normalized matrices to CSV and import into R:
```R
M <- read.csv("normalized_metabolites.csv")
R <- read.csv("normalized_rna.csv")
library(mixOmics)
# Use DIABLO, sPLS, etc.
```

Alternatively, call R from Python with `rpy2`.

**Q: Can I use METAINFORMANT for lipidomics?**

A: Yes. Lipidomics is a metabolomics subfield. Load your lipid intensity matrix and use the same functions. For lipid-specific databases (LIPID MAPS), create your own database dict mapping lipid names → m/z.

**Q: Does METAINFORMANT support metabolomics from NMR?**

A: NMR uses different preprocessing (peak alignment, bucketization). Load binned NMR spectra as an intensity matrix, then use downstream stats (normalization, differential, enrichment). Identification via accurate mass is less relevant for NMR (uses chemical shift). Future: NMR-specific utilities planned.

---

## Getting Help

- **Module docs**: See `CAPABILITIES.md` for function-by-function reference.
- **Examples**: See `EXAMPLES.md` for end-to-end workflows.
- **Source code**: `src/metainformant/metabolomics/` — functions are heavily commented.
- **Issues**: File integration bugs on GitHub with minimal reproducible example and data shapes.
