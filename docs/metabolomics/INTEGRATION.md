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

The metabolomics module is self-contained: it imports only NumPy (and SciPy when available) and does not call `metainformant.core`. It does no logging, reads no environment variables or config files, and parses CSV/MGF itself in `io.formats`. You can still use the core helpers in the scripts around it:

```python
from metainformant.core import io as core_io
from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Use core_io for JSON/config files (not MS data)
config = core_io.load_json("config.json")
```

**Key touchpoints** (all on the caller's side; none is used inside `metabolomics`):
- **`core.io`** — Reading/writing JSON, YAML, CSV and Parquet files; MS data uses the module's own `io.formats`
- **`core.utils.config`** — `load_mapping_from_file`, `load_typed_env`, `apply_env_overrides` for your own pipeline settings (metabolomics functions never read `METAINFORMANT_OUTPUT` or any other variable; see [CONFIGURATION.md](CONFIGURATION.md))
- **`core.utils.logging`** — `get_logger` for your own scripts (call `configure_logging_from_env()` to apply `CORE_LOG_LEVEL`; `get_logger` alone does not read it); metabolomics functions emit no log records

---

### With `visualization` Module

`metabolomics.visualization.plots` only **prepares plot-ready data** (volcano points, PCA scores, heatmap matrices, spectrum and chromatogram data); it draws nothing. Draw the prepared data with matplotlib, or pass it to the shared helpers in `metainformant.visualization.plots` (`pca_plot`, `volcano_plot`, `expression_heatmap`, `correlation_heatmap`), which take pandas DataFrames, return matplotlib `Axes` and save a figure when you give `output_path`.

**Integration pattern**:

```python
import pandas as pd
from metainformant.metabolomics import io, analysis
from metainformant.metabolomics.pathways.enrichment import benjamini_hochberg
from metainformant.metabolomics.visualization import plots as metab_viz
from metainformant.visualization import plots as viz  # shared matplotlib plotting helpers

# Load data
dataset = io.formats.read_csv("data.csv")
norm = analysis.identification.normalize_intensities(dataset.intensities, method="log2")

# PCA: the metabolomics helper computes the scores (metabolites x samples in),
# the shared helper draws them (PC1/PC2 columns plus the variance of each in percent)
scores, explained = metab_viz.pca_metabolomics(norm, n_components=2)
pca_df = pd.DataFrame({
    "PC1": scores[:, 0], "PC2": scores[:, 1],
    "PC1_variance": explained[0] * 100, "PC2_variance": explained[1] * 100,
})
viz.pca_plot(pca_df, output_path="pca.png")

# Heatmap: z-scored matrix from the metabolomics helper, drawn by the shared heatmap
hm = metab_viz.intensity_heatmap_data(norm, dataset.metabolites, dataset.samples)
viz.expression_heatmap(
    pd.DataFrame(hm["matrix"], index=hm["row_labels"], columns=hm["col_labels"]),
    output_path="heatmap.png",
)

# Volcano plot: the shared helper wants a log2 fold change and an adjusted p-value per metabolite
group_a, group_b = [0, 1, 2], [3, 4, 5]  # sample indices of the two groups
log2fc = analysis.identification.fold_change(dataset.intensities, group_a, group_b)  # linear-scale input
_, pvals = analysis.identification.differential_abundance(norm, group_a, group_b)
volcano_df = pd.DataFrame(
    {"log2FoldChange": log2fc, "padj": benjamini_hochberg(list(pvals))},
    index=dataset.metabolites,
)
viz.volcano_plot(volcano_df, output_path="volcano.png")
```

**Recommended**: Use `metabolomics.visualization.plots` for the metabolomics-specific data preparation (`volcano_plot_data`, `prepare_spectrum_plot`, `detect_chromatographic_peaks`, `retention_time_alignment_data`) and `metainformant.visualization` (or matplotlib directly) for drawing.

---

### With `multiomics` Module

**Primary integration point**: Metabolite-gene correlation and integrated network analysis.

```python-snippet
import numpy as np
from metainformant.multiomics.analysis.integration import MultiOmicsData, canonical_correlation

# Metabolite matrix (samples × metabolites) and gene expression matrix (samples × genes),
# as DataFrames that share sample IDs
omics = MultiOmicsData(metabolomics=metabolite_df, transcriptomics=expression_df)

# Canonical correlation between the two layers. When the layers together have more features
# than samples (typical with RNA-seq), the default regularization=0 returns some canonical
# correlations of exactly 1 (a perfect in-sample fit); pass regularization > 0 there.
X_c, Y_c, x_weights, y_weights, correlations = canonical_correlation(
    omics, layers=("metabolomics", "transcriptomics"), n_components=5
)

# Metabolites and genes contributing most to the first canonical component
top_metabolites = np.argsort(-np.abs(x_weights[:, 0]))[:10]
top_genes = np.argsort(-np.abs(y_weights[:, 0]))[:10]
```

For all-pairs metabolite-gene correlation, see the numpy implementation in Scenario A below.

**See also**: [multiomics documentation](../multiomics/README.md) for canonical correlation analysis (CCA), MOFA-style factor analysis.

---

### With `protein` Module

Similar to RNA integration but with proteomics data. `metainformant.protein` has no abundance-quantification code, so the proteomics layer goes through the same `multiomics` helpers:

```python-snippet
from metainformant.multiomics.analysis.integration import (
    MultiOmicsData,
    canonical_correlation,
    from_protein_abundance,
)

# Protein abundances (proteins × samples file -> z-scored samples × proteins DataFrame)
P = from_protein_abundance("protein_abundance.csv", transpose=True)

# Association between the metabolite and protein layers
# (metabolite_df: samples × metabolites DataFrame with the same sample IDs)
omics = MultiOmicsData(metabolomics=metabolite_df, proteomics=P)
X_c, Y_c, x_weights, y_weights, corr = canonical_correlation(
    omics, layers=("metabolomics", "proteomics"), n_components=5
)
```

**Use case**: Identify metabolite and protein features that vary together (e.g., substrate-product relationships); the canonical weights show which features drive each component.

---

### With `quality` Module

QC metrics from MS runs (e.g., total ion current, number of peaks) can be screened with the general-purpose outlier detection in `metainformant.quality`:

```python-snippet
from metainformant.quality.analysis.metrics import detect_outliers
from metainformant.metabolomics import io

spectra = io.formats.read_mgf("batch_run.mgf")

# TIC consistency: flag spectra whose total ion current is an outlier within the run.
# detect_outliers returns a dict ("outliers", "outlier_indices", "outlier_percentage", ...).
tic = detect_outliers([s.total_ion_current for s in spectra], method="iqr")
flagged = [spectra[i] for i in tic["outlier_indices"]]

# Peak count: drop spectra with fewer than min_peaks peaks (raise min_tic above 0 to also drop low-signal spectra)
usable = io.formats.filter_spectra(spectra, min_peaks=5, min_tic=0.0)
```

`metainformant.quality` has no MS-specific QC routine (for example mass accuracy drift), and there is no dedicated metabolomics QC module.

---

## Multi-Omics Integration

### Scenario A: Metabolomics + Transcriptomics

**Goal**: Identify genes whose expression correlates with metabolite levels (suggesting transcriptional regulation).

**Method**: Sample-wise correlation across matched metabolomics + RNA-seq samples.

```python
import numpy as np
import pandas as pd
from metainformant.metabolomics import io

# Load both datasets (assay matrices must be aligned by sample)
dataset = io.formats.read_csv("metabolites.csv")
metabolites = dataset.metabolites
M = dataset.intensities  # (n_metabs, n_samples); a missing cell is NaN, so impute it first (missing_value_imputation) or the correlations of that metabolite are NaN

# Merged RNA abundance table (for example the amalgkit merge output):
# tab-separated, genes as rows, samples as columns
rna = pd.read_csv("rna_merged.tsv", sep="\t", index_col=0)

# Ensure same sample order: select the RNA columns in the metabolomics sample order
# (raises KeyError if a sample is missing; use sample barcodes/metadata if the IDs differ)
rna = rna[dataset.samples]
genes = rna.index.tolist()
R = rna.to_numpy(dtype=float)  # (n_genes, n_samples)

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

```python-snippet
import numpy as np
from metainformant.multiomics.methods.factorization import mofa_simple
from metainformant.multiomics.analysis.integration import (
    from_metabolomics,
    from_protein_abundance,
    from_rna_expression,
)

# Load all omics (z-scored samples × features DataFrames with the same sample order)
M = from_metabolomics(metabolite_df)
R = from_rna_expression("expression.tsv", transpose=True)
P = from_protein_abundance("protein_abundance.csv", transpose=True)

# Joint factor analysis (MOFA+ style)
result = mofa_simple({"metabolomics": M, "transcriptomics": R, "proteomics": P}, k=5)
factors = np.asarray(result["factors"])  # shared sample factors (n_samples × 5), returned as nested lists

# Factor 1 might represent "cell proliferation" affecting all omics layers
```

**See**: `docs/multiomics/methods.md` for mathematical details.

---

### Scenario C: Metabolite Set Enrichment with Gene Sets

Combine metabolite enrichment with gene set enrichment for pathway-level cross-omics:

```python
from metainformant.metabolomics.pathways.enrichment import enrichment_with_fdr
from metainformant.ontology.pathway_enrichment.enrichment import over_representation_analysis

# Metabolite pathway enrichment (KEGG): keeps the pathways with a BH q-value <= 0.05
# (meta_sig: significant metabolite names; kegg_metab_pathways: pathway -> metabolite names)
meta_enriched = enrichment_with_fdr(meta_sig, kegg_metab_pathways, fdr_threshold=0.05)

# Gene-set over-representation (GO BP, KEGG): one result dict per gene set
# (genes_sig: significant gene symbols; gene_pathways: pathway -> gene symbols)
gene_enriched = over_representation_analysis(genes_sig, gene_pathways, correction="fdr_bh")

# Concordance analysis: are the same pathways hit in both omics?
# (both pathway databases must use the same pathway names/IDs for this comparison to mean anything)
meta_hits = {r.pathway_name for r in meta_enriched}
gene_hits = {r["term_id"] for r in gene_enriched if r["adjusted_p"] < 0.05}
shared_pathways = meta_hits & gene_hits
print(f"Concordant pathways: {shared_pathways}")
```

**Interpretation**: Concordance strengthens biological conclusion (e.g., "Glycolysis pathway enriched among both the significant genes and the significant metabolites"; over-representation does not show the direction of change).

---

## External Tool Integration

### Export to MetaboAnalyst

MetaboAnalyst is a widely-used web-based platform for metabolomics. Export METAINFORMANT results for upload:

```python
import pandas as pd
from metainformant.metabolomics import analysis, io
from metainformant.metabolomics.pathways.enrichment import benjamini_hochberg

dataset = io.formats.read_csv("intensities.csv")
group_a, group_b = [0, 1, 2], [3, 4, 5]  # sample indices of the two groups

# After normalization and differential analysis:
norm = analysis.identification.normalize_intensities(dataset.intensities, method="log2")
t, p = analysis.identification.differential_abundance(norm, group_a, group_b)

# Create MetaboAnalyst-format CSV
df = pd.DataFrame({
    "metabolite": dataset.metabolites,
    # fold_change() is the log2 ratio of group means: give it linear-scale intensities, not `norm`
    "log2FC": analysis.identification.fold_change(dataset.intensities, group_a, group_b),
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

**XCMS** feature tables exported to CSV can be arranged with columns: `name`, `mz`, `rt`, `<sample1>`, `<sample2>`, ...

```python
from metainformant.metabolomics import io

# XCMS CSV: one row per detected feature (peak)
# Column layout: name | mz | rt | sample1 | sample2 | ...

dataset = io.formats.read_csv(
    "xcms_output.csv",
    metabolite_col=0,       # name column
    data_start_col=3,       # first sample column (after name, mz, rt)
)

# The columns between the name and the first sample (mz, rt) are per-feature annotations:
# read_csv() returns them in dataset.metadata as {name: {"mz": ..., "rt": ...}}.
# Empty cells in the intensity columns (or NA) come back as NaN, and any other non-numeric
# cell raises ValueError naming the line, feature and sample.
mz_of_first_feature = dataset.metadata[dataset.metabolites[0]]["mz"]
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
# GNPS feature-based molecular networking reads a CSV/TSV feature table (feature ID, m/z, RT,
# intensity per sample; see the GNPS documentation for the exact column names) plus an MS/MS MGF.
# dataset.metadata maps feature ID -> {"mz": ..., "rt": ...} when the table was read with the
# annotation columns before data_start_col (see the XCMS example above)
feature_table = []
for i, name in enumerate(dataset.metabolites):
    feature_table.append({
        "feature_id": name,
        "mz": dataset.metadata[name]["mz"],
        "rt": dataset.metadata[name]["rt"],
        **{sample: dataset.intensities[i, j] for j, sample in enumerate(dataset.samples)},
    })

import pandas as pd
pd.DataFrame(feature_table).to_csv("gnps_feature_table.csv", index=False)
print("Upload to GNPS: https://gnps.ucsd.edu")
```

---

### Spectral Library Search (MS2)

For MS/MS data, load MGF spectra with `read_mgf()` (`write_mgf()` writes them back out for external search tools):

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

**Status**: Not implemented in METAINFORMANT. `metabolomics.io.formats` reads and writes only intensity-matrix CSV and MGF; there is no mzTab-M reader or writer.

**Workaround**: Export the intensity matrix with `io.formats.write_csv()` and convert it with an external mzTab-M tool; METAINFORMANT provides no converter.

---

### ISA-Tab (Investigation / Study / Assay)

For complete metadata + data packages (multi-omics), use ISA-Tab structure.

**Structure**:
```
i_investigation.txt
s_study.txt
a_assay.txt    # metabolomics assay
.../metabolites.txt  # data file
```

**Conversion**: Not built-in; use `isatools` Python package to generate ISA from METAINFORMANT outputs.

---

## Workflow Orchestration

### Chaining with RNA Pipeline

Real-world study: same samples underwent both metabolomics and RNA-seq. Want integrated analysis.

```python-snippet
#!/usr/bin/env python3
"""Integrated metabolomics + RNA pipeline."""

import pandas as pd
from metainformant.metabolomics import io as meta_io, analysis as meta_ana

# Step 1: RNA pipeline output (assuming samples were already downloaded/quantified and merged,
# for example with metainformant.rna.amalgkit.merge): genes as rows, samples as columns
print("=== RNA Analysis ===")
rna_matrix = pd.read_csv("output/rna/merged_abundance.tsv", sep="\t", index_col=0)

# Step 2: Load metabolomics
print("=== Metabolomics Analysis ===")
meta_dataset = meta_io.formats.read_csv("metabolomics/intensities.csv")
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

# Step 4: Integrated analysis (sample-by-sample correlation over both omics layers)
from metainformant.multiomics.analysis.integration import compute_multiomics_similarity
corr = compute_multiomics_similarity(
    {
        "metabolomics": pd.DataFrame(M_common.T, index=common),
        "transcriptomics": pd.DataFrame(R_common.T, index=common),
    },
    method="correlation",
)

# Step 5: Save
pd.DataFrame(corr, index=common, columns=common).to_csv("integrated_correlation_matrix.csv")
```

**Benefits**: Joint analysis reveals biology missed by single-omics (e.g., post-transcriptional regulation, allosteric control).

---

### Using `snakemake` or `nextflow` with METAINFORMANT

METAINFORMANT functions are **library calls**, not CLI-first. To integrate into workflow engines:

**Snakemake** rule:
```snakemake
rule metabolomics_analysis:
    input:
        csv="data/metabolites.csv"
    output:
        results="results/differential.csv",
        volcano="figures/volcano.png"
    params:
        group_a=[0,1,2],
        group_b=[3,4,5]
    script:
        "scripts/analyze_metabolomics.py"  # Your script using METAINFORMANT
```

Your `analyze_metabolomics.py` reads `snakemake.input.csv` and `snakemake.params.group_a` / `group_b`, and writes `snakemake.output.results` and `snakemake.output.volcano` (Snakemake injects the `snakemake` object into scripts run through `script:`).

**Nextflow** similar: call Python script.

---

## Cross-Module Data Flow

### Example: Spatial Metabolomics (Future)

**Potential**: Spatial transcriptomics + spatial metabolomics (imaging MS) integration via `spatial` module.

```python
from metainformant.spatial.io.visium import load_visium

# Load spatial transcriptomics (Visium): a SpatialDataset with expression, coordinates,
# barcodes, gene_names, ...
st = load_visium("visium_folder/")

# Spatial metabolomics (MALDI imaging MS) has no counterpart yet: there is no imzML reader
# in metabolomics.io.formats, and nothing in `spatial` correlates gene expression with an
# imaging-MS intensity at each spot (spatial.integration only maps scRNA-seq data onto spatial
# data: label transfer, correlation mapping and gene imputation).
```

**Not implemented**: only spatial transcriptomics loading exists (Visium above; `spatial.io` also has MERFISH and Xenium loaders). Imaging-MS (imzML) readers and spot-level gene-metabolite correlation are not in the codebase.

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
import pandas as pd

meta = pd.read_csv("sample_metadata.csv")
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

def int_normalize(x, axis=1):
    ranks = rankdata(x, axis=axis)
    return norm.ppf((ranks - 0.5) / x.shape[axis])

M_z = int_normalize(M, axis=1)  # per-metabolite (rows = metabolites, columns = samples)
R_z = int_normalize(R, axis=1)  # per-gene
# Now M_z and R_z are approximately N(0,1) per feature
```

---

## Limitations and Future Directions

### Current Gaps

1. **Limited direct integration**: No built-in `metabolomics.integration` submodule (unlike `spatial.integration`); integration handled via `multiomics`.
2. **No pathway mapping API**: Must construct `pathway_db` manually or load external JSON.
3. **No direct export to multi-omics formats** (e.g., MultiAssayExperiment R class).
4. **Batch correction has no covariates**: `metainformant.quality.batch.detection.correct_batch_combat` is parametric empirical-Bayes ComBat (Johnson, Li & Rabinovic 2007) without covariates, so biological groups that coincide with batches are removed along with them. For a design that must preserve covariates use an external implementation such as `scanpy.pp.combat`.

### Not Implemented

None of the following exists in the codebase today:

| Feature | Description |
|---------|-------------|
| `metabolomics.integration` submodule | Metabolite-gene correlation network builder |
| Built-in KEGG/Reactome pathway loader | A `load_kegg_pathways()`-style function returning a `pathway_db` dict |
| ComBat with covariates | Empirical-Bayes batch correction that preserves biological covariates (`correct_batch_combat` is ComBat without covariates) |
| Spatial metabolomics I/O | imzML and Analyze format readers |

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
# mixOmics expects samples in rows; the module's matrices are metabolites x samples, so transpose
M <- read.csv("normalized_metabolites.csv")
R <- read.csv("normalized_rna.csv")
library(mixOmics)
# Use DIABLO, sPLS, etc.
```

Alternatively, call R from Python with `rpy2`.

**Q: Can I use METAINFORMANT for lipidomics?**

A: Yes. Lipidomics is a metabolomics subfield. Load your lipid intensity matrix and use the same functions. For lipid-specific databases (LIPID MAPS), create your own database dict mapping lipid names → m/z.

**Q: Does METAINFORMANT support metabolomics from NMR?**

A: NMR uses different preprocessing (peak alignment, bucketization). Load binned NMR spectra as an intensity matrix, then use downstream stats (normalization, differential, enrichment). Identification via accurate mass is less relevant for NMR (uses chemical shift). There are no NMR-specific utilities.

---

## Getting Help

- **Module docs**: See `CAPABILITIES.md` for function-by-function reference.
- **Examples**: See `EXAMPLES.md` for end-to-end workflows.
- **Source code**: `src/metainformant/metabolomics/` — every public function has a docstring.
- **Issues**: File integration bugs on GitHub with minimal reproducible example and data shapes.
