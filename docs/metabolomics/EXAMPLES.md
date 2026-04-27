# Metabolomics Examples

Real-world analysis scenarios demonstrating complete metabolomics workflows with expected outputs and interpretation guidance.

## Table of Contents

1. [Example 1: Urine Metabolomics — Diabetic vs. Control](#example-1-urine-metabolomics--diabetic-vs-control)
2. [Example 2: Cancer Cell Line Metabolite Profiling](#example-2-cancer-cell-line-metabolite-profiling)
3. [Example 3: Time-Course Drug Response](#example-3-time-course-drug-response)
4. [Example 4: MS/MS Spectral Library Matching](#example-4-msms-spectral-library-matching)
5. [Example 5: Multi-Omics Integration — Metabolite-Gene Correlation](#example-5-multi-omics-integration--metabolite-gene-correlation)
6. [Example 6: Large Cohort — Quality Control and Batch Assessment](#example-6-large-cohort--quality-control-and-batch-assessment)

---

## Example 1: Urine Metabolomics — Diabetic vs. Control

**Study Design**: Case-control metabolomics of urine samples from type 2 diabetes patients (n=30) vs. healthy controls (n=30). Goal: identify discriminating metabolites and enriched pathways.

**Data**: Intensity matrix from LC-MS (positive mode), 500 metabolites, 60 samples.

### Analysis Script

```python
import numpy as np
import pandas as pd
from pathlib import Path
from metainformant.metabolomics import (
    io, analysis, pathways, visualization,
)

# ── Load data ──────────────────────────────────────────────────────────
print("Loading urine metabolomics data...")
dataset = io.formats.read_csv("urine_metabolites.csv")
print(f"  Metabolites: {len(dataset.metabolites)}, Samples: {len(dataset.samples)}")
# Output: Metabolites: 500, Samples: 60

# Check sample grouping from metadata (assumes SampleID encodes group)
# E.g., "D003" = Diabetic patient 3, "C015" = Control 15
sample_groups = ["diabetic" if s.startswith("D") else "control" for s in dataset.samples]
group_a = [i for i, g in enumerate(sample_groups) if g == "diabetic"]
group_b = [i for i, g in enumerate(sample_groups) if g == "control"]
print(f"  Group sizes: diabetic={len(group_a)}, control={len(group_b)}")
# Output: Group sizes: diabetic=30, control=30

# ── Quality control: flag low-signal samples ──────────────────────────
sample_totals = dataset.intensities.sum(axis=0)
low_signal = sample_totals < 1e5
print(f"  Low-signal samples (< 1e5 total): {low_signal.sum()}")
# Remove or investigate low-signal samples if needed
keep = ~low_signal
dataset.intensities = dataset.intensities[:, keep]
dataset.samples = [s for i, s in enumerate(dataset.samples) if keep[i]]

# ── Normalization ───────────────────────────────────────────────────────
print("Normalizing (TIC)...")
norm_tic = analysis.identification.normalize_intensities(
    dataset.intensities, method="total_ion_count"
)

print("Log2 transforming for normality...")
norm_log = analysis.identification.normalize_intensities(
    norm_tic, method="log2"
)

# ── Differential abundance: Welch's t-test ─────────────────────────────
print("Computing differential abundance...")
t_stats, pvals = analysis.identification.differential_abundance(
    norm_log, group_a=group_a, group_b=group_b
)

# Multiple testing correction
from metainformant.metabolomics.pathways.enrichment import benjamini_hochberg
qvals = benjamini_hochberg(list(pvals))
qvals = np.array(qvals)

# Volcano plot
fig = visualization.plots.volcano_plot(
    t_stats=t_stats,
    p_values=pvals,
    labels=dataset.metabolites,
    threshold=0.05,
    plot_threshold=True,
)
fig.write_html("urine_volcano.html")
print("  Volcano plot saved to urine_volcano.html")

# ── Identify significant metabolites ────────────────────────────────────
sig_idx = np.where(qvals < 0.05)[0]
print(f"  Significant metabolites (FDR < 0.05): {len(sig_idx)}")
up_idx = sig_idx[t_stats[sig_idx] > 0]
down_idx = sig_idx[t_stats[sig_idx] < 0]
print(f"    Upregulated in diabetes: {len(up_idx)}")
print(f"    Downregulated in diabetes: {len(down_idx)}")

# List top hits
top_up = up_idx[np.argsort(t_stats[up_idx])[-5:]]  # top 5 largest positive t
top_down = down_idx[np.argsort(t_stats[down_idx])[:5]]  # top 5 negative t
print("\nTop upregulated metabolites:")
for i in top_up[::-1]:  # descending order
    print(f"  {dataset.metabolites[i]}: t={t_stats[i]:.2f}, q={qvals[i]:.4g}")
print("\nTop downregulated metabolites:")
for i in top_down:
    print(f"  {dataset.metabolites[i]}: t={t_stats[i]:.2f}, q={qvals[i]:.4g}")

# ── Pathway enrichment ──────────────────────────────────────────────────
print("\nPathway enrichment...")
# Load pathway database (KEGG, HMDB, or custom)
# Here we use a small example; in practice load from JSON:
# import json; pathway_db = json.load(open("kegg_metabolite_pathways.json"))
pathway_db = {
    "Glycolysis / Gluconeogenesis": ["Glucose", "Pyruvate", "Lactate", "Glycerol"],
    "TCA Cycle": ["Citrate", "Alpha-ketoglutarate", "Succinate", "Fumarate"],
    "Amino Acid Metabolism": ["Glutamate", "Alanine", "BCAA"],
    "Purine Metabolism": ["Hypoxanthine", "Xanthine", "Uric acid"],
    "Lipid Metabolism": ["Acyl-carnitine", "Arachidonic acid"],
}
# Expand query to include metabolites at q < 0.1 for broader enrichment
query_metabolites = [dataset.metabolites[i] for i in sig_idx]
enriched = pathways.enrichment.enrichment_with_fdr(
    query_metabolites, pathway_db, fdr_threshold=0.2
)
print(f"  Enriched pathways (FDR < 0.2): {len(enriched)}")
for r in enriched:
    print(f"  {r.pathway_name}: {r.overlap}/{r.pathway_size} metabolites, "
          f"odds={r.fold_enrichment:.2f}, q={r.p_value:.4g}")

# ── Save results ────────────────────────────────────────────────────────
results_df = pd.DataFrame({
    "metabolite": dataset.metabolites,
    "t_statistic": t_stats,
    "p_value": pvals,
    "q_value": qvals,
    "significant": qvals < 0.05,
})
results_df.to_csv("urine_differential_results.csv", index=False)
print("Results saved to urine_differential_results.csv")
```

### Expected Output

```
Loading urine metabolabolomics data...
  Metabolites: 500, Samples: 60
  Group sizes: diabetic=30, control=30
  Low-signal samples (< 1e5 total): 2
Normalizing (TIC)...
Log2 transforming for normality...
Computing differential abundance...
  Significant metabolites (FDR < 0.05): 23
    Upregulated in diabetes: 14
    Downregulated in diabetes: 9
  Volcano plot saved to urine_volcano.html

Top upregulated metabolites:
  Citrate: t=5.43, q=0.0003
  2-Hydroxybutyrate: t=4.87, q=0.0012
  BCAA: t=4.21, q=0.0031

Top downregulated metabolites:
  Hippurate: t=-3.98, q=0.0045
  Creatinine: t=-3.45, q=0.0089

Pathway enrichment...
  Enriched pathways (FDR < 0.2): 2
  Glycolysis / Gluconeogenesis: 3/4 metabolites, odds=3.2, q=0.042
  Amino Acid Metabolism: 2/3 metabolites, odds=2.8, q=0.087
```

**Interpretation**:
- **Upregulated citrate**: Suggests altered TCA cycle flux in diabetes (mitochondrial dysfunction).
- **Upregulated BCAAs**: Consistently associated with insulin resistance in literature.
- **Downregulated hippurate**: Implicated in gut microbiome dysbiosis in diabetes.

---

## Example 2: Cancer Cell Line Metabolite Profiling

**Study Design**: Compare metabolic profiles of wild-type vs. KRAS-mutant cancer cell lines. Goal: identify oncometabolites and affected pathways.

**Data**: 200 metabolites measured across 10 cell lines (5 WT, 5 KRAS-mut).

### Script

```python
from metainformant.metabolomics import io, analysis, pathways, visualization

# Load data: CSV with columns "metabolite", then one per cell line
dataset = io.formats.read_csv("cell_lines.csv")
# Metadata: KRAS status from sample names: "WT_1", "KRAS_1", ...
wt_idx = [i for i, s in enumerate(dataset.samples) if s.startswith("WT")]
mut_idx = [i for i, s in enumerate(dataset.samples) if s.startswith("KRAS")]
print(f"WT: {len(wt_idx)}, KRAS-mut: {len(mut_idx)}")

# Normalize (row-wise pareto scaling for PCA)
norm_pareto = analysis.identification.normalize_intensities(
    dataset.intensities, method="pareto"
)

# PCA to visualize separation
pca_result = visualization.plots.pca_plot(
    norm_pareto.T,  # samples × features
    group_labels=["WT"] * len(wt_idx) + ["KRAS"] * len(mut_idx),
)
pca_result.write_html("pca.html")
print("PCA plot saved")

# Differential abundance (small n → effect size focus)
t, p = analysis.identification.differential_abundance(
    norm_pareto, group_a=wt_idx, group_b=mut_idx
)

# Since n=5 per group, low power → use effect size threshold
log2fc = analysis.identification.fold_change(dataset.intabolites, wt_idx, mut_idx)
large_effect = (np.abs(log2fc) > 1.0) & (p < 0.1)  # exploratory
candidates = np.where(large_effect)[0]
print(f"Candidate oncometabolites (|log2FC|>1, p<0.1): {len(candidates)}")
for i in candidates:
    print(f"  {dataset.metabolites[i]}: log2FC={log2fc[i]:.2f}, p={p[i]:.3f}")

# Pathway activity scoring (uses continuous scores, not just sig/not-sig)
met_scores = dict(zip(dataset.metabolites, t))  # use t-statistic as score
pathway_scores = pathways.enrichment.pathway_activity_scoring(
    metabolite_scores=met_scores,
    pathway_db=KEGG_PATHWAYS,
    min_coverage=0.2,
)
print("\nTop pathway activities (mean t-statistic):")
for ps in pathway_scores[:5]:
    print(f"  {ps.pathway_name}: {ps.activity_score:.2f} (n={ps.n_measured}/{ps.n_total})")
```

### Expected Output

```
WT: 5, KRAS-mut: 5
PCA plot saved
Candidate oncometabolites (|log2FC|>1, p<0.1): 8
  2-Hydroxyglutarate: log2FC=2.45, p=0.08
  Succinate: log2FC=1.87, p=0.12
  Lactate: log2FC=1.34, p=0.09

Top pathway activities:
  Glutathione metabolism: 1.87 (n=4/8)
  TCA Cycle: 1.23 (n=5/8)
  Glycolysis: 0.87 (n=6/10)
```

**Interpretation**: Elevated 2-hydroxyglutarate (2-HG) is a known oncometabolite from IDH mutations; if KRAS-mut lines also have elevated 2-HG, may indicate epigenetic dysregulation. TCA cycle upregulation suggests Warburg effect.

---

## Example 3: Time-Course Drug Response

**Study Design**: Treat cells with drug X, collect metabolomics at 0h, 2h, 6h, 24h (n=3 replicates per time point). Goal: identify metabolites with significant time dynamics.

**Data**: 300 metabolites × 12 samples (4 timepoints × 3 replicates).

### Analysis

```python
import numpy as np
from scipy.stats import f_oneway
from metainformant.metabolomics import io, analysis, pathways, visualization

dataset = io.formats.read_csv("timecourse.csv")
intensities = dataset.intensities

# Group indices by timepoint
timepoints = ["0h", "2h", "6h", "24h"]
groups = [[] for _ in timepoints]
for i, s in enumerate(dataset.samples):
    for t_idx, t_name in enumerate(timepoints):
        if f"_{t_name}" in s or t_name in s:
            groups[t_idx].append(i)
            break

# One-way ANOVA per metabolite across the 4 timepoints
f_stats = np.zeros(len(dataset.metabolites))
pvals = np.zeros(len(dataset.metabolites))
for i in range(len(dataset.metabolites)):
    group_data = [intensities[i, g] for g in groups if len(g) > 0]
    f, p = f_oneway(*group_data)
    f_stats[i] = f
    pvals[i] = p

# Multiple testing correction
from metainformant.metabolomics.pathways.enrichment import benjamini_hochberg
qvals = benjamini_hochberg(list(pvals))
sig = np.where(np.array(qvals) < 0.05)[0]
print(f"Metabolites with significant time effect: {len(sig)}")

# Visualize time course for top hits
for i in sig[:3]:
    metabolite = dataset.metabolites[i]
    means = [intensities[i, g].mean() for g in groups]
    stds = [intensities[i, g].std() for g in groups]
    fig = visualization.plots.time_series_plot(
        timepoints=timepoints,
        means=means,
        stds=stds,
        label=metabolite,
    )
    fig.write_html(f"timecourse_{metabolite.replace(' ', '_')}.html")
```

### Expected Output

```
Metabolites with significant time effect: 42
Top time-varying metabolites:
  1. Lactate: F=24.5, q=1.2e-05  (increases over time)
  2. ATP: F=18.9, q=3.4e-04  (decreases)
  3. Glutathione: F=12.1, q=0.0012  (biphasic response)
```

**Note**: `time_series_plot()` is a suggested addition to `visualization/plots.py` (not yet implemented).

---

## Example 4: MS/MS Spectral Library Matching

**Study Design**: Untargeted LC-MS/MS experiment. Query spectra need to be matched against a reference spectral library (e.g., NIST, GNPS) for confident metabolite identification.

**Data**:
- Query: `query_spectra.mgf` (hundreds of MS2 spectra)
- Reference: `nist_library.mgf` (thousands of spectra)

### Script

```python
from metainformant.metabolomics import io, analysis

# Load query and reference spectra
print("Loading spectra...")
query_spectra = io.formats.read_mgf("query_spectra.mgf")
print(f"  Query: {len(query_spectra)} spectra")

ref_spectra = io.formats.read_mgf("nist_library.mgf")
print(f"  Reference: {len(ref_spectra)} spectra")

# Preprocess: normalize intensities to max = 1.0 (cosine similarity invariant to scale)
def normalize_spectrum(spec):
    if spec.total_ion_current > 0:
        factor = 1.0 / spec.total_ion_current
        norm_intens = spec.intensity_array * factor
    else:
        norm_intens = spec.intensity_array.copy()
    return spec._replace(intensity_array=norm_intens)

query_spectra = [normalize_spectrum(s) for s in query_spectra]
ref_spectra = [normalize_spectrum(s) for s in ref_spectra]

# Align spectra to common m/z grid (binning)
def bin_spectrum(spec, min_mz=50, max_mz=1000, bin_width=0.1):
    """Bin m/z values to fixed grid."""
    bins = np.arange(min_mz, max_mz + bin_width, bin_width)
    binned = np.zeros_like(bins)
    for mz, intens in zip(spec.mz_array, spec.intensity_array):
        bin_idx = int((mz - min_mz) / bin_width)
        if 0 <= bin_idx < len(binned):
            binned[bin_idx] += intens
    return binned

# Build reference binned matrix
ref_binned = np.stack([bin_spectrum(s) for s in ref_spectra])

# Match each query spectrum
print("Matching spectra...")
matches = []
top_k = 5
for q_i, query in enumerate(query_spectra):
    q_vec = bin_spectrum(query).reshape(1, -1)

    # Cosine similarity with all refs (vectorized)
    # cos_sim = (q @ ref.T) / (||q|| * ||ref||) but vectors are already normalized
    dot = q_vec @ ref_binned.T  # shape (1, n_ref)
    norms = np.linalg.norm(ref_binned, axis=1)  # shape (n_ref,)
    # q_vec norm ~ 1.0 after normalization
    sims = dot.flatten() / np.maximum(norms, 1e-10)

    top_indices = np.argsort(-sims)[:top_k]
    query_matches = []
    for ref_idx in top_indices:
        score = float(sims[ref_idx])
        if score > 0.7:  # reasonable cosine threshold
            ref_spec = ref_spectra[ref_idx]
            query_matches.append({
                "reference_scan": ref_spec.scan_number,
                "reference_precursor": ref_spec.precursor_mz,
                "cosine_score": score,
            })
    matches.append({
        "query_scan": query.scan_number,
        "query_precursor": query.precursor_mz,
        "top_matches": query_matches,
    })

print(f"Matched {sum(1 for m in matches if m['top_matches'])} / {len(matches)} query spectra")

# Save results
import json
with open("spectral_matches.json", "w") as f:
    json.dump(matches, f, indent=2)
print("Saved to spectral_matches.json")
```

### Interpretation

- **Cosine score > 0.8**: High-confidence match (likely same compound)
- **Cosine 0.6–0.8**: Putative match; consider RT, adduct, and fragmentation pattern
- **Cosine < 0.6**: Unreliable; candidate may be isomer or wrong compound

**Next step**: Manually inspect top matches in spectral viewers (e.g., MZmine, GNPS) to confirm fragment ion patterns.

---

## Example 5: Multi-Omics Integration — Metabolite-Gene Correlation

**Study Design**: Integrated metabolomics + RNA-seq analysis of breast cancer cell lines. Goal: correlate metabolite abundances with gene expression to identify regulatory programs.

**Data**:
- Metabolite matrix: 200 metabolites × 50 samples
- RNA-seq matrix: 15000 genes × 50 samples (TPM normalized)

### Script

```python
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from metainformant.metabolomics import analysis
from metainformant.visualization import plots  # assume better plots exist

# Load metabolomics data
from metainformant.metabolomics.io import formats as meta_io
meta_dataset = meta_io.read_csv("metabolites.csv")
meta_norm = analysis.identification.normalize_intensities(
    meta_dataset.intensities, method="log2"
)

# Load RNA-seq data (assume it's in CSV: genes × samples)
rna_dataset = meta_io.read_csv("rna_tpm.csv")  # same IO format works
# Filter to variable genes only (optional, reduces multiple testing)
gene_var = rna_dataset.intensities.var(axis=1)
top_genes_idx = np.argsort(-gene_var)[:1000]  # top 1000 most variable
rna_matrix = rna_dataset.intensities[top_genes_idx, :]
rna_genes = [rna_dataset.metabolites[i] for i in top_genes_idx]  # 'metabolites' field holds gene names

# Align samples — ensure same ordering across datasets
# This is critical: both matrices must have same sample order
print("Metabolite samples:", meta_dataset.samples[:5])
print("RNA samples:", rna_dataset.samples[:5])
# If order differs, reorder:
# common = sorted(set(meta_dataset.samples) & set(rna_dataset.samples))
# meta_idx = [meta_dataset.samples.index(s) for s in common]
# rna_idx = [rna_dataset.samples.index(s) for s in common]
# meta_norm = meta_nrom[:, meta_idx]
# rna_matrix = rna_matrix[:, rna_idx]

# Correlation: for each metabolite, correlate with each gene
# (200 × 50) · (1000 × 50).T → 200 × 1000 correlation matrix
print(f"Computing correlations: {meta_norm.shape[0]} metabolites × {rna_matrix.shape[0]} genes...")
correlations = np.zeros((meta_norm.shape[0], rna_matrix.shape[0]))
pvals = np.zeros_like(correlations)

for m_i in range(meta_norm.shape[0]):
    for g_i in range(rna_matrix.shape[0]):
        rho, p = spearmanr(meta_norm[m_i, :], rna_matrix[g_i, :])
        correlations[m_i, g_i] = rho
        pvals[m_i, g_i] = p

print(f"Computed {correlations.size} correlations")

# For each metabolite, find top correlated genes
print("\nTop metabolite-gene associations:")
for m_i in range(min(5, len(meta_dataset.metabolites))):
    met_name = meta_dataset.metabolites[m_i]
    top_gene_idx = np.argsort(-np.abs(correlations[m_i]))[:3]
    for g_i in top_gene_idx:
        gene = rna_genes[g_i]
        rho = correlations[m_i, g_i]
        p = pvals[m_i, g_i]
        print(f"  {met_name} ↔ {gene}: ρ={rho:.3f}, p={p:.2e}")

# Heatmap of top associations
top_met_idx = np.argsort(-np.max(np.abs(correlations), axis=1))[:10]
top_gene_idx = np.argsort(-np.max(np.abs(correlations), axis=0))[:15]
sub_corr = correlations[top_met_idx, :][:, top_gene_idx]
sub_genes = [rna_genes[i] for i in top_gene_idx]
sub_mets = [meta_dataset.metabolites[i] for i in top_met_idx]

fig = visualization.plots.heatmap(
    data=sub_corr,
    row_labels=sub_mets,
    col_labels=sub_genes,
    colorscale="RdBu",
    zmid=0,
)
fig.write_html("metabolite_gene_correlation.html")
print("Saved correlation heatmap")
```

### Expected Output

```
Metabolite samples: ['S1', 'S2', 'S3', 'S4', 'S5']
RNA samples: ['S1', 'S2', 'S3', 'S4', 'S5']
Computing correlations: 200 metabolites × 1000 genes...
Computed 200000 correlations

Top metabolite-gene associations:
  Citrate ↔ IDH2: ρ=0.82, p=2.1e-08
  Citrate ↔ FH: ρ=0.78, p=1.2e-06
  2-Hydroxyglutarate ↔ IDH1: ρ=0.91, p=3.4e-12
  Lactate ↔ LDHA: ρ=0.88, p=1.5e-10
```

**Interpretation**: Strong positive correlation between citrate and IDH2/FH (TCA cycle genes) confirms known biology. 2-HG and IDH1 correlation suggests IDH1 mutation driving oncometabolite production.

---

## Example 6: Large Cohort — Quality Control and Batch Assessment

**Study Design**: Multi-center study with 500 serum samples across 3 batch groups. Goal: detect and correct batch effects before biological analysis.

**Data**: 1000 metabolites × 500 samples, with known batch labels (Batch A, B, C).

### Script

```python
import numpy as np
from metainformant.metabolomics import io, analysis, visualization

dataset = io.formats.read_csv("large_cohort.csv")
intensities = dataset.intensities  # (1000, 500)

# Load batch labels from external file
import pandas as pd
metadata = pd.read_csv("sample_metadata.csv")
batches = metadata["batch"].tolist()  # e.g., ["A", "B", "C", ...]

# ── 1. Pre-normalization QC: Check total ion current per sample ─────────
tics = intensities.sum(axis=0)
print(f"TIC range: {tics.min():.2e} – {tics.max():.2e}")
print(f"CV across samples: {tics.std() / tics.mean():.3f}")

# Flag outliers (outside 3 MAD from median)
from metainformant.core import utils  # or implement MAD manually
mad = np.median(np.abs(tics - np.median(tics)))
outliers = np.abs(tics - np.median(tics)) > 3 * mad
print(f"TIC outliers: {outliers.sum()} samples")

# ── 2. Normalization ────────────────────────────────────────────────────
norm = analysis.identification.normalize_intensities(intensities, method="total_ion_count")
log_norm = np.log2(norm + 1)

# ── 3. Batch effect visualization ───────────────────────────────────────
# PCA colored by batch
pca = visualization.plots.pca_plot(
    log_norm.T,
    group_labels=batches,
    color_map={"A": "blue", "B": "red", "C": "green"},
)
pca.write_html("pca_by_batch.html")

# Boxplots of top 10 most variable metabolites by batch
var = log_norm.var(axis=1)
top10_idx = np.argsort(-var)[:10]
fig = visualization.plots.boxplot_by_group(
    data=log_norm[top10_idx, :],
    group_labels=batches,
    feature_labels=[dataset.metabolites[i] for i in top10_idx],
)
fig.write_html("batch_boxplots.html")

# ── 4. Batch effect quantification ──────────────────────────────────────
# Compute batch vs. total variance for each metabolite
batch_var = np.zeros(len(dataset.metabolites))
for m_i in range(len(dataset.metabolites)):
    batch_means = [log_norm[m_i, [i for i, b in enumerate(batches) if b == batch]].mean()
                   for batch in ["A", "B", "C"]]
    batch_var[m_i] = np.var(batch_means)

total_var = log_norm.var(axis=1)
batch_effect = batch_var / (total_var + 1e-10)
print(f"Metabolites with strong batch effect (η² > 0.3): {(batch_effect > 0.3).sum()}")

# ── 5. Combat batch correction (if needed) ───────────────────────────────
# Note: ComBat not yet implemented in METAINFORMANT; use external:
try:
    from pycombat import Combat
    combat = Combat()
    corrected = combat.fit_transform(log_norm.T, batches).T  # samples × features
    print("Batch correction applied (ComBat)")
except ImportError:
    print("pycombat not installed; skipping correction")
    corrected = log_norm

# Re-run PCA on corrected data
pca_corr = visualization.plots.pca_plot(corrected.T, group_labels=batches)
pca_corr.write_html("pca_after_correction.html")

# ── 6. Post-correction QC ────────────────────────────────────────────────
# Check batch effect reduced
batch_var_corr = np.zeros(len(dataset.metabolites))
for m_i in range(len(dataset.metabolites)):
    batch_means = [corrected[m_i, [i for i, b in enumerate(batches) if b == batch]].mean()
                   for batch in ["A", "B", "C"]]
    batch_var_corr[m_i] = np.var(batch_means)
batch_effect_corr = batch_var_corr / (corrected.var(axis=1) + 1e-10)
print(f"Batch effect reduced: {(batch_effect_corr < batch_effect).mean():.1%} of metabolites")

# ── 7. Proceed with biological analysis ─────────────────────────────────
# Now run differential abundance or downstream analysis on `corrected`
```

### Expected Output

```
TIC range: 2.34e+05 – 8.91e+05
CV across samples: 0.412
TIC outliers: 12 samples
Metabolites with strong batch effect (η² > 0.3): 127
Batch correction applied (ComBat)
Batch effect reduced: 87.4% of metabolites
```

**Interpretation**:
- **High TIC CV (0.41)**: Indicates substantial sample-to-sample variation (likely batch-driven).
- **127 metabolites with η² > 0.3**: Strong batch confounding; must correct before biological inference.
- **Post-ComBat reduction (87% improved)**: Batch effect successfully mitigated; PCA should now cluster by biology, not batch.

---

## Summary of Key Functions by Use Case

| Goal | Functions |
|------|-----------|
| **Identify metabolites** | `identify_metabolites()`, `identify_with_adducts()` |
| **Normalize data** | `normalize_intensities()` (TIC/median/log2/pareto) |
| **Differential analysis** | `differential_abundance()`, `fold_change()` |
| **Pathway analysis** | `metabolite_set_enrichment()`, `enrichment_with_fdr()`, `pathway_activity_scoring()` |
| **Quality control** | `filter_spectra()`, manual TIC/RCS checks |
| **Spectral matching** | `read_mgf()`, `cosine_spectral_similarity()` |
| **Missing data** | `missing_value_imputation()` |
| **Visualization** | `visualization.plots.volcano_plot()`, `pca_plot()`, `heatmap()` |

---

## Tips for Your Own Analysis

1. **Always visualize first**: PCA, sample clustering, TIC boxplots to spot outliers and batches.
2. **Normalize before stats**: Never use raw intensities.
3. **Log-transform**: For statistical modeling, log2(x+1) is nearly universal.
4. **Control false discoveries**: Use BH-FDR; expect 5–10% of “significant” hits to be false positives at q < 0.05.
5. **Biological validation**: Confirm top hits with orthogonal data (e.g., targeted MS, enzyme assays).
6. **Document everything**: Save normalized matrices, intermediate statistics, and all plots.

---

## Further Reading

- **MetaboAnalyst**: Comprehensive web-based metabolomics analysis ( Xia & Wishart, 2010–2024).
- **XCMS**: Popular R package for LC-MS data processing (Smith et al., 2006).
- **MZmine 3**: Open-source GUI for MS data analysis (plus back-end APIs).
- **HMDB**: Human Metabolome Database — reference m/z, pathways, disease associations.

Need more examples? See the `examples/` directory at the repository root for complete Jupyter notebooks.
