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
import matplotlib.pyplot as plt
from pathlib import Path
from metainformant.metabolomics import (
    io, analysis, pathways, visualization,
)

# ── Load data ──────────────────────────────────────────────────────────
print("Loading urine metabolomics data...")
dataset = io.formats.read_csv("urine_metabolites.csv")
print(f"  Metabolites: {len(dataset.metabolites)}, Samples: {len(dataset.samples)}")
# Output: Metabolites: 500, Samples: 60

# ── Quality control: flag low-signal samples ──────────────────────────
sample_totals = dataset.intensities.sum(axis=0)
low_signal = sample_totals < 1e5
print(f"  Low-signal samples (< 1e5 total): {low_signal.sum()}")
# Remove or investigate low-signal samples if needed
keep = ~low_signal
dataset.intensities = dataset.intensities[:, keep]
dataset.samples = [s for i, s in enumerate(dataset.samples) if keep[i]]

# Sample grouping from metadata (assumes SampleID encodes group), computed after
# the QC filter so that the indices match the remaining columns
# E.g., "D003" = Diabetic patient 3, "C015" = Control 15
sample_groups = ["diabetic" if s.startswith("D") else "control" for s in dataset.samples]
group_a = [i for i, g in enumerate(sample_groups) if g == "diabetic"]
group_b = [i for i, g in enumerate(sample_groups) if g == "control"]
print(f"  Group sizes: diabetic={len(group_a)}, control={len(group_b)}")
# Output (with one low-signal sample in each group): Group sizes: diabetic=29, control=29

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

# Volcano plot: volcano_plot_data() prepares the points (it draws nothing); plot with matplotlib.
# fold_change() is the log2 ratio of group means, so give it the linear-scale TIC-normalized
# matrix, not the log2-transformed one.
log2fc = analysis.identification.fold_change(norm_tic, group_a, group_b)
points = visualization.plots.volcano_plot_data(
    dataset.metabolites, log2fc, qvals, fc_threshold=1.0, p_threshold=0.05
)
fig, ax = plt.subplots(figsize=(6, 5))
ax.scatter(
    [p.log2_fc for p in points],
    [p.neg_log10_p for p in points],
    c=["tab:red" if p.significant else "tab:gray" for p in points],
    s=12,
)
ax.axhline(-np.log10(0.05), linestyle="--", color="black", linewidth=0.8)
ax.set_xlabel("log2 fold change (diabetic / control)")
ax.set_ylabel("-log10 q-value")
fig.savefig("urine_volcano.png", dpi=150, bbox_inches="tight")
print("  Volcano plot saved to urine_volcano.png")

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
# Query set: the metabolites significant at q < 0.05 (use a looser cutoff for broader enrichment).
# background_size is the universe the query was drawn from, here all measured metabolites; by
# default the universe is only the pathway members plus the query, which understates enrichment
# when the measured panel is much larger.
query_metabolites = [dataset.metabolites[i] for i in sig_idx]
enriched = pathways.enrichment.enrichment_with_fdr(
    query_metabolites, pathway_db, background_size=len(dataset.metabolites), fdr_threshold=0.2
)
print(f"  Enriched pathways (FDR <= 0.2): {len(enriched)}")
for r in enriched:
    # fold_enrichment = observed / expected overlap; p_value is the raw (unadjusted) p-value
    print(f"  {r.pathway_name}: {r.overlap}/{r.pathway_size} metabolites, "
          f"fold={r.fold_enrichment:.2f}, p={r.p_value:.4g}")

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

Illustrative values for the study design above (no dataset ships with the repository, so your numbers will differ; the script prints up to five hits per direction, shortened here):

```
Loading urine metabolomics data...
  Metabolites: 500, Samples: 60
  Low-signal samples (< 1e5 total): 2
  Group sizes: diabetic=29, control=29
Normalizing (TIC)...
Log2 transforming for normality...
Computing differential abundance...
  Volcano plot saved to urine_volcano.png
  Significant metabolites (FDR < 0.05): 23
    Upregulated in diabetes: 14
    Downregulated in diabetes: 9

Top upregulated metabolites:
  Citrate: t=5.43, q=0.0003
  2-Hydroxybutyrate: t=4.87, q=0.0012
  BCAA: t=4.21, q=0.0031
  ...

Top downregulated metabolites:
  Hippurate: t=-3.98, q=0.0045
  Creatinine: t=-3.45, q=0.0089
  ...

Pathway enrichment...
  Enriched pathways (FDR <= 0.2): 2
  Glycolysis / Gluconeogenesis: 3/4 metabolites, fold=16.30, p=0.0003318
  Amino Acid Metabolism: 2/3 metabolites, fold=14.49, p=0.005913
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
import numpy as np
import matplotlib.pyplot as plt
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

# PCA to visualize separation (pca_metabolomics takes metabolites × samples and returns
# sample scores plus the explained-variance ratios; plotting is up to you)
scores, explained = visualization.plots.pca_metabolomics(norm_pareto, n_components=2)
fig, ax = plt.subplots(figsize=(5, 4))
for label, idx in (("WT", wt_idx), ("KRAS-mut", mut_idx)):
    ax.scatter(scores[idx, 0], scores[idx, 1], label=label)
ax.set_xlabel(f"PC1 ({explained[0]:.1%})")
ax.set_ylabel(f"PC2 ({explained[1]:.1%})")
ax.legend()
fig.savefig("pca.png", dpi=150, bbox_inches="tight")
print("PCA plot saved")

# Differential abundance (small n → effect size focus).
# Group A = KRAS-mut, group B = WT, so positive t / log2FC mean higher in KRAS-mut.
t, p = analysis.identification.differential_abundance(
    norm_pareto, group_a=mut_idx, group_b=wt_idx
)

# Since n=5 per group, low power → use effect size threshold
# (fold_change() is the log2 ratio of group means: pass linear-scale intensities)
log2fc = analysis.identification.fold_change(dataset.intensities, mut_idx, wt_idx)
large_effect = (np.abs(log2fc) > 1.0) & (p < 0.1)  # exploratory
candidates = np.where(large_effect)[0]
print(f"Candidate oncometabolites (|log2FC|>1, p<0.1): {len(candidates)}")
for i in candidates:
    print(f"  {dataset.metabolites[i]}: log2FC={log2fc[i]:.2f}, p={p[i]:.3f}")

# Pathway activity scoring (uses continuous scores, not just sig/not-sig)
# KEGG_PATHWAYS: dict mapping pathway name -> list of metabolite names (see CONFIGURATION.md)
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

Illustrative values (shortened; your numbers depend on your data):

```
WT: 5, KRAS-mut: 5
PCA plot saved
Candidate oncometabolites (|log2FC|>1, p<0.1): 8
  2-Hydroxyglutarate: log2FC=2.45, p=0.080
  Succinate: log2FC=1.87, p=0.070
  Lactate: log2FC=1.34, p=0.090
  ...

Top pathway activities (mean t-statistic):
  Glutathione metabolism: 1.87 (n=4/8)
  TCA Cycle: 1.23 (n=5/8)
  Glycolysis: 0.87 (n=6/10)
```

**Interpretation**: Elevated 2-hydroxyglutarate (2-HG) is a known oncometabolite from IDH mutations; if KRAS-mut lines also have elevated 2-HG, may indicate epigenetic dysregulation. Elevated lactate and glycolysis activity suggest the Warburg effect (aerobic glycolysis).

---

## Example 3: Time-Course Drug Response

**Study Design**: Treat cells with drug X, collect metabolomics at 0h, 2h, 6h, 24h (n=3 replicates per time point). Goal: identify metabolites with significant time dynamics.

**Data**: 300 metabolites × 12 samples (4 timepoints × 3 replicates).

### Analysis

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import f_oneway
from metainformant.metabolomics import io

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
qvals = np.array(benjamini_hochberg(list(pvals)))
sig = np.where(qvals < 0.05)[0]
print(f"Metabolites with significant time effect: {len(sig)}")

# Top hits by q-value
top = sig[np.argsort(qvals[sig])[:3]]
print("Top time-varying metabolites:")
for rank, i in enumerate(top, 1):
    print(f"  {rank}. {dataset.metabolites[i]}: F={f_stats[i]:.1f}, q={qvals[i]:.2g}")

# Visualize time course for the top hits
for i in top:
    metabolite = dataset.metabolites[i]
    means = [intensities[i, g].mean() for g in groups]
    stds = [intensities[i, g].std() for g in groups]
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.errorbar(range(len(timepoints)), means, yerr=stds, marker="o", capsize=3)
    ax.set_xticks(range(len(timepoints)), timepoints)
    ax.set_xlabel("Time after drug X")
    ax.set_ylabel("Intensity")
    ax.set_title(metabolite)
    fig.savefig(f"timecourse_{metabolite.replace(' ', '_')}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
```

### Expected Output

Illustrative values (your numbers depend on your data):

```
Metabolites with significant time effect: 42
Top time-varying metabolites:
  1. Lactate: F=233.2, q=1.2e-05
  2. ATP: F=83.1, q=0.00034
  3. Glutathione: F=53.8, q=0.0012
```

**Note**: `metabolomics.visualization.plots` has no time-series helper, so the loop above draws each time course with matplotlib directly.

---

## Example 4: MS/MS Spectral Library Matching

**Study Design**: Untargeted LC-MS/MS experiment. Query spectra need to be matched against a reference spectral library (e.g., NIST, GNPS) for confident metabolite identification.

**Data**:
- Query: `query_spectra.mgf` (hundreds of MS2 spectra)
- Reference: `nist_library.mgf` (thousands of spectra)

### Script

```python
import json
from dataclasses import replace

import numpy as np
from metainformant.metabolomics import io, analysis

# Load query and reference spectra
print("Loading spectra...")
query_spectra = io.formats.read_mgf("query_spectra.mgf")
print(f"  Query: {len(query_spectra)} spectra")

ref_spectra = io.formats.read_mgf("nist_library.mgf")
print(f"  Reference: {len(ref_spectra)} spectra")

# Preprocess: scale each spectrum to unit total ion current
# (cosine similarity ignores overall scale, so this only keeps intensities comparable).
# MassSpectrum is a dataclass: use dataclasses.replace() to copy it with new intensities.
def normalize_spectrum(spec):
    if spec.total_ion_current > 0:
        norm_intens = spec.intensity_array / spec.total_ion_current
    else:
        norm_intens = spec.intensity_array.copy()
    return replace(spec, intensity_array=norm_intens)

query_spectra = [normalize_spectrum(s) for s in query_spectra]
ref_spectra = [normalize_spectrum(s) for s in ref_spectra]

# Align spectra to a common m/z grid (binning): cosine_spectral_similarity()
# compares two intensity vectors defined on the same grid
def bin_spectrum(spec, min_mz=50.0, max_mz=1000.0, bin_width=0.1):
    """Sum peak intensities into fixed m/z bins."""
    edges = np.arange(min_mz, max_mz + bin_width, bin_width)
    binned, _ = np.histogram(spec.mz_array, bins=edges, weights=spec.intensity_array)
    return binned

# Build reference binned matrix
ref_binned = np.stack([bin_spectrum(s) for s in ref_spectra])

# Match each query spectrum
print("Matching spectra...")
matches = []
top_k = 5
for query in query_spectra:
    q_vec = bin_spectrum(query)

    # Cosine similarity with every reference spectrum
    sims = np.array(
        [analysis.identification.cosine_spectral_similarity(q_vec, ref_vec) for ref_vec in ref_binned]
    )

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
from metainformant.visualization import plots as viz  # shared matplotlib plotting helpers

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
# meta_norm = meta_norm[:, meta_idx]
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

# (shared heatmap helper from metainformant.visualization.plots; it returns a matplotlib Axes)
sub_df = pd.DataFrame(sub_corr, index=sub_mets, columns=sub_genes)
ax = viz.expression_heatmap(sub_df, cmap="RdBu_r", vmin=-1, vmax=1, cbar_label="Spearman rho")
ax.set_title("Metabolite-gene correlation")
ax.figure.savefig("metabolite_gene_correlation.png", dpi=150, bbox_inches="tight")
print("Saved correlation heatmap")
```

### Expected Output

Illustrative values (the script prints the top 3 genes for each of the first 5 metabolites; the list is shortened here):

```
Metabolite samples: ['S1', 'S2', 'S3', 'S4', 'S5']
RNA samples: ['S1', 'S2', 'S3', 'S4', 'S5']
Computing correlations: 200 metabolites × 1000 genes...
Computed 200000 correlations

Top metabolite-gene associations:
  Citrate ↔ IDH2: ρ=0.820, p=3.23e-13
  Citrate ↔ FH: ρ=0.780, p=2.45e-11
  2-Hydroxyglutarate ↔ IDH1: ρ=0.910, p=5.56e-20
  Lactate ↔ LDHA: ρ=0.880, p=3.91e-17
  ...
Saved correlation heatmap
```

**Interpretation**: Strong positive correlation between citrate and IDH2/FH (TCA cycle genes) confirms known biology. 2-HG and IDH1 correlation suggests IDH1 mutation driving oncometabolite production.

---

## Example 6: Large Cohort — Quality Control and Batch Assessment

**Study Design**: Multi-center study with 500 serum samples across 3 batch groups. Goal: detect and correct batch effects before biological analysis.

**Data**: 1000 metabolites × 500 samples, with known batch labels (Batch A, B, C).

### Script

```python-snippet
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from metainformant.metabolomics import io, analysis, visualization
from metainformant.quality.batch.detection import correct_batch_combat

dataset = io.formats.read_csv("large_cohort.csv")
intensities = dataset.intensities  # (1000, 500)

# Load batch labels from external file (row order must match dataset.samples)
metadata = pd.read_csv("sample_metadata.csv")
batches = metadata["batch"].tolist()  # e.g., ["A", "B", "C", ...]

# ── 1. Pre-normalization QC: Check total ion current per sample ─────────
tics = intensities.sum(axis=0)
print(f"TIC range: {tics.min():.2e} – {tics.max():.2e}")
print(f"CV across samples: {tics.std() / tics.mean():.3f}")

# Flag outliers (outside 3 MAD from median)
mad = np.median(np.abs(tics - np.median(tics)))
outliers = np.abs(tics - np.median(tics)) > 3 * mad
print(f"TIC outliers: {outliers.sum()} samples")

# ── 2. Normalization ────────────────────────────────────────────────────
norm = analysis.identification.normalize_intensities(intensities, method="total_ion_count")
log_norm = np.log2(norm + 1)

# ── 3. Batch effect visualization ───────────────────────────────────────
# PCA colored by batch. pca_metabolomics() takes metabolites × samples and returns
# sample scores and explained-variance ratios; the scatter is drawn with matplotlib.
def plot_pca_by_batch(data, labels, path):
    scores, explained = visualization.plots.pca_metabolomics(data, n_components=2)
    fig, ax = plt.subplots(figsize=(5, 4))
    for label in sorted(set(labels)):
        pick = [i for i, b in enumerate(labels) if b == label]
        ax.scatter(scores[pick, 0], scores[pick, 1], s=12, label=f"Batch {label}")
    ax.set_xlabel(f"PC1 ({explained[0]:.1%})")
    ax.set_ylabel(f"PC2 ({explained[1]:.1%})")
    ax.legend()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)

plot_pca_by_batch(log_norm, batches, "pca_by_batch.png")

# Boxplots of top 10 most variable metabolites by batch
var = log_norm.var(axis=1)
top10_idx = np.argsort(-var)[:10]
batch_names = ["A", "B", "C"]
fig, axes = plt.subplots(2, 5, figsize=(15, 6))
for ax, m_i in zip(axes.ravel(), top10_idx):
    per_batch = [log_norm[m_i, [i for i, b in enumerate(batches) if b == name]] for name in batch_names]
    ax.boxplot(per_batch)
    ax.set_xticks(range(1, len(batch_names) + 1), batch_names)
    ax.set_title(dataset.metabolites[m_i], fontsize=8)
fig.savefig("batch_boxplots.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# ── 4. Batch effect quantification ──────────────────────────────────────
# Compute batch vs. total variance for each metabolite (this ratio equals η² when the
# batches have the same number of samples)
batch_var = np.zeros(len(dataset.metabolites))
for m_i in range(len(dataset.metabolites)):
    batch_means = [log_norm[m_i, [i for i, b in enumerate(batches) if b == batch]].mean()
                   for batch in ["A", "B", "C"]]
    batch_var[m_i] = np.var(batch_means)

total_var = log_norm.var(axis=1)
batch_effect = batch_var / (total_var + 1e-10)
print(f"Metabolites with strong batch effect (η² > 0.3): {(batch_effect > 0.3).sum()}")

# ── 5. Batch correction (if needed) ──────────────────────────────────────
# correct_batch_combat() takes samples × features and returns the same shape. It is a
# simplified ComBat-style correction: each batch is centered on the grand mean and its
# standard deviation is rescaled toward the overall one (the ratio is shrunk halfway toward 1,
# a fixed shrinkage of 0.5). It is not the full parametric empirical-Bayes ComBat and has no
# covariate argument, so do not use it when batch is confounded with the biological groups
# you want to compare.
corrected = correct_batch_combat(log_norm.T, batches).T  # back to metabolites × samples
print("Batch correction applied (ComBat-style)")

# Re-run PCA on corrected data
plot_pca_by_batch(corrected, batches, "pca_after_correction.png")

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

Illustrative values for the study design above (the first four lines depend on your data):

```
TIC range: 2.34e+05 – 8.91e+05
CV across samples: 0.412
TIC outliers: 12 samples
Metabolites with strong batch effect (η² > 0.3): 127
Batch correction applied (ComBat-style)
Batch effect reduced: 100.0% of metabolites
```

**Interpretation**:
- **High TIC CV (0.41)**: Indicates substantial sample-to-sample variation (likely batch-driven).
- **127 metabolites with η² > 0.3**: Strong batch confounding; must correct before biological inference.
- **Post-correction reduction (100% of metabolites)**: `correct_batch_combat()` sets every batch mean to the grand mean, so the between-batch term drops to ~0 by construction; check that PCA now clusters by biology, not batch.

---

## Summary of Key Functions by Use Case

| Goal | Functions |
|------|-----------|
| **Identify metabolites** | `identify_metabolites()`, `identify_with_adducts()` |
| **Normalize data** | `normalize_intensities()` (TIC/median/log2/pareto) |
| **Differential analysis** | `differential_abundance()`, `fold_change()` |
| **Pathway analysis** | `metabolite_set_enrichment()`, `enrichment_with_fdr()`, `pathway_activity_scoring()` |
| **Quality control** | `filter_spectra()`, manual TIC checks (Examples 1 and 6) |
| **Spectral matching** | `read_mgf()`, `cosine_spectral_similarity()` |
| **Missing data** | `missing_value_imputation()` |
| **Batch correction** | `metainformant.quality.batch.detection.correct_batch_combat()` (simplified ComBat-style) |
| **Visualization (plot-ready data; draw with matplotlib)** | `visualization.plots.volcano_plot_data()`, `pca_metabolomics()`, `intensity_heatmap_data()`, `prepare_spectrum_plot()`, `detect_chromatographic_peaks()`, `retention_time_alignment_data()` |

---

## Tips for Your Own Analysis

1. **Always visualize first**: PCA, sample clustering, TIC boxplots to spot outliers and batches.
2. **Normalize before stats**: Raw intensities are rarely suitable input for statistical tests.
3. **Log-transform**: For statistical modeling, log2(x+1) is nearly universal.
4. **Control false discoveries**: Use BH-FDR; at q < 0.05, expect at most about 5% of the “significant” hits to be false positives.
5. **Biological validation**: Confirm top hits with orthogonal data (e.g., targeted MS, enzyme assays).
6. **Document everything**: Save normalized matrices, intermediate statistics, and all plots.

---

## Further Reading

- **MetaboAnalyst**: Comprehensive web-based metabolomics analysis ( Xia & Wishart, 2010–2024).
- **XCMS**: Popular R package for LC-MS data processing (Smith et al., 2006).
- **MZmine 3**: Open-source GUI for MS data analysis (plus back-end APIs).
- **HMDB**: Human Metabolome Database — reference m/z, pathways, disease associations.

The repository ships no Jupyter notebooks and no metabolomics-specific example scripts; the walkthroughs above are the metabolomics examples. For the cross-module pieces see [INTEGRATION.md](INTEGRATION.md).
