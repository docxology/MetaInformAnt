# Analysis

RNA-seq expression analysis including normalization, differential expression, quality control, batch effect detection, cross-species comparison, and protein integration.

## Contents

| File | Purpose |
|------|---------|
| `expression_core.py` | Count normalization (CPM, TPM, RPKM, quantile, median-ratio) and filtering |
| `expression_analysis.py` | Differential expression (DESeq2-like, t-test, Wilcoxon), PCA, and volcano data |
| `qc_metrics.py` | Sample/gene QC metrics, outlier detection, library complexity, saturation curves |
| `qc_filtering.py` | Batch effect detection, GC bias, length bias, and QC report generation |
| `cross_species.py` | Ortholog mapping, expression conservation, divergence matrices, cross-species PCA |
| `statistics_contract.py` | Predeclared analysis-provenance records, descriptive/inferential role separation, BH-FDR helper with gated inferential wrapper, and fail-closed orthology/species-tree invariants |
| `atlas_plots.py` | Atlas-style figures: species x tissue tau heatmap, per-orthogroup cross-species profile small multiples, tau-by-orthology-class strip plots (descriptive statistics only) |
| `protein_integration.py` | Translation efficiency, protein abundance prediction, ribosome profiling |
| `validation.py` | Pipeline validation: per-sample status checks and end-to-end reports |

## Key Functions

| Function | Description |
|----------|-------------|
| `normalize_counts()` | Normalize raw counts by CPM, TPM, RPKM, quantile, or median-ratio |
| `filter_low_expression()` | Remove genes below a minimum count threshold |
| `differential_expression()` | DE analysis with method selection and p-value adjustment |
| `pca_analysis()` | Principal component analysis on expression matrices |
| `compute_sample_metrics()` | Per-sample statistics: total counts, detected genes, complexity |
| `detect_outlier_samples()` | Flag outlier samples using median deviation |
| `detect_batch_effects()` | Identify confounding batch structure in expression data |
| `build_ortholog_map()` | Load one-to-one ortholog mappings between species |
| `compute_expression_conservation()` | Spearman correlation of orthologs across species |
| `compute_profile_conservation()` | Per-gene Spearman profile correlation across species over explicitly matched tissue labels |
| `summarize_profile_conservation()` | Per-gene descriptive aggregate of pairwise profile correlations |
| `compute_tpm_distribution_summary()` | Per-species/per-tissue TPM quantile summary artifacts |
| `compute_per_tissue_completeness()` | Per-tissue measured-species/gene counts (sampling-dependence record for cross-species tau) |
| `compute_tissue_overlap_summary()` | Matched-tissue counts per species pair |
| `calculate_translation_efficiency()` | Estimate translation efficiency from RNA and protein data using `ratio` or `correlation` |
| `predict_protein_abundance_from_rna()` | Predict protein abundance from RNA with the implemented `linear` method |
| `validate_all_samples()` | Check pipeline completion status for every sample |
| `AnalysisProvenance` | Frozen predeclared record: analysis id, estimand, replicate unit, seed, resampling count, null model, multiplicity family/method, tested-feature count, role, software versions |
| `validate_analysis_provenance()` | Fail-closed validation of a provenance record; role-conditional fields (family/method/feature count) must be `None`/`'not-applicable'` exactly for descriptive records and declared exactly for inferential ones |
| `render_analysis_provenance_block()` | Render additive `analysis_provenance_*: value` lines for `analysis_summary.txt`; descriptive lanes render `not-applicable` |
| `benjamini_hochberg_fdr()` | Benjamini-Hochberg FDR adjustment with fail-closed input validation |
| `declared_inferential_bh_fdr()` | GATED inferential path: requires `evidence_manifest_frozen=True` and a validated `analysis_role="inferential"` BH-FDR contract whose family size matches; returns raw plus adjusted p-values |
| `result_role()` | Return a result's declared role; refuses unlabeled output |
| `validate_orthology_profile_invariants()` | Fail-closed orthology x species presence checks: unique labels, explicit 0/1 states, declared per-orthogroup species coverage |
| `validate_species_tree_invariants()` | Fail-closed species-tree checks (Newick or nested dict): unique leaves, named taxa, declared rootedness provenance, optional bifurcating-root structural check |

## Usage

```python
from metainformant.rna.analysis.expression_core import normalize_counts, filter_low_expression
from metainformant.rna.analysis.expression_analysis import differential_expression, pca_analysis
from metainformant.rna.analysis.atlas_plots import plot_tau_heatmap, plot_orthogroup_small_multiples, plot_tau_orthology_strips

# Demo figures (synthetic data): docs/rna/figures/atlas_demo/
# Regenerate with: MPLBACKEND=Agg uv run python scripts/rna/demo_atlas_plots.py
from metainformant.rna.analysis.qc_metrics import compute_sample_metrics
from metainformant.rna.analysis.protein_integration import calculate_translation_efficiency

normalized = normalize_counts(raw_counts, method="tpm", gene_lengths=lengths)
de_results = differential_expression(counts, groups, method="deseq2_like")
pca_result = pca_analysis(normalized, n_components=3)
rna_expression = normalized.T  # samples x genes for protein integration
protein_abundance = rna_expression.copy()
translation_efficiency = calculate_translation_efficiency(rna_expression, protein_abundance, method="ratio")
```

## Validation Notes

- Cross-species fingerprint outputs (divergence matrix, feature-resampling stability) carry `attrs["role"] = "descriptive"` and never contain p-values; the only inferential path is the gated `declared_inferential_bh_fdr()`, which refuses to run without a frozen evidence manifest and a validated inferential contract. Provenance validation and the orthology/tree invariants fail closed on missing or placeholder declarations.
- Count/QC matrices must be numeric, finite, and non-negative where raw counts are expected.
- Batch labels must include every expression sample; labels for extra samples are ignored.
- GC content values for matched genes must be in `[0, 1]`; matched gene lengths must be positive.
- RNA-protein integration filters NaN measurements deterministically and rejects unsupported methods instead of returning silent empty results.
