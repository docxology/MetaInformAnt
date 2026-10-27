# Specification: gwas

## Scope

Genome-Wide Association Studies (GWAS) module for METAINFORMANT. Provides end-to-end GWAS pipelines for social insect genomics, including VCF parsing, QC, association testing (linear, logistic, mixed model), fine-mapping (SuSiE, credible sets), colocalization, eQTL analysis, heritability estimation (LDSC, GREML), compute-time benchmarking, expression data loading, and comprehensive visualization.

## Architecture

- **Dependency Level**: Domain
- **Component Type**: Source Code
- **Sub-packages**: `analysis`, `data`, `finemapping`, `heritability`, `visualization`, `workflow`

## API Exports

### analysis

| Module | Key Functions |
|--------|---------------|
| `association` | `association_test_linear()`, `association_test_logistic()` |
| `benchmarking` | `benchmark_subset_run()`, `extrapolate_full_genome_time()`, `scaling_model()` |
| `mixed_model` | `association_test_mixed()`, `run_mixed_model_gwas()` |
| `quality` | `parse_vcf_full()`, `apply_qc_filters()`, `check_haplodiploidy()` |
| `correction` | `bonferroni_correction()`, `fdr_correction()`, `genomic_control()`, `lambda_gc_from_p_values()` |
| `structure` | `compute_pca()`, `compute_kinship_matrix()`, `estimate_population_structure()` |
| `ld_pruning` | `ld_prune()` |
| `heritability` | `estimate_heritability()` |
| `calling` | `call_variants_bcftools()`, `call_variants_gatk()`, `filter_variants_bcftools()`, `merge_vcf_files()`, `index_vcf()` |
| `annotation` | `annotate_variants_with_genes()`, `classify_variant_location()` |
| `summary_stats` | `write_summary_statistics()`, `write_significant_hits()`, `create_results_summary()` |
| `eqtl` | `run_eqtl_analysis()` |

### data

| Module | Key Functions |
|--------|---------------|
| `download` | `download_reference_genome()`, `download_variant_database()`, `download_sra_run()` |
| `expression` | `ExpressionLoader.load_amalgkit_quant()` |
| `genome` | `normalize_chromosome_name()`, `parse_gff3_genes()` |
| `metadata` | `load_sample_metadata()`, `get_population_labels()` |
| `config` | `load_gwas_config()`, `validate_config_parameters()`, `estimate_runtime()` |
| `sra_download` | `download_sra_experiment()` |

### finemapping

| Module | Key Functions |
|--------|---------------|
| `credible_sets` | `compute_credible_set()`, `compute_log_bayes_factors()`, `compute_bayes_factors()`, `log_sum_exp()`, `susie_regression()`, `colocalization()`, `colocalization_from_log_bayes_factors()`, `conditional_analysis()`, `annotate_credible_set()` |
| `colocalization` | `eqtl_coloc()`, `multi_trait_coloc()`, `compute_clpp()`, `regional_coloc()` |
| `eqtl` | `cis_eqtl_scan()`, `trans_eqtl_scan()`, `conditional_eqtl()`, `eqtl_effect_sizes()`, `eqtl_summary_stats()` |

### heritability

| Module | Key Functions |
|--------|---------------|
| `estimation` | `estimate_h2_ldsc()`, `partitioned_h2()`, `genetic_correlation()`, `greml_simple()`, `haseman_elston_regression()`, `compute_liability_h2()` |

### visualization

| Module | Key Functions |
|--------|---------------|
| `general` | `manhattan_plot()`, `qq_plot()`, `regional_plot()`, `pca_plot()`, `generate_all_plots()` |
| `genomic/` | Genome-wide views, LD heatmaps, regional association |
| `population/` | PCA scatter, admixture, geographic distribution |
| `statistical/` | Effect size forest plots, method comparison |
| `interactive/` | Plotly dashboards (`interactive_manhattan()`, `interactive_pca()`, `interactive_volcano()`), fine-mapping plots, composite panels |

### workflow

| Module | Key Functions |
|--------|---------------|
| `workflow` | Re-exports from `workflow_config` and `workflow_execution` |
| `workflow_config` | `GWASWorkflowConfig`, `load_gwas_config()`, `validate_gwas_config()` |
| `workflow_execution` | `execute_gwas_workflow()`, `run_gwas()`, `run_multi_trait_gwas()` |

## Behavioral Contracts

### Association tests

- `association_test_linear()` and `association_test_logistic()` analyse complete cases. A sample is excluded when its genotype dose is negative or non-finite, when any covariate value is non-finite, and (linear test only) when its phenotype is missing (`NaN`, `None` or non-numeric). Results carry `n_samples` (samples analysed), `n_missing_dropped` and `converged`. A genotype/phenotype length mismatch raises `ValueError`.
- Linear: `status` is `"success"` or `"error"`. P-values come from Student's t distribution with `n - p` degrees of freedom (`n` complete cases, `p` fitted parameters including the intercept), with a normal approximation when SciPy is unavailable. A degenerate design (fewer than 3 complete cases, no residual degrees of freedom, a constant genotype or a singular design matrix) returns `status: "error"` with `p_value: 1.0`, `converged: False` and an `error` message.
- Logistic: `status` is `"success"`, `"failed"` (fewer than 2 cases or 2 controls among the complete cases) or `"error"` (the statsmodels fit failed). `converged` is statsmodels' convergence flag; a fit that did not converge keeps `status: "success"` with an unreliable estimate and p-value. A phenotype that is not 0 or 1, or a covariate without one value per sample, raises `ValueError`.

### Quality control

- `apply_qc_filters(vcf_input, qc_config=None)` runs the QUAL, minor allele frequency, missing-data and Hardy-Weinberg filters in that order; a variant must pass every enabled filter. `qc_config` keys: `min_maf` (default 0.01), `max_missing` (0.1), `min_hwe_p` (1e-6) and `min_qual` (optional; absent or `None` disables the QUAL filter).
- `min_qual` is inclusive (QUAL >= `min_qual`). A variant whose QUAL is missing (`.` in the VCF, `None` once parsed) is removed whenever `min_qual` is set. The result's `qual_filter` entry holds `min_qual`, `n_removed` and `n_missing_qual`. A non-numeric `min_qual` raises `ValueError`, and so does parsed data without one variant record per genotype column when `min_qual` is set.
- The result holds `status`, `num_variants_before`, `num_variants_after` and `filtered_data`, plus `qc_filters` (the thresholds applied) when any variant passes the MAF and missing-data filters, and `qual_filter` when `min_qual` is set.
- `filter_variants_bcftools(vcf_path, output_vcf, filters)` joins the conditions `QUAL>=min_qual`, `F_MISSING<=max_missing` and `DP>=min_depth` (for the keys present in `filters`) with `&&` and passes them to `bcftools filter -i`, so sites meeting every condition are kept.

### Summary statistics

- `write_summary_statistics(results, variant_info, output_path, *, genotypes=None)` and `write_significant_hits(results, variant_info, output_path, threshold=5e-8, *, genotypes=None)` write a tab-separated table with the header `CHR POS SNP REF ALT BETA SE P Q_FDR N MAF`. `write_significant_hits()` keeps results with `p_value < threshold`, sorted by p-value.
- A missing or non-finite numeric value is written as `NA` and a missing text field (`CHR`, `SNP`, `REF`, `ALT`) as `.`; a list of ALT alleles is joined with commas. BETA and SE use `.6g`, P and Q_FDR `.6e`, MAF `.4f`.
- `Q_FDR` is the result's own finite `q_value`, `q_fdr` or `Q_FDR`; otherwise it is the Benjamini-Hochberg q-value over every result with a p-value in [0, 1], computed over all of `results` before `write_significant_hits()` selects its hits.
- `MAF` is the result's own finite `maf` or `MAF`; otherwise it is computed from the matching row of `genotypes` (variants x samples; alternate-allele dosages 0, 1 and 2, any other code missing), and is `NA` when no genotype is called or `genotypes` is not given. `N` is `n_samples` or `N`.
- `ValueError` is raised when `genotypes` and `results` differ in length, and, for `write_summary_statistics()`, when `results` and `variant_info` differ in length.

### Fine-mapping

Approximate Bayes factors, credible sets and colocalization follow the model in [finemapping/SPEC.md](finemapping/SPEC.md).

### Interactive plots and plot output

- `interactive_manhattan(assoc_results, output_file=None, gene_annotations=None, significance_threshold=5e-8, title="Interactive Manhattan Plot")` reads the chromosome from the first present of `chromosome`, `chrom`, `chr`, `CHROM` and `CHR`, the position from `position`, `pos` or `POS`, and the p-value from `p_value`, `pvalue` or `P`. A result lacking a chromosome, a position or a positive finite p-value is skipped and counted in `n_skipped`; nothing is defaulted. It returns `status` (`"skipped"` when no result is usable), `output_path`, `n_variants`, `n_skipped` and `interactive` (False when Plotly is unavailable and an HTML table is written instead). `interactive_volcano()` reads the chromosome from the same keys.
- Static plot functions in `visualization/general.py` create the missing parent directories of `output_path` (`output_dir` for `generate_all_plots()`).

## Testing

- **60 GWAS test files** covering unit, integration, and end-to-end suites
- **11/11 end-to-end tests pass** (`tests/gwas/test_gwas_end_to_end.py`)
- Real-implementation policy: all tests use real functional methods (zero
  mocks; real CSV/ZIP files, real scipy-verified statistics, real bcftools
  subprocesses where the tool is installed — skipped otherwise)
- Compute-time benchmarking tests validate scaling model math

### Module test coverage (Round-4 audit 2026-09-01; Round-5 additions 2026-09-05)

Zero-mocks suites added for the previously untested modules:

| Module | Test file | Notes |
|--------|-----------|-------|
| `analysis/enrichment` | `test_gwas_enrichment.py` | Fisher p cross-checked against scipy; BH FDR known-vector tests; pure-Python hypergeometric SF verified against `scipy.stats.hypergeom.sf` (regression: log-probability sum bug fixed 2026-09-01) |
| `analysis/power` | `test_gwas_power.py` | NCP closed-form, jackknife SE structure, saturation fit (marginal-gain normalization bug fixed 2026-09-01) |
| `analysis/prs` | `test_gwas_prs.py` | C+T clumping, PRS scoring, R² against constructed causal data |
| `analysis/hwe` | `test_gwas_hwe.py` | chi-square closed-form + erfc p-value identity; sample-major flagging |
| `analysis/ld_decay` | `test_gwas_ld_decay.py` | r² vs distance, window/chromosome constraints, seed determinism |
| `analysis/ldsr` | `test_gwas_ldsr.py` | LD scores; h²/intercept recovery from synthetic χ² |
| `analysis/strain_analysis` | `test_gwas_strain_analysis.py` | Fst estimators (Weir-Cockerham degenerate-sample crash fixed 2026-09-01), private variants incl. fixed differences (detection gap fixed 2026-09-01) |
| `data/traits` + `data/extraction` | `test_gwas_traits_extraction.py` | real CSV files and real ZIP archives |
| `data/vcf_utils` | `test_gwas_vcf_utils.py` | real bcftools subprocesses (skip when absent) |
| `reporting/audit` + `validation/output_validator` | `test_gwas_audit_validator.py` | real hashes, real on-disk output trees |
| `simulation/synthetic` | `test_gwas_synthetic.py` | analytical power closed-form; deterministic synthetic VCF/phenotype generation |
| `analysis/quality` + `analysis/correction` | `test_gwas_quality.py`, `test_gwas_correction.py`, `test_gwas_workflow_execution_helpers.py`, `test_gwas_visualization_comparison.py` | Round-5 fixes (2026-09-05): λ_GC computed via χ²-median genomic-control everywhere (was `median(p)/0.456` or `-2·ln(p)/1.386` in workflow/viz helpers); `write_filtered_vcf` sample-major genotype orientation regression-tested; `analyze_genetic_architecture` inverted P-column guard fixed |
