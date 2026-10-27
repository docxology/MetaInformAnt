# Agent Directives: gwas

**Context**: Genome-Wide Association Studies (GWAS) module for METAINFORMANT.

## Capabilities

End-to-end GWAS pipelines for social insect genomics: VCF parsing, QC, association testing (linear/logistic/mixed), fine-mapping (SuSiE, credible sets), colocalization, eQTL analysis, heritability estimation (LDSC, GREML), compute-time benchmarking, and comprehensive visualization.

## Subpackages

### analysis/ (22 files)

| File | Key Functions |
|------|---------------|
| `association.py` | `association_test_linear()`, `association_test_logistic()` |
| `benchmarking.py` | `benchmark_subset_run()`, `extrapolate_full_genome_time()`, `scaling_model()` |
| `calling.py` | `call_variants_bcftools()`, `call_variants_gatk()`, `filter_variants_bcftools()`, `merge_vcf_files()` |
| `correction.py` | `bonferroni_correction()`, `fdr_correction()`, `genomic_control()`, `lambda_gc_from_p_values()` |
| `eqtl.py` | `run_eqtl_analysis()` |
| `heritability.py` | `estimate_heritability()` |
| `ld_pruning.py` | `ld_prune()` |
| `mixed_model.py` | `association_test_mixed()`, `run_mixed_model_gwas()` |
| `quality.py` | `parse_vcf_full()`, `apply_qc_filters()`, `check_haplodiploidy()` |
| `structure.py` | `compute_pca()`, `compute_kinship_matrix()`, `estimate_population_structure()` |
| `summary_stats.py` | `write_summary_statistics()`, `write_significant_hits()`, `create_results_summary()` |
| `annotation.py` | `annotate_variants_with_genes()`, `classify_variant_location()` |
| `utils.py` | `compute_r_squared()` |
| `alignment.py` | `align_reads_bwa()`, `align_and_sort_bwa()`, `sort_index_bam()` |
| `enrichment.py` | `fisher_enrichment_test()`, `bh_correction()`, `pathway_enrichment_from_annotations()` |
| `gene_annotation_api.py` | `annotate_top_hits_ncbi()`, `lookup_genes_by_region_ncbi()` |
| `hwe.py` | `hwe_chi2_test()`, `hwe_flag_variants()`, `hwe_distribution_plot()` |
| `ld_decay.py` | `compute_ld_decay()`, `ld_decay_plot()` |
| `ldsr.py` | `compute_ld_scores()`, `ldsr_regression()`, `ldsr_plot()` |
| `power.py` | `compute_power()`, `power_curve()`, `subsample_convergence()`, `jackknife_se()`, `saturation_analysis()` |
| `prs.py` | `clump_variants()`, `compute_prs()`, `prs_r2()`, `prs_full_analysis()` |
| `strain_analysis.py` | strain-level GWAS summaries and Fst utilities |

### data/ (9 files)

| File | Key Functions |
|------|---------------|
| `config.py` | `load_gwas_config()`, `validate_config_parameters()`, `estimate_runtime()` |
| `download.py` | `download_reference_genome()`, `download_variant_database()`, `download_sra_run()` |
| `expression.py` | `ExpressionLoader` — loads Amalgkit kallisto quantification |
| `genome.py` | `normalize_chromosome_name()`, `parse_gff3_genes()` |
| `metadata.py` | `load_sample_metadata()`, `get_population_labels()` |
| `sra_download.py` | `download_sra_experiment()` |
| `extraction.py` | `parse_fastq_filename()`, `build_sample_file_map()`, `extract_fastq()` |
| `traits.py` | `load_traits()`, `load_traits_with_ids()` |
| `vcf_utils.py` | `discover_sample_vcfs()`, `count_variants()`, `bgzip_and_index()`, `merge_vcfs()`, `subsample_vcf()`, `extract_sample_ids()` |

### finemapping/ (3 files)

| File | Key Functions |
|------|---------------|
| `credible_sets.py` | `compute_credible_set()`, `compute_log_bayes_factors()`, `compute_bayes_factors()`, `log_sum_exp()`, `susie_regression()`, `colocalization()`, `colocalization_from_log_bayes_factors()`, `conditional_analysis()`, `annotate_credible_set()` |
| `colocalization.py` | `eqtl_coloc()`, `multi_trait_coloc()`, `compute_clpp()`, `regional_coloc()` |
| `eqtl.py` | `cis_eqtl_scan()`, `trans_eqtl_scan()`, `conditional_eqtl()`, `eqtl_effect_sizes()`, `eqtl_summary_stats()`, `write_eqtl_results()`, `load_transcriptome_variants()` |

### heritability/ (1 file)

| File | Key Functions |
|------|---------------|
| `estimation.py` | `estimate_h2_ldsc()`, `partitioned_h2()`, `genetic_correlation()`, `greml_simple()`, `haseman_elston_regression()`, `compute_liability_h2()` |

### visualization/ (4 subdirectories + top-level files)

| Subdirectory | Purpose |
|-------------|---------|
| `genomic/` | Genome-wide views, LD heatmaps, regional association |
| `population/` | PCA scatter, admixture, geographic distribution |
| `statistical/` | Effect size forest plots, method comparison |
| `interactive/` | Plotly dashboards (`interactive_manhattan()`, `interactive_pca()`, `interactive_volcano()`), fine-mapping plots, composite panels |
| `general.py` | `manhattan_plot()`, `qq_plot()`, `regional_plot()`, `generate_all_plots()` |

### workflow/ (4 files)

| File | Key Functions |
|------|---------------|
| `workflow.py` | Re-export facade for backward compatibility |
| `workflow_config.py` | `GWASWorkflowConfig`, `load_gwas_config()`, `validate_gwas_config()` |
| `workflow_execution.py` | `execute_gwas_workflow()`, `run_gwas()`, `run_multi_trait_gwas()` |
| `rna_quantification.py` | `find_quantification_file()` — adapter to RNA quantification discovery |

### simulation/ (1 file)

| File | Key Functions |
|------|---------------|
| `synthetic.py` | `create_synthetic_data()`, `calculate_power()` |

### validation/ (1 file)

| File | Key Functions |
|------|---------------|
| `output_validator.py` | `validate_pipeline_outputs()` |

### reporting/ (2 files)

| File | Key Functions |
|------|---------------|
| `audit.py` | `sha256_file()`, `checksum_directory()`, `generate_audit_log()` |
| `html_report.py` | HTML report rendering from pipeline outputs |

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for domain data file I/O. Direct stdlib parsing is allowed in core, protocol adapters, subprocess/CLI glue, and narrow parser internals when covered by tests.
- Follow REAL IMPLEMENTATION policy — all tests must use real implementations
- Use `uv` for dependency management
- Pure Python fallbacks for all numpy/scipy operations (graceful degradation)
- All computational complexity models documented (O(n·m), O(n²·m), etc.)
- Association tests analyse complete cases: a negative or non-finite genotype dose is missing, never an allele count. Report the samples analysed as `n_samples` and the samples dropped as `n_missing_dropped`. The logistic test reports whether the fit converged in `converged`.
- A failed test is returned, not raised, and carries NaN for `beta`, `se`, `t_stat` and `p_value` (`z_stat` for the logistic test) and for `r_squared` (linear) or `odds_ratio` (logistic), plus `status`, `error` and `converged=False`. Logistic `status` is `"failed"` or `"error"`; mixed-model rows carry `status` too. NaN is missing, never p = 1 and never an effect of 0. A linear test with a constant phenotype fails with an error message, an exact fit with a nonzero genotype effect reports t of plus or minus infinity and p = 0, and zero residual variance with a zero effect fails.
- `fdr_correction()`, `bonferroni_correction()`, `qvalue_estimation()` and `adjust_p_values()` leave NaN or `None` p-values out of the family (`as_p_value()` decides what is missing). Outputs stay aligned with the input, with NaN at a missing entry, and the dictionaries of `fdr_correction()` and `bonferroni_correction()` report `n_excluded`. `create_results_summary()` counts only the tests that ran in `n_variants_tested` and reports `n_tests_failed`; `compute_comprehensive_summary()` also reports `n_tests_failed`; `per_chromosome_summary()`, `sign_test()`, `subsample_convergence()` and `jackknife_se()` skip failed tests. The writers print a failed test as `NA`. The workflows report `n_tests_failed` in `results["outputs"]` (`execute_gwas_workflow()`), `results["results"]` (`run_gwas()`) and per trait (`run_multi_trait_gwas()`).
- `apply_qc_filters()` raises `ValueError` when the genotype rows differ in length or the number of variant records differs from the number of genotype columns.
- Output tables never invent values: `write_summary_statistics()` and `write_significant_hits()` write `NA` for a missing numeric value and `.` for missing text. `Q_FDR` is the result's own q-value, or else Benjamini-Hochberg over the full set of results; `MAF` is the result's own value, or else computed from the optional `genotypes=` matrix with `summary_stats.minor_allele_frequency()`. `execute_gwas_workflow()` and `run_gwas()` give each result a `maf`, and all three workflows pass the genotypes to the writers, so the MAF column is filled.
- Functions that sample take a `seed` and use their own generator, never the global `random` state: `calculate_ld_decay()`, `compute_ld_decay()`, `ld_decay_plot()`, `subsample_convergence()` and `ldsr_plot()`.
- `run_eqtl_analysis()` is all-vs-all unless `variant_positions` and `gene_positions` are both given; then it tests only variants on the transcript's chromosome with `|position - TSS| <= cis_window` and adds a `distance` column. A `cis_window` without positions raises `ValueError`.
- `download_variant_data()` supports only `data_type="summary_stats"`; `"vcf"` (the default) and `"both"` raise `NotImplementedError` and an unknown type raises `ValueError`.
- `execute_gwas_workflow()` declares its genotype layout as `sample_major`; a helper that infers the layout from the shape assumes an undeclared square matrix is variant-major and warns.
- The QC filter `min_qual` is inclusive and removes variants whose QUAL is missing; `apply_qc_filters()` reports the counts in `qual_filter`. `filter_variants_bcftools()` keeps sites that meet every threshold (`bcftools filter -i`).
- Fine-mapping Bayes factors are Wakefield ABFs that need `standard_errors` or `sample_size`, accumulated in log space; see [finemapping/AGENTS.md](finemapping/AGENTS.md). The `ld_matrix` of `compute_credible_set()` holds LD correlations r, not r-squared (the adjustment squares it), and `regional_coloc()` returns `status: "error"` when every window fails and counts the failed windows in `n_failed_windows`.
- `interactive_manhattan()` reads the chromosome from `chromosome`, `chrom`, `chr`, `CHROM` or `CHR`. A result without a chromosome, a position and a positive p-value is skipped and counted in `n_skipped`, never plotted on a default chromosome.
- Static plot functions in `visualization/general.py` create the missing directories of their output location: the parents of `output_path`, or `output_dir` itself for `generate_all_plots()`. `generate_all_plots()` (implemented in `_general_impl.py`) reads a list of records, a DataFrame or a `.json`, `.tsv`, `.txt` or `.csv` file, leaves out results without a p-value, raises on unsupported input, and closes its figures.

## Related Documentation

- **Module guide**: [../../../docs/gwas/](../../../docs/gwas/) — In-depth usage, architecture, and examples
- **API reference**: [SPEC.md](SPEC.md) — Type signatures, data structures, error codes
- **Core infrastructure**: [../core/AGENTS.md](../core/AGENTS.md) — Shared utilities (logging, config, I/O)
- **Full module index**: [../../../docs/index.md](../../../docs/index.md) — Overview of all METAINFORMANT modules
- **Dna module**: [../dna/AGENTS.md](../dna/AGENTS.md)
- **Phenotype module**: [../phenotype/AGENTS.md](../phenotype/AGENTS.md)
- **Networks module**: [../networks/AGENTS.md](../networks/AGENTS.md)
- **Visualization module**: [../visualization/AGENTS.md](../visualization/AGENTS.md)
