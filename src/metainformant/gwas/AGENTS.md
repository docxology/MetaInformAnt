# Agent Directives: gwas

**Context**: Genome-Wide Association Studies (GWAS) module for METAINFORMANT.

## Capabilities

End-to-end GWAS pipelines for social insect genomics: VCF parsing, QC, association testing (linear/logistic/mixed), fine-mapping (SuSiE, credible sets), colocalization, eQTL analysis, heritability estimation (LDSC, GREML), compute-time benchmarking, and comprehensive visualization.

## Subpackages

### analysis/ (13 files)

| File | Key Functions |
|------|---------------|
| `association.py` | `association_test_linear()`, `association_test_logistic()` |
| `benchmarking.py` | `benchmark_subset_run()`, `extrapolate_full_genome_time()`, `scaling_model()` |
| `calling.py` | `call_variants_bcftools()`, `call_variants_gatk()`, `merge_vcf_files()` |
| `correction.py` | `bonferroni_correction()`, `fdr_correction()`, `genomic_control()` |
| `eqtl.py` | `run_eqtl_analysis()` |
| `heritability.py` | `estimate_heritability()` |
| `ld_pruning.py` | `ld_prune()` |
| `mixed_model.py` | `association_test_mixed()`, `run_mixed_model_gwas()` |
| `quality.py` | `parse_vcf_full()`, `apply_qc_filters()`, `check_haplodiploidy()` |
| `structure.py` | `compute_pca()`, `compute_kinship_matrix()`, `estimate_population_structure()` |
| `summary_stats.py` | `write_summary_statistics()`, `extract_significant_hits()` |
| `annotation.py` | `annotate_variants_with_genes()`, `classify_variant_location()` |
| `utils.py` | `compute_r_squared()` |

### data/ (7 files)

| File | Key Functions |
|------|---------------|
| `config.py` | `load_gwas_config()`, `validate_config_parameters()`, `estimate_runtime()` |
| `download.py` | `download_reference_genome()`, `download_variant_database()`, `download_sra_run()` |
| `expression.py` | `ExpressionLoader` — loads Amalgkit kallisto quantification |
| `genome.py` | `normalize_chromosome_name()`, `parse_gff3_genes()` |
| `metadata.py` | `load_sample_metadata()`, `get_population_labels()` |
| `sra_download.py` | `batch_download_sra()`, `download_sra_experiment()`, `find_sra_data_by_phenotype()` |

### finemapping/ (3 files)

| File | Key Functions |
|------|---------------|
| `credible_sets.py` | `compute_credible_set()`, `susie_regression()`, `colocalization()`, `conditional_analysis()` |
| `colocalization.py` | `eqtl_coloc()`, `multi_trait_coloc()`, `compute_clpp()`, `regional_coloc()` |
| `eqtl.py` | `cis_eqtl_scan()`, `trans_eqtl_scan()`, `conditional_eqtl()`, `eqtl_effect_sizes()` |

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
| `interactive/` | Plotly dashboards, composite panels |
| `general.py` | `manhattan_plot()`, `qq_plot()`, `regional_plot()`, `generate_all_plots()` |

### workflow/ (3 files)

| File | Key Functions |
|------|---------------|
| `workflow.py` | Re-export facade for backward compatibility |
| `workflow_config.py` | `GWASWorkflowConfig`, `load_gwas_config()`, `validate_gwas_config()` |
| `workflow_execution.py` | `execute_gwas_workflow()`, `run_gwas()`, `run_multi_trait_gwas()` |

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for file operations — never `import json` directly
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management
- Pure Python fallbacks for all numpy/scipy operations (graceful degradation)
- All computational complexity models documented (O(n·m), O(n²·m), etc.)
