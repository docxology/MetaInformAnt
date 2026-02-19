# Specification: gwas

## 🎯 Scope

Genome-Wide Association Studies (GWAS) module for METAINFORMANT. Provides end-to-end GWAS pipelines for social insect genomics, including VCF parsing, QC, association testing (linear, logistic, mixed model), fine-mapping (SuSiE, credible sets), colocalization, eQTL analysis, heritability estimation (LDSC, GREML), compute-time benchmarking, expression data loading, and comprehensive visualization.

## 🧱 Architecture

- **Dependency Level**: Domain
- **Component Type**: Source Code
- **Sub-packages**: `analysis`, `data`, `finemapping`, `heritability`, `visualization`, `workflow`

## 🔌 API Exports

### analysis

| Module | Key Functions |
|--------|---------------|
| `association` | `association_test_linear()`, `association_test_logistic()` |
| `benchmarking` | `benchmark_subset_run()`, `extrapolate_full_genome_time()`, `scaling_model()` |
| `mixed_model` | `association_test_mixed()`, `run_mixed_model_gwas()` |
| `quality` | `parse_vcf_full()`, `apply_qc_filters()`, `check_haplodiploidy()` |
| `correction` | `bonferroni_correction()`, `fdr_correction()`, `genomic_control()` |
| `structure` | `compute_pca()`, `compute_kinship_matrix()`, `estimate_population_structure()` |
| `ld_pruning` | `ld_prune()` |
| `heritability` | `estimate_heritability()` |
| `calling` | `call_variants_bcftools()`, `call_variants_gatk()`, `merge_vcf_files()`, `index_vcf()` |
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
| `sra_download` | `batch_download_sra()`, `download_sra_experiment()`, `find_sra_data_by_phenotype()` |

### finemapping

| Module | Key Functions |
|--------|---------------|
| `credible_sets` | `compute_credible_set()`, `susie_regression()`, `colocalization()`, `conditional_analysis()` |
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
| `interactive/` | Plotly dashboards, composite panels |

### workflow

| Module | Key Functions |
|--------|---------------|
| `workflow` | Re-exports from `workflow_config` and `workflow_execution` |
| `workflow_config` | `GWASWorkflowConfig`, `load_gwas_config()`, `validate_gwas_config()` |
| `workflow_execution` | `execute_gwas_workflow()`, `run_gwas()`, `run_multi_trait_gwas()` |

## 🧪 Testing

- **41 GWAS test files** covering unit, integration, and end-to-end suites
- **11/11 end-to-end tests pass** (`tests/test_gwas_end_to_end.py`)
- Zero-mock policy: all tests use real functional methods
- Compute-time benchmarking tests validate scaling model math
