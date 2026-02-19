# Specification: gwas

## 🎯 Scope

Genome-Wide Association Studies (GWAS) module for METAINFORMANT. Provides end-to-end GWAS pipelines for social insect genomics, including VCF parsing, QC, association testing (linear, logistic, mixed model), fine-mapping (SuSiE, credible sets), colocalization, eQTL analysis, heritability estimation (LDSC, GREML), and comprehensive visualization.

## 🧱 Architecture

- **Dependency Level**: Domain
- **Component Type**: Source Code
- **Sub-packages**: `analysis`, `data`, `finemapping`, `heritability`, `visualization`, `workflow`

## 🔌 API Exports

### analysis

| Module | Key Functions |
|--------|---------------|
| `association` | `association_test_linear()`, `association_test_logistic()` |
| `mixed_model` | `association_test_mixed()`, `run_mixed_model_gwas()` |
| `quality` | `parse_vcf_full()`, `apply_qc_filters()` |
| `correction` | `bonferroni_correction()`, `fdr_correction()`, `genomic_control()` |
| `structure` | `compute_pca()`, `compute_kinship_matrix()` |
| `ld_pruning` | `ld_prune()` |
| `heritability` | `estimate_heritability()` |
| `calling` | `merge_vcf_files()`, `index_vcf()` |
| `annotation` | `annotate_variants_with_genes()`, `classify_variant_location()` |
| `summary_stats` | `write_summary_stats()`, `extract_significant_hits()` |

### data

| Module | Key Functions |
|--------|---------------|
| `download` | `download_reference_genome()`, `download_variant_database()`, `download_sra_run()` |
| `genome` | `normalize_chromosome_name()`, `parse_gff3_genes()` |
| `metadata` | `load_sample_metadata()`, `get_population_labels()` |

### finemapping

| Module | Key Functions |
|--------|---------------|
| `credible_sets` | `compute_credible_set()`, `susie_regression()`, `conditional_analysis()` |
| `colocalization` | `colocalization()`, `multi_trait_coloc()`, `eqtl_coloc()`, `compute_clpp()` |
| `eqtl` | `cis_eqtl_scan()`, `trans_eqtl_scan()`, `eqtl_effect_sizes()` |

### heritability

| Module | Key Functions |
|--------|---------------|
| `estimation` | `estimate_h2_ldsc()`, `partitioned_h2()`, `genetic_correlation()`, `greml_simple()` |

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

- **537 GWAS tests** collected across unit, integration, and end-to-end suites
- **11/11 end-to-end tests pass** (`tests/test_gwas_end_to_end.py`)
- Zero-mock policy: all tests use real functional methods
