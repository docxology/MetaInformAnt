# GWAS Module Rules (`metainformant.gwas`)

## Purpose
Genome-wide association studies: variant quality control, association testing, population structure, and visualization.

## Dependencies
- **Required**: `core` (the module imports no other domain package, apart from `rna` in `workflow/rna_quantification.py`)
- **Optional**: External tools (bcftools, GATK) for variant calling

## Source Structure
```
src/metainformant/gwas/
 analysis/
 alignment.py, annotation.py, association.py, benchmarking.py, calling.py, correction.py
 enrichment.py, eqtl.py, gene_annotation_api.py, heritability.py, hwe.py
 ld_decay.py, ld_pruning.py, ldsr.py, mixed_model.py, power.py, prs.py
 quality.py, strain_analysis.py, structure.py, summary_stats.py, utils.py
 data/
 config.py, download.py, expression.py, extraction.py, genome.py
 metadata.py, sra_download.py, traits.py, vcf_utils.py
 finemapping/
 colocalization.py, credible_sets.py, eqtl.py
 heritability/
 estimation.py
 reporting/
 audit.py, html_report.py
 simulation/
 synthetic.py
 validation/
 output_validator.py
 visualization/
 config.py, general.py, _general_impl.py, utils.py
 eqtl_visualization.py, strain_plots.py, structure_plots.py
 genomic/ (genome.py, ld.py, regional.py, variants.py)
 interactive/ (composite.py, finemapping.py, interactive.py, phenotype.py, suite.py)
 population/ (geography.py, population.py, population_admixture.py, population_pca.py)
 statistical/ (comparison.py, effects.py, statistical.py)
 workflow/
 workflow.py, workflow_config.py, workflow_execution.py, rna_quantification.py
```

## Package Management
- **ALWAYS use `uv`** for all Python package management and environment operations
- Use `uv run` to execute commands: `uv run pytest`, `uv run metainformant gwas run --config config/gwas/example.yaml`
- **NEVER use `pip` directly**

## Key Submodules

### Association (`association`)
- Statistical association testing
- Linear and logistic regression models
- Effect size calculation
- P-value computation

**Patterns**:
```python
from metainformant.gwas.analysis import association

# genotype_matrix is variants x samples (0/1/2 alternate-allele counts);
# covariates (optional) is covariates x samples
results = association.run_linear_model_gwas(
    genotype_matrix=genotype_matrix,
    phenotypes=phenotypes,
    variant_info=variant_info,  # optional: [{"id": ..., "chrom": ..., "pos": ...}, ...]
    covariates=covariates,
)
# Returns: list of dicts, one per variant, with beta, se, t_stat, p_value, r_squared, ...

# Single variant: association.association_test_linear(genotypes, phenotypes, covariates=None)
# Binary traits: association.run_logistic_model_gwas / association.association_test_logistic
```

### Correction (`correction`)
- Multiple testing correction
- Bonferroni correction
- False Discovery Rate (FDR) - Benjamini-Hochberg
- Family-wise error rate control

**Patterns**:
```python
from metainformant.gwas.analysis import correction

bonferroni = correction.bonferroni_correction(p_values, alpha=0.05)
# dict: corrected_alpha, significant (one bool per test), n_significant, ...

fdr = correction.fdr_correction(p_values, alpha=0.05, method="bh")
# dict: adjusted_p_values, significant (one bool per test), n_significant, ...

adjusted = correction.adjust_p_values(p_values, method="bonferroni")
# list of adjusted p-values; method is one of "bonferroni", "fdr", "genomic_control", "qvalue"
```

### Quality (`quality`)
- Variant quality control
- Missingness filtering
- Minor allele frequency (MAF) filtering
- Hardy-Weinberg equilibrium testing

### Structure (`structure`)
- Population structure analysis
- Principal Component Analysis (`compute_pca`) and kinship matrices (`compute_kinship_matrix`)
- Population stratification correction (PCs as covariates, or a kinship matrix with `run_mixed_model_gwas`)
- Admixture proportions are only plotted (`admixture_plot` reads ADMIXTURE-format output); estimation is not implemented

### Visualization (`visualization`)
- Manhattan plots
- Q-Q plots
- Regional association plots
- Effect size plots
- Comprehensive visualization suite

**Visualization Submodules**:
- `general.py`: Public entry points (`manhattan_plot`, `qq_plot`, `regional_plot`, `pca_plot`, `kinship_heatmap`, `generate_all_plots`); implemented in `_general_impl.py`
- `config.py`: `PlotStyle`, themes and `style_from_config`
- `statistical/`: Statistical, effect-size and comparison plots
- `genomic/`: Genome-wide, regional, LD and variant-specific plots
- `population/`: Population structure, PCA, admixture and geography plots
- `interactive/`: Interactive Manhattan plots, composite panels, fine-mapping, phenotype plots and the suite runner
- `eqtl_visualization.py`, `strain_plots.py`, `structure_plots.py`: eQTL, strain and structure plots

**Patterns**:
```python
from metainformant.gwas.visualization import general as visualization

# association_results: the list of per-variant dicts from run_linear_model_gwas;
# the entry points of general.py create output/gwas/ when it is missing, but plots such as
# genomic/genome.py's manhattan_plot call plt.savefig directly and need the directory to exist
visualization.manhattan_plot(
    results=association_results,
    output_path="output/gwas/manhattan.png",
    significance_threshold=5e-8,
)
```

### Workflow (`workflow`)
- End-to-end GWAS workflow
- Configuration-based execution
- Step orchestration

### Calling (`calling`)
- Variant calling integration
- bcftools and GATK integration
- Variant calling workflows

### Download (`download`)
- Reference genome download
- Data acquisition utilities
- Genome retrieval

### SRA Download (`sra_download`)
- SRA data download utilities
- Sequence Read Archive integration
- Data retrieval workflows

### Config (`config`)
- Configuration management
- GWASWorkflowConfig class
- Config validation and loading

## Patterns

### I/O Operations
- **CRITICAL**: Always use `metainformant.core.io` for JSON/CSV/TSV operations
- Never use direct `import json` or `import csv`
- Use `io.load_json()`, `io.dump_json()` for all JSON operations

**Pattern**:
```python
from metainformant.core import io

# Save results
io.dump_json(results, "output/gwas/association/results.json")
io.write_csv(df, "output/gwas/variants/variants.csv")
```

### Path Handling
- Always use `metainformant.core.io.paths` utilities for path operations
- Always resolve and validate paths using `paths.expand_and_resolve()` and `paths.is_within()`

**Pattern**:
```python
from metainformant.core.io import paths

# Expand and resolve paths
output_path = paths.expand_and_resolve("output/gwas/analysis/")
```

### Type Hints
- Comprehensive type hints throughout
- Use `from __future__ import annotations` for forward references
- Union types: `str | Path | None`

### Configuration
- Use `GWASWorkflowConfig` for configuration
- Config files in `config/gwas/<name>.yaml`
- Environment prefix: `GWAS_` (e.g., `GWAS_THREADS`, `GWAS_WORK_DIR`)
- Use `metainformant.core.utils.config.load_mapping_from_file()` and `apply_env_overrides()`

**Config Structure**:
```yaml
work_dir: output/gwas/analysis/work
log_dir: output/gwas/analysis/logs
threads: 8
genome:
  reference: "path/to/reference.fa"
variants:
  vcf_path: "data/variants.vcf"
association:
  model: "linear"
  covariates: ["age", "sex"]
correction:
  method: "fdr"
  alpha: 0.05
```

### Output Paths
- Base: `output/gwas/<analysis_type>/`
- Association results: `output/gwas/<analysis>/association/`
- Visualizations: `output/gwas/<analysis>/plots/`
- Quality reports: `output/gwas/<analysis>/qc/`
- Variant calling: `output/gwas/<analysis>/calling/`
- Downloads: `output/gwas/<analysis>/downloads/`

### Statistical Models
- Linear, logistic and mixed-model tests live in `analysis/association.py` and `analysis/mixed_model.py`
- Support both linear and logistic regression
- Handle missing data appropriately
- Report effect sizes and confidence intervals

### Visualization
- All visualization functions accept `output_path` parameter
- Default format: PNG
- Include significance thresholds in plots
- Use consistent color schemes

### Quality Control
- Filter variants before association testing
- Report QC metrics in structured format
- Output QC reports to JSON (via `metainformant.core.io`) and text files

## Environment Variables

- `GWAS_THREADS`: Number of threads
- `GWAS_WORK_DIR`: Working directory override
- `GWAS_LOG_DIR`: Log directory override
- `GWAS_ALPHA`: Significance threshold

## Integration

- **Uses**: `core` (and `rna` for RNA quantification in `workflow/rna_quantification.py`)
- **Used by**: `multiomics` (genomics integration)
- **Integrates with**: `visualization` (plotting), external tools (bcftools, GATK)

## Testing

- **STRICTLY REAL IMPLEMENTATION**: Test real implementations only (see Global Project Rules (root `.cursorrules`) real-implementation policy)
- Test with real VCF files from `data/`
- Test statistical functions with known results
- Write test outputs to `output/gwas/test/` using `tmp_path` fixture
- Test visualization functions (check output files exist)
- Real implementations only - use real file operations and computations

## Documentation

- Document statistical methods and assumptions
- Provide examples with real GWAS data
- Document visualization customization options

## Reference

See Global Project Rules (root `.cursorrules`) for:
- Common directory structure and path handling
- Configuration patterns with env overrides
- Testing policy (real-implementation policy)
- Import patterns and code style
- Documentation guidelines


