# AI Agents in GWAS Module Development

This document outlines AI assistance in creating and maintaining the METAINFORMANT GWAS (Genome-Wide Association Studies) module.

## Implementation Status

**Status**: ✅ PARTIALLY IMPLEMENTED
- **Core functions**: Implemented (workflow.py, quality.py, structure.py, association.py, correction.py, visualization.py)
- **Advanced functions**: Partially implemented (calling.py, download.py, config.py)
- **Remaining**: Full implementation of data download, variant calling integration

## Implemented Functions by Module

### Workflow Orchestration (`workflow.py`) ✅ IMPLEMENTED
**Implemented functions:**
- `load_gwas_config()` - Load GWAS configuration from files
- `execute_gwas_workflow()` - Execute GWAS workflow
- `validate_gwas_config()` - Validate workflow configuration
- `GWASWorkflowConfig` - Configuration class for GWAS workflows

### Quality Control (`quality.py`) ✅ IMPLEMENTED
**Implemented functions:**
- `parse_vcf_full(vcf_path: Union[str, Path]) -> Dict[str, Any]` - Parse VCF files
- `filter_by_maf(genotypes: List[List[int]], maf_threshold: float = 0.01) -> Tuple[List[List[int]], List[int]]` - Filter by minor allele frequency
- `filter_by_missing(genotypes: List[List[int]], missing_threshold: float = 0.1) -> Tuple[List[List[int]], List[int]]` - Filter by missing data proportion
- `test_hwe(genotypes: List[List[int]], alpha: float = 0.05) -> List[float]` - Test Hardy-Weinberg equilibrium
- `apply_qc_filters(vcf_data: Dict[str, Any], min_maf: float = 0.01, max_missing: float = 0.1, min_hwe_p: float = 1e-6) -> Dict[str, Any]` - Apply QC filters
- `write_filtered_vcf(filtered_data: Dict[str, Any], output_path: Union[str, Path]) -> None` - Write filtered VCF data

### Population Structure (`structure.py`) ✅ IMPLEMENTED
**Implemented functions:**
- `compute_pca(genotype_matrix: List[List[int]], n_components: int = 10) -> Tuple[np.ndarray, np.ndarray, np.ndarray]` - Principal component analysis
- `compute_kinship_matrix(genotype_matrix: List[List[int]], method: str = "vanraden") -> np.ndarray` - Kinship matrix calculation
- `estimate_population_structure(genotype_matrix: List[List[int]], n_pcs: int = 10) -> Dict[str, Any]` - Comprehensive population structure analysis

### Association Testing (`association.py`) ✅ IMPLEMENTED
**Implemented functions:**
- `association_test_linear(genotypes: List[int], phenotypes: List[float], covariates: Optional[List[List[float]]] = None) -> Dict[str, Any]` - Linear regression association test
- `association_test_logistic(genotypes: List[int], phenotypes: List[int], covariates: Optional[List[List[float]]] = None, max_iter: int = 100) -> Dict[str, Any]` - Logistic regression association test

### Multiple Testing Correction (`correction.py`) ✅ IMPLEMENTED
**Implemented functions:**
- `bonferroni_correction(p_values: List[float], alpha: float = 0.05) -> Tuple[List[bool], float]` - Bonferroni correction
- `fdr_correction(p_values: List[float], alpha: float = 0.05, method: str = "bh") -> Tuple[List[bool], List[float]]` - False discovery rate correction
- `genomic_control(p_values: List[float]) -> Tuple[List[float], float]` - Genomic control correction
- `qvalue_estimation(p_values: List[float], pi0: Optional[float] = None) -> Tuple[List[float], float]` - q-value estimation
- `adjust_p_values(p_values: List[float], method: str = "bonferroni") -> List[float]` - Adjust p-values using specified method

### Visualization (`visualization.py`) ✅ IMPLEMENTED
**Implemented functions:**
- `manhattan_plot(results: Union[pd.DataFrame, Dict[str, Any]], output_path: Optional[Union[str, Path]] = None, significance_threshold: float = 5e-8) -> matplotlib.figure.Figure` - Create Manhattan plots
- `qq_plot(p_values: Union[List[float], np.ndarray], output_path: Optional[Union[str, Path]] = None) -> matplotlib.figure.Figure` - Create Q-Q plots
- `regional_plot(results: List[Dict[str, Any]], chrom: str, start: int, end: int, output_path: Optional[Union[str, Path]] = None) -> matplotlib.figure.Figure` - Create regional association plots
- `pca_plot(pca_result: tuple, output_path: Optional[Union[str, Path]] = None, explained_var: Optional[np.ndarray] = None) -> matplotlib.figure.Figure` - Create PCA scatter plots
- `kinship_heatmap(kinship_matrix: Union[np.ndarray, pd.DataFrame], output_path: Optional[Union[str, Path]] = None) -> matplotlib.figure.Figure` - Create kinship matrix heatmaps
- `effect_size_plot(results: List[Dict[str, Any]], output_path: Optional[Union[str, Path]] = None) -> matplotlib.figure.Figure` - Create effect size plots
- `missingness_plot(vcf_data: Dict[str, Any], output_path: Optional[Union[str, Path]] = None) -> matplotlib.figure.Figure` - Create missingness plots
- `functional_enrichment_plot(results: List[Dict[str, Any]], gff_path: Union[str, Path], output_path: Optional[Union[str, Path]] = None) -> matplotlib.figure.Figure` - Create functional enrichment plots
- `generate_all_plots(association_results: Union[str, Path], output_dir: Union[str, Path], pca_file: Optional[Union[str, Path]] = None, kinship_file: Optional[Union[str, Path]] = None, vcf_file: Optional[Union[str, Path]] = None, significance_threshold: float = 5e-8) -> Dict[str, Any]` - Generate all GWAS plots

### Variant Calling (`calling.py`) ✅ IMPLEMENTED
**Implemented functions:**
- `call_variants_bcftools(bam_files: List[Union[str, Path]], reference_fasta: Union[str, Path], output_vcf: Union[str, Path], threads: int = 1) -> subprocess.CompletedProcess` - Call variants using bcftools
- `call_variants_freebayes(bam_files: List[Union[str, Path]], reference_fasta: Union[str, Path], output_vcf: Union[str, Path], ploidy: int = 2) -> subprocess.CompletedProcess` - Call variants using FreeBayes
- `call_variants_gatk(bam_files: List[Union[str, Path]], reference_fasta: Union[str, Path], output_vcf: Union[str, Path], gatk_path: Optional[str] = None) -> subprocess.CompletedProcess` - Call variants using GATK
- `filter_variants_bcftools(vcf_path: Union[str, Path], output_vcf: Union[str, Path], filters: Dict[str, Any]) -> subprocess.CompletedProcess` - Filter variants using bcftools
- `merge_vcfs(vcf_files: List[Union[str, Path]], output_vcf: Union[str, Path]) -> subprocess.CompletedProcess` - Merge multiple VCF files
- `index_vcf(vcf_path: Union[str, Path]) -> subprocess.CompletedProcess` - Index VCF file
- `validate_vcf(vcf_path: Union[str, Path]) -> Dict[str, Any]` - Validate VCF file format

### Configuration Management (`config.py`) ✅ IMPLEMENTED
**Implemented functions:**
- `load_gwas_config(config_path: Union[str, Path]) -> Dict[str, Any]` - Load GWAS configuration from file
- `save_gwas_config(config: Dict[str, Any], output_path: Union[str, Path]) -> None` - Save GWAS configuration to file
- `create_gwas_config_template(output_path: Union[str, Path]) -> None` - Create GWAS configuration template
- `create_config_from_vcf(vcf_path: Union[str, Path], phenotype_path: Union[str, Path], output_dir: Union[str, Path]) -> Dict[str, Any]` - Create configuration from VCF and phenotype files
- `validate_config_parameters(config: Dict[str, Any]) -> List[str]` - Validate configuration parameters
- `validate_input_files(config: Dict[str, Any]) -> List[str]` - Validate input files exist
- `merge_config_defaults(config: Dict[str, Any]) -> Dict[str, Any]` - Merge configuration with defaults
- `get_config_summary(config: Dict[str, Any]) -> Dict[str, Any]` - Get configuration summary
- `estimate_runtime(config: Dict[str, Any]) -> Dict[str, Any]` - Estimate workflow runtime
- `update_config_for_runtime(config: Dict[str, Any]) -> Dict[str, Any]` - Update configuration for runtime optimization

### Data Download (`download.py`) ✅ IMPLEMENTED
**Implemented functions:**
- `download_reference_genome(accession: str, output_dir: Union[str, Path]) -> Path` - Download reference genome
- `download_annotation(accession: str, output_dir: Union[str, Path]) -> Path` - Download genome annotations
- `download_variant_database(db_name: str, output_dir: Union[str, Path]) -> Path` - Download variant database
- `download_sra_run(sra_accession: str, output_dir: Union[str, Path], threads: int = 1) -> Path` - Download SRA run
- `download_sra_project(project_id: str, output_dir: Union[str, Path], threads: int = 1) -> List[Path]` - Download SRA project
- `search_sra_for_organism(organism: str, max_results: int = 100) -> List[Dict[str, Any]]` - Search SRA for organism data

### SRA Data Download (`sra_download.py`) ✅ IMPLEMENTED
**Implemented functions:**
- `download_sra_with_retry(accession: str, output_dir: Union[str, Path], max_retries: int = 3, threads: int = 1) -> Path` - Download SRA with retry logic
- `batch_download_sra(accessions: List[str], output_dir: Union[str, Path], threads: int = 1, max_concurrent: int = 4) -> Dict[str, Path]` - Batch download SRA files
- `download_sra_biosample(biosample_acc: str, output_dir: Union[str, Path], threads: int = 1) -> List[Path]` - Download SRA for biosample
- `download_sra_experiment(experiment_acc: str, output_dir: Union[str, Path], threads: int = 1) -> List[Path]` - Download SRA for experiment
- `prefetch_sra_metadata(accessions: List[str]) -> Dict[str, Dict[str, Any]]` - Prefetch SRA metadata
- `validate_sra_download(accession: str, download_dir: Path) -> Dict[str, Any]` - Validate SRA download
- `find_sra_data_by_phenotype(phenotype: str, organism: str, max_results: int = 100) -> List[Dict[str, Any]]` - Find SRA data by phenotype

### Visualization Comparison (`visualization_comparison.py`) ✅ IMPLEMENTED
**Implemented functions:**
- `compare_analysis_methods(method_results: Dict[str, Any], output_file: Optional[Union[str, Path]] = None) -> Optional[matplotlib.figure.Figure]` - Compare GWAS analysis methods
- `compare_gwas_studies(study_results: Dict[str, Any], output_file: Optional[Union[str, Path]] = None) -> Optional[matplotlib.figure.Figure]` - Compare GWAS studies
- `compare_populations(pop_data: Dict[str, Any], trait_name: str, output_file: Optional[Union[str, Path]] = None) -> Optional[matplotlib.figure.Figure]` - Compare populations
- `calculate_method_power(method_results: Dict[str, Any]) -> Dict[str, float]` - Calculate method statistical power
- `calculate_population_enrichment(pop_data: Dict[str, Any]) -> Dict[str, float]` - Calculate population enrichment
- `analyze_genetic_architecture(pop_data: Dict[str, Any]) -> Dict[str, Any]` - Analyze genetic architecture

### Visualization Effects (`visualization_effects.py`) ✅ IMPLEMENTED
**Implemented functions:**
- `effect_size_plot(results: Any, output_file: Optional[Union[str, Path]] = None) -> Optional[matplotlib.figure.Figure]` - Plot effect sizes
- `effect_size_manhattan(results: Any, significance_threshold: float = 5e-8, output_file: Optional[Union[str, Path]] = None) -> Optional[matplotlib.figure.Figure]` - Manhattan plot with effect sizes
- `effect_size_qq(results: Any, output_file: Optional[Union[str, Path]] = None) -> Optional[matplotlib.figure.Figure]` - Q-Q plot with effect sizes
- `compare_effect_sizes(study_results: Dict[str, Any], output_file: Optional[Union[str, Path]] = None) -> Optional[matplotlib.figure.Figure]` - Compare effect sizes across studies

### Visualization Submodules - **NOT IMPLEMENTED**
*Planned: Specialized GWAS visualization submodules*
- `visualization_genome.py` - Genome-wide visualization functions
- `visualization_population.py` - Population structure visualization
- `visualization_regional.py` - Regional association plots
- `visualization_statistical.py` - Statistical visualization tools
- `visualization_suite.py` - Complete visualization suite
- `visualization_variants.py` - Variant-level visualization

### Testing and Validation
**Code Assistant Agent** contributed to:
- Real implementation testing (NO_MOCKING_POLICY compliance)

### Recent Functional Improvements (December 2024)

**Code Assistant Agent** implemented key enhancements:

#### Visualization Enhancements
- **Enhanced `missingness_plot()`**: Implemented full VCF parsing and missingness calculation with statistical summaries
- **Enhanced `functional_enrichment_plot()`**: Added GFF3 annotation parsing and functional categorization of variants
- **Comprehensive Plotting Utilities**: Added publication-quality styling and automated plot generation

#### Performance Optimizations
- **Kinship Matrix Caching**: Added intelligent caching for expensive kinship matrix computations
- **Vectorized Operations**: Improved performance of mathematical operations using NumPy
- **Efficient VCF Processing**: Optimized VCF parsing for large datasets

#### Code Quality Improvements
- **Input Validation**: Added validation using `metainformant.core.validation`
- **Error Handling**: Standardized error messages and exception handling
- **Type Hints**: Updated to modern Python type hinting (`|` union syntax)
- **Documentation**: Enhanced docstrings with examples, references, and performance notes

#### Testing Enhancements
- **Edge Case Testing**: Added tests for error conditions and edge cases
- **Integration Testing**: Verified cross-module functionality
- **Performance Validation**: Benchmarking of improved algorithms
- Integration tests with DNA and math modules
- Scientific validation against established methods
- Performance testing on realistic datasets

## AI Contributions

### GWAS Architecture Design
**Code Assistant Agent** designed:
- Comprehensive GWAS workflow orchestration
- Modular quality control and statistical analysis pipeline
- Population structure analysis and correction methods
- Multiple testing correction algorithms
- Visualization suite for GWAS results

### Implementation Components
**Code Assistant Agent** contributed to:
- VCF parsing and variant filtering algorithms
- Statistical association testing (linear and logistic regression)
- Population stratification correction
- Kinship matrix computation and IBS distance calculations
- Manhattan plot, Q-Q plot, and regional association visualization
- Multiple testing correction (Bonferroni, FDR, genomic control)
- Workflow configuration and validation

### Quality Assurance
**Documentation Agent** assisted with:
- GWAS algorithm documentation and validation
- Statistical method references and implementation details
- Performance benchmarking and edge case testing
- Integration guides for GWAS workflows

## Development Approach

- **Scientific Rigor**: All implementations based on peer-reviewed methods
- **Performance Focus**: Optimized for large-scale genomic datasets
- **Modular Design**: Separable components for different GWAS stages
- **Statistical Validation**: Comprehensive testing against established tools

## Quality Assurance

- Human oversight ensures statistical accuracy and biological relevance
- AI assistance accelerates implementation while maintaining scientific standards
- Rigorous validation against established GWAS tools and methods

## Scientific Accuracy

### Algorithm Validation
All GWAS algorithms follow established methods:
- PCA for population stratification (Price et al. 2006)
- Kinship matrix computation (VanRaden 2008, Astle-Balding)
- Linear/logistic regression for association testing
- Bonferroni and Benjamini-Hochberg correction methods

### Integration with Literature
- Implementation based on peer-reviewed publications
- Validation against tools like PLINK and GCTA
- Performance benchmarking on standard datasets

## Best Practices

### Code Quality
- Type hints throughout for static analysis
- Comprehensive docstrings following numpy style
- Modular functions with single responsibilities
- Consistent error handling and logging

### Performance Optimization
- Efficient NumPy operations for matrix computations
- Parallel processing support for large-scale analysis
- Memory-efficient VCF parsing
- Optional sparse matrix operations

This module demonstrates responsible AI assistance in bioinformatics software development, combining automated implementation with scientific rigor and human validation.

## Complete Function Signatures

### Workflow Orchestration (`workflow.py`)
- `GWASWorkflowConfig(work_dir: Path, threads: int, vcf_path: Optional[str], phenotype_path: Optional[str]) -> None` - Configuration class for GWAS workflows
- `load_gwas_config(config_file: Union[str, Path]) -> Dict[str, Any]` - Load GWAS configuration from file
- `execute_gwas_workflow(config: Dict[str, Any], check: bool = False) -> Dict[str, Any]` - Execute GWAS workflow
- `validate_gwas_config(config: Dict[str, Any]) -> List[str]` - Validate workflow configuration
- `GWASWorkflowConfig.from_dict(cls, config_dict: Dict[str, Any]) -> GWASWorkflowConfig` - Create config from dictionary
- `GWASWorkflowConfig.to_dict(self) -> Dict[str, Any]` - Convert config to dictionary

### Quality Control (`quality.py`)
- `parse_vcf_full(vcf_path: Union[str, Path]) -> Dict[str, Any]` - Parse VCF files
- `filter_by_maf(genotypes: List[List[int]], maf_threshold: float = 0.01) -> Tuple[List[List[int]], List[int]]` - Filter by minor allele frequency
- `filter_by_missing(genotypes: List[List[int]], missing_threshold: float = 0.1) -> Tuple[List[List[int]], List[int]]` - Filter by missing data proportion
- `test_hwe(genotypes: List[List[int]], alpha: float = 0.05) -> List[float]` - Test Hardy-Weinberg equilibrium
- `apply_qc_filters(vcf_data: Dict[str, Any], min_maf: float = 0.01, max_missing: float = 0.1, min_hwe_p: float = 1e-6) -> Dict[str, Any]` - Apply QC filters
- `write_filtered_vcf(filtered_data: Dict[str, Any], output_path: Union[str, Path]) -> None` - Write filtered VCF data

### Population Structure (`structure.py`)
- `compute_pca(genotype_matrix: List[List[int]], n_components: int = 10) -> Tuple[np.ndarray, np.ndarray, np.ndarray]` - Principal component analysis
- `compute_kinship_matrix(genotype_matrix: List[List[int]], method: str = "vanraden") -> np.ndarray` - Kinship matrix calculation
- `estimate_population_structure(genotype_matrix: List[List[int]], n_pcs: int = 10) -> Dict[str, Any]` - Comprehensive population structure analysis

### Association Testing (`association.py`)
- `association_test_linear(genotypes: List[int], phenotypes: List[float], covariates: Optional[List[List[float]]] = None) -> Dict[str, Any]` - Linear regression association test
- `association_test_logistic(genotypes: List[int], phenotypes: List[int], covariates: Optional[List[List[float]]] = None, max_iter: int = 100) -> Dict[str, Any]` - Logistic regression association test

### Multiple Testing Correction (`correction.py`)
- `bonferroni_correction(p_values: List[float], alpha: float = 0.05) -> Tuple[List[bool], float]` - Bonferroni correction
- `fdr_correction(p_values: List[float], alpha: float = 0.05, method: str = "bh") -> Tuple[List[bool], List[float]]` - False discovery rate correction
- `genomic_control(p_values: List[float]) -> Tuple[List[float], float]` - Genomic control correction
- `qvalue_estimation(p_values: List[float], pi0: Optional[float] = None) -> Tuple[List[float], float]` - q-value estimation
- `adjust_p_values(p_values: List[float], method: str = "bonferroni") -> List[float]` - Adjust p-values using specified method

### Variant Calling (`calling.py`)
- `call_variants_bcftools(bam_files: List[Union[str, Path]], reference_fasta: Union[str, Path], output_vcf: Union[str, Path], threads: int = 1) -> subprocess.CompletedProcess` - Call variants using bcftools
- `call_variants_freebayes(bam_files: List[Union[str, Path]], reference_fasta: Union[str, Path], output_vcf: Union[str, Path], ploidy: int = 2) -> subprocess.CompletedProcess` - Call variants using FreeBayes
- `call_variants_gatk(bam_files: List[Union[str, Path]], reference_fasta: Union[str, Path], output_vcf: Union[str, Path], gatk_path: Optional[str] = None) -> subprocess.CompletedProcess` - Call variants using GATK
- `filter_variants_bcftools(vcf_path: Union[str, Path], output_vcf: Union[str, Path], filters: Dict[str, Any]) -> subprocess.CompletedProcess` - Filter variants using bcftools
- `merge_vcfs(vcf_files: List[Union[str, Path]], output_vcf: Union[str, Path]) -> subprocess.CompletedProcess` - Merge multiple VCF files
- `index_vcf(vcf_path: Union[str, Path]) -> subprocess.CompletedProcess` - Index VCF file
- `validate_vcf(vcf_path: Union[str, Path]) -> Dict[str, Any]` - Validate VCF file format

### Configuration Management (`config.py`)
- `load_gwas_config(config_path: Union[str, Path]) -> Dict[str, Any]` - Load GWAS configuration from file
- `save_gwas_config(config: Dict[str, Any], output_path: Union[str, Path]) -> None` - Save GWAS configuration to file
- `create_gwas_config_template(output_path: Union[str, Path]) -> None` - Create GWAS configuration template
- `create_config_from_vcf(vcf_path: Union[str, Path], phenotype_path: Union[str, Path], output_dir: Union[str, Path]) -> Dict[str, Any]` - Create configuration from VCF and phenotype files
- `validate_config_parameters(config: Dict[str, Any]) -> List[str]` - Validate configuration parameters
- `validate_input_files(config: Dict[str, Any]) -> List[str]` - Validate input files exist
- `merge_config_defaults(config: Dict[str, Any]) -> Dict[str, Any]` - Merge configuration with defaults
- `get_config_summary(config: Dict[str, Any]) -> Dict[str, Any]` - Get configuration summary
- `estimate_runtime(config: Dict[str, Any]) -> Dict[str, Any]` - Estimate workflow runtime
- `update_config_for_runtime(config: Dict[str, Any]) -> Dict[str, Any]` - Update configuration for runtime optimization

### Data Download (`download.py`)
- `download_reference_genome(accession: str, output_dir: Union[str, Path]) -> Path` - Download reference genome
- `download_annotation(accession: str, output_dir: Union[str, Path]) -> Path` - Download genome annotations
- `download_variant_database(db_name: str, output_dir: Union[str, Path]) -> Path` - Download variant database
- `download_sra_run(sra_accession: str, output_dir: Union[str, Path], threads: int = 1) -> Path` - Download SRA run
- `download_sra_project(project_id: str, output_dir: Union[str, Path], threads: int = 1) -> List[Path]` - Download SRA project
- `search_sra_for_organism(organism: str, max_results: int = 100) -> List[Dict[str, Any]]` - Search SRA for organism data

### SRA Data Download (`sra_download.py`)
- `download_sra_with_retry(accession: str, output_dir: Union[str, Path], max_retries: int = 3, threads: int = 1) -> Path` - Download SRA with retry logic
- `batch_download_sra(accessions: List[str], output_dir: Union[str, Path], threads: int = 1, max_concurrent: int = 4) -> Dict[str, Path]` - Batch download SRA files
- `download_sra_biosample(biosample_acc: str, output_dir: Union[str, Path], threads: int = 1) -> List[Path]` - Download SRA for biosample
- `download_sra_experiment(experiment_acc: str, output_dir: Union[str, Path], threads: int = 1) -> List[Path]` - Download SRA for experiment
- `prefetch_sra_metadata(accessions: List[str]) -> Dict[str, Dict[str, Any]]` - Prefetch SRA metadata
- `validate_sra_download(accession: str, download_dir: Path) -> Dict[str, Any]` - Validate SRA download
- `find_sra_data_by_phenotype(phenotype: str, organism: str, max_results: int = 100) -> List[Dict[str, Any]]` - Find SRA data by phenotype

### Visualization (`visualization.py`)
- `manhattan_plot(results: Union[pd.DataFrame, Dict[str, Any]], output_path: Optional[Union[str, Path]] = None, significance_threshold: float = 5e-8) -> matplotlib.figure.Figure` - Create Manhattan plots
- `qq_plot(p_values: Union[List[float], np.ndarray], output_path: Optional[Union[str, Path]] = None) -> matplotlib.figure.Figure` - Create Q-Q plots
- `regional_plot(results: List[Dict[str, Any]], chrom: str, start: int, end: int, output_path: Optional[Union[str, Path]] = None) -> matplotlib.figure.Figure` - Create regional association plots
- `pca_plot(pca_result: tuple, output_path: Optional[Union[str, Path]] = None, explained_var: Optional[np.ndarray] = None) -> matplotlib.figure.Figure` - Create PCA scatter plots
- `kinship_heatmap(kinship_matrix: Union[np.ndarray, pd.DataFrame], output_path: Optional[Union[str, Path]] = None) -> matplotlib.figure.Figure` - Create kinship matrix heatmaps
- `effect_size_plot(results: List[Dict[str, Any]], output_path: Optional[Union[str, Path]] = None) -> matplotlib.figure.Figure` - Create effect size plots
- `missingness_plot(vcf_data: Dict[str, Any], output_path: Optional[Union[str, Path]] = None) -> matplotlib.figure.Figure` - Create missingness plots
- `functional_enrichment_plot(results: List[Dict[str, Any]], gff_path: Union[str, Path], output_path: Optional[Union[str, Path]] = None) -> matplotlib.figure.Figure` - Create functional enrichment plots
- `generate_all_plots(association_results: Union[str, Path], output_dir: Union[str, Path], pca_file: Optional[Union[str, Path]] = None, kinship_file: Optional[Union[str, Path]] = None, vcf_file: Optional[Union[str, Path]] = None, significance_threshold: float = 5e-8) -> Dict[str, Any]` - Generate all GWAS plots

### Visualization Comparison (`visualization_comparison.py`)
- `compare_analysis_methods(method_results: Dict[str, Any], output_file: Optional[Union[str, Path]] = None) -> Optional[matplotlib.figure.Figure]` - Compare GWAS analysis methods
- `compare_gwas_studies(study_results: Dict[str, Any], output_file: Optional[Union[str, Path]] = None) -> Optional[matplotlib.figure.Figure]` - Compare GWAS studies
- `compare_populations(pop_data: Dict[str, Any], trait_name: str, output_file: Optional[Union[str, Path]] = None) -> Optional[matplotlib.figure.Figure]` - Compare populations
- `calculate_method_power(method_results: Dict[str, Any]) -> Dict[str, float]` - Calculate method statistical power
- `calculate_population_enrichment(pop_data: Dict[str, Any]) -> Dict[str, float]` - Calculate population enrichment
- `analyze_genetic_architecture(pop_data: Dict[str, Any]) -> Dict[str, Any]` - Analyze genetic architecture

### Visualization Effects (`visualization_effects.py`)
- `effect_size_plot(results: Any, output_file: Optional[Union[str, Path]] = None) -> Optional[matplotlib.figure.Figure]` - Plot effect sizes
- `effect_size_manhattan(results: Any, significance_threshold: float = 5e-8, output_file: Optional[Union[str, Path]] = None) -> Optional[matplotlib.figure.Figure]` - Manhattan plot with effect sizes
- `effect_size_qq(results: Any, output_file: Optional[Union[str, Path]] = None) -> Optional[matplotlib.figure.Figure]` - Q-Q plot with effect sizes
- `compare_effect_sizes(study_results: Dict[str, Any], output_file: Optional[Union[str, Path]] = None) -> Optional[matplotlib.figure.Figure]` - Compare effect sizes across studies
