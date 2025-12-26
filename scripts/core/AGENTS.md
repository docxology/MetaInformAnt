# AI Agents in Script Development

This document outlines AI assistance in developing METAINFORMANT's utility scripts and automation tools.

## AI Contributions

### Script Architecture
**Code Assistant Agent** designed:
- Modular script organization with subfolder structure
- Consistent command-line interfaces
- Error handling and logging patterns
- Cross-platform compatibility approaches

### Automation Features
**Code Assistant Agent** implemented:
- Environment setup automation (`scripts/package/`)
- Test execution workflows (`scripts/package/run_tests.sh`)
- Documentation generation scripts (`scripts/package/uv_docs.sh`)
- RNA processing pipelines (`scripts/rna/`)
- GWAS analysis workflows (`scripts/gwas/`)
- Domain-specific orchestration scripts for all modules (`scripts/*/`)

### User Experience
**Documentation Agent** contributed to:
- Comprehensive help documentation
- Usage examples and best practices
- Troubleshooting guides
- Integration with development workflows

## Script Categories

### Package Management (`scripts/package/`)
- **setup_uv.sh**: Automated environment configuration
- **uv_dev_setup.sh**: Development environment optimization
- **run_tests.sh**: Comprehensive test execution with coverage
- **uv_test.sh**: Optimized testing workflows
- **uv_test_optimized.sh**: Fast test modes
- **uv_quality.sh**: Code quality assessment
- **uv_profile.sh**: Performance profiling utilities
- **uv_docs.sh**: Automated documentation generation

### RNA Processing (`scripts/rna/`)
- **batch_ena.py**: Fast ENA-based parallel downloader
- **run_multi_species_amalgkit.py**: Multi-species orchestration with cross-species analysis
- **run_multi_species_sequential.py**: Disk-space-friendly sequential processing
- **monitor_workflow.py**: Real-time monitoring dashboard
- **quant_downloaded_samples.py**: Quantification of downloaded samples
- **cleanup_quantified_sra.sh**: Safe SRA file cleanup
- **list_unquantified.sh**: Sample status reporting
- **force_fasterq.sh**: High-performance FASTQ processing
- **force_fasterq_parallel.sh**: Parallel FASTQ processing
- **process_one_srr.sh**: Individual SRR processing
- **test_*.py**: Testing and verification scripts
- **verify_skip_logic.sh**: Skip logic comprehensive verification

### Amalgkit Workflows (`scripts/rna/amalgkit/`)
- **run_amalgkit.sh**: Comprehensive RNA-seq pipeline orchestrator
- **verify_workflow.sh**: Workflow validation for any species

### GWAS Analysis (`scripts/gwas/`)
- **run_genome_scale_gwas.py**: Complete GWAS workflow orchestrator
  - SRA download through variant calling and association testing
  - Comprehensive visualization with Manhattan and Q-Q plots
  - Download-only and skip-download execution modes
  - Integration with metainformant.gwas module
- **download_genome_scale_data.sh**: Automated GWAS data acquisition
  - Genome assembly retrieval from NCBI
  - SRA sample download and organization
  - Reference genome preparation and indexing
- **download_honeybee_variants.py**: Variant data acquisition
  - Honeybee-specific variant datasets
  - Integration with NCBI VCF resources
- **query_bioproject_metadata.py**: BioProject metadata extraction
  - Automated sample information retrieval
  - Phenotype data collection and organization
  - Integration with NCBI BioProject API

### Domain-Specific Module Orchestrators
**Code Assistant Agent** created thin orchestration layers for all modules (November 2025):

#### DNA Analysis (`scripts/dna/`)
- **run_dna_analysis.py**: DNA analysis workflow orchestrator
  - Sequence processing and quality control framework
  - Variant analysis pipeline integration
  - Placeholder for metainformant.dna module integration

#### Protein Analysis (`scripts/protein/`)
- **run_protein_analysis.py**: Protein analysis workflow orchestrator
  - Structure prediction and functional annotation framework
  - Domain and motif analysis integration
  - Placeholder for metainformant.protein module integration

#### Ontology (`scripts/ontology/`)
- **run_ontology_analysis.py**: Ontology operations orchestrator
  - GO enrichment analysis framework
  - Semantic similarity calculation support
  - Pathway analysis integration placeholder

#### Phenotype (`scripts/phenotype/`)
- **run_phenotype_analysis.py**: Phenotype analysis orchestrator
  - Trait association and correlation framework
  - Heritability estimation support
  - Phenotype prediction pipeline placeholder

#### Network Analysis (`scripts/networks/`)
- **run_network_analysis.py**: Network analysis orchestrator
  - Network construction and visualization framework
  - Community detection support
  - Pathway and interaction analysis placeholder

#### Multi-Omics (`scripts/multiomics/`)
- **run_multiomics_integration.py**: Multi-omics integration orchestrator
  - Cross-omics correlation analysis framework
  - Data integration for genomics, transcriptomics, proteomics
  - Systems biology approaches placeholder

#### Mathematical Biology (`scripts/math/`)
- **run_math_modeling.py**: Mathematical modeling orchestrator
  - Differential equations and dynamical systems framework
  - Theoretical modeling support
  - Stability and bifurcation analysis placeholder

#### Machine Learning (`scripts/ml/`)
- **run_ml_pipeline.py**: Machine learning orchestrator
  - Model training and hyperparameter tuning framework
  - Cross-validation and evaluation support
  - Prediction generation placeholder

#### Single-Cell Genomics (`scripts/singlecell/`)
- **run_singlecell_analysis.py**: Single-cell analysis orchestrator
  - Quality control and normalization framework
  - Clustering and cell type annotation support
  - Trajectory and differential expression placeholder

#### Quality Control (`scripts/quality/`)
- **run_quality_control.py**: Quality control orchestrator
  - Data validation and QC metrics framework
  - Quality reporting support
  - Issue identification placeholder

#### Simulation (`scripts/simulation/`)
- **run_simulation.py**: Biological simulation orchestrator
  - Agent-based models framework
  - Population dynamics support
  - Evolutionary simulation placeholder

#### Visualization (`scripts/visualization/`)
- **run_visualization.py**: Visualization orchestrator
  - Plot generation and interactive visualization framework
  - Dashboard creation support
  - Report generation placeholder

#### Epigenome (`scripts/epigenome/`)
- **run_epigenome_analysis.py**: Epigenome analysis orchestrator
  - Methylation analysis framework
  - Histone modification support
  - Chromatin accessibility analysis placeholder

#### Ecology (`scripts/ecology/`)
- **run_ecology_analysis.py**: Ecology analysis orchestrator
  - Species diversity and community composition framework
  - Ordination and clustering support
  - Ecological modeling placeholder

### Orchestrator Design Principles
All domain-specific orchestrators follow consistent patterns:
- Standard argument parsing with --input, --output, --verbose, --dry-run
- Integration with metainformant.core utilities (logging, paths, I/O)
- Output directory defaulting to output/<module_name>/
- Placeholder TODO sections for future implementation
- Executable permissions and proper shebang lines
- Comprehensive error handling and logging

## Development Standards

### Script Quality
- Executable permissions and proper shebang lines
- Comprehensive error handling and exit codes
- Consistent option parsing and help messages
- Cross-platform shell compatibility
- No hardcoded absolute paths

### Repository Organization
- Scripts organized in logical subfolders
- `package/` for package management and testing
- `rna/` for RNA-seq processing
- `rna/amalgkit/` for amalgkit-specific scripts
- `gwas/` for GWAS analysis workflows
- `dna/`, `protein/`, `ontology/`, `phenotype/` for respective analyses
- `networks/`, `multiomics/`, `math/`, `ml/` for computational biology
- `singlecell/`, `quality/`, `simulation/` for specialized workflows
- `visualization/`, `epigenome/`, `ecology/` for domain-specific analyses
- Clear separation from output directories
- One script folder per metainformant module

### Maintenance Practices
- Regular updates with dependency changes
- Testing on multiple environments
- Clear documentation of external dependencies
- Version tracking for significant changes
- Removal of obsolete scripts with hardcoded paths

## Integration

Scripts integrate with:
- **uv** for Python package management
- **pytest** for testing framework
- **amalgkit** for RNA analysis workflows
- **SRA Toolkit** for sequence data acquisition
- **BWA, samtools, bcftools** for GWAS variant calling
- **Core utilities** for file operations and path handling

## Script Centralization Philosophy

### Methods vs Data Separation
- **Scripts** (methods) centralized in `scripts/`
- **Data** (outputs) distributed in `output/`
- Single source of truth for all scripts
- Update once, applies everywhere
- Proper version control and tracking

### Benefits
- ✅ No script duplication across species/projects
- ✅ Clean output directories (data only)
- ✅ Easy maintenance and updates
- ✅ Consistent behavior across all uses

This script collection provides essential tooling for METAINFORMANT development and deployment workflows.
