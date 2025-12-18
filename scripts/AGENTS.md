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
**Code Assistant Agent** created:
- `setup.sh`: Unified environment setup with UV and dependency management
- `test.sh`: Comprehensive test runner with multiple modes and coverage
- `verify.sh`: Environment verification with automatic fixes
- Quality assurance scripts (`uv_quality.sh`, `uv_profile.sh`, `uv_docs.sh`)
- Shared utilities (`_common.sh`) for cross-platform compatibility

### RNA Processing (`scripts/rna/`)
**Code Assistant Agent** consolidated 20+ scripts into streamlined orchestrators:
- `run_workflow.py`: Main RNA-seq workflow orchestrator
- `setup_genome.py`: Genome preparation pipeline
- `discover_species.py`: Species discovery and config generation
- `run_all_species_validation.sh`: End-to-end validation
- `triple_check_rna.py`: Documentation and code verification
- `verify_rna.py` and `verify_rna_docs.py`: Quality assurance

### GWAS Analysis (`scripts/gwas/`)
**Code Assistant Agent** implemented:
- `run_genome_scale_gwas.py`: Complete SRA→variants→GWAS pipeline
- `run_pbarbatus_gwas.py`: Species-specific GWAS analysis
- `download_genome_scale_data.sh`: Real data acquisition
- `download_real_honeybee_variants.py`: Variant data integration
- `query_bioproject_metadata.py`: NCBI metadata queries
- Enhanced visualization scripts with 32 plot types

### Core Utilities (`scripts/core/`)
**Code Assistant Agent** developed:
- `run_demo.py`: Complete workflow demonstration
- `fix_disk_space.sh`: External drive temp directory management
- `setup_temp_dirs.sh`: Repository-local temp directory setup
- Comprehensive error handling and filesystem compatibility

### Domain-Specific Orchestrators
**Code Assistant Agent** created thin orchestration layers for all modules:

#### DNA Analysis (`scripts/dna/`)
- `run_dna_analysis.py`: Sequence processing and variant analysis

#### Ecology Analysis (`scripts/ecology/`)
- `run_ecology_analysis.py`: Diversity indices and community composition

#### Epigenome Analysis (`scripts/epigenome/`)
- `run_epigenome_analysis.py`: Methylation and chromatin analysis

#### Life Events Analysis (`scripts/life_events/`)
- `run_life_events_analysis.py`: Sequence embedding and outcome prediction
- `generate_synthetic_data.py`: Synthetic life event generation
- `learn_embeddings.py`: Neural language model training
- `train_model.py`: Predictive model training
- `predict_outcomes.py`: Outcome prediction and interpretation
- `generate_all_visualizations.py`: Comprehensive visualization suite

#### Mathematical Biology (`scripts/math/`)
- `run_math_modeling.py`: Population dynamics and epidemiological models

#### Machine Learning (`scripts/ml/`)
- `run_ml_pipeline.py`: Feature selection and model evaluation

#### Multi-Omics Integration (`scripts/multiomics/`)
- `run_multiomics_integration.py`: Cross-omics correlation analysis
- Example scripts for CCA and joint analysis

#### Network Analysis (`scripts/networks/`)
- `run_network_analysis.py`: Graph theory and community detection

#### Ontology Analysis (`scripts/ontology/`)
- `run_ontology_analysis.py`: GO term queries and enrichment

#### Phenotype Analysis (`scripts/phenotype/`)
- `run_phenotype_analysis.py`: Trait correlation and statistics
- `scrape_antwiki.py`: AntWiki species data extraction
- `load_antwiki_example.py`: Data loading utilities

#### Population Genetics (`scripts/popgen/`)
- `analysis.py`: Comprehensive synthetic data generation and analysis

#### Protein Analysis (`scripts/protein/`)
- `run_protein_analysis.py`: Structure prediction and annotation

#### Quality Control (`scripts/quality/`)
- `run_quality_control.py`: Data validation and QC metrics

#### Single-Cell Analysis (`scripts/singlecell/`)
- `run_singlecell_analysis.py`: scRNA-seq processing and clustering

#### Visualization (`scripts/visualization/`)
- `run_visualization.py`: Publication-quality plot generation

## Orchestrator Design Principles

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
- Scripts organized in logical subfolders by domain
- `package/` for package management and testing
- `rna/` for RNA-seq processing
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

## Related Documentation

For detailed AI contributions to specific script domains:

### Core Scripts
- [`scripts/core/AGENTS.md`](core/AGENTS.md) - Core utility scripts development

### RNA Scripts
- [`scripts/rna/AGENTS.md`](rna/AGENTS.md) - RNA workflow scripts development
- [`scripts/rna/amalgkit/AGENTS.md`](rna/amalgkit/AGENTS.md) - Amalgkit script development

### Domain-Specific Scripts
- [`scripts/gwas/AGENTS.md`](gwas/AGENTS.md) - GWAS analysis scripts development
- [`scripts/life_events/AGENTS.md`](life_events/AGENTS.md) - Life events analysis scripts development
- Additional AGENTS.md files available for all domain modules

This script collection provides essential tooling for METAINFORMANT development and deployment workflows.

