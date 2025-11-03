# Scripts Directory

This directory contains utility scripts for METAINFORMANT development, testing, and deployment.

## Directory Structure

```
scripts/
├── package/        # Package management and testing scripts
├── rna/           # RNA-seq processing scripts
│   └── amalgkit/  # Amalgkit-specific workflow scripts
├── gwas/          # GWAS analysis scripts (genome-scale workflows)
├── dna/           # DNA analysis scripts
├── protein/       # Protein analysis scripts
├── ontology/      # Ontology operations scripts
├── phenotype/     # Phenotype analysis scripts
├── networks/      # Network analysis scripts
├── multiomics/    # Multi-omics integration scripts
├── math/          # Mathematical biology scripts
├── ml/            # Machine learning scripts
├── singlecell/    # Single-cell genomics scripts
├── quality/       # Quality control scripts
├── simulation/    # Simulation scripts
├── visualization/ # Visualization scripts
├── epigenome/     # Epigenome analysis scripts
├── ecology/       # Ecology analysis scripts
├── README.md      # This file
└── AGENTS.md      # AI agent contribution documentation
```

## Package Scripts (`scripts/package/`)

### Environment Setup
- **`setup_uv.sh`**: One-shot uv environment setup script
  - Installs uv package manager
  - Creates virtual environment
  - Installs all dependencies from `pyproject.toml`
  - Use: `bash scripts/package/setup_uv.sh`

- **`uv_dev_setup.sh`**: Development environment setup
  - Extended setup for development workflows
  - Includes additional development tools

### Testing Scripts
- **`run_tests.sh`**: Comprehensive test execution
  - Multiple execution modes (fast, network, pattern matching)
  - Coverage reporting and performance tracking
  - Use: `scripts/package/run_tests.sh --help` for options

- **`uv_test.sh`**: uv-based testing
  - Optimized test execution with uv
  - Parallel test execution

- **`uv_test_optimized.sh`**: Optimized testing workflow
  - Fastest test execution configuration

### Profiling and Quality
- **`uv_profile.sh`**: Performance profiling script
  - Memory and CPU profiling utilities

- **`uv_quality.sh`**: Code quality checking
  - Linting, formatting, and quality metrics

### Documentation
- **`uv_docs.sh`**: Documentation generation
  - Build and update documentation

## GWAS Analysis Scripts (`scripts/gwas/`)

Genome-wide association study workflows with comprehensive visualization.

### Genome-Scale Workflow
- **`run_genome_scale_gwas.py`**: Complete SRA → variants → GWAS pipeline
  - Downloads SRA data from NCBI
  - Aligns reads with BWA
  - Calls variants with bcftools
  - Runs GWAS with comprehensive visualization
  - Skip options for each step
  - **Usage**: `python scripts/gwas/run_genome_scale_gwas.py --help`
  - **Example**: `python scripts/gwas/run_genome_scale_gwas.py --config config/gwas/gwas_amellifera.yaml`

### Data Acquisition
- **`download_real_honeybee_variants.py`**: Real *Apis mellifera* data guide
  - Lists key BioProjects (PRJNA292680, PRJNA392242, PRJNA13343)
  - SRA Toolkit instructions
  - Complete workflow documentation
  - **Usage**: `python scripts/gwas/download_real_honeybee_variants.py`

- **`download_genome_scale_data.sh`**: Batch SRA download script
  - Downloads 3 samples from PRJNA292680
  - Progress reporting and validation
  - **Usage**: `bash scripts/gwas/download_genome_scale_data.sh`

### Metadata Queries
- **`query_bioproject_metadata.py`**: NCBI BioProject metadata retrieval
  - E-utilities integration
  - API queries for run counts
  - JSON output for downstream processing
  - **Usage**: `python scripts/gwas/query_bioproject_metadata.py`

### Features
- ✅ Real data integration (no synthetic data)
- ✅ 32 visualization types
- ✅ Multi-species support
- ✅ Comprehensive error handling
- ✅ Step-by-step execution with skip options
- ✅ Dependency checking

## RNA Processing Scripts (`scripts/rna/`)

### Production Workflows
- **`workflow_ena_integrated.py`**: Integrated ENA download + quantification (PRODUCTION) ⭐
- **`download_ena_robust.py`**: Robust ENA downloader with retry logic ⭐
- **`run_multi_species.py`**: Multi-species orchestration (legacy SRA-based)

### Amalgkit Workflows (`scripts/rna/amalgkit/`)
- **`run_amalgkit.sh`**: Comprehensive RNA-seq pipeline orchestrator
- **`verify_workflow.sh`**: Workflow validation for any species
- See `scripts/rna/amalgkit/README.md` for details

### Quantification and Cleanup
- **`quant_downloaded_samples.py`**: Quantify already-downloaded samples
- **`cleanup_quantified_sra.sh`**: Safe deletion of FASTQ files after quantification
- **`list_unquantified.sh`**: Report samples needing quantification
- **`manual_quant_cleanup.py`**: Manual quantification and cleanup utility
- **`quant_and_cleanup.py`**: Batch quantification and cleanup

### Monitoring
- **`monitor_comprehensive.py`**: Comprehensive real-time monitoring (recommended)
- **`monitor_workflow.py`**: Real-time workflow monitoring dashboard
- **`monitor_amalgkit_progress.sh`**: Simple progress monitor
- **`check_status.py`**: Status checking utility
- **`comprehensive_status.py`**: Detailed status reporting

### Environment and Utilities
- **`check_environment.py`**: Environment validation
- **`check_and_restart_workflows.py`**: Workflow restart automation
- **`restart_all_workflows.py`**: Restart helper

See `scripts/rna/README.md` for detailed RNA script documentation.

## GWAS Analysis Scripts (`scripts/gwas/`)

### Genome-Scale GWAS Workflows
- **`run_genome_scale_gwas.py`**: Complete GWAS workflow orchestrator
  - SRA download through variant calling
  - Association testing with comprehensive visualization
  - Supports download-only and skip-download modes
  - Use: `python3 scripts/gwas/run_genome_scale_gwas.py --config config/gwas/gwas_amellifera.yaml`

### Data Acquisition
- **`download_genome_scale_data.sh`**: Automated GWAS data download
  - Genome assembly retrieval
  - SRA sample download
  - Reference genome preparation

- **`download_real_honeybee_variants.py`**: Real variant data acquisition
  - Honeybee-specific variant retrieval
  - Integration with NCBI resources

### Metadata Management
- **`query_bioproject_metadata.py`**: BioProject metadata extraction
  - Sample information retrieval
  - Phenotype data collection
  - Automated metadata organization

## DNA Analysis Scripts (`scripts/dna/`)

### Analysis Workflow
- **`run_dna_analysis.py`**: DNA analysis workflow orchestrator
  - Sequence processing and quality control
  - Variant analysis
  - Integration with metainformant.dna module

## Protein Analysis Scripts (`scripts/protein/`)

### Analysis Workflow
- **`run_protein_analysis.py`**: Protein analysis workflow orchestrator
  - Structure prediction and functional annotation
  - Domain and motif analysis
  - Integration with metainformant.protein module

## Ontology Scripts (`scripts/ontology/`)

### Analysis Workflow
- **`run_ontology_analysis.py`**: Ontology analysis workflow orchestrator
  - GO enrichment analysis
  - Semantic similarity calculations
  - Pathway analysis

## Phenotype Scripts (`scripts/phenotype/`)

### Analysis Workflow
- **`run_phenotype_analysis.py`**: Phenotype analysis workflow orchestrator
  - Trait associations and correlations
  - Heritability estimation
  - Phenotype prediction

## Network Analysis Scripts (`scripts/networks/`)

### Analysis Workflow
- **`run_network_analysis.py`**: Network analysis workflow orchestrator
  - Network construction and visualization
  - Community detection
  - Pathway and interaction analysis

## Multi-Omics Scripts (`scripts/multiomics/`)

### Integration Workflow
- **`run_multiomics_integration.py`**: Multi-omics integration orchestrator
  - Cross-omics correlation analysis
  - Data integration across genomics, transcriptomics, proteomics
  - Systems biology approaches

## Mathematical Biology Scripts (`scripts/math/`)

### Modeling Workflow
- **`run_math_modeling.py`**: Mathematical biology orchestrator
  - Differential equations and dynamical systems
  - Theoretical modeling
  - Stability and bifurcation analysis

## Machine Learning Scripts (`scripts/ml/`)

### ML Pipeline
- **`run_ml_pipeline.py`**: Machine learning workflow orchestrator
  - Model training and hyperparameter tuning
  - Cross-validation and evaluation
  - Prediction generation

## Single-Cell Genomics Scripts (`scripts/singlecell/`)

### Analysis Workflow
- **`run_singlecell_analysis.py`**: Single-cell analysis orchestrator
  - Quality control and normalization
  - Clustering and cell type annotation
  - Trajectory and differential expression analysis

## Quality Control Scripts (`scripts/quality/`)

### QC Workflow
- **`run_quality_control.py`**: Quality control orchestrator
  - Data validation and QC metrics
  - Quality reporting
  - Issue identification

## Simulation Scripts (`scripts/simulation/`)

### Simulation Workflow
- **`run_simulation.py`**: Biological simulation orchestrator
  - Agent-based models
  - Population dynamics
  - Evolutionary simulations

## Visualization Scripts (`scripts/visualization/`)

### Visualization Workflow
- **`run_visualization.py`**: Visualization workflow orchestrator
  - Plot generation and interactive visualizations
  - Dashboard creation
  - Report generation

## Epigenome Analysis Scripts (`scripts/epigenome/`)

### Analysis Workflow
- **`run_epigenome_analysis.py`**: Epigenome analysis orchestrator
  - Methylation analysis
  - Histone modifications
  - Chromatin accessibility

## Ecology Analysis Scripts (`scripts/ecology/`)

### Analysis Workflow
- **`run_ecology_analysis.py`**: Ecology analysis orchestrator
  - Species diversity and community composition
  - Ordination and clustering
  - Ecological modeling

## Usage Guidelines

### Complete Workflow Demonstration

The repository includes a comprehensive demonstration script showing best practices:

```bash
# Run complete workflow demonstration
python3 scripts/run_complete_demo.py

# This demonstrates:
# - Configuration management
# - Input data handling
# - Data processing with statistics
# - Visualization generation with informative filenames
# - Comprehensive output organization
```

**Output Structure:**
```
output/demo/
├── workflow_configuration.json          # Workflow parameters
├── input_samples.json                   # Input data
├── processed_normalized_data.json       # Results
├── workflow_summary_report.json         # Summary
└── visualizations/                      # All plots
    ├── input_sample_values_lineplot.png
    ├── normalized_values_barplot.png
    ├── sample_correlation_heatmap.png
    └── visualization_metadata.json
```

**Key Features:**
- ✅ Informative visualization names (e.g., `input_sample_values_lineplot.png` not `plot1.png`)
- ✅ Configuration and input data saved
- ✅ Processing steps with statistics
- ✅ High-resolution visualizations (300 DPI)
- ✅ Complete workflow tracking
- ✅ Summary report generation

See `output/demo/VERIFICATION_REPORT.md` for detailed verification.

### Development Workflow
```bash
# Setup environment
bash scripts/package/setup_uv.sh

# Run tests
bash scripts/package/run_tests.sh

# Check code quality
bash scripts/package/uv_quality.sh

# Generate documentation
bash scripts/package/uv_docs.sh
```

### Production Usage
```bash
# Run RNA analysis workflow
bash scripts/rna/amalgkit/run_amalgkit.sh --config config/amalgkit/amalgkit_pbarbatus.yaml

# Multi-species RNA workflow
python3 scripts/rna/run_multi_species.py

# Monitor RNA workflow progress
python3 scripts/rna/monitor_comprehensive.py

# Run genome-scale GWAS
python3 scripts/gwas/run_genome_scale_gwas.py --config config/gwas/gwas_amellifera.yaml
```

## Script Standards

### Repository Organization
- **Scripts always belong in `scripts/`** - never in `output/`, `data/`, or other directories
- Scripts may read from `data/` and `config/` directories
- Scripts must write outputs to `output/` directory by default
- Use subdirectories for domain-specific scripts (e.g., `scripts/rna/`)

### Requirements
- All scripts must be executable (`chmod +x`)
- Include help documentation and usage examples
- Handle errors gracefully with informative messages
- Support common flags (`--help`, `--verbose`, `--dry-run`)
- Compatible with both bash and zsh

### Maintenance
- Update scripts when dependencies change
- Test scripts on multiple systems
- Document any external dependencies
- Version scripts with significant changes

## Contributing New Scripts

### Template
```bash
#!/bin/bash
# Script Name: Brief description
#
# Usage: script_name [options]
#
# Options:
#   -h, --help     Show this help message
#   -v, --verbose  Verbose output
#   -d, --dry-run  Show what would be done

set -euo pipefail

# Script implementation
```

### Review Process
1. Test script functionality thoroughly
2. Add comprehensive help documentation
3. Ensure cross-platform compatibility
4. Update this README with new script
5. **Place script in `scripts/` directory** (not `output/` or other locations)

## Integration

Scripts integrate with:
- **uv** for Python package management
- **pytest** for testing framework
- **amalgkit** for RNA analysis
- **Core utilities** for file operations

## Troubleshooting

### Common Issues
- **Permission denied**: Run `chmod +x script_name.sh`
- **Command not found**: Ensure dependencies are installed
- **Environment errors**: Source virtual environment first

### Getting Help
- Use `--help` flag on any script
- Check script source code for detailed comments
- See related documentation for context

These scripts provide essential tooling for METAINFORMANT development and usage workflows.
