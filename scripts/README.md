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
├── life_events/   # Life events analysis scripts
├── archive/       # Archived/obsolete scripts (see archive/README.md)
├── README.md      # This file
└── AGENTS.md      # AI agent contribution documentation
```

## Template Files

- **`_template_working.py`**: Reference template demonstrating proper script structure, imports, configuration, and output handling. Useful for creating new orchestrator scripts.

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

## RNA Processing Scripts (`scripts/rna/`)

### Production Workflows
- **`run_workflow.py`**: ✅ Main RNA-seq workflow orchestrator (PRODUCTION) ⭐
  - End-to-end workflow execution for single species
  - Status checking and monitoring
  - Configuration-based execution
  - Usage: `python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml`

- **`run_e2e_pbarbatus.py`**: End-to-end workflow for *Pogonomyrmex barbatus*
- **`setup_genome.py`**: Genome setup and preparation utilities
- **`discover_species.py`**: Species discovery and configuration generation
- **`convert_existing_sra.py`**: Convert existing SRA data to workflow format

### Amalgkit Workflows (`scripts/rna/amalgkit/`)
- **`run_amalgkit.sh`**: Comprehensive RNA-seq pipeline orchestrator
- **`verify_workflow.sh`**: Workflow validation for any species
- See `scripts/rna/amalgkit/README.md` for details

### Environment and Utilities
- **`check_environment.py`**: Environment validation
- **`fix_tmp_space.sh`**: Temporary space management utilities

### Testing and Development
- **`test_enhanced_heartbeat.py`**: Heartbeat monitoring tests
- **`test_genome_prep.py`**: Genome preparation tests
- **`test_getfastq_fix.py`**: FASTQ download fix tests
- **`test_heartbeat.py`**: Basic heartbeat tests

### Examples
- **`examples/discover_and_deploy_ant_species.sh`**: Example script for species discovery and deployment

See `scripts/rna/README.md` for detailed RNA script documentation.

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

## DNA Analysis Scripts (`scripts/dna/`)

### Analysis Workflow
- **`run_dna_analysis.py`**: ✅ DNA analysis workflow orchestrator (IMPLEMENTED)
  - Sequence processing and quality control
  - Variant analysis
  - Integration with metainformant.dna module
  - Usage: `python3 scripts/dna/run_dna_analysis.py --help`

## Protein Analysis Scripts (`scripts/protein/`)

### Analysis Workflow
- **`run_protein_analysis.py`**: ✅ Protein analysis workflow orchestrator (IMPLEMENTED)
  - Structure prediction and functional annotation
  - Domain and motif analysis
  - Integration with metainformant.protein module
  - Usage: `python3 scripts/protein/run_protein_analysis.py --help`

## Ontology Scripts (`scripts/ontology/`)

### Analysis Workflow
- **`run_ontology_analysis.py`**: ✅ Ontology analysis workflow orchestrator (IMPLEMENTED)
  - GO enrichment analysis
  - Semantic similarity calculations
  - Pathway analysis
  - Usage: `python3 scripts/ontology/run_ontology_analysis.py --help` or `uv run metainformant ontology run --help`

## Phenotype Scripts (`scripts/phenotype/`)

### Analysis Workflow
- **`run_phenotype_analysis.py`**: ✅ Phenotype analysis workflow orchestrator (IMPLEMENTED)
  - Trait associations and correlations
  - Statistics and correlation analysis
  - Integration with metainformant.phenotype module
  - Usage: `python3 scripts/phenotype/run_phenotype_analysis.py --help` or `uv run metainformant phenotype run --help`

## Network Analysis Scripts (`scripts/networks/`)

### Analysis Workflow
- **`run_network_analysis.py`**: ✅ Network analysis workflow orchestrator (IMPLEMENTED)
  - Network construction and visualization
  - Community detection
  - Pathway and interaction analysis
  - Network metrics and centrality calculations
  - Usage: `python3 scripts/networks/run_network_analysis.py --help` or `uv run metainformant networks run --help`

## Multi-Omics Scripts (`scripts/multiomics/`)

### Integration Workflow
- **`run_multiomics_integration.py`**: ✅ Multi-omics integration orchestrator (IMPLEMENTED)
  - Cross-omics correlation analysis
  - Data integration across genomics, transcriptomics, proteomics
  - Joint PCA, NMF, and canonical correlation analysis
  - Usage: `python3 scripts/multiomics/run_multiomics_integration.py --help` or `uv run metainformant multiomics run --help`

## Mathematical Biology Scripts (`scripts/math/`)

### Modeling Workflow
- **`run_math_modeling.py`**: ✅ Mathematical biology orchestrator (IMPLEMENTED)
  - Differential equations and dynamical systems
  - Theoretical modeling
  - Stability and bifurcation analysis
  - Usage: `python3 scripts/math/run_math_modeling.py --help`

## Machine Learning Scripts (`scripts/ml/`)

### ML Pipeline
- **`run_ml_pipeline.py`**: ✅ Machine learning workflow orchestrator (IMPLEMENTED)
  - Model training and hyperparameter tuning
  - Cross-validation and evaluation
  - Classification, regression, and feature selection
  - Usage: `python3 scripts/ml/run_ml_pipeline.py --help` or `uv run metainformant ml run --help`

## Single-Cell Genomics Scripts (`scripts/singlecell/`)

### Analysis Workflow
- **`run_singlecell_analysis.py`**: ✅ Single-cell analysis orchestrator (IMPLEMENTED)
  - Quality control and normalization
  - Clustering and cell type annotation
  - Trajectory and differential expression analysis
  - Usage: `python3 scripts/singlecell/run_singlecell_analysis.py --help` or `uv run metainformant singlecell run --help`

## Quality Control Scripts (`scripts/quality/`)

### QC Workflow
- **`run_quality_control.py`**: ✅ Quality control orchestrator (IMPLEMENTED)
  - Data validation and QC metrics
  - FASTQ quality analysis
  - Contamination detection
  - Quality reporting
  - Usage: `python3 scripts/quality/run_quality_control.py --help` or `uv run metainformant quality run --help`

## Simulation Scripts (`scripts/simulation/`)

### Simulation Workflow
- **`run_simulation.py`**: ✅ Biological simulation orchestrator (IMPLEMENTED)
  - Agent-based models
  - Sequence generation
  - Expression simulation
  - Usage: `python3 scripts/simulation/run_simulation.py --help` or `uv run metainformant simulation run --help`

## Visualization Scripts (`scripts/visualization/`)

### Visualization Workflow
- **`run_visualization.py`**: ✅ Visualization workflow orchestrator (IMPLEMENTED)
  - Plot generation and interactive visualizations
  - Lineplots, heatmaps, animations, histograms
  - Dashboard creation
  - Usage: `python3 scripts/visualization/run_visualization.py --help` or `uv run metainformant visualization run --help`

## Epigenome Analysis Scripts (`scripts/epigenome/`)

### Analysis Workflow
- **`run_epigenome_analysis.py`**: ✅ Epigenome analysis orchestrator (IMPLEMENTED)
  - Methylation analysis
  - Beta value computation
  - BedGraph track loading
  - Usage: `python3 scripts/epigenome/run_epigenome_analysis.py --help` or `uv run metainformant epigenome run --help`

## Ecology Analysis Scripts (`scripts/ecology/`)

### Analysis Workflow
- **`run_ecology_analysis.py`**: ✅ Ecology analysis orchestrator (IMPLEMENTED)
  - Species diversity and community composition
  - Diversity indices (Shannon, Simpson, etc.)
  - Beta diversity calculations
  - Usage: `python3 scripts/ecology/run_ecology_analysis.py --help` or `uv run metainformant ecology run --help`

## Life Events Analysis Scripts (`scripts/life_events/`)

### Analysis Workflow
- **`run_life_events_analysis.py`**: ✅ Life events analysis orchestrator (IMPLEMENTED)
  - Synthetic data generation for testing and demos
  - Event sequence embedding learning
  - Model training and prediction
  - Comprehensive visualization generation
  - End-to-end workflow from data to visualizations
  - Usage: `python3 scripts/life_events/run_life_events_analysis.py --help`
  - **Examples**:
    ```bash
    # Generate synthetic data and run complete analysis
    python3 scripts/life_events/run_life_events_analysis.py --synthetic --n-sequences 100 --generate-outcomes
    
    # Analyze existing sequences with config file
    python3 scripts/life_events/run_life_events_analysis.py --input data/life_events/sequences.json --config config/life_events_template.yaml
    
    # Generate all visualizations
    python3 scripts/life_events/run_life_events_analysis.py --synthetic --all-visualizations
    ```

## Usage Guidelines

### Complete Workflow Demonstration

The repository includes a comprehensive demonstration script showing best practices:

```bash
# Run complete workflow demonstration
python3 scripts/core/run_demo.py

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

See `scripts/core/run_demo.py` for the workflow demonstration. Outputs are saved to `output/demo/` with workflow configuration, processed data, visualizations, and summary reports.

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
bash scripts/rna/amalgkit/run_amalgkit.sh --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml

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
- **Core utilities** for file operations (see [Core Documentation](../docs/core/README.md))
- **Configuration management** via `core.config` (see [Config Documentation](../docs/core/config.md))
- **CLI interface** via `metainformant` command (see [CLI Documentation](../docs/cli.md))

## Testing

Orchestrator scripts are tested via `tests/test_orchestrators.py`:
- Help text verification
- CLI integration tests
- Basic functionality validation

See [Testing Documentation](../docs/testing.md) for complete testing policy.

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
