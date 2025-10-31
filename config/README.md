# Configuration Directory

This directory contains repository-level configuration files used by METAINFORMANT for various tools and workflows, with a focus on RNA-seq and GWAS analysis configurations for ant species.

## Directory Structure

```
config/
├── amalgkit/          # RNA-seq workflow configurations for ant species
├── gwas/             # GWAS workflow configurations
├── archive/          # Archived/inactive configurations
├── README.md         # This file
└── AGENTS.md        # AI agent contribution documentation
```

## Amalgkit RNA-seq Configurations

### Active Species Configurations (`amalgkit/`)

- **`amalgkit_pbarbatus.yaml`**: *Pogonomyrmex barbatus* (Red Harvester Ant) - Assembly GCF_000187915.1_Pbar_UMD_V03
- **`amalgkit_cfloridanus.yaml`**: *Camponotus floridanus* (Florida Carpenter Ant) - Assembly GCF_003227725.1_Cflo_v7.5
- **`amalgkit_mpharaonis.yaml`**: *Monomorium pharaonis* (Pharaoh Ant) - Assembly GCF_013373865.1_ASM1337386v2
- **`amalgkit_sinvicta.yaml`**: *Solenopsis invicta* (Red Fire Ant) - Assembly GCF_016802725.1_UNIL_Sinv_3.0
- **`amalgkit_template.yaml`**: Template configuration file for creating new amalgkit workflows with customizable parameters

### Archived Configurations (`archive/`)

- **`amalgkit_amellifera.yaml`**: *Apis mellifera* (Western Honey Bee) - Inactive configuration moved to archive

## GWAS Analysis Configurations

### GWAS Configurations (`gwas/`)

- **`gwas_pbarbatus.yaml`**: GWAS configuration for *Pogonomyrmex barbatus* - includes variant calling, QC, association testing
- **`gwas_template.yaml`**: Template configuration file for creating new GWAS workflows

## Usage

### Amalgkit RNA-seq Workflows

Configuration files are used by RNA analysis workflows and can be customized for different species and computational requirements.

#### Loading Configurations
```python
from metainformant.core import config
from metainformant.rna import AmalgkitWorkflowConfig

# Load configuration from file
cfg = config.load_config("config/amalgkit/amalgkit_pbarbatus.yaml")

# Create workflow configuration
workflow_cfg = AmalgkitWorkflowConfig(
    work_dir="output/amalgkit/pbarbatus/work",
    threads=cfg.threads,
    species_list=cfg.species_list,
    **cfg.steps
)
```

#### Command Line Usage
```bash
# Run amalgkit workflow with specific configuration
uv run python -m metainformant rna run \
  --config config/amalgkit/amalgkit_sinvicta.yaml \
  --work-dir output/amalgkit/sinvicta
```

### GWAS Workflows

#### Loading GWAS Configurations
```python
from metainformant.core import config

# Load GWAS configuration
gwas_cfg = config.load_config("config/gwas/gwas_pbarbatus.yaml")

# Access configuration sections
genome_cfg = gwas_cfg.genome
variants_cfg = gwas_cfg.variants
qc_cfg = gwas_cfg.qc
```

#### Command Line Usage
```bash
# Run GWAS workflow with specific configuration
uv run python -m metainformant gwas run \
  --config config/gwas/gwas_pbarbatus.yaml \
  --work-dir output/gwas/pbarbatus
```

## Configuration Parameters

### Amalgkit Common Parameters

- **`work_dir`**: Working directory for analysis outputs (e.g., `output/amalgkit/pbarbatus/work`)
- **`log_dir`**: Directory for log files (e.g., `output/amalgkit/pbarbatus/logs`)
- **`threads`**: Number of parallel threads for processing (typically 6-8)
- **`auto_install_amalgkit`**: Automatically install amalgkit if not present (true/false)
- **`species_list`**: List of species for analysis (e.g., `["Pogonomyrmex_barbatus"]`)
- **`taxon_id`**: NCBI Taxonomy ID for the species

### Amalgkit Genome Configuration

- **`genome.accession`**: NCBI assembly accession (e.g., `GCF_000187915.1`)
- **`genome.assembly_name`**: Assembly name (e.g., `Pbar_UMD_V03`)
- **`genome.annotation_release`**: NCBI annotation release number
- **`genome.dest_dir`**: Local destination for genome files
- **`genome.include`**: List of genome file types to download (genome, gff3, rna, cds, protein, etc.)
- **`genome.ftp_url`**: NCBI FTP URL for genome files

### Amalgkit Workflow Steps

- **`metadata`**: Search and download metadata from NCBI SRA
- **`integrate`**: Integrate metadata with local samples
- **`config`**: Configure workflow parameters
- **`select`**: Select samples for analysis
- **`getfastq`**: Download FASTQ files from SRA (supports aws, gcp, ncbi mirrors)
- **`quant`**: Quantify transcript abundances with kallisto/salmon
- **`merge`**: Merge quantification results across samples
- **`cstmm`**: Custom metadata management
- **`curate`**: Curate and validate results
- **`csca`**: Custom comparative analysis
- **`sanity`**: Sanity checks and validation

### GWAS Configuration Parameters

- **`genome`**: Reference genome configuration (reuses amalgkit genome structure)
- **`variants`**: Variant data sources (VCF files, calling from BAM, or downloads)
- **`qc`**: Quality control thresholds (MAF, missing rate, HWE, etc.)
- **`samples`**: Sample metadata and phenotype files
- **`structure`**: Population structure analysis (PCA, kinship)
- **`association`**: Association testing parameters (model, trait, covariates)
- **`correction`**: Multiple testing correction methods
- **`output`**: Output directories and formats

### Supported Ant Species

- **`Pogonomyrmex_barbatus`**: Red Harvester Ant (NCBI Taxonomy: 144034)
- **`Camponotus_floridanus`**: Florida Carpenter Ant (NCBI Taxonomy: 104421)
- **`Monomorium_pharaonis`**: Pharaoh Ant (NCBI Taxonomy: 307658)
- **`Solenopsis_invicta`**: Red Fire Ant (NCBI Taxonomy: 13686)

### Archived Species

- **`Apis_mellifera`**: Western Honey Bee (NCBI Taxonomy: 7460) - Configuration archived

## Integration

### Module Integration

- **RNA Module**: `src/metainformant/rna/` - Amalgkit workflow configuration and execution
- **GWAS Module**: `src/metainformant/gwas/` - GWAS workflow configuration and execution
- **Core Config**: `src/metainformant/core/config.py` - Configuration loading with environment variable overrides
- **Core Paths**: `src/metainformant/core/paths.py` - Path validation and management

### Execution Environments

- Local execution with configurable thread/memory limits
- SLURM cluster integration (planned)
- PBS cluster integration (planned)
- Cloud execution (AWS/GCP for data downloads)

### Environment Variable Overrides

All configuration parameters can be overridden via environment variables:
```bash
export METAINFORMANT_THREADS=16
export METAINFORMANT_WORK_DIR=/scratch/user/amalgkit
export NCBI_EMAIL=user@example.com
```

## Modifying Configurations

When creating or modifying these files:
1. Test changes with small datasets first
2. Validate against existing test suites
3. Update corresponding documentation
4. Consider backward compatibility
5. Document any new parameters added
6. Add appropriate test configurations

## Best Practices

### Configuration Design
- Use descriptive parameter names
- Include comments explaining complex parameters
- Provide sensible defaults for optional parameters
- Group related parameters logically
- Support environment variable overrides

### Testing Configurations
- Create test-specific configurations for validation
- Include sample datasets for configuration testing
- Document expected outputs for each configuration
- Test both success and failure scenarios

### Version Control
- Track configuration changes in version control
- Document breaking changes in release notes
- Maintain backward compatibility when possible
- Use semantic versioning for configuration formats

## Related Documentation

### Configuration Management
- `config/AGENTS.md` - AI agent contributions to configuration system
- `src/metainformant/core/config.py` - Configuration loading implementation
- `src/metainformant/core/paths.py` - Path handling and validation

### Workflow Documentation
- `docs/rna/` - RNA-seq workflow documentation (if exists)
- `docs/gwas/` - GWAS workflow documentation (if exists)
- `src/metainformant/rna/README.md` - RNA module documentation (if exists)
- `src/metainformant/gwas/README.md` - GWAS module documentation (if exists)

### Testing
- Configuration validation tests in test suite
- Integration tests for workflow execution
- Template configuration validation

---

**Last Updated**: October 31, 2025  
**Maintainer**: METAINFORMANT Development Team  
**Status**: Active development - ant species RNA-seq and GWAS analysis
