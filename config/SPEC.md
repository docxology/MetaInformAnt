# Specification: config

## Scope

YAML configuration files for all METAINFORMANT workflows and modules. Contains active configurations, templates, and archived configs. Supports environment variable overrides with module-specific prefixes.

## Architecture

- **Dependency Level**: Configuration
- **Component Type**: Configuration Files
- **Format**: YAML with environment variable override support

### Directory Structure

```
config/
├── amalgkit/           # RNA-seq workflow configs (active species)
├── archive/            # Archived/legacy amalgkit configs
├── config_base/        # Base configuration templates
├── eqtl/              # eQTL analysis configurations
├── gwas/              # GWAS pipeline configurations
├── life_events/       # Life events module configurations
├── longread/          # Long-read sequencing configurations
├── multiomics/        # Multi-omics integration configurations
├── ncbi/              # NCBI API configurations
├── networks/          # Network analysis configurations
├── phenotype/         # Phenotype analysis configurations
└── singlecell/        # Single-cell analysis configurations
```

## Data Structures

### Configuration File Types

- **{species}.yaml**: Active workflow configurations for specific species
- **{module}_template.yaml**: Template configurations with documented defaults
- **config_base/*.yaml**: Base configurations inherited by specific configs

### Environment Variable Prefixes

| Module | Prefix | Example |
|--------|--------|---------|
| RNA/Amalgkit | `AK_` | `AK_THREADS=16` |
| GWAS | `GWAS_` | `GWAS_WORK_DIR=output/gwas` |
| Core | `CORE_` | `CORE_THREADS=8` |
| DNA | `DNA_` | `DNA_WORK_DIR=output/dna` |
| Protein | `PROT_` | `PROT_CACHE_DIR=.cache` |
| Epigenome | `EPI_` | `EPI_PEAK_CALLER=macs2` |
| Ontology | `ONT_` | `ONT_GO_OBO_PATH=data/go.obo` |

### Common Config Fields

```yaml
# Standard configuration fields
work_dir: output/{module}    # Output directory
threads: 8                   # Parallel threads
cache_enabled: true          # Enable caching
log_level: INFO              # Logging verbosity
```

## Interface

### Loading Configurations

```python
from metainformant.core.config import load_mapping_from_file, apply_env_overrides

# Load with environment variable overrides
config = load_mapping_from_file("config/amalgkit/species.yaml")
config = apply_env_overrides(config, prefix="AK")
```

### Configuration Conventions

- Use YAML format exclusively
- Include comments explaining non-obvious options
- Provide sensible defaults for all fields
- Document required vs optional fields
- Use relative paths that resolve from repository root

### Creating New Configurations

1. Start from appropriate template in module subdirectory or config_base/
2. Override only necessary fields
3. Add comments for custom settings
4. Validate with module-specific config loader
