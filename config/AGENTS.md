# Agent Directives: config

## Role
Configuration file repository for all METAINFORMANT workflows.

## Directory Structure
- `amalgkit/` - RNA-seq workflow configurations (active species)
- `archive/` - Archived/legacy amalgkit configurations
- `config_base/` - Base configuration templates for amalgkit
- `gwas/` - GWAS pipeline configurations
- `phenotype/` - Phenotype analysis configurations

## Root Config Files
- `life_events_template.yaml` - Life events module template
- `multiomics_template.yaml` - Multi-omics integration template
- `networks_template.yaml` - Network analysis template
- `singlecell_template.yaml` - Single-cell analysis template
- `ncbi.yaml` - NCBI API configuration

## Rules and Constraints

### Configuration Patterns
All configs support environment variable overrides with module-specific prefixes:
- RNA/Amalgkit: `AK_` prefix
- GWAS: `GWAS_` prefix
- Core: `CORE_` prefix

### Usage
```python
from metainformant.core.config import load_mapping_from_file, apply_env_overrides

config = load_mapping_from_file("config/amalgkit/species.yaml")
config = apply_env_overrides(config, prefix="AK")
```

### File Format
- Use YAML for all configuration files
- Include comments explaining non-obvious options
- Provide sensible defaults
- Document required vs optional fields
