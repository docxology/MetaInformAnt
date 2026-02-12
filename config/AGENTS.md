# Agent Directives: config

## Role

Configuration file repository for all METAINFORMANT workflows.

## Directory Structure

- `amalgkit/` - RNA-seq workflow configurations (active species)
- `archive/` - Archived/legacy amalgkit configurations
- `config_base/` - Base configuration templates for amalgkit
- `eqtl/` - eQTL analysis configurations
- `gwas/` - GWAS pipeline configurations
- `life_events/` - Life events module configurations
- `longread/` - Long-read sequencing configurations
- `multiomics/` - Multi-omics integration configurations
- `ncbi/` - NCBI API configurations
- `networks/` - Network analysis configurations
- `phenotype/` - Phenotype analysis configurations
- `singlecell/` - Single-cell analysis configurations

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
