# Specification: Amalgkit configuration

## Scope
Configuration templates and selection-rule schemas for Amalgkit 0.16.60.

## Architecture
- **Component**: Configuration
- **Location**: `config/amalgkit/`

## Data Structures
- **Format**: YAML
- **Schema**: `metainformant.amalgkit.config`
- **Selection policy**: `select_rules.tsv`, consumed by `amalgkit select`
- **Version contract**: exact release `0.16.60`
