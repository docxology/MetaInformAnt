# Specification: Amalgkit configuration

## Scope
Configuration templates and selection-rule schemas for the pinned Amalgkit build (0.16.108).

## Architecture
- **Component**: Configuration
- **Location**: `config/amalgkit/`

## Data Structures
- **Format**: YAML
- **Schema**: `metainformant.amalgkit.config`
- **Selection policy**: `select_rules.tsv`, consumed by `amalgkit select`
- **Version contract**: exact build `0.16.108` at upstream commit `1d1566638abe32e79ee262193bb88769c3482e76` (`REQUIRED_AMALGKIT_VERSION`, `AMALGKIT_SOURCE_REVISION`)
