# Personal AI Infrastructure (PAI) - amalgkit

## 🧠 Context & Intent

- **Path**: `config/amalgkit/`
- **Purpose**: YAML configurations for amalgkit RNA-seq transcript quantification workflows
- **Domain**: config → bioinformatics → RNA-seq

## 🏗️ Virtual Hierarchy

- **Type**: Configuration Directory
- **Parent**: `config`
- **Consumers**: `scripts/rna/`, `src/metainformant/rna/`

## 📊 Production Status

| Config | Samples | Status |
|--------|---------|--------|
| `amalgkit_pbarbatus_all.yaml` | 110 | ✅ Complete (95 valid) |
| `amalgkit_pbarbatus_25sample.yaml` | 25 | Test |
| `amalgkit_pbarbatus_5sample.yaml` | 5 | Test |

## 📝 Maintenance Notes

- **Dependencies**: `amalgkit>=0.12.20`, `kallisto`, `fastp`
- **Disk Strategy**: Stream-and-clean (minimal persistent footprint)
- **Critical Settings**: `redo: no` for production runs (idempotent)
- **Shared Resources**: Genome/index in `output/amalgkit/shared/`

## 🔄 AI Workflows

- **Modification**: Test changes with 5-sample config first
- **New Species**: Copy `amalgkit_template.yaml`, adjust paths/taxon
- **Recovery**: Use `scripts/rna/recover_missing_parallel.py` for failed samples
- **Documentation**: Update this file and `README.md` when adding configs
